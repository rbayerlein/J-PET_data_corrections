/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*                (C) Copyright 2011 Department of Radiology                      *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/

/*********************************************************************************
*
*			Module Name:		calc_ad.c
*			Revision Number:	2.10
*			Date last revised:	30 November 2011
*			Programmer:			Steven Gillispie
*			Date Originated:	4 August 2011
*
*			Module Overview:	Conversion of MK contributed calc_ad.f to C
*
*			References:			
*
**********************************************************************************
*
*			Global functions defined:
*
*			Global variables defined:	
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		
*
*			Revision date:		
*
*			Revision description:
*
*********************************************************************************/


#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/* TYPES */
enum {								/* Boolean types */
	false			= 0,
	true			= 1
};
typedef int		Boolean;


/* GLOBALS */
#define				kElems		101		/* Maximum allowed elements (+1) */

static int				num_masses = 26;
static int				mass_index[] = {	/* (8 per line) */
		1, 6, 7, 8, 11, 13, 14, 15, 
		16, 17, 20, 26, 29, 32, 35, 39, 
		50, 53, 55, 57, 58, 64, 71, 74, 
		82, 83 };
static double			mass_vec[] = {	/* (8 per line) */
		1.0079, 12.011, 14.0067, 15.9994, 22.9898, 26.98154, 28.0855, 30.97376, 
		32.066, 35.453, 40.08, 55.847, 63.546, 72.59, 79.904, 88.9059, 
		118.71, 126.9045, 132.9054, 138.91, 140.116, 157.25, 174.967, 183.84, 
		207.2, 208.98 };

static double			e, A[kElems];

static double			rho, ff1;
static int				i, ie, iz, n, iztab, ientry, nentry[kElems], neg1;
static char				compound[80], out_file[30];
static Boolean			m_mask[kElems];

static double			anomr[kElems][1001], anomi[kElems][1001];
static double			x[kElems][251], ff[kElems][251], w_frac[kElems];
static double			x1, energy, lambda, lambda2, theta_step = 1.0e-6;
static double			e2wave = 1.239852e-7, t = 0.6652448;
static double			f1, f2, form, pi = 3.1415926536;
static double			theta, cos2_theta;
static long				itheta, ntheta;
static double			*sig, *cos_theta;
static double			*tcos;


/* PROTOTYPES */
void trim(	char			*str, 
			int				*n );
double interp(	double		z, 
				double*		x, 
				double*		y );
void do_calc( void );



/*********************************************************************************
*
*		Name:			trim
*
*		Summary:		
c    		calculate number of characters in str after leading and trailing
c     		s are removed; call that n.  remove the leading s from str.
*
*		Arguments:
*			char			   *str				- Character string.
*			int				   *n				- Number of valid characters.
*
*		Function return: None.
*
*********************************************************************************/

void trim(	char			*str, 
			int				*n )

{
	/*     Calculate number of characters in str after leading and trailing
	       s are removed; call that n.  Remove the leading s from str.  */
	
	int			slen;
	int			i, j, k;
	
	
	slen = strlen(str);
	if (slen == 0) {
		*n = 0;
	}
	else {
		/* find i, the index of the leading non */
		for (i=0; i<slen; i++) {
			if (! isspace(str[i])) {
				break;
			}
		}
		if (i >= slen) {
			*n = 0;
		}
		else {
			/* find j, the index of the last non */
			for (j=slen-1; j>=i; j--) {
				if (! isspace(str[j])) {
					break;
				}
			}
			
			*n = j - i + 1;		/* number of characters in between */
			
			if (i != 0) {
				/* left shift str by i */
				for (k=0; k<*n; k++) {
					str[k] = str[i];
					i++;
				}
			}
			
			str[*n] = '\0';
		}
	}
}


/*********************************************************************************
*
*		Name:			interp
*
*		Summary:		Do log-log interpolation
*
*		Arguments:
*			double			z			- True value of independent variable
*			double*			x			- List of known independent values
*			double*			y			- List of known dependent values
*
*		Function return: double -- Interpolated value for y(z) given y(x)
*
*********************************************************************************/

double interp(	double		z, 
				double*		x, 
				double*		y )

{
	double		result;						/* Result of function */
	double		zlog, x1, x2, y1, y2;
	int			i;
	
	
	if (z <= x[2]) {
		result = y[2];
	}
	else {
		i = 3;
		while (x[i] < z) {
			i++;
		}
		
		if (z == x[i-1]) {
			result = y[i-1];
		}
		else if (z == x[i]) {
			result = y[i];
		}
		else {
			x1 = log(x[i-1]);
			x2 = log(x[i]);
			y1 = log(y[i-1]);
			y2 = log(y[i]);
			zlog = log(z);
			
			result = exp( ( (zlog-x1)*y2 + (x2-zlog)*y1)/(x2-x1) );
		}
	}
	
	return ( result );
}


/*********************************************************************************
*
*		Name:			do_calc
*
*		Summary:		Do calculations
*
*		Arguments:		None.
*
*		Function return: None.
*
*********************************************************************************/

void do_calc( void )

{
	double		sig_tot, percent, *sig_cum;
	int			ipercent;
	int			iz;
	int			n;
	
	
	sig_cum = calloc( ntheta+1, sizeof(double) );
	
	trim(compound, &n);
	printf( "%d keV, material = %s\n", ie, compound );
	energy = ie;
	
	lambda = e2wave/energy;
	lambda2 = lambda * lambda;
	
	for (itheta=0; itheta<=ntheta; itheta++) {
		sig[itheta] = 0.0;
	}
	
	for (itheta=1; itheta<=ntheta; itheta++) {
		x1 = (1.0 - cos_theta[itheta])/lambda2;
		
		for (iz=1; iz<kElems; iz++) {
			if (! m_mask[iz]) {
				continue;
			}
			
			f1 = anomr[iz][ie];
			f2 = anomi[iz][ie];
			
			form = interp(x1, x[iz], ff[iz]);
			sig[itheta] += ( w_frac[iz] * ((form+f1)*(form+f1) + (f2*f2)) );
		}
		
		/*  UP TO CONSTANT FACTOR WHICH WE DONT CARE ABOUT  */
		
		sig[itheta] *= tcos[itheta];
	}
	
	sig_tot = 0.0;
	for (itheta=1; itheta<=ntheta; itheta++) {
		sig_tot += sig[itheta];
	}
	sig_cum[1] = sig[1];
	for (itheta=2; itheta<=ntheta; itheta++) {
		sig_cum[itheta] = sig_cum[itheta-1] + sig[itheta];
	}
	
	/* Renormalize */
	for (itheta=1; itheta<=ntheta; itheta++) {
		sig_cum[itheta] /= sig_tot;
	}
	
	itheta = 1;
	for (ipercent=1; ipercent<=100; ipercent++) {
		percent = ipercent/100.0;
		while (sig_cum[itheta] < percent) {
			itheta++;
		}
		
		if (itheta >= ntheta) {
			cos_theta[itheta] = -1.0;
		}
		printf( "\t%15.6E\n", cos_theta[itheta] );
	}
	
	if (sig_cum) free( sig_cum );
}



/*********************************************************************************
*
*		Name:			main
*
*		Summary:		Main program.
*
*		Arguments:		None.
*
*		Function return: None.
*
*********************************************************************************/

int main( void )

{
	int				i, j;
	double			wf;
	int				matched;
	double			sumw_frac;
	FILE*			infileFF8;
	FILE*			infileAnom9;
	
	
	ntheta = ceil(2.0/theta_step);
	sig = calloc( ntheta+1, sizeof(double) );
	cos_theta = calloc( ntheta+1, sizeof(double) );
	tcos = calloc( ntheta+1, sizeof(double) );
	if (! (sig && cos_theta && tcos)) {
		/* Memory allocation failure */
		printf( "Could not allocate main memory variables\n" );
		exit( -1 );
	}
	
	for (i=0; i<kElems; i++) {
		m_mask[i] = false;
	}
	for (i=0; i<kElems; i++) {
		for (j=0; j<250; j++) {
			ff[i][j] = 0.0;
			x[i][j] = 0.0;
		}
	}
	for (i=0; i<kElems; i++) {
		for (j=0; j<=1000; j++) {
			anomr[i][j] = 0.0;
			anomi[i][j] = 0.0;
		}
	}
	iz = 0;
	for (i=0; i<kElems; i++ ) {
		A[i] = 0.0;
		if (i == mass_index[iz]) {
			A[i] = mass_vec[iz];
			iz++;
		}
	}
	
	infileFF8 = fopen( "ff.dat", "r" );
	if ( ! infileFF8 ) {
		printf( "\tCould not open ff.dat\n" );
		exit( 8 );
	}
	infileAnom9 = fopen( "anom.dat", "r" );
	if ( ! infileAnom9 ) {
		printf( "\tCould not open anom.dat\n" );
		exit( 9 );
	}
	
	for (iz=1; iz<kElems; iz++ ) {
		ientry = 0;
		iztab = iz;
		
		while (iztab == iz) {
			fscanf( infileFF8, "%d", &iztab );
			if (iztab < 0) {
				break;
			}
			fscanf( infileFF8, "%lf %lf", &x1, &ff1 );
			ientry++;
			if (ientry > 250) {
				printf( "\tToo many values\n" );
				exit( 250 );
			}
			x[iz][ientry] = x1;
			ff[iz][ientry] = ff1;
		}
		nentry[iz] = ientry;
		
		fscanf( infileAnom9, "%d", &iztab );
		if (iztab != iz) {
			printf( "\tOut of order on iztab!\n" );
			exit( 900 );
		}
		for (ie=1; ie<=1000; ie++ ) {
			fscanf( infileAnom9, "%lf %lf %lf", 
							&energy, &(anomr[iztab][ie]), &(anomi[iztab][ie]) );
		}
		fscanf( infileAnom9, "%d", &neg1 );
		if (neg1 != -1) {
			printf( "\tOut of order on neg1!\n" );
			exit( 901 );
		}
	}
	fclose( infileFF8 );
	fclose( infileAnom9 );
	
	do {
		/*printf( "Input density (rho, g/cc), material name:\n" );*/
		scanf( "%lf,%s", &rho, compound );
		trim( compound, &n );
		strcpy( out_file, compound );
		strncat( out_file, ".ad", 20 );
		
		/*printf( "Input atomic number (-1 to stop), weight fraction:\n" );*/
		while (iz > 0) {
			scanf( "%d,%lf", &iz, &wf );
			if (iz <= 0) {
				break;
			}
			if (iz > 100) {
				printf( "\tAtomic number too large\n" );
				exit( 100 );
			}
			else {
				w_frac[iz] = wf;
			}
			
			m_mask[iz] = true;
			
			matched = 0;
			for (i=0; i<num_masses; i++) {
				if (iz == mass_index[i]) {
					matched = 1;
					break;
				}
			}
			if (matched == 0) {
				printf( "\tElement not in table\n" );
				exit( 200 );
			}
		}
		
		sumw_frac = 0.0;
		for (i=1; i<kElems; i++) {
			sumw_frac += w_frac[i];
		}
		if (fabs(1.0 - sumw_frac) > 0.01) {
			printf( "Weight fractions do not add up -- try again\n" );
			for (i=0; i<kElems; i++ ) {
				w_frac[i] = 0.0;
			}
			iz = 1;
		}
	} while (fabs(1.0 - sumw_frac) > 0.01);
	
	for (itheta=1; itheta<=ntheta; itheta++) {
		theta = 1.0 - ( itheta * theta_step );
		if (theta < -1.0) {
			theta = -1.0;
		}
		cos_theta[itheta] = theta;
		tcos[itheta] = 1.0 + (theta * theta);
	}
	
	for (ie=1; ie<=150; ie++) {
		do_calc();
	}
	
	for (ie=160; ie<=300; ie+=10) {
		do_calc();
	}
	
	for (ie=400; ie<=1000; ie+=100) {
		do_calc();
	}
	
	if (sig) free( sig );
	if (cos_theta) free( cos_theta );
	if (tcos) free( tcos );
	
	/*fprintf( stderr, "Normal exit\n" );*/
	
	return( 0 );
}

/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 2011-2012 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/

/*********************************************************************************
*
*			Module Name:		mat.c
*			Revision Number:	2.5
*			Date last revised:	17 January 2012
*			Programmer:			Steven Gillispie
*			Date Originated:	1 August 2011
*
*			Module Overview:	Conversion of MK contributed mat.f to C
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


/* PROTOTYPES */
void trim(	char			*str, 
			int				*n );



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
	double			e, sig_comp, sig_pe, sig_coh, A[kElems], factor;
	double			mu_tot[1000],prob_comp[1000],prob_pe[1000],prob_coh[1000];
	double			w_frac[kElems], avagadro = 0.602, rho, Aeff, Zeff;
	int				i, ie, iz, n;
	char			compound[20], out_file[20];
	double			wf;
	int				matched;
	double			sumw_frac;
	FILE*			infile;
	FILE*			outfile;
	
	
	sig_comp = 0;
	sig_pe = 0;
	sig_coh = 0;
	iz = 0;
	for (i=0; i<kElems; i++ ) {
		A[i] = 0.0;
		if (i == mass_index[iz]) {
			A[i] = mass_vec[iz];
			iz++;
		}
	}
	for (i=0; i<kElems; i++ ) {
		w_frac[i] = 0.0;
	}
	for (i=0; i<1000; i++ ) {
		mu_tot[i] = 0.0;
		prob_comp[i] = 0.0;
		prob_pe[i] = 0.0;
		prob_coh[i] = 0.0;
	}
	iz = 1;
	i = 1;
	Aeff = 0;
	Zeff = 0;
	
	do {
		fprintf( stderr, "Input density (rho, g/cc), material name:\n" );
		scanf( "%lf,%s", &rho, compound );
		trim( compound, &n );
		strcpy( out_file, compound );
		
		fprintf( stderr, "Input atomic number (-1 to stop), weight fraction:\n" );
		while (iz > 0) {
			scanf( "%d,%lf", &iz, &wf );
			if (iz <= 0) {
				break;
			}
			if (iz > 100) {
				fprintf( stderr, "\tAtomic number too large\n" );
				exit( 100 );
			}
			else {
				w_frac[iz] = wf;
			}
			
			matched = 0;
			for (i=0; i<num_masses; i++) {
				if (iz == mass_index[i]) {
					matched = 1;
					break;
				}
			}
			if (matched == 0) {
				fprintf( stderr, "\tElement not in table\n" );
				exit( 200 );
			}
			
			Zeff = Zeff  +  w_frac[iz] * iz;
		}
		
		for (i=1; i<kElems; i++) {
			if (w_frac[i] > 0) {
				Aeff = Aeff  +  (w_frac[i]*i/A[i]);
			}
		}
		
		Aeff = Zeff/Aeff;
		
		sumw_frac = 0.0;
		for (i=1; i<kElems; i++) {
			sumw_frac += w_frac[i];
		}
		if (fabs(1.0 - sumw_frac) > 0.01) {
			fprintf( stderr, "Weight fractions do not add up -- try again\n" );
			for (i=0; i<kElems; i++ ) {
				w_frac[i] = 0.0;
			}
			iz = 1;
		}
	} while (fabs(1.0 - sumw_frac) > 0.01);
	
	infile = fopen( "interact.dat", "r" );
	if ( ! infile ) {
		fprintf( stderr, "File interact.dat is missing from current directory\n" );
	}
	else {
		for (i=0; i<mass_index[num_masses-1]; i++) {
			fscanf( infile, "%d", &iz );
			if (w_frac[iz] > 0) {
				for (ie=0; ie<1000; ie++) {
					fscanf( infile, "%lf %lf %lf %lf", 
						&e, &sig_pe, &sig_comp, &sig_coh );
					
					factor = w_frac[iz] / A[iz];
					
					prob_pe[ie] += factor * sig_pe;
					prob_comp[ie] += factor * sig_comp;
					prob_coh[ie] += factor * sig_coh;
				}
			}
			else {
				for (ie=0; ie<1000; ie++) {
					fscanf( infile, "%lf %lf %lf %lf", 
						&e, &sig_pe, &sig_comp, &sig_coh );
				}
			}
		}
		fclose( infile );
		
		for (i=0; i<1000; i++) {
			mu_tot[i] = prob_pe[i] + prob_comp[i] + prob_coh[i];
			
			if (mu_tot[i] > 0.0) {
				prob_pe[i] /= mu_tot[i];
				prob_comp[i] /= mu_tot[i];
				prob_coh[i] /= mu_tot[i];
			}
			
			/* Normalize mu to proper units (1/cm) */
			mu_tot[i] *= rho * avagadro;
		}
		
		outfile = fopen( out_file, "w" );
		fprintf( outfile, "%s    ", compound );
		fprintf( outfile, "Density = %lf    Aeff = %lf    Zeff = %lf\n", 
							rho, Aeff, Zeff );
		for (i=0; i<1000; i++) {
			fprintf( outfile, "%8.1lf %15.5lf %15.5lf %15.5lf\n", 
				(float)(i+1), mu_tot[i], 1.0-prob_pe[i], 
				prob_comp[i]/(prob_coh[i] + prob_comp[i]) );
		}
		fclose( outfile );
	}
	
	return( 0 );
}

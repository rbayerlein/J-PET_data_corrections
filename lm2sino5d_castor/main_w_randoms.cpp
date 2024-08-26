// ************************************************************************
// Reimund Bayerlein: rbayerlein@ucdavis.edu, June, 2024
// 
// conversion from CASTOR listmode files into sinograms for jPET
// ************************************************************************

// #######################################################################################
// Parameters for 1x8 rebinning (transaxial x axial):
/*
#define SIZEOF_SINO 207*156
#define N_TX_BLK 312
#define N_TX_CRYS 312
#define N_AX_CRYS_WGAP 200
#define N_TX_CRYS_PER_BLK 1
#define N_AX_CRYS_PER_BLK 8 
#define N_AX_BLK 25
*/
// #######################################################################################
// Parameters for 4x8 rebinning (transaxial x axial):
#define SIZEOF_SINO 51*39
#define N_TX_BLK 78
#define N_TX_CRYS 312
#define N_AX_CRYS_WGAP 200
#define N_TX_CRYS_PER_BLK 4
#define N_AX_CRYS_PER_BLK 8 
#define N_AX_BLK 25

// #######################################################################################
//#define Debug_mode // comment out if do not debug the code
#define N_TX_CRYS_PER_MODULE 24
#define N_AX_CRYS_PER_MODULE 40
#define N_TOF 20 	// 300 ps per bin (as of June 18, 2024)
#define TOF_BIN_WIDTH 300 // given in ps
#define N_TIME_BIN_PER_TOF 1
#define SYMMETRIC_REBINNING 0	// 0: no, 1: yes;
#define COIN_TIME_WINDOW 3 // in ns
#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

struct Lut {
	short nv,nu;
};

struct Lm {
	uint32_t time;
	float tof;
	uint32_t c1, c2;	// castor crystal ids of crystals one and two
};

int main(int argc, char* argv[]) {

	if(SYMMETRIC_REBINNING ==  1 && float(N_TIME_BIN_PER_TOF % 2) == 0){
		cout << "Wrong SYMMETRIC_REBINNING usage: For symmetric rebinning, number of time bins per tof (N_TIME_BIN_PER_TOF) must be odd!" << endl;
		exit(1);
	}
	// check arguments
	if (argc != 5) {
		cout << "Usage: " << argv[0] << " [fname_in]  [fname_out]  [fname_lut]  [total_scan_time_seconds]" << endl;
		exit(1);
	}

	cout << "Starting " << argv[0] << " ..." << endl;

//get total sino size
	long int n_ax_squared = N_AX_BLK * N_AX_BLK;
	long int TOTAL_SINO_SIZE  = SIZEOF_SINO*n_ax_squared*N_TOF;

// read block sino lut
	string fname_lut = string(argv[3]);
	FILE* pfile_lut = fopen(fname_lut.c_str(), "rb");
	if (pfile_lut == NULL) {
		cout << fname_lut << " cannot be opened." << endl;
		exit(1);
	}
	Lut* plut = new Lut[SIZEOF_SINO]; // block sino lut
	cout << "Reading LUT of size " << SIZEOF_SINO << endl;
	fread(plut, sizeof(Lut), SIZEOF_SINO, pfile_lut);
	cout << "Done with LUT "<< endl;

// set up block sino reverse lut
	vector< vector<int> > blk_idx(N_TX_BLK, vector<int>(N_TX_BLK, -1)); // initialize values to -1 to skip the Out-of-bound sinogram index
	for (int i = 0; i < SIZEOF_SINO; i++) {
		blk_idx[plut[i].nv-1][plut[i].nu-1] = i;
	}
    cout << "BLK matrix size has  " << blk_idx.size() << "rows,and "<< blk_idx[0].size()<< "cols,\n each entry corresponds to blk# of LORA and blk#LORB in the transaxial dimension(range from 0~"<<N_TX_BLK-1<< ")" << endl;
	cout << "Done with BLK index "<< endl;

// set up module-based event rate table
	short N_AX_BLK_PER_MODULE = N_AX_CRYS_PER_MODULE/N_AX_CRYS_PER_BLK;
	short N_TX_BLK_PER_MODULE = N_TX_CRYS_PER_MODULE/N_TX_CRYS_PER_BLK;
	short N_MODULES_AX = N_AX_CRYS_WGAP/N_AX_CRYS_PER_MODULE;
	short N_MODULES_TX = N_TX_CRYS/N_TX_CRYS_PER_MODULE;
	float mod_prompt_rate[N_MODULES_TX][N_MODULES_AX] = {0};

// declare input file
	string fname_in = string(argv[1]);
	FILE* pfile_in = fopen(fname_in.c_str(), "rb");
	if (pfile_in == NULL) {
		cout << fname_in << " cannot be opened." << endl;
		exit(1);
	}else {
		cout << "done infile definition" << endl;
	}

// declare output file
	string fname_out = string(argv[2]);
	FILE* pfile_out = fopen(fname_out.c_str(), "wb");
	if (pfile_out == NULL) {
		cout << fname_out << " cannot be opened." << endl;
		exit(1);
	}else {
		cout << "done outfile definition" << endl;
	}
// total scan time (frame length)
	float frame_length = atof(argv[4]);
	cout << "total frame duration (sec): " << frame_length << endl;
// MAIN PROGRAM START 

	int read_count_total = 0;
	int read_count_outof_traxial_tof=0;  // for debugging
	Lm* plm = new Lm;	//list mode event
	unsigned short* psino_blk_out = new unsigned short[TOTAL_SINO_SIZE]; // block TOF sino for writing
	cout << "size of block sino 5d sinogram: " << TOTAL_SINO_SIZE << endl;

	int tof_offset = floor(N_TOF/2);	// = 10
	cout << "tof offset: " << tof_offset << endl;

	while (!feof(pfile_in)){
		int read_count = fread(plm, sizeof(Lm), 1, pfile_in);

		// get axial and transaxial coordinates
		int txBiA = (plm->c1 % N_TX_CRYS) / N_TX_CRYS_PER_BLK;
		int txBiB = (plm->c2 % N_TX_CRYS) / N_TX_CRYS_PER_BLK;
		int axBiA = floor(plm->c1 / N_TX_CRYS) / N_AX_CRYS_PER_BLK;
		int axBiB = floor(plm->c2 / N_TX_CRYS) / N_AX_CRYS_PER_BLK;
		int axMiA = axBiA / N_AX_BLK_PER_MODULE;
		int axMiB = axBiB / N_AX_BLK_PER_MODULE;
		int txMiA = txBiA / N_TX_BLK_PER_MODULE;
		int txMiB = txBiB / N_TX_BLK_PER_MODULE;
		int TOF_AB;
		int TOF_AB_tmp = floor(plm->tof / TOF_BIN_WIDTH);
		if (SYMMETRIC_REBINNING==0){
			TOF_AB = TOF_AB_tmp;
		}else {
			TOF_AB = TOF_AB_tmp+floor(N_TIME_BIN_PER_TOF/2);
		}
		if (TOF_AB>=0){
			TOF_AB = TOF_AB / N_TIME_BIN_PER_TOF;
		}else{
			if(N_TIME_BIN_PER_TOF > 1) TOF_AB -= N_TIME_BIN_PER_TOF;
			TOF_AB = TOF_AB / N_TIME_BIN_PER_TOF;
		}		
		if (TOF_AB+tof_offset > N_TOF-1 || TOF_AB+tof_offset < 0) {cout << "filtered event: tof out of range!" << endl; continue;}

		read_count_total++;

	// Fill Prompts Rate Array
		mod_prompt_rate[txMiA][axMiA] += 1.0;
		mod_prompt_rate[txMiB][axMiB] += 1.0;

	// calculate transaxial sinogram index
		int idx_tx_blk = blk_idx[txBiA][txBiB];
		int idx_tx_blk_reverse = blk_idx[txBiB][txBiA]; 

		int ind_blk_sino;

        if (idx_tx_blk!=-1){
         	ind_blk_sino = idx_tx_blk + SIZEOF_SINO * axBiA 
         	+ SIZEOF_SINO * N_AX_BLK * axBiB+ SIZEOF_SINO * N_AX_BLK * N_AX_BLK*(TOF_AB+tof_offset);} // 5-D sinogram 

        else if(idx_tx_blk_reverse!=-1){   
        // swap event A and event B's  position and reverse tof sign
	        TOF_AB=-TOF_AB-1;
	        idx_tx_blk=idx_tx_blk_reverse;
	        ind_blk_sino = idx_tx_blk + SIZEOF_SINO * axBiA 
	        + SIZEOF_SINO * N_AX_BLK * axBiB+ SIZEOF_SINO * N_AX_BLK * N_AX_BLK*(TOF_AB+tof_offset);
		}

		if (ind_blk_sino < 0 || ind_blk_sino >(TOTAL_SINO_SIZE-1)){
			read_count_outof_traxial_tof++;
			cout << "transaxial coordinate out of bound (" << ind_blk_sino << ") for event number " << read_count_total << endl;
			cout << "TOF: " << TOF_AB << endl;
		}else{
			psino_blk_out[ind_blk_sino] += 1;
		}

	}// end of while loop over all events

	printf("Total event is %d, %d out of our transaxial or TOF FOV\n", read_count_total,read_count_outof_traxial_tof);

// create randoms sinogram

	// divide the module array by the measurement time to get the rate
	for (int i = 0; i < N_MODULES_TX; ++i)
	{
		for (int j = 0; j < N_MODULES_AX; ++j)
		{
			mod_prompt_rate[i][j] = float(mod_prompt_rate[i][j]/frame_length);
		}
	}

	// define output file name
	string fname_out_r = string(argv[2]);
	fname_out_r.append("_r");
	FILE* pfile_out_r = fopen(fname_out_r.c_str(), "wb");
	if (pfile_out_r == NULL) {
		cout << fname_out_r << " cannot be opened." << endl;
		exit(1);
	}else {
		cout << "randoms file opened: " << fname_out_r << endl;
	}

	// define sinogram to be filled
	float* psino_blk_out_r = new float[TOTAL_SINO_SIZE/N_TOF];

	// go through every entry in that sinogram and exract the 2 blocks
	for (int i = 0; i < TOTAL_SINO_SIZE/N_TOF; ++i)
	{
		// get tof bin
		short tof = i/(N_AX_BLK*N_AX_BLK*SIZEOF_SINO);
		// extract michelogram
		int mich = i-tof*(N_AX_BLK*N_AX_BLK*SIZEOF_SINO);
		// get axial block B
		short axBiB = mich/(N_AX_BLK*SIZEOF_SINO);
		// extract remaining michelogram to get axial block A
		int mich_rem = mich - axBiB*(N_AX_BLK*SIZEOF_SINO);
		short axBiA = mich_rem/(SIZEOF_SINO);
		// extract sinogram coordinate
		int tx_blk_idx = mich_rem-axBiA*SIZEOF_SINO;
		// get transaxial blocks A and B using plut
		short txBiA = plut[tx_blk_idx].nv-1;
		short txBiB = plut[tx_blk_idx].nu-1;
		// get modules
		short txMiA = txBiA/N_TX_BLK_PER_MODULE;
		short txMiB = txBiB/N_TX_BLK_PER_MODULE;
		short axMiA = axBiA/N_AX_BLK_PER_MODULE;
		short axMiB = axBiB/N_AX_BLK_PER_MODULE;

//		cout << "i: " << i << "; tof: " << tof << "; axBiA,axBiB: " << axBiA << " " << axBiB << "; txBiA,txBiB: " << txBiA << " " << txBiB << endl;

		// calculate randoms rate per module pair
		float rr = mod_prompt_rate[txMiA][axMiA]*mod_prompt_rate[txMiB][axMiB]*2*COIN_TIME_WINDOW*0.000000001;
		if (rr<0)
		{
			cout << "Error detected: negative randoms rate!" << endl;
			cout << "i: " << i << "; tof: " << tof << "; txMiA,axMiA: " << txMiA << " " << axMiA << "; txMiB,axMiB: " << txMiB << " " << axMiB << endl;
			cout << mod_prompt_rate[txMiA][axMiA] << "\t" << mod_prompt_rate[txMiB][axMiB] << endl;
			cout << rr << endl;
			cout << "========" << endl;
		}
		// and per block pair
		float rr_blk = rr/(N_AX_BLK_PER_MODULE*N_AX_BLK_PER_MODULE*N_TX_BLK_PER_MODULE*N_TX_BLK_PER_MODULE);
		psino_blk_out_r[i]=rr_blk;
	}

	// save the output sinograms and delete everything that is not used in the cleanup
	cout << "Save sinos... "<< endl;
	fwrite(psino_blk_out, sizeof(unsigned short), TOTAL_SINO_SIZE, pfile_out);
	fwrite(psino_blk_out_r, sizeof(float), TOTAL_SINO_SIZE/N_TOF, pfile_out_r);
	cout << "Done save sinos." << endl;

	// cleanup

	cout << "cleaning up ...";
	delete[] plut;
	delete[] psino_blk_out;
	delete plm;
	fclose(pfile_lut);
	fclose(pfile_in);
  	fclose(pfile_out);
  	fclose(pfile_out_r);
    cout << "done." << endl;
	return 0;
}

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
	if (argc != 4) {
		cout << "Usage: " << argv[0] << " [fname_in]  [fname_out]  [fname_lut]" << endl;
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

	// save the output sinograms and delete everything that is not used in the cleanup
	cout << "Save sinos... "<< endl;
	fwrite(psino_blk_out, sizeof(unsigned short), TOTAL_SINO_SIZE, pfile_out);
	cout << "Done save sinos." << endl;

// cleanup

	cout << "cleaning up ...";
	delete[] plut;
	delete[] psino_blk_out;
	delete plm;
	fclose(pfile_lut);
	fclose(pfile_in);
  	fclose(pfile_out);
    cout << "done." << endl;
	return 0;
}

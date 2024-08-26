// ************************************************************************
// Reimund Bayerlein: rbayerlein@ucdavis.edu, June, 2024
// 
// conversion from listmode files into sinograms for jPET
// ************************************************************************
#define SIZEOF_SINO (51*39) // used to be: 51*39 207*156
#define N_TX_BLK 78 		// used to be: 312
#define N_AX_CRYS_WGAP 200
#define N_TX_CRYS_PER_BLK 4
#define N_AX_CRYS_PER_BLK 8 
#define N_AX_BLK 25
//#define Debug_mode // comment out if do not debug the code
#define N_TOF 20 
#define N_TIME_BIN_PER_TOF 1
#define SYMMETRIC_REBINNING 0	// 0: no, 1: yes;
#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

struct Lut {
	short nv,nu;
};

struct Lm {
	short txIDA, axIDA, txIDB, axIDB, tof;
};

int main(int argc, char* argv[]) {

	if(SYMMETRIC_REBINNING ==  1 && float(N_TIME_BIN_PER_TOF % 2) == 0){
		cout << "Wrong SYMMETRIC_REBINNING usage: For symmetric rebinning, number of time bins per tof (N_TIME_BIN_PER_TOF) must be odd!" << endl;
		exit(1);
	}
	// check arguments
	if (argc != 4) {
		cout << "Usage: " << argv[0] << " [fname_in] [fname_out] [fname_lut]" << endl;
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
	//	cout << "Starting " <<plut[i].nv-1 <<" ..."<<plut[i].nu -1 << " ..." <<i<< endl;  //Print transaxial ID (zhoujian's       'scanner.getDefaultSinogramCrystalPairs' format, remember to reverse TOF bin when swap LOR orientation)
		blk_idx[plut[i].nv-1][plut[i].nu-1] = i;
	//	blk_idx[plut[i].nu][plut[i].nv] = i;
	}
    cout << "BLK matrix size has  " << blk_idx.size() << "rows,and "<< blk_idx[0].size()<< "cols,\n each entry corresponds to blk# of LORA and blk#LORB in the transaxial dimension(range from 0~"<<N_TX_BLK-1<< ")" << endl;
	cout << "Done with BLK index "<< endl;


	// declare input file
	string fname_in = string(argv[1]);
	FILE* pfile_in = fopen(fname_in.c_str(), "rb");
	if (pfile_in == NULL) {
		cout << fname_in << "cannot be opened." << endl;
		exit(1);
	}else {
		cout << "done infile definition" << endl;
	}

	// declare output file
	string fname_out = string(argv[2]);
	FILE* pfile_out = fopen(fname_out.c_str(), "wb");
	if (pfile_out == NULL) {
		cout << fname_out << "cannot be opened." << endl;
		exit(1);
	}else {
		cout << "done outfile definition" << endl;
	}

	// main loop
        int read_count_outof_traxial_tof=1; //Because of zhoujian's  SinogramCrystalPairs has larger(22) accepted BLK Difference,which covers smaller FOV in transaxial plane, and TOF FOV is different
        int read_count_total=1; //
	Lm* plm = new Lm; // list-mode event
	float* psino_blk = new float[TOTAL_SINO_SIZE];// block TOF sino during conversion
	unsigned short* psino_blk_out = new unsigned short[TOTAL_SINO_SIZE]; // block TOF sino for writing
	cout << "size of block sino 5d sinogram: " << TOTAL_SINO_SIZE << endl;
	while (!feof(pfile_in)) {

		int read_count = fread(plm, sizeof(Lm), 1, pfile_in); 

		for (int i = 0; i < read_count; i++) {
			int txBiA = plm->txIDA / N_TX_CRYS_PER_BLK;
			int txBiB = plm->txIDB / N_TX_CRYS_PER_BLK;
			int TOF_AB;
			if (SYMMETRIC_REBINNING==0){
				TOF_AB = plm->tof;
			}else {
				TOF_AB = plm->tof+floor(N_TIME_BIN_PER_TOF/2);
			}
			if (TOF_AB>=0){
				TOF_AB = TOF_AB / N_TIME_BIN_PER_TOF;
			}else{
				TOF_AB -= N_TIME_BIN_PER_TOF;
				TOF_AB = TOF_AB / N_TIME_BIN_PER_TOF;
			}
			short tof_offset = floor(N_TOF/2);
			if (TOF_AB+tof_offset > N_TOF-1 || TOF_AB+tof_offset < 0) {/*cout << "filtered event: tof out of range!" << endl;*/ continue;}
						// 
            read_count_total++;

			
                           #ifdef Debug_mode
                           if (i/1000000==1){
                           printf("Current LOR_A w ax_gap is: %d,LOR_B w ax_gap is %d,Tof bin is %d \n", axBiA_gap,axBiB_gap,plm->tof);
                           printf("Current LOR_A minus %d gap,LOR_B minus %d gap\n", axBiA_gap/N_AX_Unit,axBiB_gap/N_AX_Unit);
                           printf("Current LOR_A w/o ax_gap is: %d,LOR_B w/o ax_gap is %d\n", axBiA_gap-axBiA_gap/N_AX_Unit,axBiB_gap-axBiB_gap/N_AX_Unit);
                           printf("Current ABLK of LOR_A is: %d,ABLK of LOR_B  is %d, TOFBLK is %d \n \n", axBiA,axBiB,TOF_AB);}
                           #endif
 			
 			int axBiA = plm->axIDA / (N_AX_CRYS_PER_BLK);
			int axBiB = plm->axIDB / (N_AX_CRYS_PER_BLK);

                        //cout << "Reading Listmode "<<txBiA<<"B "<<txBiB<< endl;
                        
			int idx_tx_blk = blk_idx[txBiA][txBiB]; // transaxial sinogram index ,78x78
			int idx_tx_blk_reverse = blk_idx[txBiB][txBiA]; // transaxial sinogram index ,78x78,reverse direction
                        int ind_blk_sino;
                        if (idx_tx_blk!=-1){
                         	ind_blk_sino = idx_tx_blk + SIZEOF_SINO * axBiA
				+ SIZEOF_SINO * N_AX_BLK * axBiB+ SIZEOF_SINO * N_AX_BLK * N_AX_BLK*(TOF_AB+tof_offset);} // 5-D sinogram 

                        else if(idx_tx_blk_reverse!=-1){   
                        // swap event A and event B's  position and reverse tof sign
                        TOF_AB=-TOF_AB;
                        idx_tx_blk=idx_tx_blk_reverse;
                        ind_blk_sino = idx_tx_blk + SIZEOF_SINO * axBiA
				+ SIZEOF_SINO * N_AX_BLK * axBiB+ SIZEOF_SINO * N_AX_BLK * N_AX_BLK*(TOF_AB+tof_offset);
						}
                        //else {read_count_outof_traxial_tof++;}
                       
                        //cout << "Reading ind_blk_sino "<<ind_blk_sino<< endl;
                        #ifdef Debug_mode
                        if (i/10==0){
                        if (ind_blk_sino<0||(TOF_AB+tof_offset)<0||(TOF_AB+tof_offset)>N_TOF-1||idx_tx_blk==(SIZEOF_SINO)||ind_blk_sino>(SIZEOF_SINO*N_AX_BLK*N_AX_BLK*N_TOF-1)){
                        	// changed to NOT exclude last sino bin (SIZEOF_SINO-1), rbayerlein, 01-12-2021
	                cout << "Index Out-of-bound sinogram index detected. Skipping the following coincidence event: "<< endl;
                          printf("Current LOR crystal ID of tranxaxial plane is: %d,crystal ID of tranxaxial plane is %d,Tof bin is %d \n", plm->txIDA,plm->txIDB,plm->tof);
                          printf("Current BLK LORA tranxaxial ID is: %d,BLK LORB tranxaxial ID is %d(R:0~119),BLK Tof bin is %d  (R:0~N_TOF-1)\n", txBiA,txBiB,TOF_AB+tof_offset);
                          printf("Current BLK LORA Axial ID  is: %d, BLK LORB  Axial ID is %d \n", axBiA,axBiB);
                          printf("Current total index of LOR_A is: %d\n\n", ind_blk_sino);
                        }
                        }
                         
                        #endif
                        

			if (ind_blk_sino <0||ind_blk_sino>(TOTAL_SINO_SIZE-1)) { // provides fault tolerance
                       //ind_blk_sino <0 represents crystal-based UIH listmode out-of-bound zhoujian's sinogram,
                       //idx_tx_blk==SIZEOF_SINO*N_AX_BLK*N_AX_BLK*N_TOF-1 represents BLK-based listmodeID covered by zhoujian's sinogram, but some of crystal UIH listmode out-of-bound because of integer division in C is equivalent to matlab 'fix', which make the one last entry of sinogram has very high value .
                        read_count_outof_traxial_tof++;
			}
			else {
				psino_blk[ind_blk_sino]+= 1.0;
		    }
		}

	}

	printf("Total event is %d, %d out of our transaxial or TOF FOV\n", read_count_total,read_count_outof_traxial_tof);

	// round 
	for (int i = 0; i < TOTAL_SINO_SIZE; ++i)
	{
		psino_blk_out[i] = round(psino_blk[i]);
		if(isnan(psino_blk_out[i])){
			cout << "Entry " << i << " is NaN" << endl;
		}
	}

	// write to file
	cout << "Save sino... "<< endl;
	fwrite(psino_blk_out, sizeof(unsigned short), TOTAL_SINO_SIZE, pfile_out);
	cout << "Done save sino." << endl;
	// cleanup


	cout << "cleaning up ...";
	delete[] plut;
	delete[] psino_blk_out;
	delete[] psino_blk;
	delete plm;
	fclose(pfile_lut);
	fclose(pfile_in);
  	fclose(pfile_out);
        // if segmentation fault, check each fclose and negative index

    cout << "done." << endl;
	return 0;
}

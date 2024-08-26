// ************************************************************************
// Reimund Bayerlein: rbayerlein@ucdavis.edu, July, 2024
// 
// creation of 4D randoms-sinograms from list-mode delayed events for jPET
// it is basically a simple list-mode to sinogram converter
// ************************************************************************
#define N_CRYSTALS 62400 	// 13*24*200; total number of crystals in the scanner
#define SIZEOF_SINO 51*39
#define N_TX_CRYS 312
#define N_AX_CRYS_WGAP 200
#define N_TX_CRYS_PER_BLK 4
#define N_AX_CRYS_PER_BLK 8 
#define N_TX_BLK 78
#define N_AX_BLK 25
#define TOTAL_SINO_SIZE (SIZEOF_SINO*N_AX_BLK*N_AX_BLK)
#define COIN_TIME_WINDOW 3 // in ns

#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <algorithm>

using namespace std;

struct Lut {
	short nv,nu;
};

struct CrysPair {
	int cID_1, cID_2;
};

int main(int argc, char* argv[]) {
	cout << "Starting..." << endl;
	// check arguments
	if (argc != 4) {
		cout << "Usage: " << argv[0] << " [fname_in]  [fname_out]  [fname_lut]" << endl;
		exit(1);
	}

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
    std::ifstream infile(fname_in.c_str());
    if (!infile.is_open()) {
        std::cerr << "Error opening file" << std::endl;
        return 1;
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

// declare the randoms block sino
	float* psino_blk_out_r = new float[TOTAL_SINO_SIZE]; // nonTOF block sinogram for writing

// MAIN PROGRAM START 
	cout << "Main program start..." << endl;
	// simple lm to block sino
	// go through data set and sort the delayed events into their corresponding block bin
	CrysPair* cp = new CrysPair; // list-mode event
	int read_count_out_of_bound  = 0;
	int read_count_total = 0;
	string line;
	while (getline(infile, line)) {
        // Find the comma position
        size_t comma_pos = line.find(',');
        // Extract the first number (before the comma)
        string num1_str = line.substr(0, comma_pos);
        // Extract the second number (after the comma)
        string num2_str = line.substr(comma_pos + 1);
        // Convert the strings to integers and fill into the CrysPairs struct
        cp->cID_1 = stoi(num1_str);
        cp->cID_2 = stoi(num2_str);
		
		read_count_total++;
		//get detector block from current crystal ID
		int txBiA = (cp->cID_1 % N_TX_CRYS) / N_TX_CRYS_PER_BLK;
		int axBiA = floor(cp->cID_1 / N_TX_CRYS) / N_AX_CRYS_PER_BLK;
		int txBiB = (cp->cID_2 % N_TX_CRYS) / N_TX_CRYS_PER_BLK;
		int axBiB = floor(cp->cID_2 / N_TX_CRYS) / N_AX_CRYS_PER_BLK;
		// get sinogram block
		int idx_tx_blk = blk_idx[txBiA][txBiB]; // transaxial sinogram index
		int idx_tx_blk_reverse = blk_idx[txBiB][txBiA]; // transaxial sinogram index
		int ind_blk_sino;
		if (idx_tx_blk!=-1){
			ind_blk_sino = idx_tx_blk + SIZEOF_SINO * axBiA + SIZEOF_SINO * N_AX_BLK * axBiB;} // 4-D sinogram 

		else if(idx_tx_blk_reverse!=-1){   
			 // swap event A and event B's  position
    	    idx_tx_blk=idx_tx_blk_reverse;
    	    ind_blk_sino = idx_tx_blk + SIZEOF_SINO * axBiA + SIZEOF_SINO * N_AX_BLK * axBiB;
		}
		if (ind_blk_sino <0||ind_blk_sino>(TOTAL_SINO_SIZE-1)) { // provides fault tolerance
			read_count_out_of_bound++;
			cout << txBiA << "\t" << axBiA << "\t" << cp->cID_1 << "\t" << txBiB << "\t" << axBiB << "\t" << cp->cID_2 << endl;

		}else {
			psino_blk_out_r[ind_blk_sino]+= 1.0;
	    }
	} // end of while

	printf("Total number of events is %d, %d out of our transaxial FOV\n", read_count_total,read_count_out_of_bound);


	// write to file
	cout << "Save sino... "<< endl;
	fwrite(psino_blk_out_r, sizeof(float), TOTAL_SINO_SIZE, pfile_out);
	cout << "Done save sino." << endl;

	// free up memory
	cout << "cleaning up ...";
	delete[] plut;
	delete[] psino_blk_out_r;
	fclose(pfile_lut);
	infile.close();
  	fclose(pfile_out);


    cout << "done." << endl;
	return 0;
}// end of main
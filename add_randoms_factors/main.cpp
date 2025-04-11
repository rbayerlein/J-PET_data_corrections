// ************************************************************************
// Reimund Bayerlein: rbayerlein@ucdavis.edu, June, 2024
// 
// calculation of randoms correction factors for CASToR-based reconstruction of jPET list-mode data
// the script adds the calculated randoms factors to the .Cdf file for CASToR reconstruction
// the script also updates the header file .Cdh
// ************************************************************************

// last updated: July 2024, rbayerlein

// #######################################################################################

// Parameters for 4x8 rebinning (transaxial x axial):
#define SIZEOF_SINO 51*39
#define N_TX_BLK 78
#define N_TX_CRYS 312
#define N_AX_CRYS_WGAP 200
#define N_TX_CRYS_PER_BLK 4
#define N_AX_CRYS_PER_BLK 8 
#define N_AX_BLK 25

#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>

using namespace std;

struct Lut {
	short nv,nu;
};

struct Lm {
	uint32_t time;
	float nf;
	float tof;
	uint32_t c1, c2;	// castor crystal ids of crystals one and two
};

struct Lm_out
{
	uint32_t time;
	float rr;	// randoms rate
	float nf;
	float tof;
	uint32_t c1, c2;	
};

void UpdateHeader(const string& inputFile, const string& outputFile, const string& oldWord, const string& newWord, bool apply_RC);
string getFileName(const string& fullPath);

int main(int argc, char const *argv[])
{
	cout << "starting..." << endl;
	// check arguments
	if (argc != 6) {
		cout << "Usage: " << argv[0] << " [fname_in_castor]  [fname_out_castor]  [rand_sino_in]  [fname_lut]  [frame_length] " << endl;
		exit(1);
	}

// read block sino lut
	string fname_lut = string(argv[4]);
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

// open input randoms sinogram
	// calculate the total sinogram size to read
	int TOTAL_SINO_SIZE = SIZEOF_SINO*N_AX_BLK*N_AX_BLK;	// 4D sinogram w/o TOF bin
	float* rand_sino = new float[TOTAL_SINO_SIZE];

	string fname_sino = string(argv[3]);
	FILE* pfile_sino = fopen(fname_sino.c_str(), "rb");
	if (pfile_sino == NULL)
	{
		cout << fname_sino << " cannot be opened." << endl;
		exit(1);
	}else{
		fread(rand_sino, sizeof(float), TOTAL_SINO_SIZE, pfile_sino);
		cout << "done reading sinogram" << endl;
	}
// total scan time (frame length)
	float frame_length = atof(argv[5]);
	cout << "total frame duration (sec): " << frame_length << endl;

// read event by event
	Lm* plm_in = new Lm;			// input list mode event
	Lm_out* plm_out = new Lm_out;	// output list mode event
	int num_LORs_per_blk_pair = N_TX_CRYS_PER_BLK*N_TX_CRYS_PER_BLK* N_AX_CRYS_PER_BLK* N_AX_CRYS_PER_BLK;
	long read_count_total = 0;
	long read_count_outof_traxial_tof = 0;
	while (!feof(pfile_in)){
		int read_count = fread(plm_in, sizeof(Lm), 1, pfile_in);
		read_count_total+=read_count;
		if(read_count==0) continue;

	// get axial and transaxial coordinates
		int txBiA = (plm_in->c1 % N_TX_CRYS) / N_TX_CRYS_PER_BLK;
		int txBiB = (plm_in->c2 % N_TX_CRYS) / N_TX_CRYS_PER_BLK;
		int axBiA = floor(plm_in->c1 / N_TX_CRYS) / N_AX_CRYS_PER_BLK;
		int axBiB = floor(plm_in->c2 / N_TX_CRYS) / N_AX_CRYS_PER_BLK;

	// calculate transaxial sinogram index
		int idx_tx_blk = blk_idx[txBiA][txBiB];
		int idx_tx_blk_reverse = blk_idx[txBiB][txBiA];

		int ind_blk_sino;

        if (idx_tx_blk!=-1){
         	ind_blk_sino = idx_tx_blk; 	// 4-D sinogram 
         } 
        else if(idx_tx_blk_reverse!=-1){   
	        ind_blk_sino=idx_tx_blk_reverse;
		}
		ind_blk_sino += SIZEOF_SINO * axBiA + SIZEOF_SINO * N_AX_BLK * axBiB;

		float rr;
		if (ind_blk_sino < 0 || ind_blk_sino >(TOTAL_SINO_SIZE-1)){
			read_count_outof_traxial_tof++;
			//cout << "transaxial coordinate out of bound (" << ind_blk_sino << ") for event number " << read_count_total << endl;
			//cout << txBiA << "\t" << txBiB << endl;
			rr=0;
		} else{
	// calculate randoms factor per LOR
			rr = rand_sino[ind_blk_sino] / num_LORs_per_blk_pair / frame_length;
		}	
	// write output lm event
		plm_out->c1=plm_in->c1;
		plm_out->c2=plm_in->c2;
		plm_out->rr=rr;
		plm_out->time=plm_in->time;
		plm_out->tof=plm_in->tof;
		plm_out->nf=plm_in->nf;

		fwrite(plm_out, sizeof(Lm_out), 1, pfile_out);
	}
	cout << "Total number of events out of bounds: " << read_count_outof_traxial_tof << " (" << (float)read_count_outof_traxial_tof/read_count_total*100 << " % of all events)" << endl;
	cout << "Total number of events processed: " << read_count_total << endl;
// save output file 

// close everything
	delete[] plut;
	fclose(pfile_lut);
	fclose(pfile_in);
  	fclose(pfile_out);
  	fclose(pfile_sino);

// ammend the header file
  	cout << "Updating the header file " << endl;
  	string inputFileName = argv[1];

    // Remove file extension and replace with ".Cdh"
    size_t dotPos = inputFileName.find_last_of('.');
    if (dotPos!= std::string::npos) {
        inputFileName.erase(dotPos);
    }
    inputFileName += ".Cdh";

    // Create a different output file name
    string outputFileName = inputFileName.substr(0, inputFileName.size() - 4) + "_w_rand.Cdh";

    // get file names without folder path
    string old_fname = getFileName(argv[1]);
    string new_fname = getFileName(argv[2]);
    cout << old_fname << "\t" << new_fname << endl;

    // make update of the header file
    UpdateHeader(inputFileName, outputFileName, old_fname, new_fname, true);

	cout << "DONE." << endl;
	return 0;
}

// ...ooo000OOO000ooo...ooo000OOO000ooo...ooo000OOO000ooo...ooo000OOO000ooo...ooo000OOO000ooo...ooo000OOO000ooo...

// Additional functions

void UpdateHeader(const string& inputFile, const string& outputFile, const string& oldWord, const string& newWord, bool apply_RC) {
    // Open the input file for reading
    ifstream inFile(inputFile);
    if (!inFile.is_open()) {
        cerr << "Unable to open input file: " << inputFile << endl;
        return;
    }

    // Read the entire content of the file into a string
    stringstream buffer;
    buffer << inFile.rdbuf();
    string fileContent = buffer.str();
    inFile.close();

    // Replace all occurrences of oldWord with newWord
    size_t pos = 0;
    while ((pos = fileContent.find(oldWord, pos)) != string::npos) {
        fileContent.replace(pos, oldWord.length(), newWord);
        pos += newWord.length(); // Move past the new word
    }

    // Append the line "Random correction flag: 1" at the end
    if(apply_RC) fileContent.append("Random correction flag: 1");

    // Open the output file for writing
    ofstream outFile(outputFile);
    if (!outFile.is_open()) {
        cerr << "Unable to open output file: " << outputFile << endl;
        return;
    }

    // Write the modified content to the output file
    outFile << fileContent;
    outFile.close();
}



string getFileName(const string& fullPath) {
    // Find the last occurrence of the directory separator
    size_t pos = fullPath.find_last_of("/\\");
    
    // If the separator is found, return the substring after it; otherwise, return the full path
    if (pos != string::npos) {
        return fullPath.substr(pos + 1);
    } else {
        return fullPath;
    }
}

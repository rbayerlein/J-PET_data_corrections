#define BUFFER_SIZE 65536
#define NUM_TX_CRYS_PER_MODULE 13
#define NUM_AX_CRYS_PER_MODULE 200
#define NUM_MODULES_PER_UNIT 24
#define TOF_RES_CM (300*0.03/2) // cm - TOF bin resolution (ps) * distance travelled by light (cm/ps) / 2

// ***** IMPORTANT *****
#define SIZEOF_PHOTON_CLASS 52 // required for unaligned data structure (non-padded size)
//#define RNG_DEBUG // comment out if a random seed is desired
//#define ADD_UNIT_GAP
//#define FORCE_NO_TOF_BLURRING // comment out if TOF blurring is desired

#include <chrono> // for RNG seed
#include <random> // for TOF blurring
#include <iostream>
#include "ConvertEvents.h"
#include <tuple>
#include <string>
#include <vector>
#include <cmath>
#include <sys/stat.h>



bool FileCheck(const std::string& NameOfFile);




using namespace std;

struct Photon {
	//double x_pos, y_pos, z_pos;
	int obj_scatters, col_scatters;
	double energy;
	double travel_distance;
    double det_x, det_y, det_z;
	int gcrys_id; // global crystal ID (see SimSET documentation)
};

struct History {
	Photon blue_photon, pink_photon;
	int8_t blue_flag, pink_flag; // flag byte
};

struct LmData { // list-mode data
	short txIDA, axIDA, txIDB, axIDB, tof;
};

void SaveHistory(History h, LmData* pbuffer, unsigned long long& buffer_index, FILE* pfile, double x_rng);
void SaveCastor(History h, LmData* pbuffer, unsigned long long& buffer_index, FILE* pfile, double x_time);
void Saveheader(int NumberOfEvents,int argc, char* argv[], string filename1);
bool readLUTFile(const std::string& path, int numberOfCrystals, std::vector<std::tuple<float, float, float>>& readedLutGeometry);
uint32_t getCastorID(float x, float y, float z, const std::vector<std::tuple<float, float, float>> &castorIDs);
std::string extractFileName(const std::string& fullPath);
double CalcAngle2D(double Hit1_x, double Hit1_y, double Hit2_x, double Hit2_y );

int main(int argc, char* argv[]) {

	if (argc != 3) {
		cout << "Usage: " << argv[0] << " [fname_in] [basename_out]" << endl;
		exit(1);
	}

	cout << "Starting " << argv[0] << "..." << endl;

#ifdef ADD_UNIT_GAP
	cout << "ADD_UNIT_GAP is turned ON." << endl;
#endif

#ifdef FORCE_NO_TOF_BLURRING
	cout << "FORCE_NO_TOF_BLURRING is turned ON." << endl;
#endif

	string fname_in = string(argv[1]);
	string basename_out = string(argv[2]);
	string fname_trues = basename_out + "_trues.lm";
	string fname_scatters = basename_out + "_scatters.lm";

    //modified
     

    static std::ios_base::openmode fOpenFlags =
      std::ios::binary |
      std::ios::out; 

    fOutputStream.open(basename_out + ".Cdf", fOpenFlags);
    string filename = basename_out;


    if (fOpenFlags == std::ios::binary | std::ios::out){
    fOpenFlags = std::ios::binary | std::ios::out | std::ios::app;}

    std::string lutFile = "/data/4/users/mdas/Softwares/simset/data_processing/hist2lm/src/Modular_rotated.lut";

    readLUTFile(lutFile, fNumberOfCrystals, fCastorIDs);
    

	// open files
	FILE* pfile_in = fopen(fname_in.c_str(), "rb");
	FILE* pfile_trues = fopen(fname_trues.c_str(), "wb");
	FILE* pfile_scatters = fopen(fname_scatters.c_str(), "wb");

	if (pfile_in == NULL) {
		cerr << fname_in << " cannot be opened." << endl;
		exit(1);
	}

	if (pfile_trues == NULL) {
		cerr << fname_trues << " cannot be opened." << endl;
		exit(1);
	}

	if (pfile_scatters == NULL) {
		cerr << fname_scatters << " cannot be opened." << endl;
		exit(1);
	}

	// declare input "buffer"
	// real buffer is not possible (see SimSET documentation)
	History h;

	// declare output buffers
	LmData* pbuffer_trues = new LmData[BUFFER_SIZE];
	LmData* pbuffer_scatters = new LmData[BUFFER_SIZE];
	unsigned long long ibuf_trues = 0; // buffer index
	unsigned long long ibuf_scatters = 0; // buffer index

	// declare RNG
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

#ifdef RNG_DEBUG
	cout << "RNG_DEBUG is turned ON. Using default seed..." << endl;
	default_random_engine generator;
#else
	default_random_engine generator(seed);
#endif
	normal_distribution<double> distribution(0,(600/2.355)*(0.03/2));
	
	// main loop
	fseek(pfile_in, 32768, SEEK_SET); // skip header
	do {
		// blue photon
		h.blue_flag = fgetc(pfile_in); // flag byte
		if (h.blue_flag == 1) {
			fread(&h.blue_photon, SIZEOF_PHOTON_CLASS, 1, pfile_in); // read photon info
		}

		// pink photon
		h.pink_flag = fgetc(pfile_in); // flag byte
		if (h.pink_flag == 1) {
			fread(&h.pink_photon, SIZEOF_PHOTON_CLASS, 1, pfile_in); // read photon info
		}

		// skip if a photon is missing
		if (h.blue_flag != 1 || h.pink_flag != 1) {
			continue;
		}
            
        
		// save history
		double x_rng = distribution(generator);
        double x_time = x_rng * (2/0.03);

        /*bool IsFileExists =  FileCheck("/data/4/users/mdas/energy.dat");
        std::ofstream res;
        if( IsFileExists )
          res.open("/data/4/users/mdas/energy.dat", std::ofstream::app);
        else
          res.open("/data/4/users/mdas/energy.dat");
        res << h.blue_photon.energy << "\t";
        res << h.pink_photon.energy;
        res << std::endl;
        res.close();*/

		if (h.blue_photon.obj_scatters + h.pink_photon.obj_scatters == 0) {
            

            //Additional cuts for J-PET
           double angle = CalcAngle2D(h.blue_photon.det_x, h.blue_photon.det_y,  h.pink_photon.det_x, h.pink_photon.det_y );
            
           if (angle>=60){
               if (h.blue_photon.det_z <= 23 &&  h.blue_photon.det_z >=(-23)&& h.pink_photon.det_z <= 23 && h.pink_photon.det_z>=(-23))
               {
            // trues
			SaveHistory(h, pbuffer_trues, ibuf_trues, pfile_trues, x_rng);
            //SaveCastor(h, pbuffer_trues, ibuf_trues, pfile_trues, x_time);


        /*bool IsFileExists =  FileCheck("/data/4/users/mdas/energy_true.dat");
        std::ofstream res1;
        if( IsFileExists )
          res1.open("/data/4/users/mdas/energy_true.dat", std::ofstream::app);
        else
          res1.open("/data/4/users/mdas/energy_true.dat");
        res1 << h.blue_photon.energy << "\t";
        res1 << h.pink_photon.energy;
        res1 << std::endl;
        res1.close();*/
               }
           }
            

		} else { // scatters

 

            double angle = CalcAngle2D(h.blue_photon.det_x, h.blue_photon.det_y,  h.pink_photon.det_x, h.pink_photon.det_y );
            
           if (angle>=60){
               if (h.blue_photon.det_z <= 23 &&  h.blue_photon.det_z >=(-23)&& h.pink_photon.det_z <= 23 && h.pink_photon.det_z>=(-23))
               { 
			SaveHistory(h, pbuffer_scatters, ibuf_scatters, pfile_scatters, x_rng);

        /*bool IsFileExists =  FileCheck("/data/4/users/mdas/energy_scatters.dat");
        std::ofstream res2;
        if( IsFileExists )
          res2.open("/data/4/users/mdas/energy_scatters.dat", std::ofstream::app);
        else
          res2.open("/data/4/users/mdas/energy_scatters.dat");
        res2 << h.blue_photon.energy << "\t";
        res2 << h.pink_photon.energy;
        res2 << std::endl;
        res2.close();*/
		} }}
	} while (!feof(pfile_in));

	// flush buffers
	fwrite(pbuffer_trues, sizeof(LmData), ibuf_trues, pfile_trues);
	fwrite(pbuffer_scatters, sizeof(LmData), ibuf_scatters, pfile_scatters);

	// cleanup
	delete[] pbuffer_trues;
	delete[] pbuffer_scatters;
	fclose(pfile_in);
	fclose(pfile_trues);
	fclose(pfile_scatters);
    
	cout << "Done." << endl;

    Saveheader(fNumberOfEvents, argc, argv, filename);

    cout<<"Castor output file has been written"<<endl;
	return 0;

}

void SaveHistory(History h, LmData* pbuffer, unsigned long long& buffer_index, FILE* pfile, double x_rng) {
	// convert to LmData format
		short UiA = h.blue_photon.gcrys_id / (NUM_TX_CRYS_PER_MODULE * NUM_AX_CRYS_PER_MODULE * NUM_MODULES_PER_UNIT);
		short UiB = h.pink_photon.gcrys_id / (NUM_TX_CRYS_PER_MODULE * NUM_AX_CRYS_PER_MODULE * NUM_MODULES_PER_UNIT);

		// module ID
		short moduleA = h.blue_photon.gcrys_id / (NUM_TX_CRYS_PER_MODULE * NUM_AX_CRYS_PER_MODULE) % NUM_MODULES_PER_UNIT;
		short moduleB = h.pink_photon.gcrys_id / (NUM_TX_CRYS_PER_MODULE * NUM_AX_CRYS_PER_MODULE) % NUM_MODULES_PER_UNIT;

		// crystal ID
		short crysA = h.blue_photon.gcrys_id % (NUM_TX_CRYS_PER_MODULE * NUM_AX_CRYS_PER_MODULE);
		short crysB = h.pink_photon.gcrys_id % (NUM_TX_CRYS_PER_MODULE * NUM_AX_CRYS_PER_MODULE);

		// put everything together
		short txIDA = (crysA % NUM_TX_CRYS_PER_MODULE) + moduleA * NUM_TX_CRYS_PER_MODULE;
		short txIDB = (crysB % NUM_TX_CRYS_PER_MODULE) + moduleB * NUM_TX_CRYS_PER_MODULE;
#ifdef ADD_UNIT_GAP
		short axIDA = (crysA / NUM_TX_CRYS_PER_MODULE) + UiA * NUM_AX_CRYS_PER_MODULE + UiA;
		short axIDB = (crysB / NUM_TX_CRYS_PER_MODULE) + UiB * NUM_AX_CRYS_PER_MODULE + UiB;
#else
		short axIDA = (crysA / NUM_TX_CRYS_PER_MODULE) + UiA * NUM_AX_CRYS_PER_MODULE;
		short axIDB = (crysB / NUM_TX_CRYS_PER_MODULE) + UiB * NUM_AX_CRYS_PER_MODULE;
#endif
		double diff_travel_distance = h.blue_photon.travel_distance - h.pink_photon.travel_distance;

#ifdef FORCE_NO_TOF_BLURRING
		short tof = round(diff_travel_distance / 2 / TOF_RES_CM);
#else
		short tof = round((diff_travel_distance / 2 + x_rng) / TOF_RES_CM);

        
#endif

		// write to buffer
		pbuffer[buffer_index].txIDA = txIDA;
		pbuffer[buffer_index].axIDA = axIDA;
		pbuffer[buffer_index].txIDB = txIDB;
		pbuffer[buffer_index].axIDB = axIDB;
		pbuffer[buffer_index].tof = tof;
   


		buffer_index++;
		if (buffer_index == BUFFER_SIZE) {
			// write buffer to file
			fwrite(pbuffer, sizeof(LmData), BUFFER_SIZE, pfile);
			buffer_index = 0;
		}
}










void SaveCastor(History h, LmData* pbuffer, unsigned long long& buffer_index, FILE* pfile, double x_time){



        short UiA = h.blue_photon.gcrys_id / (NUM_TX_CRYS_PER_MODULE * NUM_AX_CRYS_PER_MODULE * NUM_MODULES_PER_UNIT);
		short UiB = h.pink_photon.gcrys_id / (NUM_TX_CRYS_PER_MODULE * NUM_AX_CRYS_PER_MODULE * NUM_MODULES_PER_UNIT);

		// module ID
		short moduleA = h.blue_photon.gcrys_id / (NUM_TX_CRYS_PER_MODULE * NUM_AX_CRYS_PER_MODULE) % NUM_MODULES_PER_UNIT;
		short moduleB = h.pink_photon.gcrys_id / (NUM_TX_CRYS_PER_MODULE * NUM_AX_CRYS_PER_MODULE) % NUM_MODULES_PER_UNIT;

		// crystal ID
		short crysA = h.blue_photon.gcrys_id % (NUM_TX_CRYS_PER_MODULE * NUM_AX_CRYS_PER_MODULE);
		short crysB = h.pink_photon.gcrys_id % (NUM_TX_CRYS_PER_MODULE * NUM_AX_CRYS_PER_MODULE);

		// put everything together
		short txIDA = (crysA % NUM_TX_CRYS_PER_MODULE) + moduleA * NUM_TX_CRYS_PER_MODULE;
		short txIDB = (crysB % NUM_TX_CRYS_PER_MODULE) + moduleB * NUM_TX_CRYS_PER_MODULE;
#ifdef ADD_UNIT_GAP
		short axIDA = (crysA / NUM_TX_CRYS_PER_MODULE) + UiA * NUM_AX_CRYS_PER_MODULE + UiA;
		short axIDB = (crysB / NUM_TX_CRYS_PER_MODULE) + UiB * NUM_AX_CRYS_PER_MODULE + UiB;
#else
		short axIDA = (crysA / NUM_TX_CRYS_PER_MODULE) + UiA * NUM_AX_CRYS_PER_MODULE;
		short axIDB = (crysB / NUM_TX_CRYS_PER_MODULE) + UiB * NUM_AX_CRYS_PER_MODULE;
#endif
		double diff_travel_distance = h.blue_photon.travel_distance - h.pink_photon.travel_distance;
         
         //cout<<h.blue_photon.det_x <<" "<< h.blue_photon.det_y <<" "<< h.blue_photon.det_z <<endl;

      
        

//#ifdef FORCE_NO_TOF_BLURRING
		//short tof = round(diff_travel_distance / 2 / TOF_RES_CM);
//#else
		//short tof = round((diff_travel_distance / 2 + x_rng) / TOF_RES_CM);

      uint32_t time = 1;
      float tof = (diff_travel_distance * 33.3564095 )+x_time;
      uint32_t castor1 = getCastorID(h.blue_photon.det_x * 10, h.blue_photon.det_y  * 10, h.blue_photon.det_z  * 10, fCastorIDs);
      uint32_t castor2 = getCastorID(h.pink_photon.det_x * 10, h.pink_photon.det_y * 10, h.pink_photon.det_z * 10, fCastorIDs);
      fOutputStream.write(reinterpret_cast<char *>(&time), sizeof(time));
      fOutputStream.write(reinterpret_cast<char *>(&tof), sizeof(tof));
      fOutputStream.write(reinterpret_cast<char *>(&castor1), sizeof(castor1));
      fOutputStream.write(reinterpret_cast<char *>(&castor2), sizeof(castor2));
      fNumberOfEvents++;

    //cout<<castor1 <<" "<< castor2<<" "<<tof<< endl;
}










void Saveheader(int NumberOfEvents, int argc, char* argv[], string filename1)
{



	/*if (argc != 3) {
		cout << "Usage: " << argv[0] << " [fname_in] [basename_out]" << endl;
		exit(1);
	}
	string basename_out = string(argv[2]);*/
  
   string filename = filename1; \
   string filename2 = extractFileName(filename + ".Cdf");

  static uint32_t totalEvents = fNumberOfEvents;
  fOutputStream.close();
  std::string headerFile = filename1 + ".Cdh";
  fOutputStream.open(headerFile, std::ios::out);
  fOutputStream << "Scanner name: Modular" << std::endl
       << "Data filename: " << filename2  << std::endl
       << "Number of events: " << totalEvents << std::endl
       << "Data mode: list-mode" << std::endl
       << "Data type: PET" << std::endl
       // here we are setting Start time/Duration to dummy values
       // JPet Framework do not have this kind of information
       // TODO: check if this value does not influence reconstruction
       // in some way
       << "Start time (s): 0" << std::endl
       << "Duration (s): 10" << std::endl
       << "TOF information flag: 1" << std::endl
       << "TOF resolution (ps): 600 "<< std::endl
       << "List TOF measurement range (ps): 3000" << std::endl;
  fOutputStream.close();}








bool readLUTFile(const std::string& path, int numberOfCrystals, std::vector<std::tuple<float, float, float>> &readedLutGeometry) {
  FILE* LUT_file = fopen(path.c_str(), "rb");
  if (LUT_file==NULL)
  {
    cout <<"failed to open the lut file"<<endl;
    return false;
  }

  // Read data for each index
  int nb_data_read = 0;
  for (int i=0; i<numberOfCrystals; i++)
  {
    float x, y, z;
    float skip;
    // Read central crystal position X, then Y and Z
    nb_data_read += fread(&x,sizeof(x),1,LUT_file);
    nb_data_read += fread(&y,sizeof(y),1,LUT_file);
    nb_data_read += fread(&z,sizeof(z),1,LUT_file);
    // Read crystal orientation X, then Y and Z
    nb_data_read += fread(&skip,sizeof(skip),1,LUT_file);
    nb_data_read += fread(&skip,sizeof(skip),1,LUT_file);
    nb_data_read += fread(&skip,sizeof(skip),1,LUT_file);
    readedLutGeometry.push_back(std::make_tuple(x, y, z));
    
    //std::cout << "[ " << x << ", " << y << ", " << z << " ]"  << std::endl;
  }

  // Close file
  fclose(LUT_file);
  // Check reading
  if (nb_data_read!=numberOfCrystals*6) {
    cout<<"error"<<endl;
    return false;
  }
  return true;
}

uint32_t getCastorID(float x, float y, float z, const std::vector<std::tuple<float, float, float>> &castorIDs) {
  uint32_t i = 0;
  // TODO: generally doesn't metter, but maybe pick more meaningful staring lowestDistance
  float lowestDistance = 10;
  uint32_t lowestI = 0;

  for (auto t : castorIDs) {
    float diffX = std::get<0>(t) - x;
    float diffY = std::get<1>(t) - y;
    float diffZ = std::get<2>(t) - z;
    float distance = std::sqrt(diffX * diffX + diffY * diffY + diffZ * diffZ);

    if (distance < lowestDistance) {
      lowestDistance = distance;
      lowestI = i;
    }
    i++;
  }
  return lowestI;
}

std::string extractFileName(const std::string& fullPath) {
    size_t lastSlash = fullPath.find_last_of("/\\");
    if (lastSlash != std::string::npos) {
        return fullPath.substr(lastSlash + 1);
    }
    return fullPath;
}

double CalcAngle2D(double Hit1_x, double Hit1_y, double Hit2_x, double Hit2_y )
{
  double scalarProd = Hit1_x*Hit2_x + Hit1_y *Hit2_y;
  double magProd = sqrt( ( pow(Hit1_x,2) + pow(Hit1_y,2) )*
                ( pow(Hit2_x,2) + pow(Hit2_y,2) ) );
  double Angle = acos( scalarProd/magProd )*180/3.14159265;
  return Angle;
}
bool FileCheck(const std::string& NameOfFile) 
{
  struct stat buffer;   
  return (stat (NameOfFile.c_str(), &buffer) == 0); 
}

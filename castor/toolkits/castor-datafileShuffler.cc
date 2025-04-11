/*
This file is part of CASToR.

    CASToR is free software: you can redistribute it and/or modify it under the
    terms of the GNU General Public License as published by the Free Software
    Foundation, either version 3 of the License, or (at your option) any later
    version.

    CASToR is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
    details.

    You should have received a copy of the GNU General Public License along with
    CASToR (in file GNU_GPL.TXT). If not, see <http://www.gnu.org/licenses/>.

Copyright 2017-2021 all CASToR contributors listed below:

    --> Didier BENOIT, Claude COMTAT, Marina FILIPOVIC, Thibaut MERLIN, Mael MILLARDET, Simon STUTE, Valentin VIELZEUF, Zacharias CHALAMPALAKIS

This is CASToR version 3.1.1.
*/

/*!
  \file
  \ingroup utils
  \brief  This program is used to shuffle the events order of a histogram datafile.
*/

#include "gVariables.hh"
#include "gOptions.hh"
#include "oImageDimensionsAndQuantification.hh"
#include "vDataFile.hh"
#include "iDataFilePET.hh"
#include "iDataFileSPECT.hh"
#include "iDataFileCT.hh"
#include "sOutputManager.hh"
#include "sRandomNumberGenerator.hh"

// =============================================================================================================================================
// =============================================================================================================================================
// =============================================================================================================================================
//                                                        H E L P     F U N C T I O N S
// =============================================================================================================================================
// =============================================================================================================================================
// =============================================================================================================================================


/*!
  \fn      ShowHelp()
  \brief   Display main command line options for castor-datafileShuffler
*/
void ShowHelp()
{
  // Show help
  cout << endl;
  cout << "Usage:  castor-datafileShuffler  -df datafile.cdh  -(f/d)out output" << endl;
  cout << endl;
  cout << "This program can be used to shuffle the order of events of a histogram datafile." << endl;
  cout << endl;
  cout << "[Mandatory parameters]:" << endl;
  cout << "  -df datafile.cdh     : Give a CASToR list-mode datafile." << endl;
  cout << "  -fout name           : Give the root name for all output files (no default, alternative to -dout)" << endl;
  cout << "  -dout name           : Give the name of the output directory where all output files will be written (no default, alternative to -fout)" << endl;
  cout << endl;
  cout << "[Options]:" << endl;
  cout << "  -seed value          : Give a seed for the random number generator (should be >=0)" << endl;
  cout << "  -vb value            : Give the verbosity level, from 0 (no verbose) to 2 (default: 1)" << endl;
  cout << "  --help,-h,-help      : Print out this help page." << endl; // managed by main
  cout << endl;
  #ifdef BUILD_DATE
  cout << "  Build date: " << BUILD_DATE << endl;
  cout << endl;
  #endif
  #ifdef CASTOR_VERSION
  cout << "  This program is part of the CASToR release version " << CASTOR_VERSION << "." << endl;
  cout << endl;
  #endif
}

// =============================================================================================================================================
// =============================================================================================================================================
// =============================================================================================================================================
//                                                        M A I N     P R O G R A M
// =============================================================================================================================================
// =============================================================================================================================================
// =============================================================================================================================================

int main(int argc, char** argv)
{
  // ============================================================================================================
  // MPI stuff (we make all instances except the first one returning 0 directly)
  // ============================================================================================================
  int mpi_rank = 0;
  #ifdef CASTOR_MPI
  int mpi_size = 1;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  if (mpi_rank!=0) return 0;
  #endif

  // No argument, then show help
  if (argc==1)
  {
    ShowHelp();
    Exit(EXIT_SUCCESS);
  }

  // ============================================================================================================
  // Parameterized variables with their default values
  // ============================================================================================================

  // --------------------------------------------------------------------------------
  // Input settings
  // --------------------------------------------------------------------------------

  // Input datafile
  string datafile = "";

  // --------------------------------------------------------------------------------
  // Output settings
  // --------------------------------------------------------------------------------

  // Output directory name.
  string path_dout = "";
  // Or root name
  string path_fout = "";

  // --------------------------------------------------------------------------------
  // Miscellaneous settings
  // --------------------------------------------------------------------------------

  // Verbose level
  int verbose = 1;
  // Initial seed for random number generator
  int64_t random_generator_seed = -1;

  // ============================================================================================================
  // Read command-line parameters
  // ============================================================================================================

  // Must manually increment the option index when an argument is needed after an option
  for (int i=1; i<argc; i++)
  {
    // Get the option as a string
    string option = (string)argv[i];

    // --------------------------------------------------------------------------------
    // Miscellaneous settings
    // --------------------------------------------------------------------------------

    // Show help
    if (option=="-h" || option=="--help" || option=="-help")
    {
      ShowHelp();
      Exit(EXIT_SUCCESS);
    }
    // RNG seed
    else if (option=="-seed")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-datafileShuffler() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      if (ConvertFromString(argv[i+1], &random_generator_seed))
      {
        Cerr("***** castor-datafileShuffler() -> Exception when trying to read provided number '" << random_generator_seed << " for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      i++;
    }
    // General verbosity level
    else if (option=="-vb")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-datafileShuffler() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      if (ConvertFromString(argv[i+1], &verbose))
      {
        Cerr("***** castor-datafileShuffler() -> Exception when trying to read provided verbosity level '" << verbose << " for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      i++;
    }

    // --------------------------------------------------------------------------------
    // Input settings
    // --------------------------------------------------------------------------------

    // Images
    else if (option=="-df") // This is a mandatory option
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-datafileShuffler() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      datafile = (string)argv[i+1];
      i++;
    }

    // --------------------------------------------------------------------------------
    // Output settings
    // --------------------------------------------------------------------------------

    // Name of the output directory
    else if (option=="-dout") // This is a mandatory option
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-datafileShuffler() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      path_dout = argv[i+1];
      i++;
    }
    // Base name of the output files
    else if (option=="-fout") // This is a mandatory option
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-datafileShuffler() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      path_fout = argv[i+1];
      i++;
    }

    // --------------------------------------------------------------------------------
    // Unknown option!
    // --------------------------------------------------------------------------------

    else
    {
      Cerr("***** castor-datafileShuffler() -> Unknown option '" << option << "' !" << endl);
      Exit(EXIT_FAILURE);
    }
  }

  // ============================================================================================================
  // Some checks
  // ============================================================================================================

  // Data file
  if (datafile=="")
  {
    Cerr("***** castor-datafileShuffler() -> Please provide an input datafile !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Output files
  if (path_fout.empty() && path_dout.empty())
  {
    Cerr("***** castor-datafileShuffler() -> Please provide an output option for output files (-fout or -dout) !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Check that only one option has been provided
  if (!path_fout.empty() && !path_dout.empty())
  {
    Cerr("***** castor-datafileShuffler() -> Please provide either output option -fout or -dout but not both !" << endl);
    Exit(EXIT_FAILURE);
  }

  // ============================================================================================================
  // Create sOutputManager
  // ============================================================================================================
  sOutputManager* p_OutputManager = sOutputManager::GetInstance();  
  // Set verbose level
  p_OutputManager->SetVerbose(verbose);
  // Set MPI rank
  p_OutputManager->SetMPIRank(mpi_rank);
  // Initialize output directory and base name
  if (p_OutputManager->InitOutputDirectory(path_fout, path_dout))
  {
    Cerr("***** castor-datafileShuffler() -> A problem occurred while initializing output directory !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Log command line
  if (p_OutputManager->LogCommandLine(argc,argv))
  {
    Cerr("***** castor-datafileShuffler() -> A problem occurred while logging command line arguments !" << endl);
    Exit(EXIT_FAILURE);
  }

  // ============================================================================================================
  // Create sScannerManager (in order to get the datafile type)
  // ============================================================================================================
  sScannerManager* p_ScannerManager = sScannerManager::GetInstance();  
  p_ScannerManager->SetVerbose(verbose);
  // Get system name from the dataFile
  string scanner_name = "";
  if (ReadDataASCIIFile(datafile, "Scanner name", &scanner_name, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** castor-datafileShuffler() -> A problem occurred while trying to find the system name in the datafile header !" << endl);
    Exit(EXIT_FAILURE);
  } 
  if (p_ScannerManager->FindScannerSystem(scanner_name) )
  {
    Cerr("***** castor-datafileShuffler() -> A problem occurred while searching for scanner system !" << endl);
    Exit(EXIT_FAILURE);
  } 
  if (p_ScannerManager->BuildScannerObject() )
  {
    Cerr("***** castor-datafileShuffler() -> A problem occurred during scanner object construction ! !" << endl);
    Exit(EXIT_FAILURE);
  }
  if (p_ScannerManager->InstantiateScanner() )
  {
    Cerr("***** castor-datafileShuffler() -> A problem occurred while creating Scanner object !" << endl);
    Exit(EXIT_FAILURE);
  }
  if (p_ScannerManager->GetGeometricInfoFromDataFile(datafile))
  {
    Cerr("***** castor-datafileShuffler() -> A problem occurred while retrieving scanner fields from the datafile header !" << endl);
    Exit(EXIT_FAILURE);
  }
  if (p_ScannerManager->BuildLUT() )
  {
    Cerr("***** castor-datafileShuffler() -> A problem occurred while generating/reading the LUT !" << endl);
    Exit(EXIT_FAILURE);
  } 

  // ============================================================================================================
  // Create the input datafile
  // ============================================================================================================

  // Create a default image dimensions and quantification object
  oImageDimensionsAndQuantification* p_ID = new oImageDimensionsAndQuantification();
  p_ID->SetDefault();
  p_ID->SetVerbose(verbose);
  // Create the datafile based on the data type
  vDataFile* p_DataFile = NULL;
  if (p_ScannerManager->GetScannerType() == SCANNER_PET) p_DataFile = new iDataFilePET();
  else if (p_ScannerManager->GetScannerType() == SCANNER_SPECT_CONVERGENT) p_DataFile = new iDataFileSPECT(); 
  else if (p_ScannerManager->GetScannerType() == SCANNER_CT) p_DataFile = new iDataFileCT(); 
  else
  {
    Cerr("***** castor-datafileShuffler() -> Unknown scanner type (" << p_ScannerManager->GetScannerType() << ") for datafile construction ! Abort." << endl);
    Exit(EXIT_FAILURE);
  }
  p_DataFile->SetImageDimensionsAndQuantification(p_ID);
  p_DataFile->SetHeaderDataFileName(datafile);
  p_DataFile->SetVerbose(verbose);
  p_DataFile->SetBedIndex(0);
  bool do_not_affect_quantification = false;
  if (p_DataFile->ReadInfoInHeader(do_not_affect_quantification))
  {
    Cerr("***** castor-datafileShuffler() -> A problem occurred during datafile header reading ! Abort." << endl);
    Exit(EXIT_FAILURE);
  }
  if (p_DataFile->CheckParameters())
  {
    Cerr("***** castor-datafileShuffler() -> A problem occurred while checking datafile parameters ! Abort." << endl);
    Exit(EXIT_FAILURE);
  }
  if (p_DataFile->ComputeSizeEvent())
  {
    Cerr("***** castor-datafileShuffler() -> A problem occurred in datafile initialization ! Abort." << endl);
    Exit(EXIT_FAILURE);
  }
  if (p_DataFile->InitializeMappedFile())
  {
    Cerr("***** castor-datafileShuffler() -> A problem occurred in datafile initialization ! Abort." << endl);
    Exit(EXIT_FAILURE);
  }
  if (p_DataFile->PrepareDataFile())
  {
    Cerr("***** castor-datafileShuffler() -> A problem occurred in datafile preparation ! Abort." << endl);
    Exit(EXIT_FAILURE);
  }
  // Check if datafile is a histogram, otherwise, throw an error
  if (p_DataFile->GetDataMode()!=MODE_HISTOGRAM)
  {
    Cerr("***** castor-datafileShuffler() -> The input datafile is not a histogram, this program is only suitable to histogram files !" << endl);
    Exit(EXIT_FAILURE);
  }

  // ============================================================================================================
  // Create the output datafile
  // ============================================================================================================

  // Create output datafile object
  vDataFile* p_OutputDataFile = NULL;
  if (p_ScannerManager->GetScannerType() == SCANNER_PET) p_OutputDataFile = new iDataFilePET();
  else if (p_ScannerManager->GetScannerType() == SCANNER_SPECT_CONVERGENT) p_OutputDataFile = new iDataFileSPECT(); 
  else if (p_ScannerManager->GetScannerType() == SCANNER_CT) p_OutputDataFile = new iDataFileCT(); 
  // Build output data file from the input one
  if (p_OutputDataFile->SetParametersFrom(p_DataFile))
  {
    Cerr("***** castor-datafileShuffler() -> An error occurred while setting parameters of output file from input file !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Open output file
  if (p_OutputDataFile->OpenFileForWriting())
  {
    Cerr("***** castor-datafileShuffler() -> An error occurred while opening file for writing !" << endl);
    Exit(EXIT_FAILURE);
  }

  // ============================================================================================================
  // Build random event indices list
  // ============================================================================================================

  // Get the number of events
  int64_t nb_events = p_DataFile->GetSize();
  // Verbose
  if (verbose>=1) Cout("castor-datafileShuffler() -> Build a random indices list for " << nb_events << " events" << endl);
  if (verbose>=2) Cout("  --> Allocate random indices list" << endl);
  // Allocate the list of indices and a vector of random mnumbers
  vector<int64_t> p_random_indices(nb_events);
  vector<int64_t> p_random_numbers(nb_events);
  // Fill it in the correct order
  for (int64_t i=0; i<nb_events; i++) p_random_indices[i] = i;
  // Seed of the random generator
  if (random_generator_seed==-1)
  {
    // Not provided so we set the seed randomly
    mt19937 rd(time(NULL));
    random_generator_seed = rd();
  }
  else
  {
    // Check seed non negativity
    if (random_generator_seed<0)
    {
      Cerr("***** castor-datafileShuffler() -> Must provide a random generator seed that is positive !" << endl);
      Exit(EXIT_FAILURE);
    }
  }
  // Verbose
  if (verbose>=2) Cout("  --> General random seed: " << random_generator_seed << endl);
  // Create a random generator to generate as many random seeds as threads
  mt19937 random_seeder(random_generator_seed);
  // Create random generators (as many as threads)
  int nb_threads = 1;
  #ifdef CASTOR_OMP
  nb_threads = omp_get_num_procs();
  if (verbose>=2) Cout("  --> Use " << nb_threads << " threads" << endl);
  #endif
  mt19937** p_random_generators = (mt19937**)malloc(nb_threads*sizeof(mt19937*));
  for (int th=0; th<nb_threads; th++) p_random_generators[th] = new mt19937(random_seeder());
  // Create the uniform integer distribution between 0 and the size of the datafile minus 1
  uniform_int_distribution<int64_t> random_distribution(0, nb_events-1);
  // Verbose
  if (verbose>=2) Cout("  --> Shoot random numbers" << endl);
  // Shoot the random numbers and store them
  int64_t printing_index = 0;
  int64_t i = 0;
  #pragma omp parallel for private(i) schedule(static)
  for (i=0; i<nb_events; i++)
  {
    // The thread index
    int th = 0;
    #ifdef CASTOR_OMP
    th = omp_get_thread_num();
    #endif
    // Verbose
    if (th==0 && verbose>=2)
    {
      if (printing_index%5000==0)
      {
        int percent = ((int)( ((FLTNB)(i)) * 100. / ((FLTNB)(nb_events/nb_threads)) ));
        cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b      "
             << percent << " %                    " << flush;
      }
      printing_index++;
    }
    // Shoot the random number and store it
    p_random_numbers[i] = random_distribution(*p_random_generators[th]);
  }
  // Verbose
  if (verbose>=2)
  {
    cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
         << "  --> 100 %                        " << endl;
    Cout("  --> Apply random permutations" << endl);
  }
  // Shuffle the list using random permutations 
  printing_index = 0;
  for (int64_t i=0; i<nb_events; i++)
  {
    // Verbose
    if (verbose>=2)
    {
      if (printing_index%5000==0)
      {
        int percent = ((int)( ((FLTNB)(i)) * 100. / ((FLTNB)(nb_events)) ));
        cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b      "
             << percent << " %                    " << flush;
      }
      printing_index++;
    }
    // Permute the i^th event with a random event
    int64_t index = p_random_numbers[i];
    int64_t content = p_random_indices[i];
    p_random_indices[i] = p_random_indices[index];
    p_random_indices[index] = content;
  }
  // Verbose
  if (verbose>=2) cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
                       << "  --> 100 %                        " << endl;

  // ============================================================================================================
  // Write shuffled output file
  // ============================================================================================================

  // Verbose
  if (verbose>=1) Cout("castor-datafileShuffler() -> Read and write events" << endl);
  // Start the loop over events
  printing_index = 0;
  for (int64_t i=0; i<nb_events; i++)
  {
    // Verbose
    if (verbose>=2)
    {
      if (printing_index%10000==0)
      {
        int percent = ((int)( ((FLTNB)(i)) * 100. / ((FLTNB)(nb_events)) ));
        cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b      "
             << percent << " %                    " << flush;
      }
      printing_index++;
    }
    // Get the random event
    int64_t index = p_random_indices[i];
    vEvent* event = p_DataFile->GetEvent(index);
    if (event==NULL)
    {
      Cerr("***** castor-datafileShuffler() -> An error occurred while getting the event from index " << index << " !" << endl);
      if (p_OutputDataFile->CloseFile())
        Cerr("***** castor-datafileShuffler() -> An error occurred while closing file during writing !" << endl);
      Exit(EXIT_FAILURE);
    }
    // Write the event
    p_OutputDataFile->WriteEvent(event);
  }
  // Verbose
  if (verbose>=2) cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
                       << "  --> 100 %                        " << endl;
  // Close file
  if (p_OutputDataFile->CloseFile())
  {
    Cerr("***** castor-datafileShuffler() -> An error occurred while closing file after writing !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Set number of events
  p_OutputDataFile->SetNbEvents(nb_events);
  // Write header
  if (p_OutputDataFile->WriteHeader())
  {
    Cerr("***** castor-datafileShuffler() -> An error occurred while writing output header file !" << endl);
    Exit(EXIT_FAILURE);
  }

  // ============================================================================================================
  // Exit
  // ============================================================================================================

  // Ending
  #ifdef CASTOR_MPI
  MPI_Finalize();
  #endif
  return EXIT_SUCCESS;
}

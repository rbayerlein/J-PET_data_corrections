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
  \brief This program convert a GATE datafile in root format to a list-mode datafile in CASToR format
*/

#include "gVariables.hh"
#include "gOptions.hh"
#include "sOutputManager.hh"
#include <iomanip> //std::setprecision

#include "oImageSpace.hh"
#include "oDynamicModelManager.hh"
#include "sRandomNumberGenerator.hh"
#include "sChronoManager.hh"

/*!
  \fn      ShowHelp()
  \param a_returnCode
  \brief   Display main command line options for castor-imageDynamicTools
*/
void ShowHelp(int a_returnCode)
{
    cout << endl;
    cout << "Usage:  castor-imageDynamicTools      -i             path_to_ifile.root"<< endl;
    cout << "                                      -fout          name"<< endl;
    cout << "                                         -OR-         "<< endl;
    cout << "                                      -dout          name "<< endl;
    cout << "                                      -dynamic-model param" << endl;  
    cout << endl;
    cout << endl;
    cout << "[Main settings]:" << endl;
    cout << "  -i  path_to_in_image   : give an input dynamic interfile image" << endl;
    cout << "  -fout name             : Give the root name for all output files. All output files will be written as 'name_suffix.ext'." << endl;
    cout << "                           So the provided name should not end with '.' or '/' character. (no default, alternative to -dout)" << endl;
    cout << "  -dout name             : Give the name of the output directory where all output files will be written. All files will also" << endl;
    cout << "                           be prefixed by the name of the directory. The provided name should not end with '.' or '/' character." << endl;
    cout << "                          (no default, alternative to -fout)" << endl;
    cout << "  -dynamic-model param   : Dynamic model applied to either the frames of a dynamic acquisition, respiratory-gated frames, cardiac-gated frames, or simultaneously between these datasets." << endl;
    cout << "                           Select the dynamic model to be used, along with a configuration file (model:file) or the list of parameters associated to the model (model_name,param1,param2,...)." << endl;
    cout << endl;
    cout << "  -delay value           : Delay condition in seconds to exclude early frames from the application of the dynamic model, if series of images are provided as input." << endl;
    cout << endl;
    cout << "  -im param              : Provide an image-based deformation model to be used for involuntary motion correction (for testing purposes), along with a configuration file (deformation:file)" << endl;
    cout << "                           or the list of parameters associated to the projector (deformation,param1,param2,...)." << endl;
    cout << endl; 
    cout << endl;
    cout << "[Optional settings]:" << endl;
    cout << "  -it  nb                : Number of iterations for both parameters estimation and subsequent dynamic image estimation steps of the dynamic model (default = 1)." << endl;
    cout << "  -omd                   : (M)erge (D)ynamic images. Indicate if a dynamic serie of 3D images should be written on disk in one file" << endl;
    cout << "                           instead of a serie of 3D images associated with an interfile metaheader" << endl;
    cout << "  -odi                   : Output dynamic images. Generate the serie of dynamic images estimated from the parametric images of the dynamic model. By default, only the parametric images are generated" << endl;
    cout << endl;
    cout << "  -onbp prec             : By default, numbers are displayed using scientific format. This option allows to customize the format and precision" << endl;
    cout << "                         : The format is format,precision. f is the format (s=scientific, f=fixed), p is the precision" << endl;
    cout << "                           eg: -onbp f,5 --> fixed with 5 digits precision, -onbp -->  scientifici with max precision." << endl;
    cout << endl;
    #ifdef CASTOR_OMP
    cout << "  -th param              : Set the number of threads for parallel computation (default: 1). If 0 is given, then the maximum number of available threads is automatically determined." << endl;
    cout << "                           Can also give two parameters separated by a comma (e.g. 16,4), to distinguish between the number of threads for projection and image operations respectively." << endl;
    cout << endl;
    #endif
    cout << "[Miscellaneous settings]:" << endl;
    cout << "  -vb             : give the verbosity level, from 0 (no verbose) to above 5 (at the event level) (default: 1)." << endl;
    cout << endl;

  Exit(a_returnCode);
}


/*
  Main program
*/

int main(int argc, char** argv)
{

  // ============================================================================================================
  // MPI stuff (we make all instances but the first one returning 0 directly)
  // ============================================================================================================
  #ifdef CASTOR_MPI
  int mpi_rank = 0;
  int mpi_size = 1;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  if (mpi_rank!=0) return 0;
  #endif

  // No argument, then show help
  if (argc==1) ShowHelp(0);

  // ============================================================================================================
  // Parameterized variables with their default values
  // ============================================================================================================

  // String which gathers the path to the input data filename provided by the user. no default 
  string input_file = "";
  vector<string> path_to_input_image;
  
  // Path to user provided data and output files/images
  string path_fout = "";
  string path_dout = "";
  string dynamic_model_options = "";
  string motion_options = "";
  
  // Frame list descriptor
  string frame_list = "";
  // Delay for application of dynamic model
  int delayInSec = 0 ;
  // Number of images to skip based on delay criteria
  int nb_FramesToSkip =0;
  // Number of resp gates
  int nb_resp_gates = 1;
  // Number of card gates
  int nb_card_gates = 1;
  // Number of iterations for the dynamic model
  uint32_t nb_ite_dyn = 1;
  // Number of threads
  string nb_threads = "1";
  // Merge output dynamic images on disk into one file or not
  bool merge_dynamic_imgs_flag = false;
  // Output dynamic images estimated from the model
  bool output_dynamic_imgs_from_model = false;
  // Precision for output number display
  string onb_prec = "s,0";
  
  // Verbosity
  int vb = 0;
  
  // ============================================================================================================
  // Read command-line parameters
  // ============================================================================================================
  for (int i=1; i<argc; i++)
  {
    string option = (string)argv[i];
    
    if (option=="-h" || option=="--help" || option=="-help") ShowHelp(0);
    
    // Just one file is provided
    if (option=="-i")
    {
      if (path_to_input_image.size() > 0)
      {
        Cerr("***** castor-imageDynamicTools :: the following file names have already been provided (-i/-il options must be used ONCE): " << endl);    
        for (size_t i=0 ; i<path_to_input_image.size() ; i++) 
          Cout(" Selected input image file: " << path_to_input_image[i] << endl);
        Exit(EXIT_FAILURE);
      }
      else
      {
        if (argv[i+1] == NULL)
        {
          Cerr("***** castor-imageDynamicTools :: argument missing for option: " << option << endl);
          Exit(EXIT_FAILURE);
        }
        else
        {
          input_file = argv[i+1];
          path_to_input_image.push_back(input_file);
        }        
        i++;
      }
    }

    // Name of the output directory
    else if (option=="-dout") // This is a mandatory option
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-imageDynamicTools() -> Argument missing for option: " << option << endl);
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
        Cerr("***** castor-imageDynamicTools() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      path_fout = argv[i+1];
      i++;
    }

    // Flag to say that we want to save time basis functions too
    else if (option=="-omd")
    {
      merge_dynamic_imgs_flag = true;
    }
    
    else if (option=="-odi")
    {
      output_dynamic_imgs_from_model = true;
    }
    
    #ifdef CASTOR_OMP
    else if (option=="-th")
    {
    if (i>=argc-1)
    {
     Cerr("***** castor-imageDynamicTools() -> Argument missing for option: " << option << endl);
     Exit(EXIT_FAILURE);
    }
    nb_threads = (string)argv[i+1];
    i++;
    }
    #endif
    
    // Dynamic model applied to the frames/respiratory gates/cardiac gates of the dynamic acquisition
    else if (option=="-dynamic-model") 
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-imageDynamicTools() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      dynamic_model_options = (string)argv[i+1];
      i++;
    }
    // Delay in seconds
    else if (option=="-delay")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-imageDynamicTools :: Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      delayInSec = atoi(argv[i+1]);
      i++;
    }
    // Number of iterations for the dynamic model 
    else if (option=="-it") 
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-imageDynamicTools() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      
      if (ConvertFromString(argv[i+1], &nb_ite_dyn))
      {
        Cerr("***** castor-imageDynamicTools() -> Exception when trying to read provided number iterations for the dynamic model (must be strictly positive integer) '" << nb_ite_dyn << " for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      i++;
    }
    
    // Data for involuntary motion correction based on deformation
    else if (option=="-im")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      motion_options = (string)argv[i+1];
      i++;
    }
    
    // Verbosity level
    else if (option=="-vb")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-imageDynamicTools :: Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      vb = atoi(argv[i+1]);
      i++;
    }
    // Output number precision
    else if (option=="-onbp")
    {
      if (i>=argc-1)
      {
        Cerr("***** castor-recon() -> Argument missing for option: " << option << endl);
        Exit(EXIT_FAILURE);
      }
      onb_prec = argv[i+1];
      i++;
    }
    else
    {
      Cerr("***** castor-imageDynamicTools :: Unknown option '" << option << "' !" << endl);
      Exit(EXIT_FAILURE);
    }
  }


  // ============================================================================================================
  // Mandatory checks:
  // ============================================================================================================

  // Basic initialization checks (minimal initializations mandatory for the next steps)

  // data files
  if (path_to_input_image.empty() )
  {
    Cerr("***** castor-imageDynamicTools :: Please provide at least one data filename (-i)" << endl);
    ShowHelp(0);
    Exit(EXIT_FAILURE);
  }
  else
  {    
    if(vb >= 2)
    { 
      Cout(" Selected root data-file(s) to convert: " << endl);
      for (size_t i=0 ; i<path_to_input_image.size() ; i++) 
        Cout(path_to_input_image[i] << endl);
    }
  }
  
  // Output files
  if (path_fout.empty() && path_dout.empty())
  {
    Cerr("***** castor-imageDynamicTools() -> Please provide an output option for output files (-fout or -dout) !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Check that only one option has been provided
  if (!path_fout.empty() && !path_dout.empty())
  {
    Cerr("***** castor-imageDynamicTools() -> Please provide either output option -fout or -dout but not both !" << endl);
    Exit(EXIT_FAILURE);
  }

  // Check only one option has been enabled between -dynamic-model and -im
  if (dynamic_model_options != "" 
   &&        motion_options != "")
  {
    Cerr("***** castor-imageDynamicTools() -> Please provide either output option -dynamic-model or -im but not both !" << endl);
    Exit(EXIT_FAILURE);
  }
  

  // ============================================================================================================
  // Initialization:
  // ============================================================================================================

  // (Required for options using interfile I/O)
  GetUserEndianness();
    
  // ----------------------------------------------------------------------------------------
  // Create sOutputManager
  // ----------------------------------------------------------------------------------------
  sOutputManager* p_outputManager = sOutputManager::GetInstance();  
  // Set verbose level
  p_outputManager->SetVerbose(vb);
  
  // Set output dynamic image policy
  p_outputManager->SetMergeDynImagesFlag(merge_dynamic_imgs_flag);

  // Set output number precision
  p_outputManager->SetOutNbPrec(onb_prec);
  
  // Set path to the config directory
  if (p_outputManager->CheckConfigDir(""))
  {
    Cerr("***** castor-imageDynamicTools() -> A problem occurred while checking for the config directory path !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Initialize output directory and base name
  if (p_outputManager->InitOutputDirectory(path_fout, path_dout))
  {
    Cerr("***** castor-imageDynamicTools() -> A problem occurred while initializing output directory !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Log command line
  if (p_outputManager->LogCommandLine(argc,argv))
  {
    Cerr("***** castor-imageDynamicTools() -> A problem occurred while logging command line arguments !" << endl);
    Exit(EXIT_FAILURE);
  }


  // ============================================================================================================
  // Instanciate/Initialize image objects
  // ============================================================================================================
   
  Intf_fields IF;
  IntfKeyInitFields(&IF);

  int i=0; // May implement the application of the dynamic model to a series of images later
  if(IntfReadHeader(path_to_input_image[i], &IF, vb) )
  {
    Cerr("***** castor-imageDynamicTools() -> An error occurred while trying to read the interfile header of file reading file " << path_to_input_image[i] << " !" << endl);  
    Exit(EXIT_FAILURE);
  }

  // Get the frame duration from the image interfile
  vector <string> image_filenames;
  
  if ( IntfIsMHD(path_to_input_image[i], image_filenames) < 0 )
  {
    Cerr("***** castor-imageDynamicTools() -> Error : while trying to read the interfile header '"<< path_to_input_image[i] << "' !" << endl);
    Exit(EXIT_FAILURE);
  }
  
  // Generate frame list from image header metadata
  if(IF.image_duration.size()>0 || image_filenames.size() >0) // Image durations have been provided 
  {
    // We have one mhd file and a set of several 3D images
    if(image_filenames.size() > 1) 
    {
      // If image files are splitted over a series of 3D image files, we will have image durations for each image, including gates
      // As we need to recover the duration for each gate once, we skip gates if the data contains any      
      for(int fr=0 ; fr<IF.nb_time_frames ; fr++)
      {
        // First, initialize IF using the image corresponding to the 'fr' frame
        if(IntfReadHeader(image_filenames[fr], &IF, vb) )
        {
          Cerr("***** castor-imageDynamicTools() -> An error occurred while trying to read the interfile header of file reading file " << path_to_input_image[i] << " !" << endl);  
          Exit(EXIT_FAILURE);
        }
        
        // Recover image duration and dead frames (if any) in the frame list
        stringstream ss_duration, ss_startTime;
        
        // Go to the next frame (skip gate if any)
        ss_duration << IF.image_duration.at(fr * IF.nb_resp_gates* IF.nb_card_gates);
        ss_startTime << IF.image_start_time.at(fr * IF.nb_resp_gates* IF.nb_card_gates);
        string frame_duration_str = ss_duration.str();
        string frame_start_str = ss_startTime.str();

        // Building frame_list with as "StartTime:Duration,StartTime:Duration,...."
        if (stoi(frame_start_str)<delayInSec) nb_FramesToSkip ++;
        // Applied if the frame start time is greater than the provided Delay
        else
        {
          frame_list.append(frame_start_str);
          frame_list.append(":");
          frame_list.append(frame_duration_str);

          // add comma if not last frame
          if (fr + 1 < IF.nb_time_frames)
           frame_list.append(",");
        }
      }
    }
    // All dynamic images concatenated in one file 
    else
    {
      if(IF.image_duration.size() != IF.nb_time_frames)
      {
        Cerr("***** castor-imageDynamicTools() -> Interfile reading error of the input image :"<<endl);
        Cerr("      The number of provided image duration:  '"<< IF.image_duration.size() 
          << "'     does not match the actual number of time frames * respiratory gates * cardiac gates ) '"<< IF.nb_time_frames <<"' !" << endl);
        Exit(EXIT_FAILURE);
      }  

      //----------------------------------------------------------------
      for(int fr=0 ; fr<IF.nb_time_frames ; fr++)
      {
        // First frame : just add frame start
        if( fr==0  )
        {
          stringstream ss;
          // check if image_start_time key has been initialized 
          // (size == nb frames)
          // assume and init start time = 0 otherwise
          if (IF.image_start_time.size() != IF.nb_time_frames)
            IF.image_start_time.push_back(0);
    
          ss << IF.image_start_time.at(fr);
          frame_list.append(ss.str());
        }
        // If current frame starts after the previous frame
        // -> add ':' and frame duration of the previous frame
        // -> add start time of the actual frame
        else if( IF.image_start_time.size() == IF.nb_time_frames // image_start_time key provided
              && IF.image_start_time.at(fr) > ( IF.image_start_time.at(fr-1) + IF.image_duration.at(fr-1) )  )
        {
          stringstream ss;
          ss << IF.image_duration.at(fr-1);
          frame_list.append(":").append(ss.str());
          ss.str("");
          ss << IF.image_start_time.at(fr);
          frame_list.append(",").append(ss.str());
        }
        // Just add frame start otherwise
        else
        {
          stringstream ss;
          // init start time if not provided
          if (IF.image_start_time.size() != IF.nb_time_frames) 
          {
            FLTNB start_time = IF.image_start_time.at(fr-1) + IF.image_duration.at(fr-1);
            IF.image_start_time.push_back(start_time);
          }
          
          ss << IF.image_start_time.at(fr);
          frame_list.append(",").append(ss.str());
        }
        
        // If last frame : add frame duration
        if(fr+1 >= IF.nb_time_frames)
        {
          stringstream ss;
          ss << IF.image_duration.at(fr);
          frame_list.append(":").append(ss.str());
        }
      } // end of loop on frames
      
    }
    
    // Print out the skiped frame information
    if(vb >= 2) Cout ("----- Number of frames to skip from model fit: " <<  nb_FramesToSkip << endl);
    if(vb >= 2) Cout ("Frame list is : " << frame_list << endl);
  }

  // Instantiate & Initialize oImageDimensionsAndQuantification object, required for datafile generation (number of threads)
  oImageDimensionsAndQuantification* p_ID = new oImageDimensionsAndQuantification(); 

  // Check nb gating
  (IF.nb_resp_gates >1) ? nb_resp_gates = IF.nb_resp_gates : 1 ;
  (IF.nb_card_gates >1) ? nb_card_gates = IF.nb_card_gates : 1 ;
  
  // --- oImageDimensionsAndQuantification initialization ---
  p_ID->SetNbVoxX(IF.mtx_size[0]);
  p_ID->SetNbVoxY(IF.mtx_size[1]);
  p_ID->SetNbVoxZ(IF.mtx_size[2]);
  p_ID->SetNbThreads(nb_threads);
  p_ID->SetNbBeds(1);
  p_ID->SetVoxSizeX(IF.vox_size[0]);
  p_ID->SetVoxSizeY(IF.vox_size[1]);
  p_ID->SetVoxSizeZ(IF.vox_size[2]);
  p_ID->SetFOVOutMasking(0., 0);
  p_ID->SetFOVSizeX(-1.);
  p_ID->SetFOVSizeY(-1.);
  p_ID->SetFOVSizeZ(-1.);
  p_ID->SetOffsetX(0);
  p_ID->SetOffsetY(0);
  p_ID->SetOffsetZ(0);
  p_ID->SetVerbose(vb);
  p_ID->SetNbRespGates(nb_resp_gates);
  p_ID->SetNbCardGates(nb_card_gates);
  p_ID->SetFrames(frame_list);
  p_ID->SetnbFramesToSkip(nb_FramesToSkip);

  if (p_ID->CheckParameters())
  {
    Cerr("***** castor-imageDynamicTools :: A problem occurred while checking image dimensions parameters !" << endl);
    Exit(1);
  }
  if (p_ID->Initialize())
  {
    Cerr("***** castor-imageDynamicTools :: A problem occurred while initializing image dimensions !" << endl);
    Exit(1);
  }

  // Initialization of DynamicDataManager class, related 4D data splitting management 
  if (p_ID->InitDynamicData( "", 0, 0, 0, 1, 1 ) )
  {
    Cerr("***** castor-imageDynamicTools :: A problem occurred while initializing Dynamic data manager's class !" << endl);
    Exit(EXIT_FAILURE);
  }



  // ----------------------------------------------------------------------------------------
  // Create Dynamic Model Manager 
  // ----------------------------------------------------------------------------------------

  // Verbose
  if (vb>=5) Cout("----- Dynamic model initialization (if any) ... -----" << endl);
  // Create object
  oDynamicModelManager* p_DynamicModelManager = new oDynamicModelManager();
  // Set all parameters
  p_DynamicModelManager->SetImageDimensionsAndQuantification(p_ID);
  p_DynamicModelManager->SetOptions(dynamic_model_options);
  p_DynamicModelManager->SetVerbose(vb);
  // Check parameters
  if (p_DynamicModelManager->CheckParameters())
  {
    Cerr("***** castor-imageDynamicTools() -> A problem occurred while checking dynamic model manager's parameters !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Initialize optimizer manager
  if (p_DynamicModelManager->Initialize())
  {
    Cerr("***** castor-imageDynamicTools() -> A problem occurred while initializing dynamic model manager !" << endl);
    Exit(EXIT_FAILURE);
  }
  
  // Verbose
  if (vb>=5) Cout("----- Dynamic model initialization OK -----" << endl);


  // ----------------------------------------------------------------------------------------
  // Set Image Space and initialise Image
  // ----------------------------------------------------------------------------------------

  // --- Image space initialization ---
  oImageSpace* p_ImageSpace = new oImageSpace();

  p_ImageSpace->SetImageDimensionsAndQuantification(p_ID);
  p_ImageSpace->SetVerbose(vb);

  // Allocate memory for main image
  p_ImageSpace->InstantiateImage();
  p_ImageSpace->InstantiateForwardImage();
  p_ImageSpace->InstantiateBackwardImageFromDynamicBins();
  p_ImageSpace->InstantiateOutputImage();
  p_ImageSpace->InitImage(path_to_input_image[i], 0);


  // ----------------------------------------------------------------------------------------
  // Create Deformation Manager 
  // ----------------------------------------------------------------------------------------

  // Verbose
  if (vb>=5) Cout("----- Image deformation initialization (if any) ... -----" << endl);
  
  // Create object
  oDeformationManager* p_DeformationManager = new oDeformationManager();
  // Set all parameters
  p_DeformationManager->SetImageDimensionsAndQuantification(p_ID);
  p_DeformationManager->SetDataMode(0); // required to know if sensitivity image deformation should be enabled
  
  if (motion_options != "")
  {
    p_DeformationManager->SetOptions(motion_options);
    p_DeformationManager->SetNbTransformations(1); // 1 for testing
    p_DeformationManager->SetMotionType(DEF_IPAT_MOT);
  }
  else
    p_DeformationManager->SetNbTransformations(0);

  p_DeformationManager->SetVerbose(vb);
  // Check parameters
  if (p_DeformationManager->CheckParameters())
  {
    Cerr("***** castor-recon() -> A problem occurred while checking image deformation manager's parameters !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Initialize optimizer manager
  if (p_DeformationManager->Initialize())
  {
    Cerr("***** castor-recon() -> A problem occurred while initializing image deformation manager !" << endl);
    Exit(EXIT_FAILURE);
  }
  // Verbose
  if (vb >=5) Cout("----- Image deformation initialization OK -----" << endl);

  // Initial clock for execution time
  clock_t clock_start = clock();
  time_t time_start = time(NULL);
  
  // ============================================================================================================
  // Apply dynamic Model if initialized
  // ============================================================================================================ 
  if(dynamic_model_options != "")
  {
    for(uint32_t it=0; it<nb_ite_dyn ; it++)
    {
      Cout(" --- Apply Dynamic Model, iteration" <<  it+1 << " --- " << endl);
      if (p_DynamicModelManager->ApplyDynamicModel(p_ImageSpace, 0, 0))
      {
        Cerr("***** castor-imageDynamicTools -> A problem occurred while applying dynamic model to current estimate images !" << endl);
        return 1;
      }
    }
    
    // Save Parametric images
    if (p_DynamicModelManager->SaveParametricImages(-1))
    {
      Cerr("***** castor-imageDynamicTools -> A problem occurred while saving parametric images related to the dynamic model !" << endl);
      return 1;
    }
    
    if(output_dynamic_imgs_from_model)
    {
      p_ImageSpace->ComputeOutputImage();
      
      // Save output image
      if (p_ImageSpace->SaveOutputImage(-1))
      {
        Cerr("***** castor-imageDynamicTools -> A problem occurred while saving output image !" << endl);
        return 1;
      }
    }
    
    Cout(endl << " --- Parametric images generated at: " <<  ((path_fout.empty()) ? path_dout : path_fout) << " --- " << endl);
  }

  // ============================================================================================================
  // Test/Apply deformations if initialized
  // ============================================================================================================ 
  
  if(motion_options != "")
  {
    // Apply deformation on input image
    if (p_DeformationManager->TestDeformationOnImage(p_ImageSpace->m4p_image[0][0][0], p_ImageSpace->m4p_image[0][0][0], FORWARD_DEFORMATION, 0))
    {
      Cerr("***** castor-imageDynamicTools -> A problem occurred while applying dynamic model to current estimate images !" << endl);
      return 1;
    }
    
    p_ImageSpace->ComputeOutputImage();
      
    // Save output image
    if (p_ImageSpace->SaveOutputImage(-1))
    {
      Cerr("***** castor-imageDynamicTools -> A problem occurred while saving output image !" << endl);
      return 1;
    }
    
    Cout(endl << " --- Parametric images generated at: " <<  ((path_fout.empty()) ? path_dout : path_fout) << " --- " << endl);
  }
  
  
    // Final clock for execution time
  clock_t clock_stop = clock();
  time_t time_stop = time(NULL);
  
  // ============================================================================================================
  // End
  // ============================================================================================================

  if (vb>=2) Cout("oIterativeAlgorithm::IterateCPU() -> Total time spent | User: " << (FLTNB)(time_stop-time_start) 
                      << " sec | CPU: " << (clock_stop-clock_start)/((FLTNB)CLOCKS_PER_SEC) << " sec" << endl);
                      
                      
  // ============================================================================================================
  // End
  // ============================================================================================================

  if (vb>=2) Cout(" Deallocating objects ..." << endl);
      
  // Delete objects
  p_ImageSpace->DeallocateOutputImage();
  p_ImageSpace->DeallocateBackwardImageFromDynamicBins();
  p_ImageSpace->DeallocateForwardImage();
  p_ImageSpace->DeallocateImage();
  
  if(p_ImageSpace)       delete p_ImageSpace;
  if(p_ID)               delete p_ID;

  Cout(" ---          END         --- " << endl << endl);
  
  return 0;
}

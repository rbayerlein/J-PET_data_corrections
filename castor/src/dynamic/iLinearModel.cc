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
  \ingroup  dynamic
  \brief    Implementation of class iLinearModel
  \todo     Checks regarding optimization on NestedEM loops when updating parametric images
*/

#include "iLinearModel.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn iLinearModel
  \brief Constructor of iLinearModel. Simply set all data members to default values.
*/
iLinearModel::iLinearModel() : vDynamicModel() 
{
  m_nbRgateBF      = -1;
  m_nbCgateBF      = -1;
  m_nbRGModelParam = -1;
  m_nbCGModelParam = -1;

  m_OptimisationMethod = -1;
  m_nbWeightFactors=-1;

  m2p_respBasisFunctions        = NULL;
  m2p_cardBasisFunctions        = NULL;

  m_fileOptions = "";
  m_listOptions = "";

  mp_corrBasisCoeffs    = NULL;
  mp_corrBasisFunctions = NULL;
  m_nbLinearModelCycles       = 1;
  m_basisFunctionsUpdStartIte =-1;
  m_basisFunctionsUpdRatio    = 1;
  m_basisFunctionsUpdIdx      = 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ~iLinearModel
  \brief Destructor of iLinearModel
*/
iLinearModel::~iLinearModel() 
{
  if(m_initialized)
  {
    for(int rb=0 ; rb<m_nbRgateBF ; rb++)
    {
      if (m2p_respBasisFunctions[rb]) delete[] m2p_respBasisFunctions[rb];
      if (m2p_RGParametricImages[rb]) delete[] m2p_RGParametricImages[rb];
    }
  
    for(int cb=0 ; cb<m_nbCgateBF ; cb++)
    {
      if (m2p_cardBasisFunctions[cb]) delete[] m2p_cardBasisFunctions[cb];
      if (m2p_CGParametricImages[cb]) delete[] m2p_CGParametricImages[cb];
    }

    if (mp_corrBasisCoeffs != NULL) delete[] mp_corrBasisCoeffs; 
    if (mp_corrBasisFunctions != NULL) delete[]  mp_corrBasisFunctions;
    
    
    // Free NNLS variables
    if(m_OptimisationMethod == OPTIMISATION_METHOD_NNLS)
    {
      for( int th=0 ; th<mp_ID->GetNbThreadsForImageComputation() ; th++ )
      {
        for(int n=0 ; n<m_nnlsN ; n++)
          if ( m3p_nnlsA[ th ][ n ]) delete[] m3p_nnlsA[ th ][ n ];
          
        // Init 2D coefficient matrix for NNLS estimation
        if( m3p_nnlsA[ th ] ) delete[] m3p_nnlsA[ th ];

        // Init solution vector and working matrix for NNLS estimation
  
        if( m2p_nnlsB[ th ] )   delete[] m2p_nnlsB[ th ];
        if( m2p_nnlsMat[ th ] ) delete[] m2p_nnlsMat[ th ];
        if( m2p_nnlsX[ th ] )   delete[] m2p_nnlsX[ th ];
        if( m2p_nnlsWp[ th ] )  delete[] m2p_nnlsWp[ th ];
        if( m2p_nnlsIdx[ th ] ) delete[] m2p_nnlsIdx[ th ];
      }
  
      if( m3p_nnlsA )   delete[] m3p_nnlsA;
      if( m2p_nnlsB )   delete[] m2p_nnlsB;
      if( m2p_nnlsMat ) delete[] m2p_nnlsMat;
      if( m2p_nnlsX )   delete[] m2p_nnlsX;
      if( m2p_nnlsWp )  delete[] m2p_nnlsWp;
      if( m2p_nnlsIdx ) delete[] m2p_nnlsIdx;
      if( mp_w )        delete[] mp_w;
    }
  }
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ShowHelpModelSpecific
  \brief Print out specific help about the implementation of this model
         and its initialization
*/
void iLinearModel::ShowHelpModelSpecific()
{
  cout << "-- This class implements a general linear dynamic model applied between the images of a dynamic acquisition." << endl;
  cout << "-- The model is applied on a voxel-by-voxel basis between the images of the frames and/or respiratory/cardiac gates. " << endl;
  cout << "-- The main keywords 'DYNAMIC FRAMING', 'RESPIRATORY GATING' and 'CARDIAC GATING' ahead of the parameters allow to define at which level the model parameters must be applied. " << endl;
  cout << "-- Main parameters to define are:" << endl;
  cout << "   -> The number of basis functions / parametric images defined in the model" << endl;
  cout << "   -> Basis function initial values" << endl;
  cout << "   -> Parametric images initialization (optional)" << endl;
  cout << "   -> Optimisation_method : x     (mandatory) optimization method for voxelwise parameters estimation." << endl;
  cout << endl;
  cout << " It can be initialized using a configuration text file with the following keywords and information :" << endl; 
  cout << " - Mandatory keywords :" << endl;
  cout << "   The following keywords are mandatory for at least one dynamic image level (Dynamic frame, Respiratory gating or Cardiac gating) :" << endl;
  cout << "   'Number_basis_functions:'   Enter the number of basis function for each image of time frame or respiratory/cardiac gate " << endl;
  cout << "   'Basis_functions:'          Enter the basis function (bf) coefficients for each image (im) of time frame or respiratory/cardiac gate " << endl;
  cout << "                                 on successive lines, separated by ',' :" << endl;
  cout << "                                -> basis_functions: " << endl;
  cout << "                                -> coeff_bf1_im1,coeff_bf1_im2,...,coeff_bf1_imn" << endl;
  cout << "                                -> coeff_bf2_im1,coeff_bf2_im2,...,coeff_bf2_imn" << endl;
  cout << "                                -> etc..." << endl;
  cout << "   Optimisation_method : x     (mandatory) optimization method available options: " << endl;
  cout << "                                x=0: Direct ( Implementation of basis functions side by system matrix in each tomographic iterative loop " << endl;
  cout << "                                     (! Currently not compatible with motion correction)" << endl;
  cout << "                                x=1: Nested EM " << endl;
  cout << "                                x=2: Iterative non-negative Least-Square " << endl;
  cout << "                                    (C.L. Lawson and R.J. Hanson, Solving Least Squares Problems)" << endl;
  cout << endl;
  cout << " - Optional keywords :" << endl;
  cout << "   'Parametric_image_init: image_file' Set an image file to initialize the parametric images of each dynamic model. Default initialization: '1.0 for each voxel." << endl;
  cout << "   " << endl;
  cout << "   'Number_weight_values:'     Enter the number of weights to be applied for each dynamic frame for performing WLS optimisation (Optimisation method=2) " << endl;
  cout << "   'Weight_values:'            Enter the weight values to be applied for each dynamic frame ( within DYNAMIC FRAMING/ENDDF) " << endl;
  cout << "                                 on one single line, separated by ',' " << endl;
  cout << "   The previous parameters must be declared inside the couple of the following specific tags: " << endl;
  cout << "   - DYNAMIC FRAMING/ENDDF for chronological frame model  (dynamic model applied to chronological frames of a dynamic aquisition)" << endl;
  cout << "   - RESPIRATORY GATING/ENDRG for respiratory model       (dynamic model applied to respiratory gates of a dynamic aquisition)" << endl;
  cout << "   - CARDIAC GATING/ENDCG for cardiac model               (dynamic model applied to cardiac gates of a dynamic aquisition)" << endl;
  cout << "   Different levels of dynamic model can be enabled simultaneously (i.e a dynamic frame model can be used simultaneously with respiratory and/or cardiac gating model) " << endl;
  cout << "   " << endl;
  cout << "   " << endl;
  cout << "   The following keywords are optional and common to each dynamic model (Dynamic frame, Respiratory gating or Cardiac gating) :" << endl;
  cout << "   'Number_model_iterations: x'        Number of iterations of the model parameters and basis functions updates in one cycle of Nested EM" << endl;
  cout << "   (Default x ==  1)                      (one cycle consists in x iterations in which either the parametric images or the basis functions are updated) " << endl;
  cout << "                                         The ratio of parametric images / basis functions updates depends on the following parameter:" << endl;
  cout << "   'Basis_function_update_ratio: x'    Ratio of update between parametric images and basis functions updates cycle " << endl;
  cout << "   (Default x ==  0)                      Cycles consist in x iterations of the parametric images, following by x iterations of the basis functions" << endl;
  cout << "                                         If x == 0, only the parametric images are updated" << endl;
  cout << "   'Basis_function_start_ite : x'      Starting iteration for the update of basis functions " << endl;
  cout << "   (Default x == -1)                     If negative, no update of the basis functions is performed (only parametric images are updated) " << endl;
  cout << "   " << endl;
  cout << "   " << endl;
  cout << "   ---------------------------" << endl;
  cout << "   Example of initialization: " << endl;
  cout << "   DYNAMIC FRAMING" << endl;
  cout << "   Number_basis_functions     : 2" << endl;
  cout << "   Basis_functions            :" << endl;
  cout << "   23682.79, 25228.74, 26636.99, 27923.61, 29101.16" << endl;
  cout << "   5.4, 4.91, 4.48, 4.1, 3.75" << endl;
  cout << "   ENDDF" << endl;
  cout << "   " << endl;
  cout << "   RESPIRATORY GATING" << endl;
  cout << "   Number_basis_functions         : 6" << endl;
  cout << "   Basis_functions            :" << endl;  
  cout << "   1, 0.8, 0.6, 0.4, 0.2, 0.01" << endl;
  cout << "   0.7, 0.9, 0.7, 0.5, 0.3, 0.1" << endl;
  cout << "   0.4, 0.6, 0.8, 0.6, 0.4, 0.2" << endl;
  cout << "   0.2, 0.4, 0.6, 0.8, 0.6, 0.4" << endl;
  cout << "   0.1, 0.3, 0.5, 0.7, 0.9, 0.7" << endl;
  cout << "   0.01, 0.2, 0.4, 0.6, 0.8, 1" << endl;
  cout << "   ENDRG" << endl;
  cout << "   " << endl;
  cout << "   Number_model_iterations : 1" << endl;
  cout << "   Basis_function_start_ite : -1" << endl;
  cout << "   Basis_function_update_ratio : 0" << endl;
  cout << "   ---------------------------" << endl;
  cout << "   " << endl;

  // Print general help for all dynamic models
  ShowHelp();

}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ReadAndCheckConfigurationFileSpecific
  \brief This function is used to read options from a configuration file when the generic Linear Model is requested.
  \return 0 if success, other value otherwise.
*/

int iLinearModel::ReadAndCheckConfigurationFileSpecific()
{
  if(m_verbose >=3) Cout("iLinearModel::ReadAndCheckConfigurationFileSpecific ..."<< endl);

  // Apply the Generic linear Checks for all Linear Models
  if( ReadAndCheckConfigurationFileSpecificToAllLinearModels()==1)
  {
    Cerr("***** iLinearModel::ReadAndCheckConfigurationFileSpecific -> Error while trying to read and check specific configuration " << m_fileOptions << endl);
    return 1;
  }

  // Normal End
  return 0 ;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ReadAndCheckConfigurationFileSpecificToAllLinearModels
  \brief This function is used to read options from a configuration file that are generic to all linear dynamic models.
  \return 0 if success, other value otherwise.
*/
int iLinearModel::ReadAndCheckConfigurationFileSpecificToAllLinearModels()
{
  if(m_verbose >=3) Cout("iLinearModel::ReadAndCheckConfigurationFileSpecificToAllLinearModels ..."<< endl);
  
  // The file will be fully processed in the Initialize() function
  ifstream in_file(m_fileOptions.c_str(), ios::in);
  
  if(in_file)
  {
    // Check first the file contains the mandatory keyword(s)
    bool file_is_good = false;

    string line="", kword_to_search="";
    while(!in_file.eof())
    {
      getline(in_file, line);

      //remove comment
      if (line.find("#") != string::npos) line = line.substr(0, line.find_first_of("#")) ;

      if (line.find("DYNAMIC FRAMING") != string::npos)
      {
        if(kword_to_search != "")
        {
          Cerr("***** iLinearModel::ReadAndCheckConfigurationFileSpecificToAllLinearModels -> Error, found an end tag 'END**' before 'DYNAMIC FRAMING' in configuration file: " << m_fileOptions << endl);
          return 1;
        }
        else
          kword_to_search = "ENDDF";
      }

      else if (line.find("RESPIRATORY GATING") != string::npos)
      {
        if(kword_to_search != "")
        {
          Cerr("***** iLinearModel::ReadAndCheckConfigurationFileSpecificToAllLinearModels -> Error, found an end tag 'END**' before 'RESPIRATORY GATING' in configuration file: " << m_fileOptions << endl);
          return 1;
        }
        else
          kword_to_search = "ENDRG";
      }

      else if (line.find("CARDIAC GATING") != string::npos)
      {
        if(kword_to_search != "")
        {
          Cerr("***** iLinearModel::ReadAndCheckConfigurationFileSpecificToAllLinearModels -> Error, found an end tag 'END**' before 'CARDIAC GATING' in configuration file: " << m_fileOptions << endl);
          return 1;
        }
        else
          kword_to_search = "ENDCG";
      }

      else if (kword_to_search != ""
       && line.find(kword_to_search) != string::npos)
      {
        kword_to_search = "";
        file_is_good = true;
      }
    }

    if(!file_is_good)
    {
      Cerr("***** iLinearModel::ReadAndCheckConfigurationFileSpecificToAllLinearModels -> Error, no mandatory tags ('DYNAMIC FRAMING'/'ENDDF', 'RESPIRATORY GATING'/'ENDRG', 'CARDIAC GATING'/'ENDCG' detected in configuration file: "
            << m_fileOptions << ". Check help for more info" << endl);
      return 1;
    }

    if( ReadDataASCIIFile(m_fileOptions, "Number_basis_functions", &m_nbTimeBF, 1, KEYWORD_OPTIONAL, "DYNAMIC FRAMING", "ENDDF") == 1)
    {
      Cerr("***** iLinearModel::ReadAndCheckConfigurationFileSpecificToAllLinearModels -> Error while trying to read 'Number_basis_functions' flag in " << m_fileOptions << endl);
      return 1;
    }

    if( ReadDataASCIIFile(m_fileOptions, "Number_basis_functions", &m_nbRgateBF, 1, KEYWORD_OPTIONAL, "RESPIRATORY GATING", "ENDRG") == 1)
    {
      Cerr("***** iLinearModel::ReadAndCheckConfigurationFileSpecificToAllLinearModels -> Error while trying to read 'Number_basis_functions' flag in " << m_fileOptions << endl);
      return 1;
    }

    if( ReadDataASCIIFile(m_fileOptions, "Number_basis_functions", &m_nbCgateBF, 1, KEYWORD_OPTIONAL, "CARDIAC GATING", "ENDCG") == 1)
    {
      Cerr("***** iLinearModel::ReadAndCheckConfigurationFileSpecificToAllLinearModels -> Error while trying to read 'Number_basis_functions' flag in " << m_fileOptions << endl);
      return 1;
    }

    if( ReadDataASCIIFile(m_fileOptions, "Number_weight_values", &m_nbWeightFactors, 1, KEYWORD_OPTIONAL, "DYNAMIC FRAMING", "ENDDF") == 1)
    {
      Cerr("***** iLinearModel::ReadAndCheckConfigurationFileSpecificToAllLinearModels -> Error while trying to read 'Number_basis_functions' flag in " << m_fileOptions << endl);
      return 1;
    }

    if( ReadDataASCIIFile(m_fileOptions, "Number_model_iterations", &m_nbLinearModelCycles, 1, KEYWORD_OPTIONAL) == 1)
    {
      Cerr("***** iLinearModel::ReadAndCheckConfigurationFileSpecificToAllLinearModels -> Error while trying to read 'Number_model_iterations' flag in " << m_fileOptions << endl);
      return 1;
    }

    if( ReadDataASCIIFile(m_fileOptions, "Basis_function_start_ite", &m_basisFunctionsUpdStartIte, 1, KEYWORD_OPTIONAL) == 1)
    {
      Cerr("***** iLinearModel::ReadAndCheckConfigurationFileSpecificToAllLinearModels -> Error while trying to read 'Basis_function_start_ite' flag in " << m_fileOptions << endl);
      return 1;
    }

    if( ReadDataASCIIFile(m_fileOptions, "Basis_function_update_ratio", &m_basisFunctionsUpdRatio, 1, KEYWORD_OPTIONAL) == 1)
    {
      Cerr("***** iLinearModel::ReadAndCheckConfigurationFileSpecificToAllLinearModels -> Error while trying to read 'Basis_function_update_ratio' flag in " << m_fileOptions << endl);
      return 1;
    }
    // Optimisation method
    if( ReadDataASCIIFile(m_fileOptions, "Optimisation_method", &m_OptimisationMethod, 1, KEYWORD_MANDATORY))
    {
      Cerr("***** ***** iLinearModel::ReadAndCheckConfigurationFileSpecificToAllLinearModels  -> Error while trying to read 'Optimisation_method' keyword in " << m_fileOptions << endl);
      return 1;
    }

    // Print out information on optimisation method set
    if(m_verbose >=2)
    {
      if (m_OptimisationMethod == OPTIMISATION_METHOD_NNLS)
        Cout("iLinearModel::ReadAndCheckConfigurationFileSpecificToAllLinearModels() -> Selected optimization method : NNLS (2)" << endl);
      else if (m_OptimisationMethod == OPTIMISATION_METHOD_NESTEM)
        Cout("iLinearModel::ReadAndCheckConfigurationFileSpecificToAllLinearModels() -> Selected optimization method : Nested-EM (1)" << endl);
      else if (m_OptimisationMethod == OPTIMISATION_METHOD_DR)
        Cout("iLinearModel::ReadAndCheckConfigurationFileSpecificToAllLinearModels() -> Selected optimization method : Direct Dynamic (0)" << endl);
    }
  }
  else
  {
    Cerr("***** iLinearModel::ReadAndCheckConfigurationFileSpecificToAllLinearModels -> Error while trying to read configuration file at: " << m_fileOptions << endl);
    return 1;
  }


  // Initialize the number of parameters of the linear model as the number of (frame) basis functions
  // (Variable used during parametric images writing on disk)
  m_nbModelParam = m_nbTimeBF;
  m_nbRGModelParam = m_nbRgateBF;
  m_nbCGModelParam = m_nbCgateBF;
  
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ReadAndCheckOptionsList
  \brief This function is used to read parameters from a string.
  \return 0 if success, other value otherwise.
*/
int iLinearModel::ReadAndCheckOptionsList(string a_listOptions)
{
  if(m_verbose >=3) Cout("iLinearModel::ReadAndCheckOptionsList ..."<< endl); 
  
  // Just recover the string here, it will be processed in the Initialize() function
  m_listOptions = a_listOptions;

  
  // For now, just restricts the initialization using a configuration file as there are quite a lot of parameters to initialize for a use with command line options
  Cerr("***** iLinearModel::ReadAndCheckOptionsList() -> Initialization with command line options is not implemented for this class. Please use a configuration file instead" << endl);
  return 1;
  
  // Normal end
  //return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn CheckSpecificParameters
  \brief This function is used to check whether all member variables
         have been correctly initialized or not.
  \return 0 if success, positive value otherwise.
*/
int iLinearModel::CheckSpecificParameters()
{

  // Perform generic checks for the Linear Model
  if(CheckSpecificParametersForAllLinearModels())
  {
    Cerr("***** iLinearModel::CheckSpecificParameters -> A problem occurred while checking specific parameters ! " << endl);
    return 1;
  }

  // Normal end
  return 0;
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn CheckSpecificParametersForAllLinearModels
  \brief This function is used to check whether all member variables of a
  generic linear model have been correctly initialized or not.
  \return 0 if success, positive value otherwise.
*/

int iLinearModel::CheckSpecificParametersForAllLinearModels()
{

  if(m_verbose >=2) Cout("iLinearModel::CheckSpecificParametersForAllLinearModels ..."<< endl);

  // Check at least one basis function number has been initialized
  if (m_nbTimeBF<=0 && m_nbRgateBF<=0 && m_nbCgateBF<=0)
  {
    Cerr("***** iLinearModel::CheckSpecificParametersForAllLinearModels() -> Error, the variables corresponding to the number of basis function has not been initialized. There might be an error in the configuration process/file !" << endl);
    return 1;
  }
  else
  {
    // Set the other variables to 1 if not initialized
    if(m_nbTimeBF<0)  m_nbTimeBF =1;
    if(m_nbRgateBF<0) m_nbRgateBF=1;
    if(m_nbCgateBF<0) m_nbCgateBF=1;
  }

  // if weights have been provided for WLS check if their lenght is equal to number of time basis functions
  if (m_nbWeightFactors>1)
  {
    if (m_nbWeightFactors!=mp_ID->GetNbTimeFrames())
    {
      Cerr("***** iLinearModel::CheckSpecificParametersForAllLinearModels() -> Error, the number of weight factors doesn't match the number of input frames. There might be an error in the configuration process/file !" << endl);
      return 1;
    }
  }


  // Check if we have somehow both a file and a list of options for init...
  if(m_listOptions != "" && m_fileOptions != "")
  {
    Cerr("***** iLinearModel::CheckSpecificParametersForAllLinearModels -> Either a file or a list of options have to be selected to initialize the model, but not both ! " << endl);
    return 1;
  }

  // Check if we have no file not list of options for some reason...
  if(m_listOptions == "" && m_fileOptions == "")
  {
    Cerr("***** iLinearModel::CheckSpecificParametersForAllLinearModels -> Either a file or a list of options should have been provided at this point ! " << endl);
    return 1;
  }

  // In case of regular linear regression selected make sure a dynamic model with 2 basis functions is set
  if (m_OptimisationMethod==OPTIMISATION_METHOD_LS && m_nbTimeBF!=2 )
  {
    Cerr("***** iLinearModel::CheckSpecificParametersForAllLinearModels -> Optimization method Linear Regression (=3) only available when using two time basis functions ! " << endl);
    return 1;
  }

  // Check if the model has been called by imageDynamicTools and Direct optimization method has been selected
  if (m_OptimisationMethod==OPTIMISATION_METHOD_DR && !m_useModelInReconstruction )
  {
    Cerr("***** iLinearModel::CheckSpecificParametersForAllLinearModels -> Optimization method 'Direct' (=0) only available when using castor-recon ! " << endl);
    return 1;
  }

  // Normal end
  return 0;
}


// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn InitializeSpecific
  \brief This function is used to initialize the parametric images and basis functions
  \return 0 if success, other value otherwise.
*/

int iLinearModel::InitializeSpecificToAllLinearModels()
{
  if(m_verbose >=2) Cout("iLinearModel::InitializeSpecificToAllLinearModels ..."<< endl);

  // Forbid initialization without check
  if (!m_checked)
  {
    Cerr("***** iLinearModel::InitializeSpecificToAllLinearModels() -> Must call CheckParameters functions before Initialize() !" << endl);
    return 1;
  }

  // Initialize blacklisted voxels image to zero;
  mp_blackListedvoxelsImage = new FLTNB [mp_ID->GetNbVoxXYZ()];
  for (int v = 0; v < mp_ID->GetNbVoxXYZ(); v++)
  {
    mp_blackListedvoxelsImage[v] = 0.0;
  }

  // Allocate memory for Parametric images and  time basis functions of the model
  m2p_parametricImages = new FLTNB*[m_nbTimeBF];
  m2p_outputParImages  = new FLTNB*[m_nbTimeBF];
  m2p_nestedModelTimeBasisFunctions = new FLTNB*[m_nbTimeBF];

  for(int b=0 ; b<m_nbTimeBF ; b++)
  {
    m2p_nestedModelTimeBasisFunctions[b] = new FLTNB[mp_ID->GetNbTimeFrames()];
    m2p_parametricImages[b] = new FLTNB[mp_ID->GetNbVoxXYZ()];
  }

  m2p_respBasisFunctions = new FLTNB*[m_nbRgateBF];
  m2p_RGParametricImages = new FLTNB*[m_nbRgateBF];

  for(int rb=0 ; rb<m_nbRgateBF ; rb++)
  {
    m2p_respBasisFunctions[rb] = new FLTNB[mp_ID->GetNbRespGates()];
    m2p_RGParametricImages[rb] = new FLTNB[mp_ID->GetNbVoxXYZ()];
  }

  m2p_cardBasisFunctions = new FLTNB*[m_nbCgateBF];
  m2p_CGParametricImages = new FLTNB*[m_nbCgateBF];


  for(int cb=0 ; cb<m_nbCgateBF ; cb++)
  {
    m2p_cardBasisFunctions[cb] = new FLTNB[mp_ID->GetNbCardGates()];
    m2p_CGParametricImages[cb] = new FLTNB[mp_ID->GetNbVoxXYZ()];
  }

  // Memory allocation for correction images
  mp_corrBasisCoeffs = new FLTNB[mp_ID->GetNbVoxXYZ()];
  mp_corrBasisFunctions = new FLTNB[mp_ID->GetNbThreadsForProjection()];

  // Allocate output image matrices
  for(int b=0 ; b<m_nbTimeBF ; b++)
    m2p_outputParImages[b] = new FLTNB[mp_ID->GetNbVoxXYZ()];


  // TODO: Consider other cases where linear models will be applied for Cardiac / Respiratory reconstruction.

  // --- Default Initialization basis functions --- //
  for(int b=0 ; b<m_nbTimeBF ; b++)
    for (int fr=0 ; fr<mp_ID->GetNbTimeFrames() ; fr++)
      m2p_nestedModelTimeBasisFunctions[b][fr] = 1.;

  // --- Default Initialization of respiratory and cardiac basis functions  --- //
  for(int rb=0 ; rb<m_nbRgateBF ; rb++)
    for (int rg=0 ; rg<mp_ID->GetNbRespGates() ; rg++)
      m2p_respBasisFunctions[rb][rg] = 1.;

  for(int cb=0 ; cb<m_nbCgateBF ; cb++)
    for (int cg=0 ; cg<mp_ID->GetNbCardGates() ; cg++)
      m2p_cardBasisFunctions[cb][cg] = 1.;


  // Initialize gate parametric images of linear model mother class
  for(int v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
  {
    m2p_RGParametricImages[ 0 ][ v ] = 1.;
    m2p_CGParametricImages[ 0 ][ v ] = 1.;
  }

  // Allocate NNLS variables and weights.
  if(m_OptimisationMethod == OPTIMISATION_METHOD_NNLS)
  {
    m3p_nnlsA    = new FLTNB** [mp_ID->GetNbThreadsForImageComputation()];
    m2p_nnlsB    = new FLTNB*  [mp_ID->GetNbThreadsForImageComputation()];
    m2p_nnlsMat  = new FLTNB*  [mp_ID->GetNbThreadsForImageComputation()];
    m2p_nnlsX    = new FLTNB*  [mp_ID->GetNbThreadsForImageComputation()];
    m2p_nnlsWp   = new FLTNB*  [mp_ID->GetNbThreadsForImageComputation()];
    m2p_nnlsIdx  = new int*    [mp_ID->GetNbThreadsForImageComputation()];
    mp_w = new FLTNB[ mp_ID->GetNbTimeFrames() ];
  
    // Get the largest nb of parameters and samples between time frames, resp or cardiac gates
    // as all could be used as nnls samples depending on the linear model
    int nnls_max_params = (m_nbTimeBF  > m_nbRgateBF) ?
                        ( (m_nbTimeBF  > m_nbCgateBF) ? m_nbTimeBF  : m_nbCgateBF) :
                        ( (m_nbRgateBF > m_nbCgateBF) ? m_nbRgateBF : m_nbCgateBF) ;
                        
    m_nnlsN = nnls_max_params;
    
    int nnls_max_samples = (mp_ID->GetNbTimeFrames() > mp_ID->GetNbRespGates()) ?
                         ( (mp_ID->GetNbTimeFrames() > mp_ID->GetNbCardGates()) ? mp_ID->GetNbTimeFrames() : mp_ID->GetNbCardGates()) :
                         ( (mp_ID->GetNbRespGates()  > mp_ID->GetNbCardGates()) ? mp_ID->GetNbRespGates()  : mp_ID->GetNbCardGates()) ;
    
    
    for( int th=0 ; th<mp_ID->GetNbThreadsForImageComputation() ; th++ )
    {
      // Init 2D coefficient matrix for NNLS estimation
      m3p_nnlsA[ th ] = new FLTNB*[ m_nnlsN ];

      for(int n=0 ; n<m_nnlsN ; n++)
        m3p_nnlsA[ th ][n] = new FLTNB[ nnls_max_samples ];

      // Init solution vector and working matrix for NNLS estimation

      m2p_nnlsB[ th ]   = new FLTNB[ nnls_max_samples ];
      m2p_nnlsMat[ th ] = new FLTNB[ (m_nnlsN+2) * nnls_max_samples ];
      m2p_nnlsX[ th ]   = new FLTNB[ m_nnlsN ];
      m2p_nnlsWp[ th ]  = new FLTNB[ m_nnlsN ];
      m2p_nnlsIdx[ th ] = new int[ m_nnlsN ];
    }

    // Fill weights for NNLS to default values
    // todo: specific weights for gate-related models (gate duration) - Also NECR related weights or user set weights (ZC)

    // Case of no weights input - set all weights to 1
    if (m_nbWeightFactors<=1)
      for(int t=0; t<mp_ID->GetNbTimeFrames(); t++)
        mp_w[t] = 1;
  }


  // --- Data Initialization with a configuration file --- //
  if(m_fileOptions != "")
  {
    ifstream in_file(m_fileOptions.c_str(), ios::in);

    if(in_file)
    {
      // TODO: Implement a check of length of input basis functions
      // Frame basis functions Initialization
      if(m_nbTimeBF>1
          && ReadDataASCIIFile(m_fileOptions,
                               "Basis_functions",
                               m2p_nestedModelTimeBasisFunctions,
                               mp_ID->GetNbTimeFrames(),
                               m_nbTimeBF,
                               KEYWORD_OPTIONAL,
                               "DYNAMIC FRAMING",
                               "ENDDF")==1)
      {
        Cerr("***** iLinearModel::Initialize -> Error while trying to read frame basis functions coefficients !" << endl);
        Cerr("                                  'Basis_functions' keyword inside DYNAMIC FRAMING / ENDDF paragraph in " << m_fileOptions << endl);
        return 1;
      }

      // Resp gates basis functions Initialization
      if(m_nbRgateBF>1
          && ReadDataASCIIFile(m_fileOptions,
                               "Basis_functions",
                               m2p_respBasisFunctions,
                               mp_ID->GetNbRespGates(),
                               m_nbRgateBF,
                               KEYWORD_OPTIONAL,
                               "RESPIRATORY GATING",
                               "ENDRG") ==1 )
      {
        Cerr("***** iLinearModel::Initialize -> Error while trying to read respiratory gates basis functions coefficients !" << endl);
        Cerr("                                  'Basis_functions' keyword inside RESPIRATORY GATING / ENDRG paragraph in " << m_fileOptions << endl);
        return 1;
      }

      // Card gates basis functions Initialization
      if(m_nbCgateBF>1
          && ReadDataASCIIFile(m_fileOptions,
                               "Basis_functions",
                               m2p_cardBasisFunctions,
                               mp_ID->GetNbCardGates(),
                               m_nbCgateBF,
                               KEYWORD_OPTIONAL,
                               "CARDIAC GATING",
                               "ENDCG") ==1 )
      {
        Cerr("***** iLinearModel::Initialize -> Error while trying to read cardiac gates basis functions coefficients !" << endl);
        Cerr("                                  'Basis_functions' keyword inside CARDIAC GATING / ENDCG paragraph in " << m_fileOptions << endl);
        return 1;
      }

      // Optimisation method
      if( ReadDataASCIIFile(m_fileOptions,
                            "Optimisation_method",
                            &m_OptimisationMethod,
                            1,
                            KEYWORD_MANDATORY))
      {
        Cerr("***** iLinearModel::Initialize -> Error while trying to read 'Optimisation_method' keyword in " << m_fileOptions << endl);
        return 1;
      }

      //

      if(m_nbWeightFactors>1
          && ReadDataASCIIFile(m_fileOptions,
                               "Weight_values",
                               mp_w,
                               m_nbWeightFactors,
                               KEYWORD_OPTIONAL,
                               "DYNAMIC FRAMING",
                               "ENDDF")==1)
      {
        Cerr("***** iLinearModel::Initialize -> Error while trying to read frame weight factor coefficients !" << endl);
        Cerr("                                  'Weight_values' keyword inside DYNAMIC FRAMING / ENDDF paragraph in " << m_fileOptions << endl);
        return 1;
      }


      // --- Parametric image initialization --- //

      // Frame model
      string input_image = "";
      int return_value = 0;

      return_value = ReadDataASCIIFile(m_fileOptions,
                                       "Parametric_images_init",
                                       &input_image,
                                       1,
                                       KEYWORD_OPTIONAL,
                                       "DYNAMIC FRAMING",
                                       "ENDDF");

      if( return_value == 0) // Image have been provided
      {
        // Read image // INTF_LERP_DISABLED = interpolation disabled for input image reading
        if( IntfReadImgDynCoeffFile(input_image,
                                    m2p_parametricImages,
                                    mp_ID,
                                    m_nbModelParam,
                                    m_verbose,
                                    INTF_LERP_DISABLED) ) // Image have been provided
        {
          Cerr("***** iLinearModel::Initialize -> Error while trying to read the provided initialization parametric images : " << input_image << endl);
          return 1;
        }
      }
      else if( return_value == 1) // Error during reading
      {
        Cerr("***** iLinearModel::Initialize -> Error while trying to read dynamic frame model parametric images !" << endl);
        Cerr("                                  'Parametric_image_init' keyword in " << m_fileOptions << endl);
        return 1;
      }
      else //(return_value >= 1 ) // Keyword not found : no initialization provided
      {
        // Standard initialization
        for(int b=0 ; b<m_nbTimeBF ; b++)
          for(int v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
            m2p_parametricImages[b][v] = 1.;
      }


      // Resp gate model
      input_image = "";
      return_value = 0;

      return_value = ReadDataASCIIFile(m_fileOptions,
                                       "Parametric_images_init",
                                       &input_image,
                                       1,
                                       KEYWORD_OPTIONAL,
                                       "RESPIRATORY GATING",
                                       "ENDRG");

      if( return_value == 0) // Image have been provided
      {
        // Read image // INTF_LERP_DISABLED = interpolation disabled for input image reading
        if( IntfReadImgDynCoeffFile(input_image,
                                    m2p_RGParametricImages,
                                    mp_ID,
                                    m_nbRGModelParam,
                                    m_verbose,
                                    INTF_LERP_DISABLED) ) // Image have been provided
        {
          Cerr("***** iLinearModel::Initialize -> Error while trying to read the provided initialization parametric images : " << input_image << endl);
          return 1;
        }
      }
      else if( return_value == 1) // Error during reading
      {
        Cerr("***** iLinearModel::Initialize -> Error while trying to read respiratory gate model parametric images !" << endl);
        Cerr("                                  'Parametric_image_init' keyword in " << m_fileOptions << endl);
        return 1;
      }
      else //(return_value >= 1 ) // Keyword not found : no initialization provided
      {
        // Standard initialization
        for(int rb=0 ; rb<m_nbRgateBF ; rb++)
          for (int v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
            m2p_RGParametricImages[rb][v] = 1.;
      }



      // Card gate model
      input_image = "";
      return_value = 0;

      return_value = ReadDataASCIIFile(m_fileOptions,
                                       "Parametric_images_init",
                                       &input_image,
                                       1,
                                       KEYWORD_OPTIONAL,
                                       "CARDIAC GATING",
                                       "ENDCG");

      if( return_value == 0) // Image have been provided
      {
        // Read image // INTF_LERP_DISABLED = interpolation disabled for input image reading
        if( IntfReadImgDynCoeffFile(input_image,
                                    m2p_CGParametricImages,
                                    mp_ID,
                                    m_nbCGModelParam,
                                    m_verbose,
                                    INTF_LERP_DISABLED) ) // Image have been provided
        {
          Cerr("***** iLinearModel::Initialize -> Error while trying to read the provided initialization parametric images : " << input_image << endl);
          return 1;
        }
      }
      else if( return_value == 1) // Error during reading
      {
        Cerr("***** iLinearModel::Initialize -> Error while trying to read cardiac gate model parametric images !" << endl);
        Cerr("                                  'Parametric_image_init' keyword in " << m_fileOptions << endl);
        return 1;
      }
      else //(return_value >= 1 ) // Keyword not found : no initialization provided
      {
        // Standard initialization
        for(int cb=0 ; cb<m_nbCgateBF ; cb++)
          for (int v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
            m2p_CGParametricImages[cb][v] = 1.;
      }
    }

    else
    {
      Cerr("***** iLinearModel::Initialize() -> Error while trying to read configuration file at: " << m_fileOptions << endl);
      return 1;
    }
  }

  // If weight values have been provided check they are not negative
  if (m_nbWeightFactors>1)
    for(int t=0; t<mp_ID->GetNbTimeFrames(); t++)
    {
     if (mp_w[t]<=0)
     {
      Cerr("***** iLinearModel::CheckSpecificParametersForAllLinearModels() -> Error, negative weight factors found. Only positive non-zero factors allowed !" << endl);
      return 1;
     }
     if(m_verbose >=4) Cout("iLinearModel::NNLS optimisation weight for frame: "<< t << " set at "<< mp_w[t]<< endl);
    }



  // --- Data Initialization with a list of options --- //

  if(m_listOptions != "")
  {
    // TODO
  }

  // Allocate output image matrices
  for(int b=0 ; b<m_nbTimeBF ; b++)
    m2p_outputParImages[b] = new FLTNB[mp_ID->GetNbVoxXYZ()];

  // If optimisation method requires an input of basis functions set m_ModelSpecificBasisFunctionsRequired flat
  if (m_OptimisationMethod==OPTIMISATION_METHOD_DR)
  {
    m_ModelSpecificBasisFunctionsRequired = true;
    m_noImageUpdateFlag = true;
  }

  // TODO : print output parametric images rg et cg


  // Normal End
  return 0;
}


// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ShowBasisFunctions
  \brief This function is used to initialize the parametric images and basis functions
  \return 0 if success, other value otherwise.
*/

void iLinearModel::ShowBasisFunctions()
{

  // Display Basis Functions set for the model
  if ((m_verbose>=2) & (m_nbTimeBF>1))
  {
    Cout("***** iLinearModel::InitializeSpecificToAllLinearModels() -> Time Basis Function coefficients :" << endl);
    for(int b=0 ; b<m_nbTimeBF ; b++)
    {
      Cout("                              ");
      Cout("Basis function["<<b+1<<"] : ");
      for(int fr=0 ; fr<mp_ID->GetNbTimeFrames()-1 ; fr++)
      {
        Cout(m2p_nestedModelTimeBasisFunctions[b][fr] << ", ");
      }
      // Print last value
      Cout(m2p_nestedModelTimeBasisFunctions[b][mp_ID->GetNbTimeFrames()-1] << endl);
    }
  }


}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn InitializeSpecific
  \brief This function is used to initialize the parametric images and basis functions
  \return 0 if success, other value otherwise.
*/
int iLinearModel::InitializeSpecific()
{
  if(m_verbose >=2) Cout("iLinearModel::InitializeSpecific ..."<< endl); 

  // Forbid initialization without check
  if (!m_checked)
  {
    Cerr("***** oDynamicModelManager::InitializeSpecific() -> Must call CheckParameters functions before Initialize() !" << endl);
    return 1;
  }

  // Run generic Initialization for all Linear Models
  if (InitializeSpecificToAllLinearModels())
  {
    Cerr("***** iLinearModel::InitializeSpecific() -> Error while performing generic initialisations for linear models !" << endl);
    return 1;
  }

  
  // Normal end
  m_initialized = true;
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn EstimateModelParameters
  \param ap_ImageS : pointer to the ImageSpace
  \param a_ite : index of the actual iteration (not used)
  \param a_sset : index of the actual subset (not used)
  \brief Estimate model parameters (parametric images and basis functions)
  \return 0 if success, other value otherwise.
*/
int iLinearModel::EstimateModelParameters(oImageSpace* ap_ImageS, int a_ite, int a_sset) 
{

  #ifdef CASTOR_DEBUG
    if (!m_initialized)
    {
      Cerr("***** iLinearModel::EstimateModelParameters() -> Called while not initialized !" << endl);
      Exit(EXIT_DEBUG);
    }
  #endif

  if(m_verbose >=3) Cout(" iLinearModel::EstimateModelParameters -> Optimisation method to use: " << m_OptimisationMethod << endl) ;

  // Use the nested EM parametric image estimation method
  if(m_OptimisationMethod == OPTIMISATION_METHOD_NESTEM)
  {
    if(NestedEM(ap_ImageS, a_ite) )
    {
      Cerr("***** iLinearModel::EstimateModelParameters() -> An error occured while using the nested EM parametric image estimation method !" << endl);
      return 1;
    }
  }
  // Direct reconstruction ( non nested ) - no update required
  else if(m_OptimisationMethod == OPTIMISATION_METHOD_DR)
  {
    Cout("iLinearModel::EstimateModelParameters() -> Skipping nested calculation, not required for direct method " << endl);
  }

  // Least Square for Post Reconsutrction Patlak Analysis
  else if(m_OptimisationMethod == OPTIMISATION_METHOD_NNLS )
  {
    if(EstimateParametersWithNNLS(ap_ImageS, a_ite) )
    {
      Cerr("***** iLinearModel::EstimateModelParameters() -> An error occured while using the Patlak LS Estimation Method !" << endl);
      return 1;
    }
  }
    // least-square linear regression
  else if(m_OptimisationMethod == OPTIMISATION_METHOD_LS )
  {
    if(Patlak_LS(ap_ImageS, a_ite) )
    {
      Cerr("***** iLinearModel::EstimateModelParameters() -> An error occured while using the Least-Squares linear regression Estimation Method !" << endl);
      return 1;
    }
  }

  else
  {
    Cerr("***** iLinearModel::EstimateModelParameters() -> Error : unknown method to estimate images ! !" << endl);
    return 1;
  }

  if(m_verbose >=3)
  {
    for (int fb=0 ; fb<m_nbTimeBF ; fb++)
    {
      FLTNB avg = 0;
      for (int v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
        if(mp_maskModel != NULL && mp_maskModel[ v ] == 1)
          avg += m2p_parametricImages[fb][v];

      avg /= m_nbVoxelsMask;
      Cout( "iLinearModel::EstimateModelParameters() -> Frame parametric image["<<fb<<"] avg value:" << avg << endl);
    }
  }

  return 0;
}
  

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn NestedEM
  \param ap_ImageS : pointer to the ImageSpace
  \param a_ite : index of the actual iteration (not used)
  \brief Estimate parametric images and basis functions (if enabled) using the nested EM method
  \return 0 if success, other value otherwise.
*/
int iLinearModel::NestedEM(oImageSpace* ap_ImageS, int a_ite) 
{
  if(m_verbose >=3) Cout("iLinearModel::NestedEM() ..." <<endl);
  
  for (uint32_t it=0 ; it<m_nbLinearModelCycles ; it++)
  {
    if(m_verbose >=3) Cout("iLinearModel::NestedEM() cycle "<< it+1 << "/" << m_nbLinearModelCycles <<endl);
    
    // Step 1 : Generate a difference image from voxel-by-voxel division of the current estimation of the image and the model image generated from the current parametric images / basis functions of the model
    // The backward image matrix (which is useless at this point of the reconstruction) is used as provisional image to gather the model image
    
    int v;
    #pragma omp parallel for private(v) schedule(static, 1)
    for (int fr=0 ; fr<mp_ID->GetNbTimeFrames() ; fr++)
      for (int rg=0 ; rg<mp_ID->GetNbRespGates() ; rg++)
        for (int cg=0 ; cg<mp_ID->GetNbCardGates() ; cg++)
          for (v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
          {
            // If a mask has been provided, check if the model applies in this voxel
            if(mp_maskModel != NULL && mp_maskModel[ v ] == 0)
              continue;
      
            // Reset this voxel to 0
            ap_ImageS->m6p_backwardImage[0][0][fr][rg][cg][v] = 0;
  
            for (int fb=0 ; fb<m_nbTimeBF ; fb++)
              for (int rb=0 ; rb<m_nbRgateBF ; rb++)
                for (int cb=0 ; cb<m_nbCgateBF ; cb++)
                {
                  // Retrieve current estimation of image according to coeffs/basis functions
                  ap_ImageS->m6p_backwardImage[0][0][fr][rg][cg][v] += m2p_parametricImages[fb][v] 
                                                                     * m2p_nestedModelTimeBasisFunctions[fb][fr]
                                                                     * m2p_RGParametricImages[rb][v] 
                                                                     * m2p_respBasisFunctions[rb][rg]
                                                                     * m2p_CGParametricImages[cb][v] 
                                                                     * m2p_cardBasisFunctions[cb][cg] ;
                }
              
            // Recover correction images using the ratio of the current estimation of the image (m4p_image) and the image generated with coffs/basis functions
            // (Again, use backward image as provisional image to recover the result)

            // TODO: Find where NaN is generated
            if (ap_ImageS->m6p_backwardImage[0][0][fr][rg][cg][v] == 0 )
            {
              mp_blackListedvoxelsImage[v] = 1.0;
              ap_ImageS->m6p_backwardImage[0][0][fr][rg][cg][v] = 1;
            }
            else
            {
              ap_ImageS->m6p_backwardImage[0][0][fr][rg][cg][v] =
                  ap_ImageS->m4p_image[fr][rg][cg][v] / ap_ImageS->m6p_backwardImage[0][0][fr][rg][cg][v];
            }
          }

    // Step 2 : Estimate either basis functions or parametric images
    
    // Step 2a : Basis functions estimation // TODO: shouldn't this condition be a single & as we want all conditions ??
    if( m_basisFunctionsUpdStartIte >= 0      // if (m_basisFunctionsUpdStartIte < 0) --> no update of basis functions, just update the parametric images
    &&  m_basisFunctionsUpdStartIte <= a_ite // if (m_basisFunctionsUpdStartIte >= 0) --> Check condition of minimal iteration before starting to update the basis functions 
    && !(int((m_basisFunctionsUpdIdx-1)/m_basisFunctionsUpdRatio)&1) ) // Check if we must estimate basis functions or coefficients at this stage (depends on m_basisFunctionsUpdRatio)
    {

      FLTNB parametric_image_norm = 0;
    
      if( m_nbTimeBF>1 )
      {
        if(m_verbose >=3) Cout("iLinearModel::NestedEM() -> Estimate Basis Functions - Frame basis functions estimation step" <<endl);
          
        for (int fb=0 ; fb<m_nbTimeBF ; fb++) 
        { 
          // Compute normalization related to the vowelwise time basis functions coefficients.
          parametric_image_norm = 0;
    
          // Compute normalization related to the parametric images.
          for (int rb=0 ; rb<m_nbRgateBF ; rb++)
            for (int cb=0 ; cb<m_nbRgateBF ; cb++)
                for (int v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
                  if(mp_maskModel != NULL && mp_maskModel[ v ] == 1)
                    parametric_image_norm += m2p_parametricImages[fb][v]
                                           * m2p_RGParametricImages[rb][v]
                                           * m2p_CGParametricImages[cb][v];
                  
          parametric_image_norm *= mp_ID->GetNbRespGates()
                                 * mp_ID->GetNbCardGates();
                                 
          for (int fr=0 ; fr<mp_ID->GetNbTimeFrames() ; fr++)
          {
            // Initialization of corrections factor for time basis functions 
            for (int th=0 ; th<mp_ID->GetNbThreadsForImageComputation() ; th++)
              mp_corrBasisFunctions[th]=0;
    
            for (int rg=0 ; rg<mp_ID->GetNbRespGates() ; rg++)
              for (int cg=0 ; cg<mp_ID->GetNbCardGates() ; cg++)
                for (int rb=0 ; rb<m_nbRgateBF ; rb++)
                  for (int cb=0 ; cb<m_nbCgateBF ; cb++)
                  {
                    int v;
                    #pragma omp parallel for private(v) schedule(static, 1)
                    for (v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
                    {
                      // If a mask has been provided, check if the model applies in this voxel
                      if(mp_maskModel != NULL && mp_maskModel[ v ] == 0)
                        continue;
              
                      int th = 0;
                      #ifdef CASTOR_OMP 
                      th = omp_get_thread_num();
                      #endif
                      
                      mp_corrBasisFunctions[th] += ap_ImageS->m6p_backwardImage[0][0][fr][rg][cg][v] 
                                                 * m2p_parametricImages[fb][v] 
                                                 * m2p_RGParametricImages[rb][v]
                                                 * m2p_CGParametricImages[cb][v];
                    }
                  }
    
            // Reduce
            for (int th=1 ; th<mp_ID->GetNbThreadsForImageComputation() ; th++)
              mp_corrBasisFunctions[0] += mp_corrBasisFunctions[th];
              
            // Apply corrections and normalization to the temporal basis functions values
            if (mp_corrBasisFunctions[0] > 0.)
              m2p_nestedModelTimeBasisFunctions[fb][fr] *= mp_corrBasisFunctions[0]/parametric_image_norm;
          }
          
        }
      }
      
      if( m_nbRgateBF>1 )
      {
        if(m_verbose >=3) Cout("iLinearModel::NestedEM() -> Estimate Basis Functions - Respiratory gate basis functions estimation step" <<endl); 
          
        for (int rb=0 ; rb<m_nbRgateBF ; rb++) 
        {
          // Compute normalization related to the vowelwise time basis functions coefficients.
          parametric_image_norm = 0;
          
          for (int fb=0 ; fb<m_nbTimeBF ; fb++)
            for (int cb=0 ; cb<m_nbRgateBF ; cb++)
                for (int v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
                  if(mp_maskModel != NULL && mp_maskModel[ v ] == 1)
                    parametric_image_norm += m2p_parametricImages[fb][v]
                                           * m2p_RGParametricImages[rb][v]
                                           * m2p_CGParametricImages[cb][v];
                    
          parametric_image_norm *= mp_ID->GetNbTimeFrames()
                                 * mp_ID->GetNbCardGates();
                                 
          for (int rg=0 ; rg<mp_ID->GetNbRespGates() ; rg++)
          {
            // Initialization of  corrections for time basis functions 
            for (int th=0 ; th<mp_ID->GetNbThreadsForImageComputation() ; th++)
              mp_corrBasisFunctions[th]=0;
              
            for (int fr=0 ; fr<mp_ID->GetNbTimeFrames() ; fr++)
              for (int cg=0 ; cg<mp_ID->GetNbCardGates() ; cg++)
                for (int fb=0 ; fb<m_nbTimeBF ; fb++)
                  for (int cb=0 ; cb<m_nbCgateBF ; cb++)
                  {
                    int v;
                    #pragma omp parallel for private(v) schedule(static, 1)
                    for (v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
                    {
                      // If a mask has been provided, check if the model applies in this voxel
                      if(mp_maskModel != NULL && mp_maskModel[ v ] == 0)
                        continue;
                        
                      int th = 0;
                      #ifdef CASTOR_OMP
                      th = omp_get_thread_num();
                      #endif
                      mp_corrBasisFunctions[th] += ap_ImageS->m6p_backwardImage[0][0][fr][rg][cg][v]
                                                 * m2p_parametricImages[fb][v] 
                                                 * m2p_RGParametricImages[rb][v]
                                                 * m2p_CGParametricImages[cb][v];
                    }
                  }
            
            // Reduce
            for (int th=1 ; th<mp_ID->GetNbThreadsForImageComputation() ; th++)
              mp_corrBasisFunctions[0] += mp_corrBasisFunctions[th];
              
            // Apply corrections and normalization to the temporal basis functions values 
            if (mp_corrBasisFunctions[0] > 0.)
              m2p_respBasisFunctions[rb][rg] *= mp_corrBasisFunctions[0]/parametric_image_norm;
            
          }
        }
      }
          
          
      if( m_nbCgateBF>1 )
      {
        if(m_verbose >=3) Cout("iLinearModel::NestedEM() -> Estimate Basis Functions - Cardiac gate basis functions estimation step" <<endl);
          
        for (int cb=0 ; cb<m_nbCgateBF ; cb++) 
        {
          // Compute normalization related to the vowelwise time basis functions coefficients.
          parametric_image_norm = 0;
          
          for (int fb=0 ; fb<m_nbTimeBF ; fb++)
            for (int rb=0 ; rb<m_nbRgateBF ; rb++) 
                for (int v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
                  if(mp_maskModel != NULL && mp_maskModel[ v ] == 1)
                    parametric_image_norm += m2p_parametricImages[fb][v]
                                           * m2p_RGParametricImages[rb][v]
                                           * m2p_CGParametricImages[cb][v];
                                         
          parametric_image_norm *= mp_ID->GetNbTimeFrames()
                                 * mp_ID->GetNbRespGates();
                    
          for (int cg=0 ; cg<mp_ID->GetNbCardGates() ; cg++)
          {
            // Initialization of  corrections for time basis functions 
            for (int th=0 ; th<mp_ID->GetNbThreadsForImageComputation() ; th++)
              mp_corrBasisFunctions[th]=0;
              
            for (int fr=0 ; fr<mp_ID->GetNbTimeFrames() ; fr++)
              for (int rg=0 ; rg<mp_ID->GetNbRespGates() ; rg++)
                for (int fb=0 ; fb<m_nbTimeBF ; fb++)
                  for (int rb=0 ; rb<m_nbRgateBF ; rb++) 
                  {
                    int v;
                    #pragma omp parallel for private(v) schedule(static, 1)
                    for (v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
                    {
                      // If a mask has been provided, check if the model applies in this voxel
                      if(mp_maskModel != NULL && mp_maskModel[ v ] == 0)
                        continue;
              
                      int th = 0;
                      #ifdef CASTOR_OMP
                      th = omp_get_thread_num();
                      #endif
                      mp_corrBasisFunctions[th] += ap_ImageS->m6p_backwardImage[0][0][fr][rg][cg][v]
                                                 * m2p_parametricImages[fb][v] 
                                                 * m2p_RGParametricImages[rb][v]
                                                 * m2p_CGParametricImages[cb][v];
                    }
                  }
            
            // Reduce
            for (int th=1 ; th<mp_ID->GetNbThreadsForImageComputation() ; th++)
              mp_corrBasisFunctions[0] += mp_corrBasisFunctions[th];
              
            // Apply corrections and normalization to the temporal basis functions values 
            if (mp_corrBasisFunctions[0] > 0.)
              m2p_cardBasisFunctions[cb][cg] *= mp_corrBasisFunctions[0]/parametric_image_norm;
            
          }
        }
      }
      
      m_basisFunctionsUpdIdx++;
  
      // Some feedback :
      if(m_verbose >=3)
      {
        for (int fb=0 ; fb<m_nbTimeBF ; fb++) 
          for (int fr=0 ; fr<mp_ID->GetNbTimeFrames() ; fr++)
            Cout( "iLinearModel::NestedEM() -> Basis function ["<<fb<<"] coefficients for frame ["<<fr<<"] :" << m2p_nestedModelTimeBasisFunctions[fb][fr] << endl);
            
        for (int rb=0 ; rb<m_nbRgateBF ; rb++)
          for (int rg=0 ; rg<mp_ID->GetNbRespGates() ; rg++)
            Cout( "iLinearModel::NestedEM() -> Basis function ["<<rb<<"] coefficients for resp gate ["<<rg<<"] :" << m2p_respBasisFunctions[rb][rg] << endl);
      
        for (int cb=0 ; cb<m_nbCgateBF ; cb++)
          for (int cg=0 ; cg<mp_ID->GetNbCardGates() ; cg++)
            Cout( "iLinearModel::NestedEM() -> Basis function ["<<cb<<"] coefficients for card gate ["<<cg<<"] :" << m2p_cardBasisFunctions[cb][cg] << endl);
      }
    } // end of if loop (basis functions update)



    // Step 2b : Image Coefficients estimation
    else 
    {
      // Normalization general factor
      FLTNB basis_functions_norm = 0.;
  
      // Regularisation according to frame basis functions
      if( m_nbModelParam>1 )
      {
        if(m_verbose >=3) Cout("iLinearModel::NestedEM() -> Estimate Coefficients - Frame coeffs estimation step" <<endl);
          
        for (int fb=0 ; fb<m_nbTimeBF ; fb++)
        {
          // Initialization of voxelwise corrections coefficients for time basis functions coefficients
          for (int v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
            mp_corrBasisCoeffs[v] = 0;
          
          basis_functions_norm = 0.;
          
          
          // Compute normalization related to the temporal basis functions.
          for (int fr=0 ; fr<mp_ID->GetNbTimeFrames() ; fr++)
            for (int rg=0 ; rg<mp_ID->GetNbRespGates() ; rg++)
              for (int cg=0 ; cg<mp_ID->GetNbCardGates() ; cg++)
                for (int rb=0 ; rb<m_nbRgateBF ; rb++)
                  for (int cb=0 ; cb<m_nbCgateBF ; cb++)
                    basis_functions_norm += m2p_nestedModelTimeBasisFunctions[fb][fr]
                                          * m2p_respBasisFunctions[rb][rg]
                                          * m2p_cardBasisFunctions[cb][cg];
          
          int v;
          #pragma omp parallel for private(v) schedule(static, 1)
          for (v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
          {
            // If a mask has been provided, check if the model applies in this voxel
            if(mp_maskModel != NULL && mp_maskModel[ v ] == 0)
              continue;

            for (int fr=0 ; fr<mp_ID->GetNbTimeFrames() ; fr++)
              for (int rg=0 ; rg<mp_ID->GetNbRespGates() ; rg++)
                for (int cg=0 ; cg<mp_ID->GetNbCardGates() ; cg++)
                  for (int rb=0 ; rb<m_nbRgateBF ; rb++)
                    for (int cb=0 ; cb<m_nbCgateBF ; cb++)
                    {
                      mp_corrBasisCoeffs[v] += ap_ImageS->m6p_backwardImage[0][0][fr][rg][cg][v] 
                                             * m2p_nestedModelTimeBasisFunctions[fb][fr]
                                             * m2p_respBasisFunctions[rb][rg]
                                             * m2p_cardBasisFunctions[cb][cg];
                    }
                    
            // Apply corrections and normalization to the frame model parametric image.
            if (mp_corrBasisCoeffs[v] >= 0.) // could be == 0 
              m2p_parametricImages[fb][v] *= mp_corrBasisCoeffs[v]/basis_functions_norm; 
          }
        }
      }
         
      // Regularisation according to respiratory basis functions
      
      if( m_nbRGModelParam>1 )
      {
        if(m_verbose >=3) Cout("iLinearModel::NestedEM() -> Estimate Coefficients - Respiratory Gate coeffs estimation step" <<endl);
          
        for (int rb=0 ; rb<m_nbRgateBF ; rb++)
        {
          // Initialization of voxelwise corrections for time basis functions coefficients
          for (int v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
            mp_corrBasisCoeffs[v] = 0;
          
          basis_functions_norm = 0.; 
          
          for (int fr=0 ; fr<mp_ID->GetNbTimeFrames() ; fr++)
            for (int rg=0 ; rg<mp_ID->GetNbRespGates() ; rg++) 
              for (int cg=0 ; cg<mp_ID->GetNbCardGates() ; cg++)
                for (int fb=0 ; fb<m_nbTimeBF ; fb++)
                  for (int cb=0 ; cb<m_nbCgateBF ; cb++)
                    basis_functions_norm += m2p_nestedModelTimeBasisFunctions[fb][fr]
                                          * m2p_respBasisFunctions[rb][rg]
                                          * m2p_cardBasisFunctions[cb][cg];
                
          // Compute corrections for the time vowelwise basis functions coefficients  
          int v;
          #pragma omp parallel for private(v) schedule(static, 1)
          for (v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++) 
          {
            // If a mask has been provided, check if the model applies in this voxel
            if(mp_maskModel != NULL && mp_maskModel[ v ] == 0)
              continue;
              
            for (int fr=0 ; fr<mp_ID->GetNbTimeFrames() ; fr++)
              for (int rg=0 ; rg<mp_ID->GetNbRespGates() ; rg++)
                for (int cg=0 ; cg<mp_ID->GetNbCardGates() ; cg++)
                  for (int fb=0 ; fb<m_nbTimeBF ; fb++)
                    for (int cb=0 ; cb<m_nbCgateBF ; cb++)
                    {
                      mp_corrBasisCoeffs[v] += ap_ImageS->m6p_backwardImage[0][0][fr][rg][cg][v] 
                                             * m2p_nestedModelTimeBasisFunctions[fb][fr]
                                             * m2p_respBasisFunctions[rb][rg]
                                             * m2p_cardBasisFunctions[cb][cg];
                    }
                    
            // Apply corrections and normalization to the respiratory model parametric image.
            if (mp_corrBasisCoeffs[v] >= 0.)
              m2p_RGParametricImages[rb][v] *= mp_corrBasisCoeffs[v]/basis_functions_norm;
          }
        }
      }
  
  
      // Regularisation according to cardiac basis functions
      if( m_nbCGModelParam>1 )
      {
        if(m_verbose >=3) Cout("iLinearModel::NestedEM() -> Estimate Coefficients - Cardiac Gate coeffs estimation step" <<endl);
          
        for (int cb=0 ; cb<m_nbCgateBF ; cb++)
        {
          // Initialization of voxelwise corrections for time basis functions coefficients
          for (int v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
            mp_corrBasisCoeffs[v] = 0;
          
          basis_functions_norm = 0.; 
          
          for (int fr=0 ; fr<mp_ID->GetNbTimeFrames() ; fr++)
            for (int rg=0 ; rg<mp_ID->GetNbRespGates() ; rg++) 
              for (int cg=0 ; cg<mp_ID->GetNbCardGates() ; cg++)
                for (int fb=0 ; fb<m_nbTimeBF ; fb++)
                  for (int rb=0 ; rb<m_nbRgateBF ; rb++)
                    basis_functions_norm += m2p_nestedModelTimeBasisFunctions[fb][fr]
                                          * m2p_respBasisFunctions[rb][rg]
                                          * m2p_cardBasisFunctions[cb][cg];

          // Compute corrections for the time vowelwise basis functions coefficients
          int v;
          #pragma omp parallel for private(v) schedule(static, 1)
          for (v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
          {
            // If a mask has been provided, check if the model applies in this voxel
            if(mp_maskModel != NULL && mp_maskModel[ v ] == 0)
              continue;
    
            for (int fr=0 ; fr<mp_ID->GetNbTimeFrames() ; fr++)
              for (int rg=0 ; rg<mp_ID->GetNbRespGates() ; rg++)
                for (int cg=0 ; cg<mp_ID->GetNbCardGates() ; cg++)
                  for (int fb=0 ; fb<m_nbTimeBF ; fb++)
                    for (int rb=0 ; rb<m_nbRgateBF ; rb++)
                    {     
                      mp_corrBasisCoeffs[v] += ap_ImageS->m6p_backwardImage[0][0][fr][rg][cg][v] 
                                             * m2p_nestedModelTimeBasisFunctions[fb][fr]
                                             * m2p_respBasisFunctions[rb][rg]
                                             * m2p_cardBasisFunctions[cb][cg];
                    }
                    
            // Apply corrections and normalization to the cardiac model parametric image.
            if (mp_corrBasisCoeffs[v] >= 0.)
              m2p_CGParametricImages[cb][v] *= mp_corrBasisCoeffs[v]/basis_functions_norm;
          }
        }
      }
      
      m_basisFunctionsUpdIdx++;
      
      // Some feedback :
      if(m_verbose >=3)
      {
        if( m_nbModelParam>1 )
          for (int fb=0 ; fb<m_nbTimeBF ; fb++) 
          { 
            FLTNB avg = 0;
            for (int v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
              if(mp_maskModel != NULL && mp_maskModel[ v ] == 1)
                avg += m2p_parametricImages[fb][v];
            
            avg /= m_nbVoxelsMask;
            Cout( "iLinearModel::NestedEM() -> Frame parametric image["<<fb<<"] avg value:" << avg << endl);
          }
  
        if( m_nbRGModelParam>1 )
          for (int rb=0 ; rb<m_nbRgateBF ; rb++) 
          { 
            FLTNB avg = 0;
            
            for (int v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
              if(mp_maskModel != NULL && mp_maskModel[ v ] == 1)
                avg += m2p_RGParametricImages[rb][v];
            
            avg /= m_nbVoxelsMask;
            Cout( "iLinearModel::NestedEM() -> Resp gate parametric image["<<rb<<"] avg value:" << avg << endl);
          }
        
        if( m_nbCGModelParam>1 )
          for (int cb=0 ; cb<m_nbCgateBF ; cb++) 
          { 
            FLTNB avg = 0;
            
            for (int v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
              if(mp_maskModel != NULL && mp_maskModel[ v ] == 1)
                avg += m2p_CGParametricImages[cb][v];
            
            avg /= m_nbVoxelsMask;
            Cout( "iLinearModel::NestedEM() -> Card gate parametric image["<<cb<<"] avg value:" << avg << endl);
          }
      }
      
    } // end of else loop (parametric images update)
    
  } // end of loop on cycles (linear model iterations)

  return 0;
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn EstimateParametersWithNNLS
  \param ap_ImageS : pointer to the ImageSpace
  \param a_ite : index of the actual iteration (not used)
  \brief Estimate parametric images using the NNLS method
  \return 0 if success, other value otherwise.
*/
int iLinearModel::EstimateParametersWithNNLS(oImageSpace* ap_ImageS, int a_ite) 
{
  if(m_verbose >=3) Cout("iLinearModel::EstimateParametersWithNNLS() ..." <<endl);

  // Regularisation according to frame basis functions
  if( m_nbModelParam>1 )
  {
    if(m_verbose >=3) Cout("iLinearModel::EstimateParametersWithNNLS() -> Estimate Coefficients - Frame coeffs estimation step" <<endl);
          
    for (int rg=0 ; rg<mp_ID->GetNbRespGates() ; rg++)
      for (int cg=0 ; cg<mp_ID->GetNbCardGates() ; cg++)
      {
        int v;
        #pragma omp parallel for private(v) schedule(guided)
        for (v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
        {
          // If a mask has been provided, check if the model applies in this voxel
          if(mp_maskModel != NULL && mp_maskModel[ v ] == 0)
            continue;
          
          int th=0;
          #ifdef CASTOR_OMP
          th = omp_get_thread_num();
          #endif
          
          FLTNB** pp_nnls_A = m3p_nnlsA[ th ];
          FLTNB*  p_nnls_B = m2p_nnlsB[ th ];
          FLTNB   nnls_rnorm;
          FLTNB  *p_nnls_zz=NULL;
          
          // Fill NNLS matrices
          for(int t=0; t<mp_ID->GetNbTimeFrames(); t++)
          {
            // Fill NNLS A matrix:
            for (int fb=0 ; fb<m_nbTimeBF ; fb++)
              pp_nnls_A[ fb ][ t ] = m2p_nestedModelTimeBasisFunctions[ fb ][ t ] * mp_w[ t ];
              
            // Fill NNLS B array: tissue
            p_nnls_B[ t ]=ap_ImageS->m4p_image[ t ][ rg ][ cg ][ v ] * mp_w[ t ];
          }
      
          // Launch NNLS
          if(NNLS(pp_nnls_A,
                  mp_ID->GetNbTimeFrames(),
                  m_nbTimeBF,
                  p_nnls_B,
                  m2p_nnlsX[ th ],
                  &nnls_rnorm,
                  m2p_nnlsWp[ th ],
                  p_nnls_zz,
                  m2p_nnlsIdx[ th ]) )
            continue;
          
          // Recover parameters
          for (int fb=0 ; fb<m_nbTimeBF ; fb++)
            m2p_parametricImages[fb][v] = m2p_nnlsX[ th ][ fb ];
            
        } // end of voxel loop
        
      }
  }
     
  // Regularisation according to respiratory basis functions
  
  if( m_nbRGModelParam>1 )
  {
    if(m_verbose >=3) Cout("iLinearModel::EstimateParametersWithNNLS() -> Estimate Coefficients - Respiratory Gate coeffs estimation step" <<endl);    
    
    for (int fr=0 ; fr<mp_ID->GetNbTimeFrames() ; fr++)
      for (int cg=0 ; cg<mp_ID->GetNbCardGates() ; cg++)
      {
        int v;
        #pragma omp parallel for private(v) schedule(guided)
        for (v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
        {
          // If a mask has been provided, check if the model applies in this voxel
          if(mp_maskModel != NULL && mp_maskModel[ v ] == 0)
            continue;
          
          int th=0;
          #ifdef CASTOR_OMP
          th = omp_get_thread_num();
          #endif
          
          FLTNB** pp_nnls_A = m3p_nnlsA[ th ];
          FLTNB*  p_nnls_B = m2p_nnlsB[ th ];
          FLTNB   nnls_rnorm;
          FLTNB  *p_nnls_zz=NULL;
          
          // Fill NNLS matrices
          for (int rg=0 ; rg<mp_ID->GetNbRespGates() ; rg++)
          {
            // Fill NNLS A matrix:
            for (int rb=0 ; rb<m_nbRgateBF ; rb++)
              pp_nnls_A[ rb ][ rg ] = m2p_respBasisFunctions[ rb ][ rg ];
              
            // Fill NNLS B array: tissue
            p_nnls_B[ rg ]=ap_ImageS->m4p_image[ fr ][ rg ][ cg ][ v ];
          }
      
          // Launch NNLS
          if(NNLS(pp_nnls_A,
                  mp_ID->GetNbRespGates(),
                  m_nbRgateBF,
                  p_nnls_B,
                  m2p_nnlsX[ th ],
                  &nnls_rnorm,
                  m2p_nnlsWp[ th ],
                  p_nnls_zz,
                  m2p_nnlsIdx[ th ]) )
            continue;
          
          // Recover parameters
          for (int rb=0 ; rb<m_nbRgateBF ; rb++)
            m2p_RGParametricImages[rb][v] = m2p_nnlsX[ th ][ rb ];
          
        }  // end of loop on voxels
        
      }
  }


  // Regularisation according to cardiac basis functions
  if( m_nbCGModelParam>1 )
  {
    if(m_verbose >=3) Cout("iLinearModel::EstimateParametersWithNNLS() -> Estimate Coefficients - Cardiac Gate coeffs estimation step" <<endl);
    
    for (int fr=0 ; fr<mp_ID->GetNbTimeFrames() ; fr++)
      for (int rg=0 ; rg<mp_ID->GetNbRespGates() ; rg++)
      {
        int v;
        #pragma omp parallel for private(v) schedule(guided)
        for (v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
        {
          // If a mask has been provided, check if the model applies in this voxel
          if(mp_maskModel != NULL && mp_maskModel[ v ] == 0)
            continue;
          
          int th=0;
          #ifdef CASTOR_OMP
          th = omp_get_thread_num();
          #endif
          
          FLTNB** pp_nnls_A = m3p_nnlsA[ th ];
          FLTNB*  p_nnls_B = m2p_nnlsB[ th ];
          FLTNB   nnls_rnorm;
          FLTNB  *p_nnls_zz=NULL;
          
          // Fill NNLS matrices
          for (int cg=0 ; cg<mp_ID->GetNbCardGates() ; cg++)
          {
            // Fill NNLS A matrix:
            for (int cb=0 ; cb<m_nbCgateBF ; cb++)
              pp_nnls_A[ cb ][ cg ] = m2p_cardBasisFunctions[ cb ][ cg ];
              
            // Fill NNLS B array: tissue
            p_nnls_B[ cg ]=ap_ImageS->m4p_image[ fr ][ rg ][ cg ][ v ];
          }
      
          // Launch NNLS
          if(NNLS(pp_nnls_A,
                  mp_ID->GetNbCardGates(),
                  m_nbCgateBF,
                  p_nnls_B,
                  m2p_nnlsX[ th ],
                  &nnls_rnorm,
                  m2p_nnlsWp[ th ],
                  p_nnls_zz,
                  m2p_nnlsIdx[ th ]) )
            continue;
          
          // Recover parameters
          for (int cb=0 ; cb<m_nbCgateBF ; cb++)
            m2p_CGParametricImages[cb][v] = m2p_nnlsX[ th ][ cb ];
            
        }  // end of loop on voxels
        
      }
  }
  
  
  // Some feedback :
  if(m_verbose >=3)
  {
    if( m_nbModelParam>1 )
      for (int fb=0 ; fb<m_nbTimeBF ; fb++) 
      { 
        FLTNB avg = 0;
        for (int v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
          if(mp_maskModel != NULL && mp_maskModel[ v ] == 1)
            avg += m2p_parametricImages[fb][v];
        
        avg /= m_nbVoxelsMask;
        Cout( "iLinearModel::EstimateParametersWithNNLS() -> Frame parametric image["<<fb<<"] avg value:" << avg << endl);
      }

    if( m_nbRGModelParam>1 )
      for (int rb=0 ; rb<m_nbRgateBF ; rb++) 
      { 
        FLTNB avg = 0;
        
        for (int v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
          if(mp_maskModel != NULL && mp_maskModel[ v ] == 1)
            avg += m2p_RGParametricImages[rb][v];
        
        avg /= m_nbVoxelsMask;
        Cout( "iLinearModel::EstimateParametersWithNNLS() -> Resp gate parametric image["<<rb<<"] avg value:" << avg << endl);
      }
    
    if( m_nbCGModelParam>1 )
      for (int cb=0 ; cb<m_nbCgateBF ; cb++) 
      { 
        FLTNB avg = 0;
        
        for (int v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
          if(mp_maskModel != NULL && mp_maskModel[ v ] == 1)
            avg += m2p_CGParametricImages[cb][v];
        
        avg /= m_nbVoxelsMask;
        Cout( "iLinearModel::EstimateParametersWithNNLS() -> Card gate parametric image["<<cb<<"] avg value:" << avg << endl);
      }
  }

  return 0;
}


// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn EstimateParametersWithLS (linear regresion)
  \param ap_ImageS : pointer to the ImageSpace
  \param a_ite : index of the actual iteration (not used)
  \brief Estimate parametric images using LS linear regression
  \return 0 if success, other value otherwise.
*/

int iLinearModel::Patlak_LS(oImageSpace* ap_ImageS, int a_ite)
{
  if(m_verbose >=3) Cout("iLinearPatlakModel::Patlak_LS ..." <<endl);

  int v;
  #pragma omp parallel for private(v) schedule(guided)
  for (v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
  {
    // If a mask has been provided, check if the model applies in this voxel
    if(mp_maskModel != NULL && mp_maskModel[ v ] == 0)
      continue;

    FLTNB xmean = 0.,
            ymean = 0.,
            K = 0.,
            Vd = 0.,
            Cov= 0.,
            Var = 0.;


    // Compute means for this voxel
    for (int fr=0 ; fr<mp_ID->GetNbTimeFrames() ; fr++)
    {
      xmean += m2p_nestedModelTimeBasisFunctions[0][fr]/m2p_nestedModelTimeBasisFunctions[1][fr];
      ymean += ap_ImageS->m4p_image[fr][0][0][v]/m2p_nestedModelTimeBasisFunctions[1][fr];
    }

    xmean /= mp_ID->GetNbTimeFrames();
    ymean /= mp_ID->GetNbTimeFrames();

    FLTNB Yt,Xt ;

    // Compute covariance & variance
    for (int fr=0 ; fr<mp_ID->GetNbTimeFrames() ; fr++)
    {
      Yt = ap_ImageS->m4p_image[fr][0][0][v]/m2p_nestedModelTimeBasisFunctions[1][fr];
      Xt = m2p_nestedModelTimeBasisFunctions[0][fr]/m2p_nestedModelTimeBasisFunctions[1][fr];
      Cov += (Xt-xmean) * (Yt-ymean);
      Var += (Xt-xmean) * (Xt-xmean);
    }

    // Slope
    K = (Var != 0) ? Cov/Var : 0.;

    // Non-negativity constraint
    if(K<0) K = 0.;

    // Intercept
    Vd = ymean - K*xmean;

    // Non-negativity constraint
    if(Vd<0) Vd = 0.;

    m2p_parametricImages[0][v] = K;
    m2p_parametricImages[1][v] = Vd;
  }

  return 0;
}


// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn EstimateImageWithModel
  \param ap_ImageS : pointer to the ImageSpace
  \param a_ite : index of the actual iteration (not used)
  \param a_sset : index of the actual subset (not used)
  \brief Re-estimate image using the linear parametric images and basis functions
  \return 0 if success, other value otherwise.
*/
int iLinearModel::EstimateImageWithModel(oImageSpace* ap_ImageS, int a_ite, int a_sset) 
{
  if(m_verbose >= 3) Cout("iLinearModel::EstimateImageWithModel ... " <<endl);

  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** iLinearModel::EstimateImageWithModel() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif

  int v;
  #pragma omp parallel for private(v) schedule(static, 1)
  for (int fr=0 ; fr<mp_ID->GetNbTimeFrames() ; fr++)
    for (int rg=0 ; rg<mp_ID->GetNbRespGates() ; rg++)
      for (int cg=0 ; cg<mp_ID->GetNbCardGates() ; cg++)
        for (v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
        {
          // If a mask has been provided, check if the model applies in this voxel
          if(mp_maskModel != NULL && mp_maskModel[ v ] == 0)
            continue;
    
          // Reset current estimated value
          ap_ImageS->m4p_image[fr][rg][cg][v] = 0.;

          for (int fb=0 ; fb<m_nbTimeBF ; fb++)
            for (int rb=0 ; rb<m_nbRgateBF ; rb++)
              for (int cb=0 ; cb<m_nbCgateBF ; cb++)
              {
                // Retrieve current estimation of image according to coeffs/basis functions
                ap_ImageS->m4p_image[fr][rg][cg][v] += m2p_parametricImages[fb][v] 
                                                     * m2p_nestedModelTimeBasisFunctions[fb][fr]
                                                     * m2p_RGParametricImages[rb][v] 
                                                     * m2p_respBasisFunctions[rb][rg]
                                                     * m2p_CGParametricImages[cb][v] 
                                                     * m2p_cardBasisFunctions[cb][cg] ;
              }
        }
  
  return 0;
}

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

/*
  Implementation of class i1TCModel

  - separators: X
  - doxygen: X
  - default initialization: none  (require user inputs (basis functions), no sense to provide 'standard' configuration file as for optimizer/projector)
  - CASTOR_DEBUG: X
  - CASTOR_VERBOSE: X
*/

/*!
  \file
  \ingroup  dynamic
  \brief    Implementation of class i1TCModel
*/


#include "i1TCModel.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn i1TCModel
  \brief Constructor of i1TCModel. Simply set all data members to default values.
*/
i1TCModel::i1TCModel() : vDynamicModel() 
{
  m_nbTimeBF = 2; // Two basis functions (Plasma curve (Cp) and integral of plasma curve (Cpi))
  m_nbModelParam = 3; // K1, k2, Va
  m2p_parametricImages = NULL;
  m2p_nestedModelTimeBasisFunctions = NULL;

  mp_parUpperBounds = NULL;
  mp_parLowerBounds = NULL;
  m_RRcst = -1.;
  m_ridgeRegressionFlag = false;
  
  m2p_ct = NULL;
  m2p_cti = NULL;
  mp_w = NULL;
  m3p_nnlsA = NULL;
  m2p_nnlsB = NULL;
  m2p_nnlsMat = NULL;
  mp_VaImage = NULL;
  
  // Matrices for LS optimization method
  mp_Y = NULL;
  mp_X = NULL;
  mp_Xt = NULL;
  mp_XtX = NULL;
  mp_LSnum = NULL;
  mp_LSden = NULL;
  mp_Theta = NULL;
  // RRLS
  mp_RRm = NULL;
  mp_RRw = NULL;
  mp_RRnum = NULL;
    
  m_fileOptions = "";
  m_listOptions = "";
  m_OptimisationMethodFlag = 0;
    
  mp_DT2 = NULL;
  
  m_intMethodFlag = METHOD_INT_WPO;
  mp_wpoQ = NULL;
  mp_wpoA = NULL;
  m2p_wpoP = NULL;
  m2p_wpoFD = NULL;
  m2p_wpoBD = NULL;

}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ~i1TCModel
  \brief Destructor of i1TCModel
*/
i1TCModel::~i1TCModel() 
{
  if(m_initialized)
  {
    if(m2p_ct && m2p_cti)
    {
      for( int th=0 ; th<mp_ID->GetNbThreadsForImageComputation() ; th++ )
      {
        if(m2p_ct[ th ])  delete[] m2p_ct[ th ];
        if(m2p_cti[ th ]) delete[] m2p_cti[ th ];
      }
      
      delete[] m2p_ct;
      delete[] m2p_cti;
    }

    // NNLS variables
    if(m3p_nnlsA && m2p_nnlsB && m2p_nnlsMat)
    {      
      for( int th=0 ; th<mp_ID->GetNbThreadsForImageComputation() ; th++ )
      {
        //for( int n=0 ; n<m_nnlsN ; n++ )
        //  if(m3p_nnlsA[ th ][ n ]) delete[] m3p_nnlsA[ th ][ n ];

        if(m3p_nnlsA[ th ] )   delete[] m3p_nnlsA[ th ];
        if(m2p_nnlsB[ th ] )   delete[] m2p_nnlsB[ th ];
        if(m2p_nnlsMat[ th ] ) delete[] m2p_nnlsMat[ th ];
      }

      delete[] m3p_nnlsA;
      delete[] m2p_nnlsB;
      delete[] m2p_nnlsMat;
    }
    
    if( mp_w ) delete[] mp_w;
       
    // Matrices for LS method
    if( m_OptimisationMethodFlag== METHOD_1CPT_LS)
    {  
      for( int th=0 ; th<mp_ID->GetNbThreadsForImageComputation() ; th++ )
      {
        if( mp_Y[ th ] )      delete mp_Y[ th ]      ;
        if( mp_X[ th ] )      delete mp_X[ th ]      ;
        if( mp_Xt[ th ] )     delete mp_Xt[ th ]    ;
        if( mp_XtX[ th ] )    delete mp_XtX[ th ]    ;
        if( mp_LSnum[ th ] )  delete mp_LSnum[ th ]  ;
        if( mp_LSden[ th ] )  delete mp_LSden[ th ]  ;
        if( mp_Theta[ th ] )  delete mp_Theta[ th ]  ;
        if( mp_RRnum[ th ] ) delete mp_RRnum[ th ] ;
      }
      
      if( mp_Y ) delete[] mp_Y ;
      if( mp_X ) delete[] mp_X ;
      if( mp_Xt ) delete[] mp_Xt;
      if( mp_XtX ) delete[] mp_XtX;
      if( mp_LSnum ) delete[] mp_LSnum;
      if( mp_LSden ) delete[] mp_LSden;
      if( mp_Theta ) delete[] mp_Theta;
      if( mp_RRnum ) delete[] mp_RRnum;
      if( mp_RRm ) delete[] mp_RRm;
      if( mp_RRw ) delete[] mp_RRw;
    }

    if( mp_parUpperBounds ) delete[]mp_parUpperBounds;
    if( mp_parLowerBounds ) delete[]mp_parLowerBounds;
    
    if( mp_VaImage ) delete[] mp_VaImage;
    if( mp_blackListedvoxelsImage ) delete[] mp_blackListedvoxelsImage;

    if( mp_DT2 ) delete[] mp_DT2;
    if(mp_wpoQ) delete[] mp_wpoQ;
    if(mp_wpoA) delete[] mp_wpoA;
    if(m2p_wpoP) delete[] m2p_wpoP;
    if(m2p_wpoFD) delete[] m2p_wpoFD;
    if(m2p_wpoBD) delete[] m2p_wpoBD;
  }
  
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ShowHelpModelSpecific
  \brief Print out specific help about the implementation of the model
   and its initialization
*/
void i1TCModel::ShowHelpModelSpecific()
{
  cout << "-- This class implements a 2 compartments kinetic model " << endl;
  cout << "-- or 1 Tissue Compartment Model (1TCM) for perfusion quantitation " << endl;
  cout << "-- for radiotracers such as [15-O]H2O. " << endl;
  cout << endl;
  cout << "-- *--------*   K1    *--------*  " << endl;
  cout << "-- |        | ----->  |        |  " << endl;
  cout << "-- |   Cp   |         |   Ct   |  " << endl;
  cout << "-- |        | <-----  |        |  " << endl;
  cout << "-- *--------*   k2    *--------*  " << endl;
  cout << endl;
  cout << "-- This model contains 3 parameters, including 2 rate constants, " << endl;
  cout << "-- K1 (v/s/v, where v is a volume unit), k2 (s-1) and the arterial volume fraction in tissue (Va), " << endl;
  cout << "-- as described by the following equations : " << endl;
  cout << "-- (1) Cpet ( t )    = Ct( t )    + Va*Cp( t ) " << endl;
  cout << "-- (2) dCt( t ) / dt = K1*Cp( t ) - k2*Ct( t ) " << endl;
  cout << "-- where Cpet( t ) is the image value in voxel for frame t, and" << endl;
  cout << "-- Ct and Cp are activity concentration in tissue and plasma," << endl;
  cout << "-- The model estimates K1, k2 and Va parametric images." << endl;
  cout << "-- The input function (Cp) must be provided by the user." << endl;
  cout << endl;
  cout << "--------------------------------------" << endl;
  cout << "   The model can be initialized using either a configuration text file, or a list of options :" << endl; 
  cout << endl;
  cout << " - CONFIGURATION FILE : (command-line options : _1TCM:path/to/conf/file.txt:)" << endl;
  cout << "   " << endl;
  cout << "   'Input_function:'      (mandatory) " << endl;
  cout << "                            Enter the activity concentration values (kBq/v/s)  of the plasma function (cp) for each time frame (tf)," << endl;
  cout << "                            separated by ','. Time unit is seconds : " << endl;
  cout << "                            -> Input functions: " << endl;
  cout << "                            -> val_cp_tf1 , val_cp_tf2 , ..., val_cp_tfn" << endl;
  cout << endl;
  cout << "   'Integral_input_function:' (optional) " << endl; 
  cout << "                            Enter the activity concentration values of the integral of the plasma function (Icp) for each time frame (tf)," << endl;
  cout << "                            separated by ','. Time unit is seconds :" << endl;
  cout << "                            -> Integral_input_function: " << endl;
  cout << "                            -> val_Icp_tf1, val_Icp_tf2, ..., val_Icp_tfn" << endl;
  cout << "                            NOTE: If not provided, the integral will be estimated from the input function" << endl;
  cout << "                                  Set the 'Integral method' flag below to select the method for TAC integration (default: WPO))" << endl;
  cout << endl;
  cout << "   'Optimisation_method:'  (optional) Define the method to use for parameters estimation. Only least-squares methods " << endl;
  cout << "    (default: NNLS)                   are currently available to estimate O = (K1, k2, Va)" << endl;
  cout << "                                           Ô = [ X'WX ]-1 X'Wy, with " << endl;
  cout << "                                           with  y = data vector                      " << endl; 
  cout << "                                                 X = model matrix (Cp, Icp, Cpet)     " << endl; 
  cout << "                                                 W = weights vector (frame duration)  " << endl; 
  cout << "                                      If the estimated parameter values are negative, the activity value for the related voxels will be kept to their original number" << endl;
  cout << "                                      0 (default) = Non-negative least-square (NNLS)" << endl;
  cout << "                                                    Derived from Turku PET center libraries, authors: Vesa Oikonen and Kaisa Sederholm" << endl;
  cout << "                                                    (http://www.turkupetcentre.net/petanalysis/index.html)" << endl;
  cout << "                                                    Routine based on the text and fortran code in C.L. Lawson and R.J. Hanson," << endl;
  cout << "                                                    Solving Least Squares Problems, Prentice-Hall, Englewood Cliffs, New Jersey, 1974." << endl;
  cout << "                                                    Note: K1 estimated values could still be negative as they are computed from a substraction of the estimated NNLS parameters." << endl;
  cout << "                                                          Activity values will be kept to their original values for the voxels involved" << endl;
  cout << "                                      1           = Standard Least Square (LS) optimization  " << endl; 
  cout << endl;
  cout << "   'Integration_method:'   (optional) Define the method to use for TAC integration over the time samples" << endl;
  cout << "                                      0 (default) = Weighed parabola overlapping (WPO) (Z.Wang, D.Feng, Int. Sys. Sci 23 (1992), pp.1361-69)" << endl;
  cout << "                                      1           = Trapezoidal " << endl;
  cout << "   'Ridge_parameter:'      (optional) Constant for Ridge Regression during Least-Square optimization (only available with Least-Square algorithm and not NNLS )   " << endl;
  cout << "    (default: 0)                      Bounds must be provided with the eponym options below in order to compute ridge weights and means for the new cost function:  " << endl;
  cout << "                                      Ô = [ X'*W*X + t.Rw ]-1 [ X'*W*y  ]        " << endl;
  cout << "                                        + [ X'*W*X + t.Rw ]-1 [ t.Rw*Rm ]    " << endl;
  cout << "                                      with  y = data vector                      " << endl; 
  cout << "                                            X = model matrix (Cp, Icp, Cpet)     " << endl; 
  cout << "                                            W = weights vector (frame duration)  " << endl; 
  cout << "                                            t = Ridge constant                   " << endl; 
  cout << "                                            Rw= Ridge weights                    " << endl; 
  cout << "                                            Rm= Ridge means                      " << endl;   
  cout << "   'Bounds:'               (optional) Upper / Lower Bounds for each 3 parameters, to define ridge means mr0, and weights wr0," << endl;
  cout << "                                      such as mr0 = (Max+Min)/2 and wr0 = 1 / (Max-Min)^2 " << endl;
  cout << "                                      They must be provided as in the following syntax: " << endl;
  cout << "                                      Bounds: K1Max, K1Min, k2Max, k2Min, VaMax, VaMin" << endl;
  cout << "                                      Default: K1(Max,Min) = 0.1, 0. " << endl;
  cout << "                                               K2(Max,Min) = 0.1, 0. " << endl;
  cout << "                                               Va(Max,Min) = 1. , 0. " << endl;
  cout << "   'VA_image:'             (optional) Path to an interfile image containing the arterial volume fraction value in tissue for each voxel," << endl;
  cout << "                                      only K1 and k2 rate constants will be estimated (Default: All parameters are estimated) " << endl;
  cout << endl;
  cout << " - LINE OPTIONS : (command-line options : _1TCM,options:)" << endl;
  cout << "   " << endl;
  cout << "   Command-line options just require the samples of the plasma input curve, separated by commas. " << endl;
  cout << "   All other options will be set to default. The optimization algorithm will be NNLS, the integration method for TAC will be WPO " << endl;
  cout << "   -> 1cptModel,val_cp_tf1 , val_cp_tf2 , ..., val_cp_tfn" << endl;
  cout << endl;

  // Print general help for all dynamic models
  ShowHelp();
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ReadAndCheckConfigurationFileSpecific
  \brief This function is used to read options from a configuration file.
  \return 0 if success, other value otherwise.
*/
int i1TCModel::ReadAndCheckConfigurationFileSpecific()
{
  if(m_verbose >=3) Cout("i1TCModel::ReadAndCheckConfigurationFileSpecific ..."<< endl); 
  
  // The file will be processed in the Initialize() function
  
  ifstream in_file(m_fileOptions.c_str(), ios::in);
  
  if(in_file)
  {
    // Method to be used to compute integral of TACs
    // (Check this first in case we must estimate plasma integral
    // so we need to know which method to use)
    // Default : WPO
    
    string dtag = "Integration_method";
    if( ReadDataASCIIFile(m_fileOptions,
                          dtag,
                          &m_intMethodFlag,
                          1,
                          KEYWORD_OPTIONAL) == 1)
    {
      Cerr("***** i1TCModel::ReadAndCheckConfigurationFileSpecific() -> Error while trying to read '"<< dtag <<"' flag in "<< m_fileOptions << endl);
      return 1;
    }
    
    // Optimisation method for the estimation of parameters
    dtag = "Optimisation_method";
    if( ReadDataASCIIFile(m_fileOptions,
                          dtag,
                          &m_OptimisationMethodFlag,
                          1,
                          KEYWORD_OPTIONAL) == 1)
    {
      Cerr("***** i1TCModel::ReadAndCheckConfigurationFileSpecific() -> Error while trying to read '"<< dtag <<"' flag in "<< m_fileOptions << endl);
      return 1;
    }
    
    // Save blacklisted voxels images on disk ?
    dtag = "Save_blacklisted_voxels_images";
    if( ReadDataASCIIFile(m_fileOptions,
                          dtag,
                          &m_saveBlacklistedImageMaskFlag,
                          1,
                          KEYWORD_OPTIONAL) == 1)
    {
      Cerr("***** vDynamicModel::ReadAndCheckConfigurationFile() -> Error while trying to read '"<< dtag <<"' flag in " << m_fileOptions << endl);
      return 1;
    }
      
    // --- Constraints for ridge regression ---
    const int nb_params = 3;
    // Default upper/lower bounds
    HPFLTNB p_bounds_params[ 6 ] = {0.1, 0., 0.1, 0., 1., 0.};
    
    // Ridge Regression parameter
    dtag = "Ridge_Parameter";
    if( ReadDataASCIIFile(m_fileOptions,
                          dtag,
                          &m_RRcst,
                          1,
                          KEYWORD_OPTIONAL) == 1)
    {
      Cerr("***** i1TCModel::ReadAndCheckConfigurationFileSpecific() -> Error while trying to read '"<< dtag <<"' flag in "<< m_fileOptions << endl);
      return 1;
    }

    // If a ridge regression parameter constant has been provided,
    // set up the related flag and increase the number of LS fonction models
    if (m_RRcst>=0) m_ridgeRegressionFlag = true;
    
    dtag = "Bounds";
    if( ReadDataASCIIFile(m_fileOptions,
                          dtag,
                          p_bounds_params,
                          2*nb_params,
                          KEYWORD_OPTIONAL) == 1)
    {
      Cerr("***** i1TCModel::ReadAndCheckConfigurationFileSpecific() -> Error while trying to read '"<< dtag <<"' flag (upper/lower bounds) in "<< m_fileOptions << endl);
      return 1;
    }

    // init parameters bounds
    mp_parUpperBounds = new FLTNB[ nb_params ];
    mp_parLowerBounds = new FLTNB[ nb_params ];
    
    for (int p=0 ; p<nb_params ; p++ )
    {
      mp_parUpperBounds[ p ] = p_bounds_params[ p*2 ];
      mp_parLowerBounds[ p ] = p_bounds_params[ p*2 + 1 ];
    }
  }
  else
  {
    Cerr("***** i1TCModel::ReadAndCheckConfigurationFileSpecific() -> Error while trying to read configuration file at: " << m_fileOptions << endl);
    return 1;
  }
    
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ReadAndCheckOptionsList
  \param const string& a_optionsList : a list of parameters separated by commas
  \brief This function is used to read parameters from a string.
  \return 0 if success, other value otherwise.
*/
int i1TCModel::ReadAndCheckOptionsList(string a_listOptions)
{
  if(m_verbose >=2) Cout("i1TCModel::ReadAndCheckOptionsList ..."<< endl); 
  
  // Just get the options here and perform initialization in InitializeSpecific()
  // function, as some memory allocations must be performed first (images & TACs)
  m_listOptions = a_listOptions;
  
  // Normal end
  return 0;
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
int i1TCModel::CheckSpecificParameters()
{
  if(m_verbose >=2) Cout("i1TCModel::CheckSpecificParameters ..."<< endl); 
  
  // Check image dimensions
  if (mp_ID==NULL)
  {
    Cerr("***** i1TCModel::CheckSpecificParameters() -> ImageDimensions object has not been provided !" << endl);
    return 1;
  }
  
  // Check number time basis functions
  if (m_nbTimeBF<0)
  {
    Cerr("***** i1TCModel::CheckSpecificParameters() -> Wrong number of time frame basis functions !" << endl);
    return 1;
  }
  
  // Check if we have somehow both a file and a list of options for init...
  if(m_listOptions != "" && m_fileOptions != "")
  {
    Cerr("***** i1TCModel::CheckSpecificParameters() -> Either a file or a list of options have to be selected to initialize the model, but not both ! " << endl);
    return 1;
  }
  
  // Check if we have no file not list of options for some reason...
  if(m_listOptions == "" && m_fileOptions == "")
  {
    Cerr("***** i1TCModel::CheckSpecificParameters -> Either a file or a list of options should have been provided at this point ! " << endl);
    return 1;
  }
  
  // Check if we reconstruct gated data. Throw warning if it is the case
  if(mp_ID->GetNbRespGates()>1 || mp_ID->GetNbCardGates()>1)
  {
    Cerr("***** i1TCModel::CheckSpecificParameters -> The implemented model is not compatible yet with gated reconstruction (parametric images will be the same for each gate)! " << endl);
    return 1;
  }
  
  // Check if ridge regression and NNLS has been enabled (RR only compatible with LS)
  if ( m_OptimisationMethodFlag == METHOD_1CPT_NNLS
    && m_ridgeRegressionFlag )
  {
    Cerr("***** i1TCModel::CheckSpecificParameters() -> Error, ridge regression is not compatible with NNLS algorithm. Switch 'Optimisation method' flag to LS in order to use ridge-regression !" << endl);
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
  \brief This function is used to initialize the model parametric images and basis functions
  \return 0 if success, other value otherwise.
*/
int i1TCModel::InitializeSpecific()
{
  if(m_verbose >=2) Cout("i1TCModel::InitializeSpecific() ..."<< endl); 

  // Forbid initialization without check
  if (!m_checked)
  {
    Cerr("***** oDynamicModelManager::InitializeSpecific() -> Must call CheckParameters functions before Initialize() !" << endl);
    return 1;
  }
  
  // --- Memory Allocation --- //
  
  // Allocate memory for Parametric images and functions
  m2p_nestedModelTimeBasisFunctions = new FLTNB*[m_nbTimeBF];
  m2p_parametricImages = new FLTNB*[m_nbModelParam];

  for(int b=0 ; b<m_nbTimeBF ; b++)
    m2p_nestedModelTimeBasisFunctions[b] = new FLTNB[mp_ID->GetNbTimeFrames()];

  // Init with negative values
  for(int b=0 ; b<m_nbTimeBF ; b++)
    for(int t=0 ; t<mp_ID->GetNbTimeFrames() ; t++)
      m2p_nestedModelTimeBasisFunctions[b][t] = -1.;
    
  for(int p=0 ; p<m_nbModelParam ; p++)
    m2p_parametricImages[p] = new FLTNB[mp_ID->GetNbVoxXYZ()];


  // Initialize blacklisted voxels image to zero;
  mp_blackListedvoxelsImage = new FLTNB[mp_ID->GetNbVoxXYZ()];

  for (int v = 0; v < mp_ID->GetNbVoxXYZ(); v++)
    mp_blackListedvoxelsImage[v] = 0.0;
  
  
  // --- Data Initialization with a configuration file --- //
  // Data requiring memory allocation (Image and TACs initializations)
  
  if(m_fileOptions != "")
  {
    ifstream in_file(m_fileOptions.c_str(), ios::in);
    
    if(in_file)
    {      
      // --- Plasma input functions --- 
      
      // Recover Plasma input function TAC (mandatory)
      string dtag = "Input_function";
      if( ReadDataASCIIFile(m_fileOptions,
                            dtag,
                            m2p_nestedModelTimeBasisFunctions[ 1 ],
                            mp_ID->GetNbTimeFrames(),
                            KEYWORD_MANDATORY) )
      {
        Cerr("***** i1TCModel::InitializeSpecific() -> Error while trying to read 1TCPM input functions !" << endl);
        Cerr("                                         '"<< dtag <<"' keyword in " << m_fileOptions << endl);
        return 1;
      }

      // Recover Integral of plasma input function TAC (optional)
      dtag = "Integral_input_function";
      if( ReadDataASCIIFile(m_fileOptions,
                            dtag,
                            m2p_nestedModelTimeBasisFunctions[ 0 ],
                            mp_ID->GetNbTimeFrames(),
                            KEYWORD_OPTIONAL) == 1 )
      {
        Cerr("***** i1TCModel::InitializeSpecific() -> Error while trying to read 1TCPM input functions !" << endl);
        Cerr("                                         '"<< dtag <<"' keyword in " << m_fileOptions << endl);
        return 1;
      }
      
      
      // Parametric images initialization
      string input_image = "";
      int return_value = 0;
      
      dtag = "Parametric_images_init";
      return_value = ReadDataASCIIFile(m_fileOptions,
                                       dtag,
                                       &input_image,
                                       1,
                                       KEYWORD_OPTIONAL);
      
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
          Cerr("***** i1TCModel::InitializeSpecific() -> Error while trying to read the provided initialization parametric images : " << input_image << endl);
          return 1;
        }
      }
      else if( return_value == 1) // Error during reading
      {
        Cerr("***** i1TCModel::InitializeSpecific() -> Error while trying to read input functions coefficients !" << endl);
        Cerr("                                          '"<< dtag <<"' keyword in " << m_fileOptions << endl);
        return 1;
      }
      else //(return_value >= 1 ) // Keyword not found : no initialization provided
      {
        // Standard initialization
        for(int p=0 ; p<m_nbModelParam ; p++)
          for(int v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
            m2p_parametricImages[p][v] = 0.;
      }
      

      string path_to_VAimage;
      dtag = "VA_image";
      if( ReadDataASCIIFile(m_fileOptions, 
                                     dtag, 
                         &path_to_VAimage, 
                                        1, 
                         KEYWORD_OPTIONAL) == 1)
      {
        Cerr("***** i1TCModel::InitializeSpecific() -> Error while trying to read '"<< dtag <<"' keyword in " << m_fileOptions << endl);
        return 1;
      }
      
      if(!path_to_VAimage.empty())
      {
        mp_VaImage = new FLTNB[ mp_ID->GetNbVoxXYZ() ];
        
        if(IntfReadImage(path_to_VAimage, mp_VaImage, mp_ID, m_verbose, INTF_LERP_DISABLED))
        {
          Cerr("***** i1TCModel::InitializeSpecific()-> Error reading Interfile : " << path_to_VAimage << " !" << endl);  
          return 1;
        }
        
        // Check all voxels in VaImage = [0 ,1]
        for(int v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
        {
          if(mp_VaImage[ v ] <0. 
           ||mp_VaImage[ v ] >1.)
          {
            Cerr("***** i1TCModel::InitializeSpecific() -> Voxels values in provided Va image ;"<< path_to_VAimage << " must be between [0,1]. One voxel value is "<< m_fileOptions << endl);
            Cerr("                                          Erroneous value found:"<< mp_VaImage[ v ] << endl);
            return 1;
          }
        }
        
        // Set the number of parameters to estimate to 2
        m_nnlsN=2;
      }
    }
    
    else
    {
      Cerr("***** i1TCModel::InitializeSpecific() -> Error while trying to read configuration file at: " << m_fileOptions << endl);
      return 1;
    }
  }
  // --- Data Initialization with a list of options  --- //
  else if( m_listOptions != "" )
  {
    // We expect here the plasma tac only
    
    // Read it
    if (ReadStringOption(m_listOptions,
                         m2p_nestedModelTimeBasisFunctions[ 1 ],
                         mp_ID->GetNbTimeFrames(),
                         ",",
                         "1cpt model configuration"))
    {
      Cerr("***** i1TCModel::InitializeSpecific() -> Failed to correctly read the list of options !" << endl);
      return 1;
    }
    
    // Standard initialization for the parametric images
    for(int p=0 ; p<m_nbModelParam ; p++)
      for(int v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
        m2p_parametricImages[p][v] = 0.;
      
  }
  
  
  // WPO for integrals
  // Half-frame time duration
  mp_DT2 = new FLTNB[mp_ID->GetNbTimeFrames()];
  
  for(int t=0; t<mp_ID->GetNbTimeFrames(); t++) 
    mp_DT2[t] = mp_ID->GetFrameDurationInSec(0,t) * 0.5 ;
    
  if( m_intMethodFlag == METHOD_INT_WPO )
  {
    mp_wpoQ  = new FLTNB[ mp_ID->GetNbTimeFrames() ];
    mp_wpoA  = new FLTNB[ mp_ID->GetNbTimeFrames() ];
    m2p_wpoP  = new FLTNB*[ mp_ID->GetNbThreadsForImageComputation() ];
    m2p_wpoFD = new FLTNB*[ mp_ID->GetNbThreadsForImageComputation() ];
    m2p_wpoBD = new FLTNB*[ mp_ID->GetNbThreadsForImageComputation() ];
    
    for( int th=0 ; th<mp_ID->GetNbThreadsForImageComputation() ; th++ )
    {
      m2p_wpoP[ th ]  = new FLTNB[ mp_ID->GetNbTimeFrames() ];
      m2p_wpoFD[ th ] = new FLTNB[ mp_ID->GetNbTimeFrames() ];
      m2p_wpoBD[ th ] = new FLTNB[ mp_ID->GetNbTimeFrames() ];
    }
    
    // Compute WPO_q and WPO_a
    for(int32_t t=0 ; t<mp_ID->GetNbTimeFrames() ; t++)
    {
      if( t == mp_ID->GetNbTimeFrames() -1)
        mp_wpoQ[t] = 1;
      else
        mp_wpoQ[t] = mp_DT2[t] / mp_DT2[t+1];
    }
    
    
    for(int32_t t=0 ; t<mp_ID->GetNbTimeFrames() ; t++)
    {
      if( t==0 )
        mp_wpoA[t] = 1;
      else if( t==(mp_ID->GetNbTimeFrames() -2) )
        mp_wpoA[t] = 0.5;
      else if( t==(mp_ID->GetNbTimeFrames() -1) )
        mp_wpoA[t] = 0.;
      else
        mp_wpoA[t] = (mp_DT2[ t+1 ] + 2*mp_DT2[ t ] ) / (2 * (mp_DT2[ t+2 ] + mp_DT2[ t+1 ] + mp_DT2[ t ]) );
      //mp_wpoA[t] = (mp_DT2[ t ] + 2*mp_DT2[ t-1 ] ) / (2 * (mp_DT2[ t+1 ] + mp_DT2[ t ] + mp_DT2[ t-1 ]) );
    }
  }
     
    
  // Variables for Tissue TAC
  m2p_ct  = new FLTNB*[mp_ID->GetNbThreadsForImageComputation()];
  m2p_cti = new FLTNB*[mp_ID->GetNbThreadsForImageComputation()];

  for( int th=0 ; th<mp_ID->GetNbThreadsForImageComputation() ; th++ )
  {
    m2p_ct[ th ]  = new FLTNB[mp_ID->GetNbTimeFrames()];
    m2p_cti[ th ] = new FLTNB[mp_ID->GetNbTimeFrames()];
  }

  // --- Input functions ---

  // Check IF has been provided
  for(int t=0; t<mp_ID->GetNbTimeFrames(); t++)
    if(m2p_nestedModelTimeBasisFunctions[ 1 ][ t ] < 0)
    {
      Cerr("***** i1TCModel::InitializeSpecific() -> Error, plasma input curve has not been initialized !" << endl);
      return 1;
    }
    
  // Compute integration of IF if not provided
  if( m2p_nestedModelTimeBasisFunctions[ 0 ][ 0 ] < 0) 
  {
    // --- Compute integral ---
    FLTNB* Icp = new FLTNB[ mp_ID->GetNbTimeFrames() ];
    
    // Compute plasma tac integral from plasma tac in m2p_nestedModelTimeBasisFunctions[ 1 ]
    IntegrateTAC(m2p_nestedModelTimeBasisFunctions[ 1 ], Icp, 0);
    
    // Recover in model function vectors
    for(int t=0; t<mp_ID->GetNbTimeFrames(); t++) 
      m2p_nestedModelTimeBasisFunctions[ 0 ][ t ] = Icp[ t ];
    
    delete[] Icp;
  }


  // NNLS variables
  m3p_nnlsA    = new FLTNB** [mp_ID->GetNbThreadsForImageComputation()];
  m2p_nnlsB    = new FLTNB*  [mp_ID->GetNbThreadsForImageComputation()];
  m2p_nnlsMat  = new FLTNB*  [mp_ID->GetNbThreadsForImageComputation()];
  
  for( int th=0 ; th<mp_ID->GetNbThreadsForImageComputation() ; th++ )
  {
    // Init 2D coefficient matrix for NNLS estimation
    m3p_nnlsA[ th ] = new FLTNB*[ m_nnlsN ];

    for(int n=0 ; n<m_nnlsN ; n++)
      m3p_nnlsA[ th ][n] = new FLTNB[ mp_ID->GetNbTimeFrames() ];
  
    // Init solution vector and working matrix for NNLS estimation
    m2p_nnlsB[ th ]   = new FLTNB[ mp_ID->GetNbTimeFrames() ];
    m2p_nnlsMat[ th ] = new FLTNB[ (m_nnlsN+2) * mp_ID->GetNbTimeFrames() ];
  }
  
  
  // Compute weights for NNLS
  mp_w = new FLTNB[ mp_ID->GetNbTimeFrames() ];

  for(int t=0; t<mp_ID->GetNbTimeFrames(); t++) 
    mp_w[t] = sqrt(mp_ID->GetFrameDurationInSec(0,t)) ;  // (convert time in min);

  // Init matrix for LS method
  if( m_OptimisationMethodFlag== METHOD_1CPT_LS)
  {
    mp_Y = new oMatrix*[mp_ID->GetNbThreadsForImageComputation()];
    mp_X = new oMatrix*[mp_ID->GetNbThreadsForImageComputation()];
    mp_Xt = new oMatrix*[mp_ID->GetNbThreadsForImageComputation()];
    mp_XtX = new oMatrix*[mp_ID->GetNbThreadsForImageComputation()];
    mp_LSnum = new oMatrix*[mp_ID->GetNbThreadsForImageComputation()];
    mp_LSden = new oMatrix*[mp_ID->GetNbThreadsForImageComputation()];
    mp_RRnum = new oMatrix*[mp_ID->GetNbThreadsForImageComputation()];
    mp_Theta = new oMatrix*[mp_ID->GetNbThreadsForImageComputation()];
    mp_RRm = new oMatrix(m_nnlsN, 1);
    mp_RRw = new oMatrix(m_nnlsN, m_nnlsN);

    for( int th=0 ; th<mp_ID->GetNbThreadsForImageComputation() ; th++ )
    {
      mp_Y[ th ]      = new oMatrix(mp_ID->GetNbTimeFrames(),1);
      mp_X[ th ]      = new oMatrix(mp_ID->GetNbTimeFrames(), m_nnlsN);
      mp_Xt[ th ]     = new oMatrix(m_nnlsN, mp_ID->GetNbTimeFrames());
      mp_XtX[ th ]    = new oMatrix(m_nnlsN, m_nnlsN);
      mp_LSnum[ th ]  = new oMatrix(m_nnlsN, 1);
      mp_LSden[ th ]  = new oMatrix(m_nnlsN, m_nnlsN);
      mp_Theta[ th ]  = new oMatrix(m_nnlsN, 1);
      mp_RRnum[ th ] = new oMatrix(m_nnlsN, 1);
    }

    // Init ridge regression matrices
    if( m_ridgeRegressionFlag )
    {
      // Init RR diagonal  matrices
      for (int pl=0 ; pl<m_nnlsN ; pl++)
        for (int pc=0 ; pc<m_nnlsN ; pc++)
        {
          if( pl==pc)
          {
            mp_RRw->SetMatriceElt( pl , pc , 1 / ( (mp_parUpperBounds[ pl ]-mp_parLowerBounds[ pl ])*(mp_parUpperBounds[ pl ]-mp_parLowerBounds[ pl ]) ) ) ;
            mp_RRm->SetMatriceElt( pl ,  0 , (mp_parUpperBounds[ pl ] + mp_parLowerBounds[ pl ]) / 2. ) ;
            
            for(int th=0 ; th<mp_ID->GetNbThreadsForImageComputation()  ; th++)
              mp_RRnum[ th ]->SetMatriceElt( pl ,  0 , 0. ) ;
          }
          else
          {
            mp_RRw->SetMatriceElt( pl , pc , 0 ) ;
          }
        }
  
    }
  }
 
    
  // Display TACs and other info
  if(m_verbose >=2)
  {
    Cout("i1TCModel::InitializeSpecific() -> Input Plasma TAC coefficients :" << endl);
    Cout("                              ");
    for(int fr=0 ; fr<mp_ID->GetNbTimeFrames() ; fr++)
      Cout(m2p_nestedModelTimeBasisFunctions[0][fr] << ", ");
    Cout(endl);
    Cout("i1TCModel::InitializeSpecific() -> Plasma TAC coefficients :" << endl);
    Cout("                              ");
    for(int fr=0 ; fr<mp_ID->GetNbTimeFrames() ; fr++)
      Cout(m2p_nestedModelTimeBasisFunctions[1][fr] << ", ");
    Cout(endl);
    
    if(m_noImageUpdateFlag)
    Cout("i1TCModel::InitializeSpecific() -> Image update from estimated parametric images is disabled" << endl);
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
  \brief Estimate K1, k2, Va parametric images
  \return 0 if success, other value otherwise.
*/
int i1TCModel::EstimateModelParameters(oImageSpace* ap_ImageS, int a_ite, int a_sset) 
{
  if(m_verbose >=2) Cout("i1TCModel::EstimateModelParameters ..." <<endl);

  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** i1TCModel::EstimateModelParameters() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif

  switch (m_OptimisationMethodFlag)
  {
    case METHOD_1CPT_LS:
      if( EstimateModelParametersWithLS(ap_ImageS) )
      {
        Cerr("***** i1TCModel::EstimateModelParameters() -> A problem occured when estimating 1cpt model parameters with LS method !" << endl);
        return 1;
      }
      break;
      
    case METHOD_1CPT_BF   :
      if ( EstimateModelParametersWithBF(ap_ImageS) )
      {
        Cerr("***** i1TCModel::EstimateModelParameters() -> A problem occured when estimating 1cpt model parameters with Basis functions method !" << endl);
        return 1;
      } 
      break;
      
    case METHOD_1CPT_NNLS    :
      if ( EstimateModelParametersWithNNLS(ap_ImageS) )
      {
        Cerr("***** i1TCModel::EstimateModelParameters() -> A problem occured when estimating 1cpt model parameters with NNLS method !" << endl);
        return 1;
      }
      break; 
  }


  // --- Voxel-based Constraints ---

  // Loop on all voxels
  int v;
  #pragma omp parallel for private(v) schedule(guided)
  for(v=0; v<mp_ID->GetNbVoxXYZ(); v++) 
  {    
    // Set to 0. if value is Nan or < 0.
    if (!isfinite(m2p_parametricImages[0][v])
    ||  !isfinite(m2p_parametricImages[1][v])
    ||  !isfinite(m2p_parametricImages[2][v]) 
      )
    {
      m2p_parametricImages[2][v] = 0.;
      m2p_parametricImages[1][v] = 0.;
      m2p_parametricImages[0][v] = 0.;
    }

    // Biological constraints on Va (positivity and max == 1)
    if( m2p_parametricImages[2][v] < 0.) m2p_parametricImages[2][v]=0.;
    if( m2p_parametricImages[2][v] > 1.) m2p_parametricImages[2][v]=1.;
 
 
    // Biological constaints on k2
    // (set to 0 if:
    // -> Va ~= 1
    // -> Va << && K1 <<
    if( m2p_parametricImages[2][v] > 0.95 || 
      ( m2p_parametricImages[2][v] < 0.1 && m2p_parametricImages[0][v]<0.00001) )
    m2p_parametricImages[1][v] = 0;

    // Set debug image
    if( // At least one parameter is non nil
        ( m2p_parametricImages[ 0 ][ v ]>0. 
       || m2p_parametricImages[ 1 ][ v ]>0. 
       || m2p_parametricImages[ 2 ][ v ]>0.
         )
         // No negative parameter
    &&  ( m2p_parametricImages[ 0 ][ v ]>=0. 
       && m2p_parametricImages[ 1 ][ v ]>=0. 
       && m2p_parametricImages[ 2 ][ v ]>=0. ) 
        ) 
      mp_blackListedvoxelsImage[ v ] = (mp_blackListedvoxelsImage[v] == 2) ? 2 : 0.;
    else
      mp_blackListedvoxelsImage[ v ] = (mp_blackListedvoxelsImage[v] == 2) ? 2 : 1.;

  }
  
  
  return 0;
}
  


// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn EstimateModelParametersWithNNLS
  \param ap_ImageS : pointer to the ImageSpace
  \brief Estimate K1, k2, Va parametric images using NLLS
  \return 0 if success, other value otherwise.
*/
int i1TCModel::EstimateModelParametersWithNNLS(oImageSpace* ap_ImageS) 
{
  if(m_verbose >=2) Cout("i1TCModel::EstimateModelParametersWithNNLS ..." <<endl);

  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** i1TCModel::EstimateModelParametersWithNNLS() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif
  
  int rg=0, cg=0;
    
  
  // Loop on all voxels
  int v;
  #pragma omp parallel for private(v) schedule(guided)
  for(v=0; v<mp_ID->GetNbVoxXYZ(); v++) 
  {
    // If a mask has been provided, check if the model applies in this voxel
    if(mp_maskModel != NULL && mp_maskModel[ v ] == 0)
      continue;
      
    int th=0;
    #ifdef CASTOR_OMP
    th = omp_get_thread_num();
    #endif
    
    uint32_t nb_tf = mp_ID->GetNbTimeFrames();
        
    FLTNB** _2p_nnlsA  = m3p_nnlsA[ th ];
    FLTNB*  _p_nnlsB   = m2p_nnlsB[ th ];
    FLTNB*  _p_nnlsMat = m2p_nnlsMat[ th ];
    
    // m-array of working space
    FLTNB *nnls_zz;
    
    // Solution vectors
    FLTNB  *nnls_x = new FLTNB[ m_nnlsN ];
  
    //n-array of working space for the dual solution y 
    FLTNB  *nnls_wp = new FLTNB[ m_nnlsN ];
    
    // n-array of working space
    int   *nnls_index = new int[ m_nnlsN ];
    
    // Euclidean norm of the residual vector 
    FLTNB  nnls_rnorm;

    FLTNB *dptr=_p_nnlsMat;
    dptr=_p_nnlsMat;
    for(int n=0 ; n<m_nnlsN; n++) 
    {
      _2p_nnlsA[n]=dptr;
      dptr+=nb_tf;
    }
  
    _p_nnlsB=dptr;
    dptr+=nb_tf;
    nnls_zz=dptr;

    // Integral of plasma TAC
    FLTNB* cpi = m2p_nestedModelTimeBasisFunctions[0];
    // Plasma TAC
    FLTNB* cp  = m2p_nestedModelTimeBasisFunctions[1];
    
    // Pointers for tissue TAC and integral of tissue TAC
    FLTNB* ct  = m2p_ct[th];
    FLTNB* cti = m2p_cti[th];
        
    //------------------------
    //  Estimate K1, k2 and Va
    //
    // Recover CPET value
    for(size_t t=0; t<nb_tf; t++) 
      ct[t] = ap_ImageS->m4p_image[t][rg][cg][v];

    // If a VaImage has been provided, subtract Va*Ca from tissue TAC
    if( mp_VaImage != NULL)
    {
      for(size_t t=0; t<nb_tf; t++) 
      {
        ct[t] -= mp_VaImage[v]*cp[t] ;
        if (ct[t]<0) ct[t]=0.;
      }
    }

    // Compute integral of tissue TAC
    IntegrateTAC(ct, cti, th);
    
    // Estimate K1, k2, Va parameters with NNLS 
    if( mp_VaImage == NULL)
    {
      // Set up NNLS variables
      for(size_t t=0; t<nb_tf; t++)
      {
        if(cti[t]<0) cti[t] = 0.;
        // Fill NNLS A matrix:
        // function #1: tissue integral x -1
        _2p_nnlsA[0][t]=-cti[t]*mp_w[t];
        // function #2: integral of input
        _2p_nnlsA[1][t]=cpi[t]*mp_w[t];
        // function #3: input curve
        _2p_nnlsA[2][t]=cp[t]*mp_w[t];
          
        // Fill NNLS B array: tissue
        _p_nnlsB[t]=ct[t]*mp_w[t];
      }
    
      if(NNLS(_2p_nnlsA, nb_tf, m_nnlsN, _p_nnlsB, nnls_x, &nnls_rnorm,
               nnls_wp, nnls_zz, nnls_index) )
        continue; // no solution is possible

      // Computer 1cpt model K1, k2, Va from NNLS estimated parameters
      FLTNB Va=nnls_x[2]; 
    
      // Recover K1
      m2p_parametricImages[0][v] = nnls_x[1];
      // If Va>0, we fitted to CPET instead of CT
      // In this case we didn't estimate K1, but K1 + Va*k2 as there is a Va*k2 contribution to CPi
      // We have to substract this contribution to get K1
      m2p_parametricImages[0][v] -= Va*nnls_x[0];
      
      m2p_parametricImages[1][v] = nnls_x[0];
      m2p_parametricImages[2][v] = Va;
    }
    
    else // Va value from user provided image
    {
      FLTNB Va = mp_VaImage[v];
            
      // Fill matrix without blood volume contribution
      for(size_t t=0; t<nb_tf; t++)
      {
        // Fill NNLS A matrix:
        // function #1: tissue integral x -1
        _2p_nnlsA[0][t] = -cti[t]*mp_w[t];
        // function #2: integral of input
        _2p_nnlsA[1][t] = cpi[t]*mp_w[t];
          
        // Fill NNLS B array: tissue
        _p_nnlsB[t] = ct[t]*mp_w[t];
      }
      
      // Estimate with NNLS (nnls_n - 1 as Vb is already corrected)
      if(NNLS(_2p_nnlsA, nb_tf, m_nnlsN, _p_nnlsB, nnls_x, &nnls_rnorm,
               nnls_wp, nnls_zz, nnls_index) )  
         continue; // no solution is possible
  
      // Recover K1
      m2p_parametricImages[0][v] = nnls_x[1];
      m2p_parametricImages[1][v] = nnls_x[0];
      m2p_parametricImages[2][v] = Va;
    }
    
    // Free memory
    if(nnls_x)     delete[]nnls_x;
    if(nnls_wp)    delete[]nnls_wp;
    if(nnls_index) delete[]nnls_index;
    
  } // multithreaded loop on voxels
  
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn EstimateModelParametersWithBF
  \param ap_ImageS : pointer to the ImageSpace
  \brief Estimate K1, k2, Va parametric images using basis functions approach
  \return 0 if success, other value otherwise.
*/
int i1TCModel::EstimateModelParametersWithBF(oImageSpace* ap_ImageS) 
{
  if(m_verbose >=2) Cout("i1TCModel::EstimateModelParametersWithBF ..." <<endl);

  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** i1TCModel::EstimateModelParametersWithBF() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif
  
  Cerr("i1TCModel::EstimateModelParametersWithBF -> Not yet implemented !!!" <<endl);
  return 1;
}










// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn EstimateModelParametersWithLS
  \param ap_ImageS : pointer to the ImageSpace
  \brief Estimate K1, k2, Va parametric images using LS
  \return 0 if success, other value otherwise.
*/
int i1TCModel::EstimateModelParametersWithLS(oImageSpace* ap_ImageS) 
{
  if(m_verbose >=2) Cout("i1TCModel::EstimateModelParametersWithLS ..." <<endl);

  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** i1TCModel::EstimateModelParametersWithLS() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif
  
  int rg=0, cg=0;

  uint32_t nb_tf = mp_ID->GetNbTimeFrames();

  bool error_in_th_loop = false;
  string error_msg ="";
  
  
  // Loop on all voxels
  int v;
  #pragma omp parallel for private(v) schedule(guided)
  for(v=0; v<mp_ID->GetNbVoxXYZ(); v++) 
  {
    // If a mask has been provided, check if the model applies in this voxel
    if(mp_maskModel != NULL && mp_maskModel[ v ] == 0)
      continue;
    
    int th=0;
    #ifdef CASTOR_OMP
    th = omp_get_thread_num();
    #endif

    // Pointers for tissue TAC and integral of tissue TAC
    FLTNB* ct  = m2p_ct[th];
    FLTNB* cti = m2p_cti[th];
    
    // Integral of plasma TAC
    FLTNB* cpi = m2p_nestedModelTimeBasisFunctions[0];
    // Plasma TAC
    FLTNB* cp  = m2p_nestedModelTimeBasisFunctions[1];
    
    // Recover CPET value 
    for(size_t t=0; t<nb_tf; t++) 
      ct[t] = ap_ImageS->m4p_image[t][rg][cg][v];
    
    // If a VaImage has been provided, subtract Va*Ca from tissue TAC
    
    if( mp_VaImage != NULL)
    {
      for(size_t t=0; t<nb_tf; t++) 
      {
        ct[t] -= mp_VaImage[v]*cp[t] ;
        if (ct[t]<0) ct[t]=0.;
      }
    }
    
    // Compute integral of tissue TAC
    IntegrateTAC(ct, cti, th);
    
    // Solution vector
    FLTNB  *LS_res = new FLTNB[ m_nnlsN ];
    
    // Model matrix
    FLTNB** LS_Mod  = m3p_nnlsA[ th ];
        
    for(int16_t t=0; t<mp_ID->GetNbTimeFrames(); t++)
    {
      // function #1: tissue integral x -1
      LS_Mod[0][t] = -cti[t];
      // function #2: integral of input
      LS_Mod[1][t] = cpi[t];
      // function #3: input curve
      if( m_nnlsN>2) // Va image not provided
        LS_Mod[2][t] = cp[t];
    }
    
    // LS optimization using ridge-regression or not
    if (m_ridgeRegressionFlag )
    {
      // Estimate parameters with LS with ridge-regression
      if(RRLS(m_nnlsN,
              mp_ID->GetNbTimeFrames(),
              LS_Mod,
              ct,
              mp_w,
              LS_res) )
      {
        error_in_th_loop = true;
        error_msg = "***** i1TCModel::EstimateModelParametersWithLS() -> An error occurred when minimizing WRSS function with ridge-regression!"; 
      }
    }
    else
    {
      // Estimate parameters with standard LS
      if(LS(m_nnlsN,
            mp_ID->GetNbTimeFrames(),
            LS_Mod,
            ct,
            mp_w,
            LS_res) )
      {
        error_in_th_loop = true;
        error_msg = "***** i1TCModel::EstimateModelParametersWithLS() -> An error occurred when minimizing WRSS function with standard LS!"; 
      }    
    }
    
    // Recover Va from LS parameter, or Va image if provided
    m2p_parametricImages[2][v] = (mp_VaImage == NULL) ? 
                                          LS_res[ 2 ] : 
                                        mp_VaImage[v] ;
                                        
    // k2
    m2p_parametricImages[1][v] = LS_res[ 0 ];
    //K1
    m2p_parametricImages[0][v] = LS_res[ 1 ] 
                               - m2p_parametricImages[2][v]*LS_res[ 0 ];
    
    if( LS_res ) delete[] LS_res;
  }
  
  if ( error_in_th_loop == true )
  {
    Cerr(error_msg << endl);
    return 1;
  }
  
  return 0;
}



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn RRLS
  \param a_nP      : number of parameters (P)
  \param a_nT      : number of data samples (T)
  \param a2p_model : matrix containing the model functions (TxP elts)
  \param ap_data   : data vector (T elts)
  \param ap_w      : vector containing the weights (T elts)
  \param ap_result : vector recovering the estimated parameters (P elts)
  \brief   Non-linear least square estimator with Ridge-Regression
  \details This function will estimate a set of parameters (O) by 
           minimizing the weighted residual sum of squares (WRSS) 
           difference between the data and the model function
           Ô = [ X'*W*X + t.Rw ]-1 [ X'*W*y ]  
             + [ X'*W*X + t.Rw ]-1 [ t.Rw*Rm ]
           y = data vector
           X = model matrix
           W = weights vector
           t = Ridge constant
           Rw= Ridge weights
           Rm= Ridge means
           
  \return 0 if success, other value otherwise.
*/
int i1TCModel::RRLS(uint16_t  a_nP,
                     uint16_t  a_nT,
                       FLTNB **a2p_model,
                       FLTNB  *ap_data,
                       FLTNB  *ap_w,
                       FLTNB  *ap_result
                    )
{
  #ifdef CASTOR_DEBUG
  if(m_verbose >=4) Cout("i1TCModel::RRLS ..." <<endl);
  
  if (!m_initialized)
  {
    Cerr("***** i1TCModel::RRLS() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif
  
  // Check that the LS flag is on, so that we are sure that all variables 
  // have been correctly allocated
  if(m_OptimisationMethodFlag != METHOD_1CPT_LS)
  {
    Cerr("***** i1TCModel::RRLS() -> Function is called while the optimization method is not Least-square !" << endl);
    return 1;
  }
  
  int th=0;
  #ifdef CASTOR_OMP
  th = omp_get_thread_num();
  #endif

  // Get the working matrices
  oMatrix*  Y       = mp_Y[ th ];
  oMatrix*  X       = mp_X[ th ];
  oMatrix*  Xt      = mp_Xt[ th ];
  oMatrix*  XtX     = mp_XtX[ th ];
  oMatrix*  Theta   = mp_Theta[ th ];
  oMatrix*  LSnum   = mp_LSnum[ th ];
  oMatrix*  LSden   = mp_LSden[ th ];
  
  oMatrix *RRw   = mp_RRw; 
  oMatrix *RRm   = mp_RRm; 
  oMatrix *RRnum = mp_RRnum[ th ];
  
  // Init data vector
  for( uint16_t t=0 ; t<a_nT ; t++ )
    Y->SetMatriceElt( t , 0, ap_data[ t ]*ap_w[ t ] );

  // Init model matrix
  for( uint16_t p=0 ; p<a_nP ; p++ )
    for( uint16_t t=0 ; t<a_nT ; t++ )
      X->SetMatriceElt( t , p , a2p_model[ p ][ t ] );


  #ifdef CASTOR_DEBUG
    if(m_verbose >=4)
    {
      cout << " ----------------------------------------------- " << endl;
      cout << " X->Describe() " << endl;
      cout << " ----------------------------------------------- " << endl;
    }
  #endif
  
  // Compute X'
  X->Transpose(Xt);

  // Once X' is computed, add weights to X
  for( uint16_t p=0 ; p<a_nP ; p++ )
    for( uint16_t t=0 ; t<a_nT ; t++ )
      X->SetMatriceElt( t , p , a2p_model[ p ][ t ] * ap_w[ t ] );
      
      
  #ifdef CASTOR_DEBUG
    if(m_verbose >=4)
    {
      cout << " ----------------------------------------------- " << endl;
      cout << " Xt->Describe() " << endl;
      Xt->Describe();
      cout << " ----------------------------------------------- " << endl;
    }
  #endif
  
  // Compute (X'*W*X), recover the result in XtX matrix
  Xt->Multiplication( X, XtX );

  // Add Ridge weights, recover result in XtX
  for (int pl=0 ; pl<a_nP ; pl++)
    for (int pc=0 ; pc<a_nP ; pc++)
      XtX->SetMatriceElt( pl , 
                          pc , 
                          XtX->GetMatriceElt(pl,pc) + m_RRcst*RRw->GetMatriceElt(pl,pc) );
      
      
  #ifdef CASTOR_DEBUG
    if(m_verbose >=4)
    {
      cout << " ----------------------------------------------- " << endl;
      cout << " X'WX->Describe() " << endl;
      XtX->Describe();
      cout << " ----------------------------------------------- " << endl;
    }
  #endif

  // Compute LS denominator [X'*W*X + t.*Rw]-1, recover the result in LSden matrix
  // If nil on diagonal, set all elts of result vector to 0. and stop here 
  if ( XtX->Inverse( LSden ) == 1 )
  {
    for( uint16_t p=0 ; p<a_nP ; p++ )
      ap_result[ p ] = 0.;
      
    return 0.;
  }

  #ifdef CASTOR_DEBUG
  if(m_verbose >=4)
    { 
      cout << " ----------------------------------------------- " << endl;
      cout << " RRLS denominator ->Describe() " << endl;
      LSden->Describe();
      cout << " ----------------------------------------------- " << endl;
    }
  #endif
  
  // Compute LS numerator (X'*W*y), recover the result in LS_num matrix
  Xt->Multiplication( Y, LSnum );

  #ifdef CASTOR_DEBUG
    if(m_verbose >=4)
    { 
      cout << " ----------------------------------------------- " << endl;
      cout << " LS numerator ->Describe() " << endl;
      LSnum->Describe();
      cout << " ----------------------------------------------- " << endl;
    }
  #endif

  // Compute t.Rw*Rm and write result in rrtWM
  RRw->Multiplication( RRm, RRnum );

  // Compute complete numerator [ X'*W*y ] [ t.Rw*Rm ], recover in LSnum
  for (int p=0 ; p<a_nP ; p++)
    LSnum->SetMatriceElt( p , 0 , LSnum->GetMatriceElt(p,0) 
                                + m_RRcst * RRnum->GetMatriceElt(p,0) ) ;
    
  #ifdef CASTOR_DEBUG
    if(m_verbose >=4)
    { 
      cout << " ----------------------------------------------- " << endl;
      cout << " RRLS numerator ->Describe() " << endl;
      LSnum->Describe();
      cout << " ----------------------------------------------- " << endl;
    }
  #endif
  
  
  // Compute final parameters, recover in 'Theta'
  // Ô = [ X'*W*X + t.Rw ]-1 [ X'*W*y ]  
  //   + [ X'*W*X + t.Rw ]-1 [ t.Rw*Rm ],    
  LSden->Multiplication( LSnum, Theta );
;
  #ifdef CASTOR_DEBUG
    if(m_verbose >=4)
    { 
      cout << " ----------------------------------------------- " << endl;
      cout << " Theta->Describe() " << endl;
      Theta->Describe();
    }
  #endif
  
  // Recover estimated parameters in return vector
  for( uint16_t p=0 ; p<a_nP ; p++ )
    ap_result[ p ] = Theta->GetMatriceElt(p,0);

  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn LS
  \param a_nP      : number of parameters (P)
  \param a_nT      : number of data samples (T)
  \param a2p_model : matrix containing the model functions (TxP elts)
  \param ap_data   : data vector (T elts)
  \param ap_w      : vector containing the weights (T elts)
  \param ap_result : vector recovering the estimated parameters (P elts)
  \brief   Non-linear least square estimator 
  \details This function will estimate a set of parameters (O) by 
           minimizing the weighted residual sum of squares (WRSS) 
           difference between the data and the model function
           Ô = [ X'WX ]-1 X'Wy 
           y = data vector
           X = model matrix
           W = weights vector
           
  \return 0 if success, other value otherwise.
*/
int i1TCModel::LS(uint16_t  a_nP,
                   uint16_t  a_nT,
                     FLTNB **a2p_model,
                     FLTNB  *ap_data,
                     FLTNB  *ap_w,
                     FLTNB  *ap_result
                  )
{
  #ifdef CASTOR_DEBUG
  if(m_verbose >=4) Cout("i1TCModel::LS ..." <<endl);
      
  if (!m_initialized)
  {
    Cerr("***** i1TCModel::LS() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif
  
  // Check that the LS flag is on, so that we are sure that all variables 
  // have been correctly allocated
  if(m_OptimisationMethodFlag != METHOD_1CPT_LS)
  {
    Cerr("***** i1TCModel::LS() -> Function is called while the optimization method is not Least-square !" << endl);
    return 1;
  }
  
  int th=0;
  #ifdef CASTOR_OMP
  th = omp_get_thread_num();
  #endif

  // Get the working matrices
  oMatrix*  Y       = mp_Y[ th ];
  oMatrix*  X       = mp_X[ th ];
  oMatrix*  Xt      = mp_Xt[ th ];
  oMatrix*  XtX     = mp_XtX[ th ];
  oMatrix*  Theta   = mp_Theta[ th ];
  oMatrix*  LSnum   = mp_LSnum[ th ];
  oMatrix*  LSden   = mp_LSden[ th ];
  
  // Init data vector
  for( uint16_t t=0 ; t<a_nT ; t++ )
    Y->SetMatriceElt( t , 0, ap_data[ t ]*ap_w[ t ] );

  // Init model matrix
  for( uint16_t p=0 ; p<a_nP ; p++ )
    for( uint16_t t=0 ; t<a_nT ; t++ )
      X->SetMatriceElt( t , p , a2p_model[ p ][ t ] );


  #ifdef CASTOR_DEBUG
    if(m_verbose >=4)
    {
      cout << " ----------------------------------------------- " << endl;
      cout << " X->Describe() " << endl;
      cout << " ----------------------------------------------- " << endl;
    }
  #endif
  
  // Compute X'
  X->Transpose(Xt);

  // Once X' is computed, add weights to X
  for( uint16_t p=0 ; p<a_nP ; p++ )
    for( uint16_t t=0 ; t<a_nT ; t++ )
      X->SetMatriceElt( t , p , a2p_model[ p ][ t ] * ap_w[ t ] );
      
      
  #ifdef CASTOR_DEBUG
    if(m_verbose >=4)
    {
      cout << " ----------------------------------------------- " << endl;
      cout << " Xt->Describe() " << endl;
      Xt->Describe();
      cout << " ----------------------------------------------- " << endl;
    }
  #endif
  
  // Compute (X'WX), recover the result in XtX matrix
  Xt->Multiplication( X, XtX );

  #ifdef CASTOR_DEBUG
    if(m_verbose >=4)
    {
      cout << " ----------------------------------------------- " << endl;
      cout << " X'WX->Describe() " << endl;
      XtX->Describe();
      cout << " ----------------------------------------------- " << endl;
    }
  #endif
  
  // Compute LS numerator (X'Wy), recover the result in LS_num matrix
  Xt->Multiplication( Y, LSnum );
 
  #ifdef CASTOR_DEBUG
    if(m_verbose >=4)
    { 
      cout << " ----------------------------------------------- " << endl;
      cout << " LS numerator ->Describe() " << endl;
      LSnum->Describe();
      cout << " ----------------------------------------------- " << endl;
    }
  #endif

  // Compute LS denominator [X'WX]-1, recover the result in LSden matrix
  // If nil on diagonal, set all elts of result vector to 0. and stop here 
  if ( XtX->Inverse( LSden ) == 1 )
  {
    for( uint16_t p=0 ; p<a_nP ; p++ )
      ap_result[ p ] = 0.;
    
    return 0.;
  }

  #ifdef CASTOR_DEBUG
  if(m_verbose >=4)
    { 
      cout << " ----------------------------------------------- " << endl;
      cout << " LS denominator ->Describe() " << endl;
      LSden->Describe();
      cout << " ----------------------------------------------- " << endl;
    }
  #endif
  
  // Get final parameter Ô = [ X'WX ]-1 X'Wy 
  LSden->Multiplication( LSnum, Theta );

  #ifdef CASTOR_DEBUG
  if(m_verbose >=4)
  { 
    cout << " ----------------------------------------------- " << endl;
    cout << " Theta->Describe() " << endl;
    Theta->Describe();
  }
  #endif
  
  // Recover estimated parameters in return vector
  for( uint16_t p=0 ; p<a_nP ; p++ )
    ap_result[ p ] = Theta->GetMatriceElt(p,0);
  
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
  \brief Estimate image using model parametric images and basis functions
  \todo  Maybe add the possibility to use specific thresholds
         ex: if Va~1 (eg >0.8) && K1~0 && k2~0 --> Va=1, K1=k2=0
  \return 0 if success, other value otherwise.
*/
int i1TCModel::EstimateImageWithModel(oImageSpace* ap_ImageS, int a_ite, int a_sset) 
{
  if(m_verbose >= 2) Cout("i1TCModel::EstimateImageWithModel ... " <<endl);

  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** i1TCModel::EstimateImageWithModel() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif

  if( !m_noImageUpdateFlag  
  &&  a_ite >= m_startIteUpdateFlag)
  {
    // Compute average of time shift between two samples in TACs formulas
    // Required for TACs computation
    HPFLTNB* DTavg = new HPFLTNB[mp_ID->GetNbTimeFrames() ]; 
    DTavg[ 0 ] = mp_DT2[ 0 ];
    for (int32_t t=1 ; t<mp_ID->GetNbTimeFrames() ; t++)
      DTavg[ t ] = (mp_DT2[ t ] + mp_DT2[ t-1 ] ) / 2;
      
    for (int rg=0 ; rg<mp_ID->GetNbRespGates() ; rg++) 
      for (int cg=0 ; cg<mp_ID->GetNbCardGates() ; cg++)
      {
        int v;
        #pragma omp parallel for private(v) schedule(guided)
        for (v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
        {
          FLTNB K1 = m2p_parametricImages[0][v];
          FLTNB k2 = m2p_parametricImages[1][v];
          FLTNB Va = m2p_parametricImages[2][v];
          
          // Keep reconstructed image value if one parameter value <0
          // or if all parameters == 0
          if( ( K1>0. || k2>0. ||  Va>0.) // At least one parameter is non nil
          &&  ( K1>=0. && k2>=0. && Va>=0. ) ) // No negative parameter
          {
            FLTNB ct=0., // tissue tac
                b_ct=0., // tissue tac previous sample
                 cti=0., // cumulative integral of tissue tac
                   b=0., // output tissue concentration
               bb_ct=0., // tissue tact penultinum previous sample
                n_ct=0.; // tissue_tac next sample 
                    
            for (int32_t t=0 ; t<mp_ID->GetNbTimeFrames() ; t++)
            {
              bb_ct = b_ct;
              b_ct  = ct;
              b = cti + DTavg[t]*b_ct ;
              
              ct = ( K1*m2p_nestedModelTimeBasisFunctions[ 0 ][ t ] - k2*b ) 
                 / ( 1 + k2*DTavg[ t ] );
                              

              // Approximating next cti with weighed parabola method
              if( m_intMethodFlag == METHOD_INT_WPO ) 
              {
                // Pre-estimate ct+1 (n_ct) with Trapz
                if(t>0)
                {
                  FLTNB n_cti = cti + mp_DT2[t]*(ct + b_ct);
                  FLTNB n_b = n_cti + mp_DT2[t]*ct ; // output tissue concentration
                  n_ct = (K1*m2p_nestedModelTimeBasisFunctions[0][t] - k2*n_b) 
                       / (1 + k2*DTavg[t]);
                }
                
                cti += WPOinc( t, ct, b_ct, bb_ct, n_ct);
              }
              else // Approximating cti with trapezoidal rule
                cti += mp_DT2[t]*(ct + b_ct);
              
              // CPET (add blood contribution)
              FLTNB new_value  = ct + m2p_nestedModelTimeBasisFunctions[1][t]*Va;
              ap_ImageS->m4p_image[t][rg][cg][v] = (new_value > 0) ? 
                                                         new_value : 
                                ap_ImageS->m4p_image[t][rg][cg][v] ;
            }
          } // end of condition on kinetic parameters
          else
            for (int32_t t=0 ; t<mp_ID->GetNbTimeFrames() ; t++)
            ap_ImageS->m4p_image[t][rg][cg][v] = 0;
          
        } // end of loop on voxels
        
      }  // end of loop on rg, cg
      
    delete[] DTavg; 
    
  } // update image condition

  return 0;
}


// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn      WPOinc
  \param   a_time : time point for which the integral will be computed
  \param   tac : tac value at time t
  \param   b_tac : tac value at time t-1
  \param   bb_tac : tac value at time t-2
  \param   n_tac : tac value at time t+1
  \brief   Estimate the next integration value for a specific time point of a tac using WPO
  \return  The estimated value of the integral for this time point
*/
FLTNB i1TCModel::WPOinc(uint32_t a_time, FLTNB tac, FLTNB b_tac, FLTNB bb_tac, FLTNB n_tac)
{
  FLTNB Itac=0.;
  FLTNB wpo_p=0., wpo_p_n=0, wpo_fd_n=0., wpo_bd=0.;
  uint32_t t = a_time;

  if(t==0) // first sample (b_tac, bb_tac not known)
    Itac = mp_DT2[t]*tac*0.5;
  else if( t==1 ) // second sample (bb_tac not known)
  {
    wpo_p_n  = ( tac  -  ( mp_DT2[t+1]*b_tac + mp_DT2[t]*n_tac ) / ( mp_DT2[t+1]+mp_DT2[t] ) )  /   6;
    wpo_fd_n = 2*mp_DT2[1] * ( 0.5*( tac +  b_tac )  + wpo_p_n * mp_wpoQ[t+1] );
    wpo_bd = mp_DT2[t-1] * ( b_tac + tac  );
    
    Itac += (1-mp_wpoA[ t ]) * wpo_bd  +   mp_wpoA[ t ]*wpo_fd_n; // WPO rule
  }
  else if( t < ((uint32_t)mp_ID->GetNbTimeFrames()-1) )
  {
    // todo : try to get tac+1 with trapZ ?           
    wpo_p_n  = (   tac  -  ( mp_DT2[t+1]*b_tac   +  mp_DT2[t]    *n_tac )  / ( mp_DT2[t+1]+ mp_DT2[t]   )   )  /   6;
    wpo_p    = ( b_tac  -  ( mp_DT2[t]  *bb_tac  +  mp_DT2[t-1]  *tac )    / ( mp_DT2[t]  + mp_DT2[t-1] )   )  /   6;
    wpo_fd_n = 2*mp_DT2[t]   * ( 0.5*( tac +  b_tac )  + wpo_p_n * mp_wpoQ[t+1] );
    wpo_bd   = 2*mp_DT2[t-1] * ( 0.5*( b_tac +  tac )  + wpo_p   / mp_wpoQ[t] );
    
    Itac += (1-mp_wpoA[ t ]) * wpo_bd  +   mp_wpoA[ t ]*wpo_fd_n; // WPO rule
  }
  else
  {
    wpo_p    = ( b_tac  -  (  ( ( mp_DT2[t]  *bb_tac  )  + mp_DT2[t-1]  *tac )  ) / ( mp_DT2[t]+mp_DT2[t-1]   )   )  /   6;
    wpo_bd = 2*mp_DT2[t-1] * ( 0.5*( b_tac +  tac   )  + wpo_p / mp_wpoQ[t] );
    
    Itac += wpo_bd; // WPO rule
  }
  return Itac;
  
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*!
  \fn      IntegrateTAC
  \param   ap_tac : vector with nb time frames elements containing the time activity curve from which integral will be computed
  \param   ap_citac : (returned) vector with nb time frames elements recovering the cumulative integral for each time point t 
  \param   a_th : thread index
  \brief   Call one of the TAC integration method
  \return  0 if success, positive value otherwise
*/
int i1TCModel::IntegrateTAC(FLTNB* ap_tac, FLTNB* ap_citac, int a_th)
{
  if( m_intMethodFlag == METHOD_INT_WPO )
    return WPO(ap_tac, ap_citac, a_th);
  else // ( m_intMethodFlag == METHOD_INT_TRAP )
    return Trapz(ap_tac, ap_citac);
          
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn      WPO
  \param   ap_tac : vector with nb time frames elements containing the time activity curve from which integral will be computed
  \param   ap_citac : (returned) vector with nb time frames elements recovering the cumulative integral for each time point t 
  \param   a_th : thread index
  \brief   return integral for each time point t (WPO_S), and cumulative integral (cumWPO_S )
  \return  0 if success, positive value otherwise
*/
int i1TCModel::WPO(FLTNB* ap_tac, FLTNB* ap_citac, int a_th)
{
  uint32_t nb_tf = mp_ID->GetNbTimeFrames(); 

  FLTNB* DT2 = mp_DT2; // Frame duration / 2
  FLTNB itac = 0.; // integral at time point t
  
  // Init vectors 
  FLTNB* WPO_P  = m2p_wpoP[ a_th ];
  FLTNB* WPO_FD = m2p_wpoFD[ a_th ];
  FLTNB* WPO_BD = m2p_wpoBD[ a_th ];
  
  //  WPO_P,  WPO_BD, WPO_FD
  for(uint32_t t=0 ; t<nb_tf ; t++)
  {
    if( t==0 )
    {
      WPO_P[t]  = 0. ;
      WPO_BD[t] = ap_tac[0] ;
      WPO_FD[t] = ap_tac[0] ;
    }
    else if( t==1 )
    {
      WPO_P[t]  = 0. ;
      WPO_FD[t] = DT2[t-1] *( ( ap_tac[ t-1 ] )  );
      WPO_BD[t] = DT2[t-1] *( ( ap_tac[ t-1 ]  +  ap_tac[ t ] ));  
    }
    else
    {
      WPO_P[t]  = ( ap_tac[ t-1 ]  -  ( DT2[t]*ap_tac[ t-2 ]  +  DT2[t-1]*ap_tac[ t ] ) / ( DT2[t] + DT2[t-1] )   )  /   6;
      WPO_FD[t] = 2*DT2[t-1] * ( 0.5 *( ap_tac[ t-1 ] +  ap_tac[ t-2 ] )  + WPO_P[t] * mp_wpoQ[t] );
      WPO_BD[t] = 2*DT2[t-1] * ( 0.5 *( ap_tac[ t-1 ] +  ap_tac[ t ]   )  + WPO_P[t] / mp_wpoQ[t] );
    }
  }
  
  // TODO : adapt Simpson rule to the frame duration (only the isotropic frame in the peak)
  for(uint32_t t=0 ; t<nb_tf ; t++)
  {
    // Compute tac integral at time point t
    if( t==0 )
      itac = ap_tac[ t ] * DT2[ t ] ;
    else if( t == (nb_tf-1) )
      itac = WPO_BD[ t ];
    else
      itac =  (1-mp_wpoA[ t ]) * WPO_BD[ t ]  +   mp_wpoA[ t ]*WPO_FD[ t+1 ]; // WPO rule
    
    // Compute cumumative integral
    if( t==0 )
      ap_citac[ t ] = itac  ;
    else
      ap_citac[ t ] = ap_citac[ t-1] + itac;
  }

  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn      Trapz
  \param   ap_tac : vector with nb time frames elements containing the time activity curve from which integral will be computed 
  \param   ap_citac : (returned) vector with nb time frames elements recovering the cumulative integral for each time point t 
  \brief   return integral for each time point t (WPO_S), and cumulative integral (cumWPO_S )
  \return  0 if success, positive value otherwise
*/
int i1TCModel::Trapz(FLTNB* ap_tac, FLTNB* ap_citac)
{
  uint32_t nb_tf = mp_ID->GetNbTimeFrames(); 
  
  FLTNB* DT2 = mp_DT2; // Frame duration / 2
  FLTNB itac = 0.; // integral at time point t
  
  for(uint32_t t=0 ; t<nb_tf ; t++)
  {
    // Init
    ap_citac[t] = 0;
    
    if( t==0 )
    {
      itac  = ap_tac[0] * DT2[0] *0.5;
      ap_citac[0] = itac;
    }
    else
    {
      itac = ( ap_tac[t]+ap_tac[t-1] ) * DT2[t];
      ap_citac[t] +=   ap_citac[t-1] + itac;
    }
  }

  return 0;
}




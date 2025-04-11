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
  \brief    Implementation of class iLinearPatlakModel
  \todo     Checks regarding optimization on NestedEM loops when updating parametric images
*/

#include "iLinearPatlakModel.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn iLinearPatlakModel
  \brief Constructor of iLinearPatlakModel. Simply set all data members to default values.
*/
iLinearPatlakModel::iLinearPatlakModel() : iLinearModel() 
{
  // Hard coded number as Patlak has always only two basis functions and hence nbModelParameters
  m_nbTimeBF       = 2;
  m_nbModelParam   = 2;
  m_nnlsN = 2;
  m_nbRgateBF = 1; m_nbRGModelParam=1;
  m_nbCgateBF = 1; m_nbCGModelParam=1;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ~iLinearPatlakModel
  \brief Destructor of iLinearPatlakModel
*/
iLinearPatlakModel::~iLinearPatlakModel() 
{
  if(m_initialized)
  {

  }
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ShowHelpModelSpecific
  \brief Print out specific help about the implementation of the Patlak
         model and its initialization
*/
void iLinearPatlakModel::ShowHelpModelSpecific()
{
  cout << "-- This class implements the Patlak Reference Tissue Model : " << endl;
  cout << "-- Patlak CS, Blasberg RG: Graphical evaluation of blood-to-brain transfer constants from multiple-time uptake data" << endl;
  cout << "-- J Cereb Blood Flow Metab 1985, 5(4):5 84-590." << endl;
  cout << "-- DOI http://dx.doi.org/10.1038/jcbfm.1985.87" << endl;
  cout << "-- It is used to model radiotracers which follows as 2-tissue compartment model with irreversible trapping  " << endl;
  cout << "-- The Patlak temporal basis functions are composed of the Patlak slope (integral of the reference TAC from the injection time " << endl;
  cout << "   divided by the instantaneous reference activity), and intercept (reference tissue TAC)  " << endl;
  cout << endl;
  cout << "  It can be initialized using either an ASCII file or a list of option with the following keywords and information :" << endl;
  cout << endl;
  cout << "   For ASCII file options:  " << endl;
  cout << "      As this class inherits from the iLinearModel class, the following parameters must be declared inside the couple of the following specific tags: " << endl;
  cout << "      - DYNAMIC FRAMING/ENDDF " << endl;
  cout << "      - The ASCII file must contain the following keywords :" << endl;
  cout << "      'Basis_functions:'             (optional) Enter the coefficients of the Patlak basis functions for each time frame (tf) " << endl;
  cout << "                                                  on two successive lines, separated by ',' :" << endl;
  cout << "                                      -> Patlak_functions: " << endl;
  cout << "                                      -> coeff_Pplot_tf1,coeff_Pplot_tf2,...,coeff_Pplot_tfn" << endl;
  cout << "                                      -> coeff_Pintc_tf1,coeff_Pintc_tf2,...,coeff_Pintc_tfn" << endl;
  cout << endl;
  cout << "      'AIC_input_file: path/to/file' (optional) As an alternative to direct input of the basis functions, the sampled Arterial Input Function can be given for  " << endl;
  cout << "                                                estimation of the Patlak basis " << endl;
  cout << "                                         The file must contain the following information in successive lines, separated by ',' " << endl;
  cout << "                                      -> AIC_number_of_points: " << endl;
  cout << "                                      -> AIC_time_points: " << endl;
  cout << "                                      -> AIC_data_points: " << endl;
  cout << "                                      -> AIC_units: 'seconds' or 'minutes'  " << endl;
  cout << endl;
  cout << "     'Parametric_image_init: path'   (optional) path to an interfile image to be used as initialization for the parametric images." << endl;
  cout << endl;
  cout << "     'Optimisation_method:' x        (mandatory) optimization method available options: " << endl;
  cout << "                                        x=0: Direct ( Implementation of basis functions side by system matrix within each tomographic iterative loop )  " << endl;
  cout << "                                        x=1: Nested EM  " << endl;
  cout << "                                        x=2: Iterative non-negative Least-Square " << endl;
  cout << "                                             (C.L. Lawson and R.J. Hanson, Solving Least Squares Problems)" << endl;
  cout << "                                        x=3: Least-Squares linear regression " << endl;
  cout << endl;
  cout << "   For command line list of options:  " << endl;
  cout << "      The list of options must contain the coefficients of both Patlak functions separated by commas, with the following template :" << endl;
  cout << "      coeff_Pplot_tf1,coeff_Pplot_tf2,...,coeff_Pplot_tfn,";
  cout << "      coeff_Pintc_tf1,coeff_Pintc_tf2,...,coeff_Pintc_tfn "<< endl;
  cout << "      Default optimization method is Nested EM "<< endl;
  cout << "      Parametric images will be uniformly initialized with 0.001 and 1. for Patlak slope and intercept by default " << endl;
  cout << "      The parametric images estimations will be written on disk for each iteration" << endl;
  cout << "       " << endl;

  // Print general help for all dynamic models
  ShowHelp();

}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ReadAndCheckConfigurationFileSpecific
  \param const string& a_configurationFile : ASCII file containing information about a dynamic model
  \brief This function is used to read options from a configuration file.
  \return 0 if success, other value otherwise.
*/
int iLinearPatlakModel::ReadAndCheckConfigurationFileSpecific()
{
  if(m_verbose >=3) Cout("iLinearPatlakModel::ReadAndCheckConfigurationFileSpecific ..."<< endl); 
  
  // Apply the generic linear parameter read for all Linear Models
  if(ReadAndCheckConfigurationFileSpecificToAllLinearModels())
  {
    Cerr("***** iLinearPatlakModel::ReadAndCheckConfigurationFileSpecific -> Error while trying to read configuration file for generic options of all linear models" << endl);
    return 1;
  }
  // Normal End
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
int iLinearPatlakModel::ReadAndCheckOptionsList(string a_listOptions)
{
  if(m_verbose >=3) Cout("iLinearPatlakModel::ReadAndCheckOptionsList ..."<< endl); 
  
  // Just recover the string here, it will be processed in the Initialize() function
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
int iLinearPatlakModel::CheckSpecificParameters()
{

  // Perform generic checks that apply for the Linear Models
  if(CheckSpecificParametersForAllLinearModels())
  {
    Cerr("***** iLinearPatlakModel::CheckSpecificParameters -> A problem occurred while checking specific parameters ! " << endl);
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
  \brief This function is used to initialize Patlak parametric images and basis functions
  \return 0 if success, other value otherwise.
*/
int iLinearPatlakModel::InitializeSpecific()
{
  if(m_verbose >=2) Cout("iLinearPatlakModel::InitializeSpecific ..."<< endl);

  // Run generic Initialization for all Linear Models
  if (InitializeSpecificToAllLinearModels())
  {
    Cerr("***** iLinearPatlakModel::InitializeSpecific() -> Error while performing generic initialisations for linear models !" << endl);
    return 1;
  }

  // Initialization with a list of options
  if(m_listOptions != "")
  {
    // Default initialization : Nested-EM
    m_OptimisationMethod = OPTIMISATION_METHOD_NESTEM;
    
    // We expect here the coefficients of the functions for each time point
    // Patlak slope before Patlak intercept
  
    // Allocate memory to recover the elements in one tmp vector
    FLTNB *p_coeffs = new FLTNB[m_nbTimeBF * mp_ID->GetNbTimeFrames()];
  
    // Read them
    if (ReadStringOption(m_listOptions,
                         p_coeffs,
                         m_nbTimeBF * mp_ID->GetNbTimeFrames(),
                         ",",
                         "Patlak model configuration"))
    {
      Cerr("***** iPatlakModel::Initialize() -> Failed to correctly read the list of parameters in command-line options  !" << endl);
      return 1;
    }
  
    // Affect coeffs
    for(int c=0 ; c<m_nbTimeBF * mp_ID->GetNbTimeFrames() ; c++)
    {
      int bf = int(c/mp_ID->GetNbTimeFrames()); // Patlak basis function index
      int fr = int(c%mp_ID->GetNbTimeFrames()); // Frame index
      m2p_nestedModelTimeBasisFunctions[ bf ][ fr ] = p_coeffs[ c ];
    }
  
    // Delete the tmp vector
    delete[] p_coeffs;
  }

  // --- Default Initialization of time basis functions --- //
  for(int v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
  {
    m2p_parametricImages[0][v] = 1.0;
    m2p_parametricImages[1][v] = 1.0;
  }
  
  if (m_AICfileProvided)
  {
    if(m_verbose >=2) Cout("iLinearPatlakModel::Estimating Patlak basis functions ..."<< endl);
    // Calculate the Patlak Basis functions and set the m2p_nestedModelTimeBasisFunctions --------------------------------------/
    uint32_t* a_frameTimeStopInMs = mp_ID->GetFramesTimeStopArray(0);
    uint32_t* a_frameTimeStartInMs = mp_ID->GetFramesTimeStartsArray(0);
    int a_nbTimeFrames = mp_ID->GetNbTimeFrames();

    // Get pointer to interpolated Arterial Input Curve
    HPFLTNB* a_AICIntrpY = mp_ArterialInputCurve ->GetInterpolatedAIC();
    // Calculation of the 'Patlak Time' integral in the discrete level with the trapezoidal method
    // For uniform spacing that is $ \frac{\Delta T}{2} ( f(x_0) + \sum_{k=0}^{N-1} f(x_k) + f(x_N)) $
    HPFLTNB* a_PatlakTAC = new HPFLTNB[(a_frameTimeStopInMs[a_nbTimeFrames-1])];
    HPFLTNB RunningSum = 0 ;
    HPFLTNB halfDeltaT = (0.001/2);
    // First value of integration over the AIC will be the first value of the AIC * DeltaT
    a_PatlakTAC[0] = a_AICIntrpY[0] * 0.001;
    // Then iterate though all the other points using the trapezoidal function for uniform spacing
    for (uint32_t i=1 ;i<(a_frameTimeStopInMs[a_nbTimeFrames - 1]);i++)
    {
      RunningSum += a_AICIntrpY[i-1] ;
      a_PatlakTAC[i] = (a_AICIntrpY[0] +  2 * RunningSum + a_AICIntrpY[i])*halfDeltaT ;
    }

    // Allocate memory for the Patlak Basis functions
    m2p_nestedModelTimeBasisFunctions = new FLTNB* [2];
    // Setting initial values to 0
    for (int tbf=0; tbf<2; tbf++) m2p_nestedModelTimeBasisFunctions[tbf] = new FLTNB [a_nbTimeFrames];

    // Calculate first basis function - Integration of Patlak time over each frame and division by duration-//
    for (int fr = 0; fr < a_nbTimeFrames; fr++)
    {
      RunningSum  = 0 ;
      // iteration over all the datapoints within the frame
      for (uint32_t i = a_frameTimeStartInMs[fr]+1 ;i <(a_frameTimeStopInMs[fr]);i++)
      {
        RunningSum+=a_PatlakTAC[i];
      }
      m2p_nestedModelTimeBasisFunctions[0][fr] = (FLTNB) ((a_PatlakTAC[a_frameTimeStartInMs[fr]] +
                                                2 * RunningSum + a_PatlakTAC[a_frameTimeStopInMs[fr]])*halfDeltaT) ;
      m2p_nestedModelTimeBasisFunctions[0][fr] /=  ((FLTNB)(a_frameTimeStopInMs[fr] - a_frameTimeStartInMs[fr]))/1000;
    }

    // Calculate second basis function - Average of interpolated AIC over each frame and division by duration-//
    for (int fr = 0; fr < a_nbTimeFrames; fr++)
    {
      RunningSum = 0;
      // iteration over all the datapoints within the frame
      for (uint32_t i = a_frameTimeStartInMs[fr]+1; i < (a_frameTimeStopInMs[fr]); i++)
      {
        RunningSum += a_AICIntrpY[i];
      }
      m2p_nestedModelTimeBasisFunctions[1][fr] = (FLTNB) ((a_AICIntrpY[a_frameTimeStartInMs[fr]] +
                                                2 * RunningSum + a_AICIntrpY[a_frameTimeStopInMs[fr]]) * halfDeltaT);
      m2p_nestedModelTimeBasisFunctions[1][fr] /=  ((FLTNB)(a_frameTimeStopInMs[fr] - a_frameTimeStartInMs[fr]))/1000;
    }

    // Clear memory for Patlak time
    if ( a_PatlakTAC ) delete[] a_PatlakTAC;
    // --------------------------------------------------------------------------------------------------------------/

  }

  // Print the basis functions
  ShowBasisFunctions();

  // Normal end
  m_initialized = true;
  return 0;
}
  

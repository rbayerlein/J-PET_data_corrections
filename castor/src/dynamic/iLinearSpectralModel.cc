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
  \brief    Implementation of class iDynamicModelTemplate
*/

#include "iLinearSpectralModel.hh"
#include "cmath"
#include "math.h"



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn iDynamicModelTemplate
  \brief Constructor of iDynamicModelTemplate. Simply set all data members to default values.
*/
iLinearSpectralModel::iLinearSpectralModel() : iLinearModel()
{
  // --- Parameters inherited from vDynamicModel class --- //
  
  m_nbTimeBF = -1; // Number of basis functions in the model un-initialised

  // Negative initialization of specific model values
  m_fast_exp = -1;
  m_slow_exp = -1;
  m_full_trapping_basis_flag = false;
  m_blood_fraction_fasis_flag = false;
  m_spectral_bank_size = -1;
  mp_spectral_bank=NULL;

  // Initialize number of additional basis functions to 0;
  m_additionalBF_size = 0  ;

}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ~iDynamicModelTemplate
  \brief Destructor of iDynamicModelTemplate
*/
iLinearSpectralModel::~iLinearSpectralModel()
{

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
void iLinearSpectralModel::ShowHelpModelSpecific()
{
  cout << "-- This class implements the Spectral Model : " << endl;
  cout << "-- first introduced by J. Cunningham et al. then used in PET reconstruction for temporal regularisation by Andrew Reader et al " << endl;
  cout << "-- It is used to model radiotracer dynamic behaviour with a set of decaying exponential functions ( exp(-β*t) )  " << endl;
  cout << "-- The decaying exponential coefficients β are logarithmically equally spaced within the selected range of values   " << endl;
  cout << "   All decaying exponentials are convolved with the interpolated Arterial Input Curve, before being discretised to the duration of the reconstruction frames  " << endl;
  cout << endl;
  cout << "  The model must be initialized using a ASCII file with the following keywords and information :" << endl;
  cout << "   As this class inherits from the iLinearModel class, the following parameters must be declared inside the couple of the following specific tags: " << endl;
  cout << "   - DYNAMIC FRAMING/ENDDF " << endl;
  cout << "   - The ASCII file must contain the following keywords :" << endl;
  cout << "    'AIC_input_file: path/to/file'  (mandatory) The file containing the sampled Arterial Input Function " << endl;
  cout << "                                         The file must contain the following information in successive lines, separated by ',' " << endl;
  cout << "                                      -> AIC_number_of_points: " << endl;
  cout << "                                      -> AIC_time_points: " << endl;
  cout << "                                      -> AIC_data_points: " << endl;
  cout << "                                      -> AIC_units: 'seconds' or 'minutes'  " << endl;
  cout << "            The following options are required for the specification of the spectral functions. "<< endl;
  cout << "             As in the spectral analysis method, the spectral functions are convolved with the interpolated Input Curve. "<< endl;
  cout << "             By default a function which is the interpolated Input Curve is added in the dataset (equivalent to convolution of the IC with a dirac function)" << endl;
  cout << "             The spectral functions created will be logarithmically spaced within the selected range " << endl;
  cout << "    Number_spectral_functions: x    (mandatory) The number of spectral functions" << endl;
  cout << "    Fastest_rate: x                 (mandatory) The rate (1/min) for the fastest decaying exponential " << endl;
  cout << "    Slowest_rate: x                 (mandatory) The rate (1/min) for the slowest decaying exponential" << endl;
  cout << "    Full_trapping_basis_function:  1 or 0 (optional) Option to say whether we want to include a basis function to model full trapping of tracer (1) or not (0 by default) " << endl;
  cout << "    Blood_fraction_basis_function: 1 or 0 (optional) Option to say whether we want to include a basis function to model the blood fraction of the tracer (1) or not (0 by default) " << endl;
  cout << "            In this implementation and for the blood fraction it is assumed that the whole blood input curve and plasma input curve are the same (no change due to metabolism)  " << endl;
  cout << endl;
  cout << "   'Parametric_image_init: path ' (optional) path to an interfile image to be used as initialization for the parametric images." << endl;
  cout << endl;
  cout << "   -> Optimisation_method : x       (mandatory) optimization method available options: " << endl;
  cout << "                                      x=0: Direct ( Implementation of basis functions side by system matrix in each tomographic iterative loop " << endl;
  cout << "                                      x=1: Nested EM  " << endl;
  cout << "                                      x=2: Iterative non-negative Least-Square " << endl;
  cout << "                                           (C.L. Lawson and R.J. Hanson, Solving Least Squares Problems)" << endl;
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
int iLinearSpectralModel::ReadAndCheckConfigurationFileSpecific()
{
  if(m_verbose >=3) Cout("iLinearSpectralModel::ReadAndCheckConfigurationFileSpecific ..."<< endl);
  
  // ===================================================================
  // Implement here the reading of any options specific to this dynamic model 
  // (i.e : parameters or path to deformation files), through a configuration file
  // The ReadDataASCIIFile() functions could be helpful to recover data from a file
  // (check other dynamicModels for examples)
  // ===================================================================

  // Apply the Generic linear Checks for all Linear Models
  if( ReadAndCheckConfigurationFileSpecificToAllLinearModels()==1)
  {
    Cerr("***** iLinearSpectralModel::ReadAndCheckConfigurationFileSpecific -> Error while trying to read and check generic configuration for linear models " << m_fileOptions << endl);
    return 1;
  }

  // The file will be fully processed in the Initialize() function
  ifstream in_file(m_fileOptions.c_str(), ios::in);

  if(in_file)
  {
    // Number of requested basis/spectral functions
    if( ReadDataASCIIFile(m_fileOptions, "Number_spectral_functions", &m_spectral_bank_size, 1, KEYWORD_MANDATORY) == 1)
    {
      Cerr("***** iLinearModel::ReadAndCheckConfigurationFileSpecific -> Error while trying to read 'Number_spectral_functions' flag in " << m_fileOptions << endl);
      return 1;
    }


    // Is trapping basis function requested ?
    int full_trapping_basis_input =-1 ;
    if( ReadDataASCIIFile(m_fileOptions, "Full_trapping_basis_function", &full_trapping_basis_input, 1, KEYWORD_OPTIONAL) == 1)
    {
      Cerr("***** iLinearModel::ReadAndCheckConfigurationFileSpecific -> Error while trying to read 'Full_trapping_basis_function' flag in " << m_fileOptions << endl);
      return 1;
    }
    if (full_trapping_basis_input>0) m_full_trapping_basis_flag = true;

    // Is blood fraction basis function requested ?
    int blood_fraction_basis_input =-1 ;
    if( ReadDataASCIIFile(m_fileOptions, "Blood_fraction_basis_function", &blood_fraction_basis_input, 1, KEYWORD_OPTIONAL) == 1)
    {
      Cerr("***** iLinearModel::ReadAndCheckConfigurationFileSpecific -> Error while trying to read 'Blood_fraction_basis_function' flag in " << m_fileOptions << endl);
      return 1;
    }
    if (blood_fraction_basis_input>0) m_blood_fraction_fasis_flag = true;

    // Set model parameters and number of basis functions equal to the number of basis functions (including the additional)
    // Add another basis function if the trapping basis flag has been provided
    if (m_full_trapping_basis_flag) m_additionalBF_size++;
    // Add another basis function if the blood fraction flag has been provided
    if (m_blood_fraction_fasis_flag) m_additionalBF_size++;

    m_nbModelParam = m_spectral_bank_size + m_additionalBF_size ;
    m_nbTimeBF = m_spectral_bank_size + m_additionalBF_size ;

    // Parameters for spectral Model -> Lower, Upper decay rate constants
    if( ReadDataASCIIFile(m_fileOptions, "Fastest_rate", &m_fast_exp, 1, KEYWORD_MANDATORY) == 1)
    {
      Cerr("***** iLinearSpectralModel::ReadAndCheckConfigurationFileSpecific -> Error while trying to read 'Fastest_rate' flag in " << m_fileOptions << endl);
      return 1;
    }
    // Make rate units to 1/sec
    m_fast_exp/=60. ;

    if( ReadDataASCIIFile(m_fileOptions, "Slowest_rate", &m_slow_exp, 1, KEYWORD_MANDATORY) == 1)
    {
      Cerr("***** iLinearSpectralModel::ReadAndCheckConfigurationFileSpecific -> Error while trying to read 'Slowest_rate' flag in " << m_fileOptions << endl);
      return 1;
    }
    // Make rate units to 1/sec
    m_slow_exp/=60. ;
  }
  else
  {
    Cerr("***** iLinearPatlakModel::ReadAndCheckConfigurationFileSpecific -> Error while trying to read configuration file at: " << m_fileOptions << endl);
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
  \param a_optionsList : a list of parameters separated by commas
  \brief This function is used to read parameters from a string.
  \return 0 if success, other value otherwise.
*/
int iLinearSpectralModel::ReadAndCheckOptionsList(string a_listOptions)
{
  if(m_verbose >=3) Cout("iLinearSpectralModel::ReadAndCheckConfigurationFileSpecific ..."<< endl);

  // Must be initialized with configuration file
  Cerr("***** iLinearSpectralModel::ReadAndCheckOptionsList -> This model needs to use a file for its initialization. Use help for more info !" << endl);
  return 1;
  
  // ===================================================================
  // Implement here the reading of any options specific to this deformation model,
  // through a list of options separated by commas
  // The ReadStringOption() function could be helpful to parse the list of parameters in an array
  // ===================================================================

  // Just recover the string here, it will be processed in the Initialize() function
  //m_listOptions = a_listOptions;

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
int iLinearSpectralModel::CheckSpecificParameters()
{
  // ===================================================================
  // Implement here checks over parameters which should be read using either
  // ReadAndCheckConfigurationFile() and ReadAndCheckOptionsList() functions
  // ===================================================================

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
  \brief This function is used to initialize the model parametric images and basis functions
  \return 0 if success, other value otherwise.
*/
int iLinearSpectralModel::InitializeSpecific()
{
  if (!m_checked)
  {
    Cerr("***** iLinearSpectralModel::InitializeSpecific() -> Must call CheckParameters functions before Initialize() !" << endl);
    return 1;
  }


  if(m_verbose >=2) Cout("iLinearSpectralModel::InitializeSpecific ..."<< endl);
  // Run generic Initialization for all Linear Models
  if (InitializeSpecificToAllLinearModels())
  {
    Cerr("***** iLinearSpectralModel::InitializeSpecific() -> Error while performing generic initialisations for linear models !" << endl);
    return 1;
  }

  // Assumption of framing set equally to all bed position at the moment !!
  int bed = 0;
  // Get pointer to interpolated Arterial Input Curve
  HPFLTNB* a_AICIntrpY = mp_ArterialInputCurve->GetInterpolatedAIC();
  // Get Downsampled AIF and get it into a_AICIntrpYDown
  mp_ArterialInputCurve->Downsample();
  HPFLTNB* a_AICIntrpYDown = mp_ArterialInputCurve->GetDownsampledAIC();

  if (m_full_trapping_basis_flag)
  {
    // --- Calculation of φ0 function for estimation of the tracer trapping ( if included ) -------------------------- //

    // Calculation of the 'Intagrated Time' integral in the discrete level with the trapezoidal method
    // For uniform spacing that is $ \frac{\Delta T}{2} ( f(x_0) + \sum_{k=0}^{N-1} f(x_k) + f(x_N)) $
    HPFLTNB* a_IntegAICY = new HPFLTNB[(mp_ID->GetFrameTimeStopInMs(bed, mp_ID->GetNbTimeFrames() - 1))];

    HPFLTNB RunningSum = 0;
    HPFLTNB halfDeltaT = (0.001 / 2);
    // First value of integration over the AIC will be the first value of the AIC * DeltaT
    a_IntegAICY[0] = a_AICIntrpY[0] * 0.001;
    // Then iterate though all the other points using the trapezoidal function for uniform spacing
    for (uint32_t i = 1; i < (mp_ID->GetFrameTimeStopInMs(bed, mp_ID->GetNbTimeFrames() - 1)); i++)
    {
      RunningSum += a_AICIntrpY[i - 1];
      a_IntegAICY[i] = (a_AICIntrpY[0] + 2 * RunningSum + a_AICIntrpY[i]) * halfDeltaT;
    }

    // Calculate first basis function - Integration of sampled integral of AIC over each frame and division by duration-//
    for (int fr = 0; fr < mp_ID->GetNbTimeFrames(); fr++)
    {
      RunningSum = 0;
      // iteration over all the datapoints within the frame
      for (uint32_t i = mp_ID->GetFrameTimeStartInMs(bed, fr) + 1; i < (mp_ID->GetFrameTimeStopInMs(bed, fr)); i++)
      {
        RunningSum += a_IntegAICY[i];
      }
      m2p_nestedModelTimeBasisFunctions[0][fr] = (FLTNB) ((a_IntegAICY[mp_ID->GetFrameTimeStartInMs(bed, fr)] +
                                                           2 * RunningSum +
                                                           a_IntegAICY[mp_ID->GetFrameTimeStopInMs(bed, fr)]) *
                                                          halfDeltaT);
      m2p_nestedModelTimeBasisFunctions[0][fr] /=
          ((FLTNB) (mp_ID->GetFrameTimeStopInMs(bed, fr) - mp_ID->GetFrameTimeStartInMs(bed, fr))) / 1000;
    }
    // Clear Integrated TAC from memory
    if (a_IntegAICY) delete [] (a_IntegAICY);
  }

  // --- Calculation of spectral functions φ1 .. φN-1-------------------------------------------------------- //
  // Allocate spectral bank
  HPFLTNB* mp_spectral_bank = new HPFLTNB[m_spectral_bank_size];
  // Calculate the spectral coefficients to be equally spaced in a log scale
  HPFLTNB ln_slow_rate = log(m_slow_exp);
  HPFLTNB ln_fast_rate = log(m_fast_exp);
  HPFLTNB step = (ln_fast_rate - ln_slow_rate) / (HPFLTNB) (m_spectral_bank_size - 1);
  for (int b = 0; b < m_spectral_bank_size; b++) mp_spectral_bank[b] = ln_slow_rate + step * b;
  // Get spectral coefficients back to scale from log scale
  for (int b = 0; b < m_spectral_bank_size; b++) mp_spectral_bank[b] = exp(mp_spectral_bank[b]);

  // --- Sample the exponential over the requested frames duration --------------------------------------------- //
  if(m_verbose >=3) Cout( " iLinearSpectralModel::InitializeSpecific() -> Sampling exponential functions" << endl);

  // Calculated total samples for interval of 0.1 seconds
  int total_samples =
      (int) ((((float) mp_ID->GetFrameTimeStopInMs(bed, mp_ID->GetNbTimeFrames() - 1)) / 1000.) * 10);

  HPFLTNB** spectral_exp = new HPFLTNB*[m_spectral_bank_size];
  for (int ex = 0; ex < m_spectral_bank_size; ex++)
  {
    spectral_exp[ex] = new HPFLTNB[total_samples];
  }

  // Evaluating exponents over total time. Scaling time for 0.1 second intervals
  for (int ex = 0; ex < m_spectral_bank_size; ex++)
  {
    for (int i = 0; i < total_samples; i++)
    {
      spectral_exp[ex][i] = (FLTNB) exp(-0.1 * mp_spectral_bank[ex] * i);
    }
  }
  // --- Get Interpolated AIC and convolve with every spectral function ----------------------------------------- //

  // Allocate memory for the convolved spectral functions.
  HPFLTNB** ConvSpectralFunctions = new HPFLTNB*[m_spectral_bank_size];
  for (int ex = 0; ex < m_spectral_bank_size; ex++)
  {
    ConvSpectralFunctions[ex] = new HPFLTNB[total_samples];
  }

  // For every function perform the convolution ! use of multithreading for each spectral function convolution
  int ex;
  #pragma omp parallel for private(ex) schedule(static, 1)
  for (ex = 0; ex < m_spectral_bank_size; ex++)
  {
    // Double loop for convolution of AICIntrpY with a sampled exponential function
    for (int i = 0; i < total_samples; i++)
    {
      for (int j = i; j < total_samples; j++)
      {
        ConvSpectralFunctions[ex][j] += (a_AICIntrpYDown[i] * spectral_exp[ex][j - i]);
      }
      // Normalise each value with the size of each step
      ConvSpectralFunctions[ex][i] *= 0.1;
      // Check and remove negative values - might occur if multiplication results to lower/higher values than double !!
      if ( ConvSpectralFunctions[ex][i] <0 )
      {
        Cout(" Negative spectral function value detected at " << i << ", with value: "<< ConvSpectralFunctions[ex][i] <<
         "  -> Setting value to zero "<< endl);
        ConvSpectralFunctions[ex][i] = 0 ;
      }
      if (std::isnan(ConvSpectralFunctions[ex][i]))
      {
        Cout(" NaN spectral function value detected at " << i << "--> Setting value to zero "<< endl);
        ConvSpectralFunctions[ex][i] = 0 ;
      }
    }
  }
  // --- Average convolved basis functions into the study frames
  // If full trapping BF set, then start spectral functions index 1 as φ1 ; else from index 0 as φ0.
  int bf_index_start = (m_full_trapping_basis_flag==true) ? 1 : 0 ;
  // Average functions over frames to calculate Basis Functions
  HPFLTNB RunningSum;
  HPFLTNB halfDeltaT = (0.1/2);
  // Loop over all spectral functions
  for (int ex = 0; ex < m_spectral_bank_size; ex++)
  {
    // Loop over frames to get the sum of the values within each frame
    for (int fr = 0; fr < mp_ID->GetNbTimeFrames(); fr++)
    {
      RunningSum  = 0 ;
      // iteration over all the datapoints within the frame
      for (uint32_t i = (mp_ID->GetFrameTimeStartInMs(bed, fr)/100)+1 ;i <(mp_ID->GetFrameTimeStopInMs(bed, fr)/100);i++)
      {
        RunningSum+=ConvSpectralFunctions[ex][i];
      }
      // Calculate integral using the trapezoidal method
      m2p_nestedModelTimeBasisFunctions[ex+bf_index_start][fr] = (FLTNB) ((ConvSpectralFunctions[ex][(mp_ID->GetFrameTimeStartInMs(bed, fr) / 100)] +
                  2 * RunningSum + ConvSpectralFunctions[ex][(mp_ID->GetFrameTimeStopInMs(bed, fr) / 100)]) * halfDeltaT);
      // Average for the frame duration
      m2p_nestedModelTimeBasisFunctions[ex+bf_index_start][fr] /=  (FLTNB)((mp_ID->GetFrameTimeStopInMs(bed, fr)) - (mp_ID->GetFrameTimeStartInMs(bed, fr)))/1000 ;
    }
  }

  // Clear memory from ConvSpectralFunctions
  if (mp_spectral_bank) delete [] (mp_spectral_bank);
  for (ex = 0; ex < m_spectral_bank_size; ex++)
    if (spectral_exp[ex]) delete[] (spectral_exp[ex]);
  if (spectral_exp) delete[] (spectral_exp);

  for (ex = 0; ex < m_spectral_bank_size; ex++)
  {
    if (ConvSpectralFunctions[ex]) delete[] (ConvSpectralFunctions[ex]);
  }
  if (ConvSpectralFunctions) delete[] (ConvSpectralFunctions);

  // If blood fraction BF set, set it to the φN coefficient
  if (m_blood_fraction_fasis_flag)
  {
    // --- Calculation of functions φN for estimation of blood fraction ------------------------------------------- //
    halfDeltaT = (0.001 / 2);
    for (int fr = 0; fr < mp_ID->GetNbTimeFrames(); fr++)
    {
      RunningSum = 0;
      // iteration over all the datapoints within the frame
      for (uint32_t i = (mp_ID->GetFrameTimeStartInMs(bed, fr)) + 1; i < (mp_ID->GetFrameTimeStopInMs(bed, fr)); i++)
      {
        RunningSum += a_AICIntrpY[i];
      }
      // Calculate integral using the trapezoidal method
      m2p_nestedModelTimeBasisFunctions[m_nbTimeBF - 1][fr] = (FLTNB) (
          (a_AICIntrpY[(mp_ID->GetFrameTimeStartInMs(bed, fr))] +
           2 * RunningSum + a_AICIntrpY[(mp_ID->GetFrameTimeStopInMs(bed, fr))]) * halfDeltaT);
      // Average for the frame duration
      m2p_nestedModelTimeBasisFunctions[m_nbTimeBF - 1][fr] /=
          (FLTNB) ((mp_ID->GetFrameTimeStopInMs(bed, fr)) - (mp_ID->GetFrameTimeStartInMs(bed, fr))) / 1000;
    }
  }

  // Print the basis functions
  ShowBasisFunctions();

  // Normal end
  m_initialized = true;
  return 0;
}

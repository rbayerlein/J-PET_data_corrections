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
  \ingroup  optimizer
  \brief    Implementation of class iOptimizerOneStepLate
*/

#include "iOptimizerOneStepLate.hh"
#include "sOutputManager.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iOptimizerOneStepLate::iOptimizerOneStepLate() : vOptimizer()
{
  // ---------------------------
  // Mandatory member parameters
  // ---------------------------

  // Initial value at 1
  m_initialValue = 1.;
  // Only one backward image
  m_nbBackwardImages = 1;
  // OSL accepts penalties, which must have a derivative order of 1 minimum
  m_requiredPenaltyDerivativesOrder = 1;
  // OSL is compatible with listmode and histogram data
  m_listmodeCompatibility = true;
  m_histogramCompatibility = true;
  // OSL is only compatible with emission data
  m_emissionCompatibility = true;
  m_transmissionCompatibility = true;

  // --------------------------
  // Specific member parameters
  // --------------------------

  m_dataSpaceDenominatorThreshold = -1.;
  m_minimumImageUpdateFactor = -1.;
  m_maximumImageUpdateFactor = -1.;
  m4p_firstDerivativePenaltyImage = NULL;
  m_displayWarningFlag = true;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iOptimizerOneStepLate::~iOptimizerOneStepLate()
{
  // Delete the penalty image
  // Note: there is no need to deallocate the images themselves as they are allocate using the
  //       miscellaneous image function from the image space, which automatically deals with
  //       memory deallocations.
  if (m4p_firstDerivativePenaltyImage)
  {
    // Loop over time basis functions
    for (int tbf=0; tbf<mp_ImageDimensionsAndQuantification->GetNbTimeBasisFunctions(); tbf++)
    {
      if (m4p_firstDerivativePenaltyImage[tbf])
      {
        // Loop over respiratory basis functions
        for (int rbf=0; rbf<mp_ImageDimensionsAndQuantification->GetNbRespBasisFunctions(); rbf++)
        {
          if (m4p_firstDerivativePenaltyImage[tbf][rbf])
          {
            free(m4p_firstDerivativePenaltyImage[tbf][rbf]);
          }
        }
        free(m4p_firstDerivativePenaltyImage[tbf]);
      }
    }
    free(m4p_firstDerivativePenaltyImage);
  }
  // If the warning about negative sensitivity+penalty has been printed out, we print it again
  // here at the end to be sure that the user is paying attention to it.
  if (!m_displayWarningFlag)
  {
    Cerr("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl);
    Cerr("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl);
    Cerr("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl);
    Cerr("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl);
    Cerr("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl);
    Cerr("!!!!! iOptimizerOneStepLate::Destructor() -> Several times, the update factor was not positive and was set to 1. !!!!!" << endl);
    Cerr("!!!!! For convenience, it was advertised only the first time. But this may be a sign that the penalty strength   !!!!!" << endl);
    Cerr("!!!!! is too high for this approximative algorithm. It may also be due to voxels close to FOV extremities. Be    !!!!!" << endl);
    Cerr("!!!!! sure to double check your images !                                                                         !!!!!" << endl);
    Cerr("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl);
    Cerr("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl);
    Cerr("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl);
    Cerr("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl);
    Cerr("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl);
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iOptimizerOneStepLate::ShowHelpSpecific()
{
  cout << "This optimizer is the One-Step-Late algorithm from P. J. Green, IEEE TMI, Mar 1990, vol. 9, pp. 84-93." << endl;
  cout << "Subsets can be used as for OSEM. It accepts penalty terms that have a derivative order of at least one." << endl;
  cout << "Without penalty, it is stricly equivalent to the MLEM algorithm." << endl;
  cout << "It is numerically implemented in the multiplicative form (as opposed to the gradient form)." << endl;
  cout << "The following options can be used (in this particular order when provided as a list):" << endl;
  cout << "  initial image value: to set the uniform voxel value for the initial image" << endl;
  cout << "  denominator threshold: to set the threshold of the data space denominator under which the ratio is set to 1" << endl;
  cout << "  minimum image update: to set the minimum of the image update factor under which it stays constant (0 or a negative value" << endl;
  cout << "                        means no minimum thus allowing a 0 update)" << endl;
  cout << "  maximum image update: to set the maximum of the image update factor over which it stays constant (0 or a negative value means" << endl;
  cout << "                        no maximum)" << endl;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerOneStepLate::ReadConfigurationFile(const string& a_configurationFile)
{
  string key_word = "";
  // Read the initial image value option
  key_word = "initial image value";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_initialValue, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iOptimizerOneStepLate::ReadAndCheckConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Read the denominator threshold option
  key_word = "denominator threshold";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_dataSpaceDenominatorThreshold, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iOptimizerOneStepLate::ReadAndCheckConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Read the minimum image update option
  key_word = "minimum image update";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_minimumImageUpdateFactor, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iOptimizerOneStepLate::ReadAndCheckConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Read the maximum image update option
  key_word = "maximum image update";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_maximumImageUpdateFactor, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iOptimizerOneStepLate::ReadAndCheckConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerOneStepLate::ReadOptionsList(const string& a_optionsList)
{
  // There are 4 floating point variables as options
  const int nb_options = 4;
  FLTNB options[nb_options];
  // Read them
  if (ReadStringOption(a_optionsList, options, nb_options, ",", "OneStepLate configuration"))
  {
    Cerr("***** iOptimizerOneStepLate::ReadAndCheckConfigurationFile() -> Failed to correctly read the list of options !" << endl);
    return 1;
  }
  // Affect options
  m_initialValue = options[0];
  m_dataSpaceDenominatorThreshold = options[1];
  m_minimumImageUpdateFactor = options[2];
  m_maximumImageUpdateFactor = options[3];
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerOneStepLate::CheckSpecificParameters()
{
  // Check that initial image value is strictly positive
  if (m_initialValue<=0.)
  {
    Cerr("***** iOptimizerOneStepLate->Initialize() -> Provided initial image value (" << m_initialValue << ") must be strictly positive !" << endl);
    return 1;
  }
  // Check that denominator threshold value is strictly positive
  if (m_dataSpaceDenominatorThreshold<=0.)
  {
    Cerr("***** iOptimizerOneStepLate->Initialize() -> Provided data space denominator threshold (" << m_dataSpaceDenominatorThreshold << ") must be strictly positive !" << endl);
    return 1;
  }
  // Check that maximum image update factor is higher than the minimum
  if (m_minimumImageUpdateFactor>0. && m_maximumImageUpdateFactor>0. && m_maximumImageUpdateFactor<m_minimumImageUpdateFactor)
  {
    Cerr("***** iOptimizerOneStepLate->Initialize() -> Provided minimum/maximum (" << m_minimumImageUpdateFactor << "/" << m_maximumImageUpdateFactor << " are inconsistent !" << endl);
    return 1;
  }
  // Cannot deal with list-mode transmission data
  if (m_dataSpec==SPEC_TRANSMISSION && m_dataMode==MODE_LIST)
  {
    Cerr("***** iOptimizerMLEM->CheckSpecificParameters() -> Cannot reconstruct list-mode transmission data !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerOneStepLate::InitializeSpecific()
{
  // Allocate and create the penalty image
  m4p_firstDerivativePenaltyImage = (FLTNB****)malloc(mp_ImageDimensionsAndQuantification->GetNbTimeBasisFunctions()*sizeof(FLTNB***));
  // Loop over time basis functions
  for (int tbf=0; tbf<mp_ImageDimensionsAndQuantification->GetNbTimeBasisFunctions(); tbf++)
  { 
    m4p_firstDerivativePenaltyImage[tbf] = (FLTNB***)malloc(mp_ImageDimensionsAndQuantification->GetNbRespBasisFunctions()*sizeof(FLTNB**));
    // Loop over respiratory basis functions
    for (int rbf=0; rbf<mp_ImageDimensionsAndQuantification->GetNbRespBasisFunctions(); rbf++)
    {
      m4p_firstDerivativePenaltyImage[tbf][rbf] = (FLTNB**)malloc(mp_ImageDimensionsAndQuantification->GetNbRespBasisFunctions()*sizeof(FLTNB*));
      // Loop over cardiac basis functions
      for (int cbf=0; cbf<mp_ImageDimensionsAndQuantification->GetNbCardBasisFunctions(); cbf++)
      {
        // Get a pointer to a newly allocated image
        m4p_firstDerivativePenaltyImage[tbf][rbf][cbf] = mp_ImageSpace -> AllocateMiscellaneousImage();
        // Initialize to 0, in case the penalty is not used
        for (INTNB v=0; v<mp_ImageDimensionsAndQuantification->GetNbVoxXYZ(); v++) m4p_firstDerivativePenaltyImage[tbf][rbf][cbf][v] = 0.;
      }
    }
  }
  // Force the display of the OSL warning (when sensitivity+penalty is negative or null in the image update step)
  m_displayWarningFlag = true;
  // Verbose
  if (m_verbose>=2)
  {
    Cout("iOptimizerOneStepLate::Initialize() -> Use the OneStepLate optimizer" << endl);
    if (m_verbose>=3)
    {
      Cout("  --> Initial image value: " << m_initialValue << endl);
      Cout("  --> Data space denominator threshold: " << m_dataSpaceDenominatorThreshold << endl);
      if (m_minimumImageUpdateFactor>0.) Cout("  --> Minimum image update factor: " << m_minimumImageUpdateFactor << endl);
      else Cerr("!!!!! The minimum update value is not set, if using subsets, voxels could be trapped in 0 value causing some negative bias !" << endl);
      if (m_maximumImageUpdateFactor>0.) Cout("  --> Maximum image update factor: " << m_maximumImageUpdateFactor << endl);
    }
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerOneStepLate::PreImageUpdateSpecificStep()
{
  // ==========================================================================================
  // If no penalty, then exit (the penalty image term has been initialized to 0)
  if (mp_Penalty==NULL) return 0;
  // Set the number of threads
  #ifdef CASTOR_OMP
  omp_set_num_threads(mp_ImageDimensionsAndQuantification->GetNbThreadsForImageComputation());
  #endif
  // Verbose
  if (m_verbose>=1) Cout("iOptimizerOneStepLate::PreImageUpdateSpecificStep() -> Compute penalty term" << endl);
  // ==========================================================================================
  // Global precomputation step if needed by the penalty
  if (mp_Penalty->GlobalPreProcessingStep())
  {
    Cerr("***** iOptimizerOneStepLate::PreImageUpdateSpecificStep() -> A problem occurred while computing the penalty pre-processing step !" << endl);
    return 1;
  }
  // ==========================================================================================
  // Loop over time basis functions
  for (int tbf=0; tbf<mp_ImageDimensionsAndQuantification->GetNbTimeBasisFunctions(); tbf++)
  {
    // Loop over respiratory basis functions
    for (int rbf=0; rbf<mp_ImageDimensionsAndQuantification->GetNbRespBasisFunctions(); rbf++)
    {
      // Loop over cardiac basis functions
      for (int cbf=0; cbf<mp_ImageDimensionsAndQuantification->GetNbCardBasisFunctions(); cbf++)
      {
        // In order to detect problems in the multi-threaded loop
        bool problem = false;
        // Voxel index
        INTNB v;
        // Multi-threading over voxels
        #pragma omp parallel for private(v) schedule(guided)
        for (v=0; v<mp_ImageDimensionsAndQuantification->GetNbVoxXYZ(); v++)
        {
          // Get the thread index
          int th = 0;
          #ifdef CASTOR_OMP
          th = omp_get_thread_num(); 
          #endif
          // Local precomputation step if needed by the penalty
          if (mp_Penalty->LocalPreProcessingStep(tbf,rbf,cbf,v,th))
          {
            Cerr("***** iOptimizerOneStepLate::PreImageUpdateSpecificStep() -> A problem occurred while computing the penalty local pre-processing step for voxel " << v << " !" << endl);
            problem = true;
          }
          // Compute the penalty term at first order
          m4p_firstDerivativePenaltyImage[tbf][rbf][cbf][v] = mp_Penalty->ComputeFirstDerivative(tbf,rbf,cbf,v,th);
        }
        // Check for problems
        if (problem)
        {
          Cerr("***** iOptimizerOneStepLate::PreImageUpdateSpecificStep() -> A problem occurred inside the multi-threaded loop, stop now !" << endl);
          return 1;
        }
      }
    }
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
int iOptimizerOneStepLate::SensitivitySpecificOperations( FLTNB a_data, FLTNB a_forwardModel, FLTNB* ap_weight,
                                                          FLTNB a_multiplicativeCorrections, FLTNB a_additiveCorrections, FLTNB a_blankValue,
                                                          FLTNB a_quantificationFactor, oProjectionLine* ap_Line )
{
  // Line weight here is simply 1
  *ap_weight = 1.;
  // That's all
  return 0;
}
// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerOneStepLate::DataSpaceSpecificOperations( FLTNB a_data, FLTNB a_forwardModel, FLTNB* ap_backwardValues,
                                                        FLTNB a_multiplicativeCorrections, FLTNB a_additiveCorrections, FLTNB a_blankValue,
                                                        FLTNB a_quantificationFactor, oProjectionLine* ap_Line )
{
  // Case for emission tomography
  if (m_dataSpec==SPEC_EMISSION)
  {
    // Truncate data to 0 if negative
    if (a_data<0.) a_data = 0.;
    // If the foward model is too close to zero, then ignore this data (set to 1 and not 0 because this line is taken into account in the sensitivity)
    if (a_forwardModel>m_dataSpaceDenominatorThreshold) *ap_backwardValues = a_data / a_forwardModel;
    else *ap_backwardValues = 1.;
  }
  // Case for transmission tomography (equivalent of log-converted MLEM from Nuyts et al 1998)
  else if (m_dataSpec==SPEC_TRANSMISSION)
  {
    // Subtract scatters
    a_data -= a_additiveCorrections;
    a_forwardModel -= a_additiveCorrections;
    // Safely ignore this data (set backward value to 1 and return) in 2 different cases:
    // - if data or model is inferior to 1
    // - if data or model is higher than blank value
    if (a_data<1. || a_forwardModel<1. || a_data>a_blankValue || a_forwardModel>a_blankValue)
    {
      *ap_backwardValues = 1.;
      return 0;
    }
    // Log-convert the data and the model
    a_data = log(a_blankValue/a_data);
    a_forwardModel = log(a_blankValue/a_forwardModel);
    // If the foward model is to close to zero, then ignore this data (set to 1 and not 0 because this line is taken into account in the sensitivity)
    if (a_forwardModel>m_dataSpaceDenominatorThreshold) *ap_backwardValues = a_data / a_forwardModel;
    else *ap_backwardValues = 1.;
  }
  // That's all
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerOneStepLate::ImageSpaceSpecificOperations( FLTNB a_currentImageValue, FLTNB* ap_newImageValue,
                                                         FLTNB a_sensitivity, FLTNB* ap_correctionValues,
                                                         INTNB a_voxel, int a_tbf, int a_rbf, int a_cbf )
{
  // The update factor
  FLTNB image_update_factor = 0.;
  // Compute the sensitivity + penalty term (the penalty term is divided by the current number of subsets, for balance)
  FLTNB sensitivity_plus_penalty = a_sensitivity + m4p_firstDerivativePenaltyImage[a_tbf][a_rbf][a_cbf][a_voxel]/(FLTNB)mp_nbSubsets[m_currentIteration];
  // If negative update (which can be caused by the limitation of the OSL algorithm itself)
  if (sensitivity_plus_penalty <= 0.)
  {
    // Do not update the image, so set the update factor to 1
    image_update_factor = 1.;
    if (m_displayWarningFlag)
    {
      Cerr("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl);
      Cerr("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl);
      Cerr("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl);
      Cerr("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl);
      Cerr("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl);
      Cerr("!!!!! iOptimizerOneStepLate::ImageSpaceSpecificOperations() -> The update factor was not positive ! It has been set to 1.  !!!!!" << endl);
      Cerr("!!!!! This may be a sign that the penalty strength is too high, but it may only be due to voxels close to FOV extremities. !!!!!" << endl);
      Cerr("!!!!! For convenience, this warning will appear only for this update. Be careful about your results !                      !!!!!" << endl);
      Cerr("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl);
      Cerr("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl);
      Cerr("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl);
      Cerr("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl);
      Cerr("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl);
      // Display this warning only once, otherwise, it can be invasive (many voxels next to FOV extremitites can be affected, even with reasonable beta values)
      m_displayWarningFlag = false;
    }
  }
  // Otherwise, normal update
  else
  {
    // Compute image update factor
    image_update_factor = *ap_correctionValues / sensitivity_plus_penalty;
    // Apply minimum image update factor
    if ( m_minimumImageUpdateFactor > 0. && image_update_factor < m_minimumImageUpdateFactor ) image_update_factor = m_minimumImageUpdateFactor;
    // Apply maximum image update factor
    if ( m_maximumImageUpdateFactor > 0. && image_update_factor > m_maximumImageUpdateFactor ) image_update_factor = m_maximumImageUpdateFactor;
  }
  // Update image
  *ap_newImageValue = a_currentImageValue * image_update_factor;
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

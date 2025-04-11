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
  \brief    Implementation of class iOptimizerPenalizedPreconditionedGradientML
*/

#include "iOptimizerPenalizedPreconditionedGradientML.hh"
#include "sOutputManager.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iOptimizerPenalizedPreconditionedGradientML::iOptimizerPenalizedPreconditionedGradientML() : vOptimizer()
{
  // ---------------------------
  // Mandatory member parameters
  // ---------------------------

  // Initial value at 1
  m_initialValue = 1.;
  // Only one backward image
  m_nbBackwardImages = 1;
  // This algorithm accepts penalty and requires the first and second derivatives
  m_requiredPenaltyDerivativesOrder = 2;
  // PreconditionedGradientMAP is only compatible with histogram data
  m_listmodeCompatibility = false;
  m_histogramCompatibility = true;
  // PreconditionedGradientMAP is only compatible with emission data
  m_emissionCompatibility = true;
  m_transmissionCompatibility = false;

  // --------------------------
  // Specific member parameters
  // --------------------------

  m_dataSpaceDenominatorThreshold = -1.;
  m_minimumImageUpdateFactor = -1.;
  m_maximumImageUpdateFactor = -1.;
  m4p_firstDerivativePenaltyImage = NULL;
  m4p_secondDerivativePenaltyImage = NULL;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iOptimizerPenalizedPreconditionedGradientML::~iOptimizerPenalizedPreconditionedGradientML()
{
  // Note: there is no need to deallocate the images themselves as they are allocate using the
  //       miscellaneous image function from the image space, which automatically deals with
  //       memory deallocations.
  // Delete the first order derivative penalty image
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
  // Delete the second order derivative penalty image
  if (m4p_secondDerivativePenaltyImage)
  {
    // Loop over time basis functions
    for (int tbf=0; tbf<mp_ImageDimensionsAndQuantification->GetNbTimeBasisFunctions(); tbf++)
    {
      if (m4p_secondDerivativePenaltyImage[tbf])
      {
        // Loop over respiratory basis functions
        for (int rbf=0; rbf<mp_ImageDimensionsAndQuantification->GetNbRespBasisFunctions(); rbf++)
        {
          if (m4p_secondDerivativePenaltyImage[tbf][rbf])
          {
            free(m4p_secondDerivativePenaltyImage[tbf][rbf]);
          }
        }
        free(m4p_secondDerivativePenaltyImage[tbf]);
      }
    }
    free(m4p_secondDerivativePenaltyImage);
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iOptimizerPenalizedPreconditionedGradientML::ShowHelpSpecific()
{
  cout << "This optimizer is the Penalized Preconditioned Gradient algorithm from J. Nuyts et al, IEEE TNS, Feb 2002, vol. 49, pp. 56-60." << endl;
  cout << "As usually described by its inventor, it is a heuristic but effective gradient ascent algorithm" << endl;
  cout << "for penalized maximum-likelihood reconstruction. It addresses the shortcoming of One-Step-Late when large" << endl;
  cout << "penalty strengths can create numerical problems. Penalty terms must have a derivative order of at least two." << endl;
  cout << "Subsets can be used as for OSEM. Without penalty, it is equivalent to the gradient ascent form of the MLEM algorithm." << endl;
  cout << "Based on likelihood gradient and penalty, a multiplicative update factor is computed and its range is limited by provided parameters." << endl;
  cout << "Thus, negative values cannot occur and voxels cannot be trapped into 0 values, providing a first positive estimate." << endl;
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

int iOptimizerPenalizedPreconditionedGradientML::ReadConfigurationFile(const string& a_configurationFile)
{
  string key_word = "";
  string buffer = "";
  // Read the initial image value option
  key_word = "initial image value";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_initialValue, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iOptimizerPenalizedPreconditionedGradientML::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Read the denominator threshold option
  key_word = "denominator threshold";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_dataSpaceDenominatorThreshold, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iOptimizerPenalizedPreconditionedGradientML::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Read the minimum image update option
  key_word = "minimum image update";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_minimumImageUpdateFactor, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iOptimizerPenalizedPreconditionedGradientML::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Read the maximum image update option
  key_word = "maximum image update";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_maximumImageUpdateFactor, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iOptimizerPenalizedPreconditionedGradientML::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerPenalizedPreconditionedGradientML::ReadOptionsList(const string& a_optionsList)
{
  // There are 4 floating point variables as options
  const int nb_options = 4;
  FLTNB options[nb_options];
  // Read them
  if (ReadStringOption(a_optionsList, options, nb_options, ",", "PenalizedPreconditionedGradient configuration"))
  {
    Cerr("***** iOptimizerPenalizedPreconditionedGradientML::ReadOptionsList() -> Failed to correctly read the list of options !" << endl);
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

int iOptimizerPenalizedPreconditionedGradientML::CheckSpecificParameters()
{
  // Check that initial image value is strictly positive
  if (m_initialValue<=0.)
  {
    Cerr("***** iOptimizerPenalizedPreconditionedGradientML->CheckSpecificParameters() -> Provided initial image value (" << m_initialValue << ") must be strictly positive !" << endl);
    return 1;
  }
  // Check that denominator threshold value is strictly positive
  if (m_dataSpaceDenominatorThreshold<=0.)
  {
    Cerr("***** iOptimizerPenalizedPreconditionedGradientML->CheckSpecificParameters() -> Provided data space denominator threshold (" << m_dataSpaceDenominatorThreshold << ") must be strictly positive !" << endl);
    return 1;
  }
  // Check that maximum image update factor is higher than the minimum
  if (m_minimumImageUpdateFactor>0. && m_maximumImageUpdateFactor>0. && m_maximumImageUpdateFactor<m_minimumImageUpdateFactor)
  {
    Cerr("***** iOptimizerPenalizedPreconditionedGradientML->CheckSpecificParameters() -> Provided minimum/maximum (" << m_minimumImageUpdateFactor << "/" << m_maximumImageUpdateFactor << " are inconsistent !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerPenalizedPreconditionedGradientML::InitializeSpecific()
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
  // Allocate and create the second order derivative penalty image
  m4p_secondDerivativePenaltyImage = (FLTNB****)malloc(mp_ImageDimensionsAndQuantification->GetNbTimeBasisFunctions()*sizeof(FLTNB***));
  // Loop over time basis functions
  for (int tbf=0; tbf<mp_ImageDimensionsAndQuantification->GetNbTimeBasisFunctions(); tbf++)
  { 
    m4p_secondDerivativePenaltyImage[tbf] = (FLTNB***)malloc(mp_ImageDimensionsAndQuantification->GetNbRespBasisFunctions()*sizeof(FLTNB**));
    // Loop over respiratory basis functions
    for (int rbf=0; rbf<mp_ImageDimensionsAndQuantification->GetNbRespBasisFunctions(); rbf++)
    {
      m4p_secondDerivativePenaltyImage[tbf][rbf] = (FLTNB**)malloc(mp_ImageDimensionsAndQuantification->GetNbRespBasisFunctions()*sizeof(FLTNB*));
      // Loop over cardiac basis functions
      for (int cbf=0; cbf<mp_ImageDimensionsAndQuantification->GetNbCardBasisFunctions(); cbf++)
      {
        // Get a pointer to a newly allocated image
        m4p_secondDerivativePenaltyImage[tbf][rbf][cbf] = mp_ImageSpace -> AllocateMiscellaneousImage();
        // Initialize to 0, in case the penalty is not used
        for (INTNB v=0; v<mp_ImageDimensionsAndQuantification->GetNbVoxXYZ(); v++) m4p_secondDerivativePenaltyImage[tbf][rbf][cbf][v] = 0.;
      }
    }
  }  
  // Verbose
  if (m_verbose>=2)
  {
    Cout("iOptimizerPenalizedPreconditionedGradientML::InitializeSpecific() -> Use the Penalized Preconditioned Gradient optimizer" << endl);
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

int iOptimizerPenalizedPreconditionedGradientML::PreImageUpdateSpecificStep()
{
  // ==========================================================================================
  // If no penalty, then exit (the penalty image term has been initialized to 0)
  if (mp_Penalty==NULL) return 0;
  // Set the number of threads
  #ifdef CASTOR_OMP
  omp_set_num_threads(mp_ImageDimensionsAndQuantification->GetNbThreadsForImageComputation());
  #endif
  // Verbose
  if (m_verbose>=1) Cout("iOptimizerPenalizedPreconditionedGradientML::PreImageUpdateSpecificStep() -> Compute penalty term" << endl);
  // ==========================================================================================
  // Global precomputation step if needed by the penalty
  if (mp_Penalty->GlobalPreProcessingStep())
  {
    Cerr("***** iOptimizerPenalizedPreconditionedGradientML::PreImageUpdateSpecificStep() -> A problem occurred while computing the penalty pre-processing step !" << endl);
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
            Cerr("***** iOptimizerPenalizedPreconditionedGradientML::PreImageUpdateSpecificStep() -> A problem occurred while computing the penalty local pre-processing step for voxel " << v << " !" << endl);
            problem = true;
          }
          // Compute first and second derivative order penalty terms
          m4p_firstDerivativePenaltyImage[tbf][rbf][cbf][v] = mp_Penalty->ComputeFirstDerivative(tbf,rbf,cbf,v,th);
          m4p_secondDerivativePenaltyImage[tbf][rbf][cbf][v] = mp_Penalty->ComputeSecondDerivative(tbf,rbf,cbf,v,th);
        }
        // Check for problems
        if (problem)
        {
          Cerr("***** iOptimizerPenalizedPreconditionedGradientML::PreImageUpdateSpecificStep() -> A problem occurred inside the multi-threaded loop, stop now !" << endl);
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

int iOptimizerPenalizedPreconditionedGradientML::SensitivitySpecificOperations( FLTNB a_data, FLTNB a_forwardModel, FLTNB* ap_weight,
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

int iOptimizerPenalizedPreconditionedGradientML::DataSpaceSpecificOperations( FLTNB a_data, FLTNB a_forwardModel, FLTNB* ap_backwardValues,
                                                                              FLTNB a_multiplicativeCorrections, FLTNB a_additiveCorrections, FLTNB a_blankValue,
                                                                              FLTNB a_quantificationFactor, oProjectionLine* ap_Line )
{
  // Truncate data to 0 if negative
  if (a_data<0.) a_data = 0.;
  // Compute numerator
  FLTNB numerator = a_data - a_forwardModel;

  // If the foward model is too close to zero, then ignore this data (set to 0 as the update is additive)
  if (a_forwardModel>m_dataSpaceDenominatorThreshold) *ap_backwardValues = numerator / a_forwardModel;
  else *ap_backwardValues = 0.;
/*
  // Compute denominator that will be strictly positive
  FLTNB denominator = max(a_forwardModel,m_dataSpaceDenominatorThreshold);
  // Update backward values
  *ap_backwardValues = numerator / denominator;
*/
  // That's all
  return 0;

}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerPenalizedPreconditionedGradientML::ImageSpaceSpecificOperations( FLTNB a_currentImageValue, FLTNB* ap_newImageValue,
                                                                               FLTNB a_sensitivity, FLTNB* ap_correctionValues,
                                                                               INTNB a_voxel, int a_tbf, int a_rbf, int a_cbf )
{
  // Note 1: the penalty terms are divided by the current number of subsets, for balance.
  // Note 2: we did not deal with negative values as suggest by original Johan's paper, because it traps voxels
  //         into 0 value. We rather compute a multiplicative update factor and limit its range as for MLEM.
  // Compute numerator
  FLTNB numerator = *ap_correctionValues - m4p_firstDerivativePenaltyImage[a_tbf][a_rbf][a_cbf][a_voxel] / (FLTNB)mp_nbSubsets[m_currentIteration];
  // Compute denominator
  FLTNB denominator = a_sensitivity + a_currentImageValue * m4p_secondDerivativePenaltyImage[a_tbf][a_rbf][a_cbf][a_voxel] / (FLTNB)mp_nbSubsets[m_currentIteration];
  // Compute multiplicative image update factor
  FLTNB image_update_factor = 1. + numerator / denominator;
  // Apply minimum image update factor
  if ( m_minimumImageUpdateFactor > 0. && image_update_factor < m_minimumImageUpdateFactor ) image_update_factor = m_minimumImageUpdateFactor;
  // Apply maximum image update factor
  if ( m_maximumImageUpdateFactor > 0. && image_update_factor > m_maximumImageUpdateFactor ) image_update_factor = m_maximumImageUpdateFactor;
  // Update image
  *ap_newImageValue = a_currentImageValue * image_update_factor;
  // Check if it is a number, if not, keep the value unchanged
  if (!isfinite(*ap_newImageValue)) *ap_newImageValue = a_currentImageValue;
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

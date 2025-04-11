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
  \brief    Implementation of class iOptimizerModifiedEMForMRFPenalizedML
*/

#include "iOptimizerModifiedEMForMRFPenalizedML.hh"
#include "sOutputManager.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iOptimizerModifiedEMForMRFPenalizedML::iOptimizerModifiedEMForMRFPenalizedML() : vOptimizer()
{
  // ---------------------------
  // Mandatory member parameters
  // ---------------------------

  // Initial value at 1
  m_initialValue = 1.;
  // Only one backward image
  m_nbBackwardImages = 1;
  // OTPML accepts penalties, which must have a derivative order of 1 minimum
  m_requiredPenaltyDerivativesOrder = 1;
  // OTPML is compatible with listmode and histogram data
  m_listmodeCompatibility = true;
  m_histogramCompatibility = true;
  // OTPML is only compatible with emission data
  m_emissionCompatibility = true;
  m_transmissionCompatibility = false;

  // --------------------------
  // Specific member parameters
  // --------------------------

  m_dataSpaceDenominatorThreshold = -1.;
  m_minimumImageUpdateFactor = -1.;
  m_maximumImageUpdateFactor = -1.;
  m4p_smoothedImage = NULL;
  m4p_sumOfWeightsImage = NULL;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iOptimizerModifiedEMForMRFPenalizedML::~iOptimizerModifiedEMForMRFPenalizedML()
{
  // Delete the penalty images
  // Note: there is no need to deallocate the images themselves as they are allocate using the
  //       miscellaneous image function from the image space, which automatically deals with
  //       memory deallocations.
  if (m4p_smoothedImage)
  {
    // Loop over time basis functions
    for (int tbf=0; tbf<mp_ImageDimensionsAndQuantification->GetNbTimeBasisFunctions(); tbf++)
    {
      if (m4p_smoothedImage[tbf])
      {
        // Loop over respiratory basis functions
        for (int rbf=0; rbf<mp_ImageDimensionsAndQuantification->GetNbRespBasisFunctions(); rbf++)
        {
          if (m4p_smoothedImage[tbf][rbf])
          {
            free(m4p_smoothedImage[tbf][rbf]);
          }
        }
        free(m4p_smoothedImage[tbf]);
      }
    }
    free(m4p_smoothedImage);
  }
  if (m4p_sumOfWeightsImage)
  {
    // Loop over time basis functions
    for (int tbf=0; tbf<mp_ImageDimensionsAndQuantification->GetNbTimeBasisFunctions(); tbf++)
    {
      if (m4p_sumOfWeightsImage[tbf])
      {
        // Loop over respiratory basis functions
        for (int rbf=0; rbf<mp_ImageDimensionsAndQuantification->GetNbRespBasisFunctions(); rbf++)
        {
          if (m4p_sumOfWeightsImage[tbf][rbf])
          {
            free(m4p_sumOfWeightsImage[tbf][rbf]);
          }
        }
        free(m4p_sumOfWeightsImage[tbf]);
      }
    }
    free(m4p_sumOfWeightsImage);
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iOptimizerModifiedEMForMRFPenalizedML::ShowHelpSpecific()
{
  cout << "This optimizer is based on the algorithm from A. De Pierro, IEEE TMI, vol. 14, pp. 132-137, 1995." << endl;
  cout << "This algorithm uses optimization transfer techniques to derive an exact and convergent algorithm" << endl;
  cout << "for maximum likelihood reconstruction including a MRF penalty with different potential functions." << endl;
  cout << "The algorithm is convergent and is numerically robust to high penalty strength." << endl;
  cout << "It is stricly equivalent to MLEM without penalty, but can be unstable with extremly low penalty strength." << endl;
  cout << "Currently, it only implements the quadratic penalty." << endl;
  cout << "To be used, a MRF penalty still need to be defined accordingly (at least to define the neighborhood)." << endl;
  cout << "Subsets can be used as for OSEM, without proof of convergence though." << endl;
  cout << "The algorithm is compatible with list-mode or histogram data." << endl;
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

int iOptimizerModifiedEMForMRFPenalizedML::ReadConfigurationFile(const string& a_configurationFile)
{
  string key_word = "";
  // Read the initial image value option
  key_word = "initial image value";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_initialValue, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iOptimizerModifiedEMForMRFPenalizedML::ReadAndCheckConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Read the denominator threshold option
  key_word = "denominator threshold";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_dataSpaceDenominatorThreshold, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iOptimizerModifiedEMForMRFPenalizedML::ReadAndCheckConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Read the minimum image update option
  key_word = "minimum image update";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_minimumImageUpdateFactor, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iOptimizerModifiedEMForMRFPenalizedML::ReadAndCheckConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Read the maximum image update option
  key_word = "maximum image update";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_maximumImageUpdateFactor, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iOptimizerModifiedEMForMRFPenalizedML::ReadAndCheckConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerModifiedEMForMRFPenalizedML::ReadOptionsList(const string& a_optionsList)
{
  // There are 4 floating point variables as options
  const int nb_options = 4;
  FLTNB options[nb_options];
  // Read them
  if (ReadStringOption(a_optionsList, options, nb_options, ",", "Penalized ML with optimization transfer configuration"))
  {
    Cerr("***** iOptimizerModifiedEMForMRFPenalizedML::ReadAndCheckConfigurationFile() -> Failed to correctly read the list of options !" << endl);
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

int iOptimizerModifiedEMForMRFPenalizedML::CheckSpecificParameters()
{
  // Check that initial image value is strictly positive
  if (m_initialValue<=0.)
  {
    Cerr("***** iOptimizerModifiedEMForMRFPenalizedML->Initialize() -> Provided initial image value (" << m_initialValue << ") must be strictly positive !" << endl);
    return 1;
  }
  // Check that denominator threshold value is strictly positive
  if (m_dataSpaceDenominatorThreshold<=0.)
  {
    Cerr("***** iOptimizerModifiedEMForMRFPenalizedML->Initialize() -> Provided data space denominator threshold (" << m_dataSpaceDenominatorThreshold << ") must be strictly positive !" << endl);
    return 1;
  }
  // Check that maximum image update factor is higher than the minimum
  if (m_minimumImageUpdateFactor>0. && m_maximumImageUpdateFactor>0. && m_maximumImageUpdateFactor<m_minimumImageUpdateFactor)
  {
    Cerr("***** iOptimizerModifiedEMForMRFPenalizedML->Initialize() -> Provided minimum/maximum (" << m_minimumImageUpdateFactor << "/" << m_maximumImageUpdateFactor << " are inconsistent !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerModifiedEMForMRFPenalizedML::InitializeSpecific()
{
  // -------------------------------------------------------------------
  // First, check if a Penalty has been initialized
  if ( mp_Penalty != NULL )
  {
    // Check that Penalty is MRF
    if ( mp_Penalty->GetPenaltyID() != "MRF" )
    {
      Cerr("***** iOptimizerModifiedEMForMRFPenalizedML->Initialize() -> This optimizer is only compatible with the Markov Random Field penalty (MRF) !" << endl);
      return 1;
    }
    // If Penalty is MRF, now check potential function is quadratic
    else if ((dynamic_cast<iPenaltyMarkovRandomField*>(mp_Penalty))->GetPotentialType() != MRF_POTENTIAL_QUADRATIC )
    {
      Cerr("***** iOptimizerModifiedEMForMRFPenalizedML->Initialize() -> This optimizer is only compatible with Quadratic potential function for the penalty !" << endl);
      Cerr("                                                                     Please use 'potential function: quadratic' in the penalty initialization file!" << endl);
      return 1;
    }
  }

  // -------------------------------------------------------------------  
  // Allocate and create the smoothed image and the sum_of_weights image
  m4p_smoothedImage = (FLTNB****)malloc(mp_ImageDimensionsAndQuantification->GetNbTimeBasisFunctions()*sizeof(FLTNB***));
  m4p_sumOfWeightsImage = (FLTNB****)malloc(mp_ImageDimensionsAndQuantification->GetNbTimeBasisFunctions()*sizeof(FLTNB***));
  // Loop over time basis functions
  for (int tbf=0; tbf<mp_ImageDimensionsAndQuantification->GetNbTimeBasisFunctions(); tbf++)
  {
    m4p_smoothedImage[tbf] = (FLTNB***)malloc(mp_ImageDimensionsAndQuantification->GetNbRespBasisFunctions()*sizeof(FLTNB**));
    m4p_sumOfWeightsImage[tbf] = (FLTNB***)malloc(mp_ImageDimensionsAndQuantification->GetNbRespBasisFunctions()*sizeof(FLTNB**));
    // Loop over respiratory basis functions
    for (int rbf=0; rbf<mp_ImageDimensionsAndQuantification->GetNbRespBasisFunctions(); rbf++)
    {
      m4p_smoothedImage[tbf][rbf] = (FLTNB**)malloc(mp_ImageDimensionsAndQuantification->GetNbRespBasisFunctions()*sizeof(FLTNB*));
      m4p_sumOfWeightsImage[tbf][rbf] = (FLTNB**)malloc(mp_ImageDimensionsAndQuantification->GetNbRespBasisFunctions()*sizeof(FLTNB*));
      // Loop over cardiac basis functions
      for (int cbf=0; cbf<mp_ImageDimensionsAndQuantification->GetNbCardBasisFunctions(); cbf++)
      {
        // Get a pointer to a newly allocated image
        m4p_smoothedImage[tbf][rbf][cbf] = mp_ImageSpace->AllocateMiscellaneousImage();
        m4p_sumOfWeightsImage[tbf][rbf][cbf] = mp_ImageSpace->AllocateMiscellaneousImage();
        // Initialize to 0, in case the penalty is not used
        for (INTNB v=0; v<mp_ImageDimensionsAndQuantification->GetNbVoxXYZ(); v++)
        {
          m4p_smoothedImage[tbf][rbf][cbf][v] = 0.;
          m4p_sumOfWeightsImage[tbf][rbf][cbf][v] = 0.;
        }
      }
    }
  }

  // -------------------------------------------------------------------
  // Verbose
  if (m_verbose>=2)
  {
    Cout("iOptimizerModifiedEMForMRFPenalizedML::Initialize() -> Use the 1995 De Pierro's algorithm for convergent penalized ML" << endl);
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

int iOptimizerModifiedEMForMRFPenalizedML::PreImageUpdateSpecificStep()
{
  // If no penalty, then exit (the penalty image term has been initialized to 0)
  if (mp_Penalty==NULL) return 0;
  // Set the number of threads
  #ifdef CASTOR_OMP
  omp_set_num_threads(mp_ImageDimensionsAndQuantification->GetNbThreadsForImageComputation());
  #endif  
  
  // Verbose
  if (m_verbose>=1) Cout("iOptimizerModifiedEMForMRFPenalizedML::PreImageUpdateSpecificStep() -> Compute penalty term" << endl);

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
          // Cast the penalty to an MRF
          iPenaltyMarkovRandomField* p_mrf_penalty = (dynamic_cast<iPenaltyMarkovRandomField*>(mp_Penalty));
          // Build specific neighborhood of the MRF penalty (do not check the return value as this function cannot fail, and we are in a threaded loop)
          p_mrf_penalty->LocalPreProcessingStep(tbf, rbf, cbf, v, th);
          // Compute the weighted sum and sum of weights
          HPFLTNB weighted_sum = 0.;
          HPFLTNB sum_of_weights = 0.;
          // Sum on the contribution of the different neighbors
          for (INTNB n = 0; n<p_mrf_penalty->GetNeighborhoodMaxNbVoxels(); n++)
          {
            // Get the neighbor index
            INTNB neighbor = p_mrf_penalty->GetNeighborhoodIndices()[th][n];
            // If the voxel is not in the neighborhood, skip it
            if (neighbor==-1) continue;
            else
            {
              // Get the proximity and similarity factors as the weight
              HPFLTNB weight = p_mrf_penalty->GetProximityKernel()[n] * p_mrf_penalty->GetSimilarityFactors(th)[n];
              // Add the contribution of this neighbor to the sums
              weighted_sum += weight 
                           * ( mp_ImageSpace-> m4p_image[tbf][rbf][cbf][v] 
                             + mp_ImageSpace-> m4p_image[tbf][rbf][cbf][neighbor] );
              sum_of_weights += weight;
            }
          }
          // Affect penalty image and sum of weights image
          m4p_smoothedImage[tbf][rbf][cbf][v] = (FLTNB) weighted_sum;
          m4p_sumOfWeightsImage[tbf][rbf][cbf][v] = (FLTNB) sum_of_weights;
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

int iOptimizerModifiedEMForMRFPenalizedML::SensitivitySpecificOperations( FLTNB a_data, FLTNB a_forwardModel, FLTNB* ap_weight,
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

int iOptimizerModifiedEMForMRFPenalizedML::DataSpaceSpecificOperations( FLTNB a_data, FLTNB a_forwardModel, FLTNB* ap_backwardValues,
                                                                                FLTNB a_multiplicativeCorrections, FLTNB a_additiveCorrections, FLTNB a_blankValue,
                                                                                FLTNB a_quantificationFactor, oProjectionLine* ap_Line )
{
  // Truncate data to 0 if negative
  if (a_data<0.) a_data = 0.;
  // If the foward model is to close to zero, then ignore this data (set to 1 and not 0 because this line is taken into account in the sensitivity)
  if (a_forwardModel>m_dataSpaceDenominatorThreshold) *ap_backwardValues = a_data / a_forwardModel;
  else *ap_backwardValues = 1.;
  // That's all
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerModifiedEMForMRFPenalizedML::ImageSpaceSpecificOperations( FLTNB a_currentImageValue, FLTNB* ap_newImageValue,
                                                                                 FLTNB a_sensitivity, FLTNB* ap_correctionValues,
                                                                                 INTNB a_voxel, int a_tbf, int a_rbf, int a_cbf )
{
  // Note: HPFLTNB is used here, which should be double precision at least, otherwise numerical problems happen (Inf and NaN)

  // Bit by bit, image value can become zero, so we exit here to avoid numerical problems
  if (fpclassify(a_currentImageValue)==FP_ZERO) return 0;

  // ===== Start by computing the MLEM update
  HPFLTNB mlem_image_update_factor = *ap_correctionValues / a_sensitivity;
  // Apply minimum image update factor
  if ( m_minimumImageUpdateFactor > 0. && mlem_image_update_factor < m_minimumImageUpdateFactor ) mlem_image_update_factor = m_minimumImageUpdateFactor;
  // Apply maximum image update factor
  if ( m_maximumImageUpdateFactor > 0. && mlem_image_update_factor > m_maximumImageUpdateFactor ) mlem_image_update_factor = m_maximumImageUpdateFactor;
  // Update image
  HPFLTNB mlem_update = a_currentImageValue * mlem_image_update_factor;
  // Set the new value to the MLEM value by default
  *ap_newImageValue = mlem_update;

  // ===== If a penalty exists, then include it using De Pierro's derivation for quadratic penalty
  if (mp_Penalty!=NULL)
  {
    // Get penalty strength and divide it by the number of subsets to get correct balance between likelihood and penalty
    HPFLTNB beta = mp_Penalty->GetPenaltyStrength() / ((HPFLTNB)(mp_nbSubsets[m_currentIteration]));
    // Do the rest only if beta is superior to 0, otherwise the formula is not valid
    if (beta>0.)
    {
      // Compute regularized version of the current estimate, minus the sensitivity (note here that the neighbors coefficients
      // only used the proximity factors that are intrisically normalized to 1, explaining why they do not appear here).
      HPFLTNB term = beta * m4p_smoothedImage[a_tbf][a_rbf][a_cbf][a_voxel] - a_sensitivity;
      // Compute the new image (note that we multiply the current MLEM update by the sensitivity because we already applied it)
      *ap_newImageValue = (term + sqrt( term*term + 8.*((HPFLTNB)(m4p_sumOfWeightsImage[a_tbf][a_rbf][a_cbf][a_voxel]))*beta*((HPFLTNB)(a_sensitivity))*mlem_update ))
                        / (4.*((HPFLTNB)(m4p_sumOfWeightsImage[a_tbf][a_rbf][a_cbf][a_voxel]))*beta);
    }
    // Apply the rule that limits multiplicative updates to reasonable steps
    // (this operation works because we checked for negative current image value at the beginning of the function)
    FLTNB image_update_factor = *ap_newImageValue / a_currentImageValue;
    // Apply minimum image update factor
    if ( m_minimumImageUpdateFactor > 0. && image_update_factor < m_minimumImageUpdateFactor ) image_update_factor = m_minimumImageUpdateFactor;
    // Apply maximum image update factor
    if ( m_maximumImageUpdateFactor > 0. && image_update_factor > m_maximumImageUpdateFactor ) image_update_factor = m_maximumImageUpdateFactor;
    // Update image
    *ap_newImageValue = a_currentImageValue * image_update_factor;
    // Stabilize float numbers (this situation has to be controlled, otherwise the image update factor at the next iteration will be Inf, then Nan)
    if (fpclassify(*ap_newImageValue)!=FP_NORMAL) *ap_newImageValue = 0.;
  }

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

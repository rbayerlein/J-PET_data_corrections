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
  \brief    Implementation of class iOptimizerBSREM
*/

#include "iOptimizerBSREM.hh"
#include "sOutputManager.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iOptimizerBSREM::iOptimizerBSREM() : vOptimizer()
{
  // ---------------------------
  // Mandatory member parameters
  // ---------------------------

  // Initial value at 1
  m_initialValue = 1.;
  // Only one backward image for BSREM
  m_nbBackwardImages = 1;
  // BSREM accepts penalties, which must have a derivative order of 1 minimum
  m_requiredPenaltyDerivativesOrder = 1;
  // BSREM is only compatible with histogram data
  m_listmodeCompatibility = false;
  m_histogramCompatibility = true;
  // BSREM is only compatible with emission data
  m_emissionCompatibility = true;
  m_transmissionCompatibility = false;
  // BSREM needs the global sensitivity which is computed at the beginning of the reconstruction process.
  // In Ahn & Fessler 2003, it is clearly explained that the algorithm cannot use a scale function that
  // changes over subsets, otherwise the convergence to "the desired optimum point" is not guaranteed.
  m_needGlobalSensitivity = true;

  // --------------------------
  // Specific member parameters
  // --------------------------

  m_minimumImageValue = -1.;
  m_maximumImageValue = -1.;
  m4p_firstDerivativePenaltyImage = NULL;
  m_relaxationFactorType = BSREM_NOT_DEFINED;
  m_relaxationFactorInitialValue = -1.;
  m_relaxationFactorStepSize = -1.;
  m_relaxationFactor = -1.;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iOptimizerBSREM::~iOptimizerBSREM()
{
  // Delete the penalty image
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
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iOptimizerBSREM::ShowHelpSpecific()
{
  cout << "This optimizer is the Block Sequential Regularized Expectation Maximization (BSREM) algorithm from S. Ahn and" << endl;
  cout << "J. Fessler, IEEE TMI, May 2003, vol. 22, pp. 613-626. Its abbreviated name in this paper is BSREM-II." << endl;
  cout << "This algorithm is the only one to have proven convergence using subsets. Its implementation is entirely based" << endl;
  cout << "on the reference paper. It may have numerical problems when a full field-of-view is used, because of the sharp" << endl;
  cout << "sensitivity loss at the edges of the cylindrical field-of-view. As it is simply based on the gradient, penalty" << endl;
  cout << "terms must have a derivative order of at least one. Without penalty, it reduces to OSEM but where the sensitivity" << endl;
  cout << "is not dependent on the current subset. This is a requirement of the algorithm, explaining why it starts by" << endl;
  cout << "computing the global sensitivity before going through iterations. The algorithm is restricted to histograms." << endl;
  cout << "The following options can be used:" << endl;
  cout << "  initial image value: to set the uniform voxel value for the initial image" << endl;
  cout << "  minimum image value: to set the minimum allowed image value (parameter 't' in the reference paper)" << endl;
  cout << "  maximum image value: to set the maximum allowed image value (parameter 'U' in the reference paper)" << endl;
  cout << "  relaxation factor type: type of relaxation factors (can be one of the following: 'classic')" << endl;
  cout << "Relaxation factors of type 'classic' correspond to what was proposed in the reference paper in equation (31)." << endl;
  cout << "This equation gives: alpha_n = alpha_0 / (gamma * iter_num + 1)" << endl;
  cout << "The iteration number 'iter_num' is supposed to start at 0 so that for the first iteration, alpha_0 is used." << endl;
  cout << "This parameter can be provided using the following keyword: 'relaxation factor classic initial value'." << endl;
  cout << "The 'gamma' parameter can be provided using the following keyword: 'relaxation factor classic step size'." << endl;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerBSREM::ReadConfigurationFile(const string& a_configurationFile)
{
  string key_word = "";
  string buffer = "";
  
  // Read the initial image value option
  key_word = "initial image value";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_initialValue, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iOptimizerBSREM::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Read the minimum image value option
  key_word = "minimum image value";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_minimumImageValue, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iOptimizerBSREM::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Read the maximum image value option
  key_word = "maximum image value";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_maximumImageValue, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iOptimizerBSREM::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  if (m_maximumImageValue<0.) m_maximumImageValue = numeric_limits<FLTNB>::infinity();
  // Read the relaxation factors type option
  key_word = "relaxation factor type";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &buffer, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iOptimizerBSREM::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  if (buffer == "classic")
  {
    m_relaxationFactorType = BSREM_RELAXATION_CLASSIC;
    // Read the relaxation factors initial value option
    key_word = "relaxation factor classic initial value";
    if (ReadDataASCIIFile(a_configurationFile, key_word, &m_relaxationFactorInitialValue, 1, KEYWORD_MANDATORY))
    {
      Cerr("***** iOptimizerBSREM::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
      return 1;
    }
    // Read the relaxation factors step size
    key_word = "relaxation factor classic step size";
    if (ReadDataASCIIFile(a_configurationFile, key_word, &m_relaxationFactorStepSize, 1, KEYWORD_MANDATORY))
    {
      Cerr("***** iOptimizerBSREM::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
      return 1;
    }
  }
  else
  {
    Cerr("***** iOptimizerBSREM::ReadConfigurationFile() -> Provided type of the relaxation factors (" << buffer << ") is not recognized" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerBSREM::ReadOptionsList(const string& a_optionsList)
{
  // For the moment, as there is only the classic relaxation parameter type, let the possibility to read the parameters from a list.

  // There are 5 floating point variables as options
  const int nb_options = 5;
  FLTNB options[nb_options];
  // Read them
  if (ReadStringOption(a_optionsList, options, nb_options, ",", "BSREM configuration"))
  {
    Cerr("***** iOptimizerBSREM::ReadOptionsList() -> Failed to correctly read the list of options !" << endl
      << "      Note: if you are not using the classic relaxation, you must for now use the configuration file" << endl);
    return 1;
  }
  // Affect options
  m_initialValue = options[0];
  m_minimumImageValue = options[1];
  m_maximumImageValue = options[2];
  m_relaxationFactorType = BSREM_RELAXATION_CLASSIC;
  m_relaxationFactorInitialValue = options[3];
  m_relaxationFactorStepSize = options[4];
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerBSREM::CheckSpecificParameters()
{
  // Check relaxation factors type
  if (m_relaxationFactorType==BSREM_NOT_DEFINED)
  {
    Cerr("***** iOptimizerBSREM::CheckSpecificParameters() -> Should provide a relaxation factor type !" << endl);
    return 1;
  }
  // Check that initial image value is strictly positive
  if (m_initialValue<=0.)
  {
    Cerr("***** iOptimizerBSREM::CheckSpecificParameters() -> Provided initial image value (" << m_initialValue << ") must be strictly positive !" << endl);
    return 1;
  }
  // Check that the minimum image value is strictly positive
  if (m_minimumImageValue<=0.)
  {
    Cerr("***** iOptimizerBSREM::CheckSpecificParameters() -> Provided minimum image value (" << m_minimumImageValue << ") must be strictly positive ! !" << endl);
    return 1;
  }
  // Check that the maximum image value is superior to the min image value
  if (m_maximumImageValue<=2*m_minimumImageValue || m_maximumImageValue<=0.)
  {
    Cerr("***** iOptimizerBSREM::CheckSpecificParameters() -> Provided maximum image value (" << m_maximumImageValue << ") must be positive and superior to twice the minimum image value !" << endl);
    return 1;
  }
  // Check parameters for the relaxation factor
  if (m_relaxationFactorType == BSREM_RELAXATION_CLASSIC)
  {
    // Check that the relaxation factors initial value is strictly positive
    if (m_relaxationFactorInitialValue<=0.)
    {
      Cerr("***** iOptimizerBSREM::CheckSpecificParameters() -> Provided relaxation factors initial value (" << m_relaxationFactorInitialValue << ") must be strictly positive ! !" << endl);
      return 1;
    }
    // Check that the relaxation factors step size is strictly positive
    if (m_relaxationFactorStepSize<=0.)
    {
      Cerr("***** iOptimizerBSREM::CheckSpecificParameters() -> Provided relaxation factors step size (" << m_relaxationFactorInitialValue << ") must be strictly positive ! !" << endl);
      return 1;
    }
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerBSREM::InitializeSpecific()
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
        //Get a pointer to a newly allocated image
        m4p_firstDerivativePenaltyImage[tbf][rbf][cbf] = mp_ImageSpace -> AllocateMiscellaneousImage();
        // Initialize to 0, in case the penalty is not used
        for (INTNB v=0; v<mp_ImageDimensionsAndQuantification->GetNbVoxXYZ(); v++) m4p_firstDerivativePenaltyImage[tbf][rbf][cbf][v] = 0.;
      }
    }
  }
  // Verbose
  if (m_verbose>=VERBOSE_NORMAL)
  {
    Cout("iOptimizerBSREM::InitializeSpecific() -> Use the BSREM optimizer" << endl);
    if (m_verbose>=VERBOSE_DETAIL)
    {
      Cout("  --> Initial image value: " << m_initialValue << endl);
      Cout("  --> Minimum image value: " << m_minimumImageValue << endl);
      Cout("  --> Maximum image value: " << m_maximumImageValue << endl);
      if (m_relaxationFactorType == BSREM_RELAXATION_CLASSIC)
      {
        Cout("  --> Relaxation factors type: classic" << endl);
        Cout("  --> Relaxation factors initial value: " << m_relaxationFactorInitialValue << endl);
        Cout("  --> Relaxation factors step size: " << m_relaxationFactorStepSize << endl);
      }
      else
      {
        Cerr("***** iOptimizerBSREM::InitializeSpecific() -> Provided relaxation factors type (" << m_relaxationFactorType << ") is not recognized !" << endl);
        return 1;
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

int iOptimizerBSREM::PreImageUpdateSpecificStep()
{
  // ==========================================================================================
  // If no penalty, then exit (the penalty image term has been initialized to 0)
  if (mp_Penalty==NULL) return 0;
  // Set the number of threads
  #ifdef CASTOR_OMP
  omp_set_num_threads(mp_ImageDimensionsAndQuantification->GetNbThreadsForImageComputation());
  #endif
  // Verbose
  if (m_verbose>=1) Cout("iOptimizerBSREM::PreImageUpdateSpecificStep() -> Compute penalty term" << endl);
  // ==========================================================================================
  // Global precomputation step if needed by the penalty
  if (mp_Penalty->GlobalPreProcessingStep())
  {
    Cerr("***** iOptimizerBSREM::PreImageUpdateSpecificStep() -> A problem occurred while computing the penalty pre-processing step !" << endl);
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
            Cerr("***** iOptimizerBSREM::PreImageUpdateSpecificStep() -> A problem occurred while computing the penalty local pre-processing step for voxel " << v << " !" << endl);
            problem = true;
          }
          // Compute the penalty term at first order
          m4p_firstDerivativePenaltyImage[tbf][rbf][cbf][v] = mp_Penalty->ComputeFirstDerivative(tbf,rbf,cbf,v,th);
        }
        // Check for problems
        if (problem)
        {
          Cerr("***** iOptimizerBSREM::PreImageUpdateSpecificStep() -> A problem occurred inside the multi-threaded loop, stop now !" << endl);
          return 1;
        }
      }
    }
  }
  // ==========================================================================================
  // Update the relaxation factor and the number of subsets if we are in a new loop
  if (m_currentSubset == 0)
  {
    if (m_relaxationFactorType == BSREM_RELAXATION_CLASSIC)
    {
      m_relaxationFactor = m_relaxationFactorInitialValue / (m_relaxationFactorStepSize * ((HPFLTNB)m_currentIteration) + 1.);
      if (m_verbose>=2) Cout("  --> Relaxation factor for this iteration: " << m_relaxationFactor << endl);
    }
    else
    {
      Cerr("***** iOptimizerBSREM::PostDataUpdateSpecificStep -> Unknown relaxation factors type (" << m_relaxationFactorType << ") !" << endl);
      return 1;
    }
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerBSREM::SensitivitySpecificOperations( FLTNB a_data, FLTNB a_forwardModel, FLTNB* ap_weight,
                                                   FLTNB a_multiplicativeCorrections, FLTNB a_additiveCorrections, FLTNB a_blankValue,
                                                   FLTNB a_quantificationFactor, oProjectionLine* ap_Line )
{
  // Line weight here is simply 1 or any other value. It will actually not be used
  *ap_weight = 1.;
  // That's all
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerBSREM::DataSpaceSpecificOperations( FLTNB a_data, FLTNB a_forwardModel, FLTNB* ap_backwardValues,
                                                 FLTNB a_multiplicativeCorrections, FLTNB a_additiveCorrections, FLTNB a_blankValue,
                                                 FLTNB a_quantificationFactor, oProjectionLine* ap_Line )
{
  // Compute numerator
  FLTNB numerator = a_data - a_forwardModel;
  // Check if the denominator is zero
  if (a_forwardModel==0.) *ap_backwardValues = 0.;
  // Update backward values
  else *ap_backwardValues = numerator / a_forwardModel;
  // That's all
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerBSREM::ImageSpaceSpecificOperations( FLTNB a_currentImageValue, FLTNB* ap_newImageValue,
                                                   FLTNB a_sensitivity, FLTNB* ap_correctionValues,
                                                   INTNB a_voxel, int a_tbf, int a_rbf, int a_cbf )
{
  // Note: This algorithm is sensitive to precision of computations, so everything is performed in HPFLTNB.

  // Check if correction values (backward matrix) is zero, this means that no data are contributing to the update of this voxel.
  // In this case, we do not update the voxel value only based on penalty, because it quickly causes numerical problems.
  if (*ap_correctionValues==0.) return 0;

  // ===== Implementation of BSREM-I algorithm

  // Scale penalty with respect to the number of subsets to get correct balance between likelihood and penalty
  HPFLTNB scaled_penalty = ((HPFLTNB)(m4p_firstDerivativePenaltyImage[a_tbf][a_rbf][a_cbf][a_voxel])) / ((HPFLTNB)(mp_nbSubsets[m_currentIteration]));
  // Additive update factor without the image scaling (equation (21) in the reference paper)
  HPFLTNB image_update_factor = m_relaxationFactor * (((HPFLTNB)(*ap_correctionValues))-scaled_penalty) / ((HPFLTNB)a_sensitivity);
  // Image scaling of the update factor depends on image value (equation (22) in the reference paper)
  if (a_currentImageValue<=m_maximumImageValue/2.) image_update_factor *= ((HPFLTNB)a_currentImageValue);
  else image_update_factor *= (m_maximumImageValue-((HPFLTNB)a_currentImageValue));
  // Update the image
  *ap_newImageValue = ((HPFLTNB)a_currentImageValue) + image_update_factor;

  // ===== Up to here, this is the BSREM-I algorithm.
  //       Then, the two following lines implement the BSREM-II algorithm (equation (23) and footnote (9) in the reference paper).

  // Check if the new image value is lower than 0
  if (*ap_newImageValue < 0.) *ap_newImageValue = m_minimumImageValue;
  // Check if the new image value is higher than the max image value
  if (*ap_newImageValue > m_maximumImageValue) *ap_newImageValue = m_maximumImageValue - m_minimumImageValue;

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

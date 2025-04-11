
/*!
  \file
  \ingroup  optimizer
  \brief    Implementation of class iOptimizerSens
*/

#include "iOptimizerSens.hh"
#include "sOutputManager.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iOptimizerSens::iOptimizerSens() : vOptimizer()
{
  // ---------------------------
  // Mandatory member parameters
  // ---------------------------

  // Initial value at 1
  m_initialValue = 1.;
  // Only one backward image for MLEM
  m_nbBackwardImages = 1;
  // MLEM does not accept penalties
  //m_penaltyEnergyFunctionDerivativesOrder = 0;
  // Compatible with listmode and histogram data
  m_listmodeCompatibility = true;
  m_histogramCompatibility = true;
  // Compatible with both emission and log-converted transmission data
  m_emissionCompatibility = true;
  m_transmissionCompatibility = true;

  // --------------------------
  // Specific member parameters
  // --------------------------

  m_dataSpaceDenominatorThreshold = -1.;
  m_minimumImageUpdateFactor = -1.;
  m_maximumImageUpdateFactor = -1.;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iOptimizerSens::~iOptimizerSens()
{
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iOptimizerSens::ShowHelpSpecific()
{
  cout << "This is a specific optimizer to compute a sensitivity image from a uniform acquisition, such as a uniform cylinder covering the FOV." << endl;
  cout << "To use it, just add '-opti SENS' in the command-line option." << endl;
  cout << "The data update step treats the input image as a mumap in cm-1 (data update step = 1 / exp(forwardProj*0.1)" << endl;
  cout << "The image update step just take the image update factor as the actual image value" << endl;
  cout << "This is just a quick/cheap implementation of such process, so here is a few remarks to make it work:" << endl;
  cout << "- Iterations/subsets MUST both be 1 (-it 1:1)" << endl;
  cout << "- For attenuation, an input mumap (cm-1) could be incorporated with the option -img. Make sure however it has the same numbers of voxel and voxel sizes than reconstruction" << endl;
  cout << "- Calibration factor in the datafile header MUST be set to 1" << endl;
  cout << "- The command-line options MUST disable frame duration correction, use the following command: 'ignore-corr fdur' " << endl;
  cout << "- The algorithm will still try to compute a sensitivity image before the reconstruction process, whereas we just want the reconstruction to compute our sensitivity image" << endl;
  cout << "  Just input any image with the -sens option to bypass that step (the values from this image will not be taken into account anyway)." << endl;
  cout << "  Make sure it is the same dimensions as reconstruction dimensions and voxel sizes" << endl;  
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerSens::ReadConfigurationFile(const string& a_configurationFile)
{
  string key_word = "";
  // Read the initial image value option
  key_word = "initial image value";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_initialValue, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iOptimizerSens::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Read the denominator threshold option
  key_word = "denominator threshold";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_dataSpaceDenominatorThreshold, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iOptimizerSens::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Read the minimum image update option
  key_word = "minimum image update";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_minimumImageUpdateFactor, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iOptimizerSens::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Read the maximum image update option
  key_word = "maximum image update";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_maximumImageUpdateFactor, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iOptimizerSens::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerSens::ReadOptionsList(const string& a_optionsList)
{
  // There are 4 floating point variables as options
  const int nb_options = 4;
  FLTNB options[nb_options];
  
  // Read them
  if (ReadStringOption(a_optionsList, options, nb_options, ",", "MLEM configuration"))
  {
    Cerr("***** iOptimizerSens::ReadOptionsList() -> Failed to correctly read the list of options !" << endl);
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

int iOptimizerSens::CheckSpecificParameters()
{
  // Check that initial image value is strictly positive
  if (m_initialValue<0.)
  {
    Cerr("***** iOptimizerSens->CheckSpecificParameters() -> Provided initial image value (" << m_initialValue << ") must be positive !" << endl);
    return 1;
  }
  // Check that denominator threshold value is strictly positive
  if (m_dataSpaceDenominatorThreshold<=0.)
  {
    Cerr("***** iOptimizerSens->CheckSpecificParameters() -> Provided data space denominator threshold (" << m_dataSpaceDenominatorThreshold << ") must be strictly positive !" << endl);
    return 1;
  }
  // Check that maximum image update factor is higher than the minimum
  if (m_minimumImageUpdateFactor>0. && m_maximumImageUpdateFactor>0. && m_maximumImageUpdateFactor<m_minimumImageUpdateFactor)
  {
    Cerr("***** iOptimizerSens->CheckSpecificParameters() -> Provided minimum/maximum (" << m_minimumImageUpdateFactor << "/" << m_maximumImageUpdateFactor << " are inconsistent !" << endl);
    return 1;
  }

  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iOptimizerSens::InitializeSpecific()
{
  // Verbose
  if (m_verbose>=VERBOSE_NORMAL)
  {
    Cout("iOptimizerSens::InitializeSpecific() -> Use the SENS optimizer" << endl);
    if (m_verbose>=VERBOSE_DETAIL)
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

int iOptimizerSens::SensitivitySpecificOperations( FLTNB a_data, FLTNB a_forwardModel, FLTNB* ap_weight,
                                                   FLTNB a_multiplicativeCorrections, FLTNB a_additiveCorrections, FLTNB a_blankValue,
                                                   FLTNB a_quantificationFactor, oProjectionLine* ap_Line )

//int iOptimizerSens::SensitivitySpecificOperations( FLTNB a_data, FLTNB a_forwardModel, FLTNB* ap_weight,
//                                                   FLTNB a_multiplicativeCorrections, FLTNB a_additiveCorrections,
//                                                   FLTNB a_quantificationFactor, oProjectionLine* ap_Line )
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


int iOptimizerSens::DataSpaceSpecificOperations( FLTNB a_data, FLTNB a_forwardModel, FLTNB* ap_backwardValues,
                                                 FLTNB a_multiplicativeCorrections, FLTNB a_additiveCorrections, FLTNB a_blankValue,
                                                 FLTNB a_quantificationFactor, oProjectionLine* ap_Line )
{
  // ----------------------------------------------------------------------------------------------
  // Part 0: check the multiplication and quantification factors
  // ----------------------------------------------------------------------------------------------

  // If multiplicative correction factor is null, then skip this event
  if (a_multiplicativeCorrections<=0.) return 0;

  // If quantification factor is null, then skip this event
  if (a_quantificationFactor<=0.) return 0;


  // First set to 1 each voxel of the sensitivitty image crossed by the line
  // (unless no voxel will be visited by the restriction on sensitivity image value in vOptimizer::ImageUpdate
  for (int b=0; b<ap_Line->GetNbTOFBins(); b++)
    for (int vl=0 ; vl<ap_Line->GetCurrentNbVoxels(FORWARD, b) ; vl++)
      mp_ImageSpace->m5p_sensitivity[0][0][0][0][ap_Line->GetVoxelIndex(FORWARD, b,vl)] = 1.;


  // compute mu*x correcting by the 'a_quantificationFactor' applied
  // in 'FLTNB oProjectionLine::ForwardProject()' and convert from cm-1 to mm-1
  
  FLTNB mux = a_forwardModel * (FLTNB)0.1; 

  // The Attenuation Correction Factor
  FLTNB acf = exp(-mux);

  // '*ap_backwardValues' is taken into account as the sensitivity
  // not applying 'a_quantificationFactor' because is already included in
  // 'm_multiplicativeCorrection' in method 'oProjectionLine::BackwardProject()'
  *ap_backwardValues = acf;

  // That's all
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
int iOptimizerSens::ImageSpaceSpecificOperations( FLTNB a_currentImageValue, FLTNB* ap_newImageValue,
                                                  FLTNB a_sensitivity, FLTNB* ap_correctionValues,
                                                  INTNB a_voxel, int a_tbf, int a_rbf, int a_cbf )

{
  //cout << "*ap_correctionValues " << *ap_correctionValues << endl;
  // Compute image update factor
  FLTNB image_update_factor = *ap_correctionValues ;
  // Apply minimum image update factor
  if ( m_minimumImageUpdateFactor > 0. && image_update_factor < m_minimumImageUpdateFactor ) image_update_factor = m_minimumImageUpdateFactor;
  // Apply maximum image update factor
  if ( m_maximumImageUpdateFactor > 0. && image_update_factor > m_maximumImageUpdateFactor ) image_update_factor = m_maximumImageUpdateFactor;
  // Update image
  *ap_newImageValue = image_update_factor;

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================


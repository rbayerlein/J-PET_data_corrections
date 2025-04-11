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
  \brief    Implementation of class oOptimizerManager
*/

#include "oOptimizerManager.hh"
#include "sOutputManager.hh"
#include "sAddonManager.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

oOptimizerManager::oOptimizerManager()
{
  // Image dimensions
  mp_ImageDimensionsAndQuantification = NULL;
  // Image space
  mp_ImageSpace = NULL;
  // Data mode
  m_dataMode = MODE_UNKNOWN;
  // Data type
  m_dataType = TYPE_UNKNOWN;
  // Data spec
  m_dataSpec = SPEC_UNKNOWN;
  // Number of TOF bins
  m_nbTOFBins = 0;
  // Optimizer and penalty options
  m_optionsOptimizer = "";
  m_optionsPenalty = "";
  m_penaltyStrength = -1.;
  // Optimizer and penalty
  mp_Optimizer = NULL;
  mp_Penalty = NULL;
  // Verbosity
  m_verbose = 0;
  // Optimizer FOM computation, and image update stat flags
  m_optimizerFOMFlag = false;
  m_optimizerImageStatFlag = false;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

oOptimizerManager::~oOptimizerManager() 
{
  if (mp_Optimizer) delete(mp_Optimizer);
  if (mp_Penalty) delete(mp_Penalty);
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oOptimizerManager::CheckParameters()
{
  // Check image dimensions
  if (mp_ImageDimensionsAndQuantification==NULL)
  {
    Cerr("***** oOptimizerManager::CheckParameters() -> No image dimensions provided !" << endl);
    return 1;
  }
  // Check image space
  if (mp_ImageSpace==NULL)
  {
    Cerr("***** oOptimizerManager::CheckParameters() -> No image space provided !" << endl);
    return 1;
  }
  // Check data mode
  if (m_dataMode==MODE_UNKNOWN || (m_dataMode!=MODE_LIST && m_dataMode!=MODE_HISTOGRAM))
  {
    Cerr("***** oOptimizerManager::CheckParameters() -> No or meaningless data mode provided !" << endl);
    return 1;
  }
  // Check data type
  if (m_dataType==TYPE_UNKNOWN || (m_dataType!=TYPE_PET && m_dataType!=TYPE_SPECT && m_dataType!=TYPE_CT))
  {
    Cerr("***** oOptimizerManager::CheckParameters() -> No or meaningless data type provided !" << endl);
    return 1;
  }
  // Check data spec
  if (m_dataSpec==SPEC_UNKNOWN || (m_dataSpec!=SPEC_EMISSION && m_dataSpec!=SPEC_TRANSMISSION))
  {
    Cerr("***** oOptimizerManager::CheckParameters() -> No or meaningless data specificity provided (emission or transmission) !" << endl);
    return 1;
  }
  // Check optimizer options
  if (m_optionsOptimizer=="")
  {
    Cerr("***** oOptimizerManager::CheckParameters() -> No optimizer options provided !" << endl);
    return 1;
  }
  // Check verbosity
  if (m_verbose<0)
  {
    Cerr("***** oOptimizerManager::CheckParameters() -> Wrong verbosity level provided !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oOptimizerManager::Initialize()
{
  // Verbose
  if (m_verbose>=1) Cout("oOptimizerManager::Initialize() -> Initialize optimizer and penalty" << endl);

  // Parse projector options and initialize them
  if (ParseOptionsAndInitializeOptimizerAndPenalty())
  {
    Cerr("***** oOptimizerManager::Initialize() -> A problem occurred while parsing optimizer options and initializing it !" << endl);
    return 1;
  }

  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oOptimizerManager::ParseOptionsAndInitializeOptimizerAndPenalty()
{
  // ---------------------------------------------------------------------------------------------------
  // Manage optimizer options
  // ---------------------------------------------------------------------------------------------------

  string name_optimizer = "";
  string list_options_optimizer = "";
  string file_options_optimizer = "";

  // This is for the automatic initialization of the optimizer and penalty
  typedef vOptimizer *(*maker_optimizer) ();

  // ______________________________________________________________________________
  // Get the optimizer name in the options and isolate the real optimizer's options

  // First check emptyness of the options
  if (m_optionsOptimizer=="")
  {
    Cerr("***** oOptimizerManager::ParseOptionsAndInitializeOptimizerAndPenalty() -> No optimizer provided !" << endl);
    return 1;
  }

  // Search for a colon ":", this indicates that a configuration file is provided after the optimizer name
  size_t colon = m_optionsOptimizer.find_first_of(":");
  size_t comma = m_optionsOptimizer.find_first_of(",");

  // Case 1: we have a colon
  if (colon!=string::npos)
  {
    // Get the optimizer name before the colon
    name_optimizer = m_optionsOptimizer.substr(0,colon);
    // Get the configuration file after the colon
    file_options_optimizer = m_optionsOptimizer.substr(colon+1);
    // List of options is empty
    list_options_optimizer = "";
  }
  // Case 2: we have a comma
  else if (comma!=string::npos)
  {
    // Get the optimizer name before the first comma
    name_optimizer = m_optionsOptimizer.substr(0,comma);
    // Get the list of options after the first comma
    list_options_optimizer = m_optionsOptimizer.substr(comma+1);
    // Configuration file is empty
    file_options_optimizer = "";
  }
  // Case 3: no colon and no comma (a single optimizer name)
  else
  {
    // Get the optimizer name
    name_optimizer = m_optionsOptimizer;
    // List of options is empty
    list_options_optimizer = "";
    // Build the default configuration file
    sOutputManager* p_output_manager = sOutputManager::GetInstance();
    file_options_optimizer = p_output_manager->GetPathToConfigDir() + "/optimizer/" + name_optimizer + ".conf";
  }

  // ______________________________________________________________________________
  // Construct and initialize the optimizer

  // Get optimizer's listfrom addon manager
  std::map <string,maker_optimizer> list_optimizer = sAddonManager::GetInstance()->mp_listOfOptimizers;
  // Create the optimizer
  if (list_optimizer[name_optimizer]) mp_Optimizer = list_optimizer[name_optimizer]();
  else
  {
    Cerr("***** oOptimizerManager::ParseOptionsAndInitializeOptimizerAndPenalty() -> Optimizer '" << name_optimizer << "' does not exist !" << endl);
    return 1;
  }
  // Set parameters
  mp_Optimizer->SetOptimizerID(name_optimizer);
  mp_Optimizer->SetImageDimensionsAndQuantification(mp_ImageDimensionsAndQuantification);
  mp_Optimizer->SetImageSpace(mp_ImageSpace);
  mp_Optimizer->SetDataMode(m_dataMode);
  mp_Optimizer->SetDataType(m_dataType);
  mp_Optimizer->SetDataSpec(m_dataSpec);
  mp_Optimizer->SetNbTOFBins(m_nbTOFBins);
  mp_Optimizer->SetFOMFlag(m_optimizerFOMFlag);
  mp_Optimizer->SetImageStatFlag(m_optimizerImageStatFlag);
  mp_Optimizer->SetVerbose(m_verbose);
  // Provide configuration file if any (child specific function)
  if (file_options_optimizer!="" && mp_Optimizer->ReadConfigurationFile(file_options_optimizer))
  {
    Cerr("***** oOptimizerManager::ParseOptionsAndInitializeOptimizerAndPenalty() -> A problem occurred while reading and checking optimizer's configuration file !" << endl);
    return 1;
  }
  // Provide options if any (child specific function)
  if (list_options_optimizer!="" && mp_Optimizer->ReadOptionsList(list_options_optimizer))
  {
    Cerr("***** oOptimizerManager::ParseOptionsAndInitializeOptimizerAndPenalty() -> A problem occurred while parsing and reading optimizer's options !" << endl);
    return 1;
  }
  // Check parameters (mother generic function that will call the child specific function at the end)
  if (mp_Optimizer->CheckParameters())
  {
    Cerr("***** oOptimizerManager::ParseOptionsAndInitializeOptimizerAndPenalty() -> A problem occurred while checking optimizer parameters !" << endl);
    return 1;
  }

  // ---------------------------------------------------------------------------------------------------
  // Manage penalty options
  // ---------------------------------------------------------------------------------------------------

  string name_penalty = "";
  string list_options_penalty = "";
  string file_options_penalty = "";

  // This is for the automatic initialization of the penalty
  typedef vPenalty *(*maker_penalty) ();

  // ______________________________________________________________________________
  // Get the penalty name in the options and isolate the real penalty's options

  // Do we have some penalty provided ?
  if (m_optionsPenalty!="")
  {

    // Otherwise, check if the optimizer accepts penalties
    if (!mp_Optimizer->GetAcceptPenalty())
    {
      Cerr("***** oOptimizerManager::ParseOptionsAndInitializeOptimizerAndPenalty() -> Penalty provided while the selected optimizer does not accept penalties !" << endl);
      Cerr("                                                                           Remove penalty or change for another optimizer that accepts penalties." << endl);
      return 1;
    }
    // Search for a colon ":", this indicates that a configuration file is provided after the penalty name
    colon = m_optionsPenalty.find_first_of(":");
    comma = m_optionsPenalty.find_first_of(",");
    // Case 1: we have a colon
    if (colon!=string::npos)
    {
      // Get the penalty name before the colon
      name_penalty = m_optionsPenalty.substr(0,colon);
      // Get the configuration file after the colon
      file_options_penalty = m_optionsPenalty.substr(colon+1);
      // List of options is empty
      list_options_penalty = "";
    }
    // Case 2: we have a comma
    else if (comma!=string::npos)
    {
      // Get the penalty name before the first comma
      name_penalty = m_optionsPenalty.substr(0,comma);
      // Get the list of options after the first comma
      list_options_penalty = m_optionsPenalty.substr(comma+1);
      // Configuration file is empty
      file_options_penalty = "";
    }
    // Case 3: no colon and no comma (a single penalty name)
    else
    {
      // Get the penalty name
      name_penalty = m_optionsPenalty;
      // List of options is empty
      list_options_penalty = "";
      // Build the default configuration file
      sOutputManager* p_output_manager = sOutputManager::GetInstance();
      file_options_penalty = p_output_manager->GetPathToConfigDir() + "/optimizer/" + name_penalty + ".conf";
    }
    // ______________________________________________________________________________
    // Construct and initialize the penalty
    // Get penalty's list from addon manager
    std::map <string,maker_penalty> list_penalty = sAddonManager::GetInstance()->mp_listOfPenalties;
    // Create the penalty
    if (list_penalty[name_penalty]) mp_Penalty = list_penalty[name_penalty]();
    else
    {
      Cerr("***** oOptimizerManager::ParseOptionsAndInitializeOptimizerAndPenalty() -> Penalty '" << name_penalty << "' does not exist !" << endl);
      return 1;
    }
    // Set parameters
    mp_Penalty->SetPenaltyID(name_penalty);
    mp_Penalty->SetPenaltyStrength(m_penaltyStrength);
    mp_Penalty->SetImageDimensionsAndQuantification(mp_ImageDimensionsAndQuantification);
    mp_Penalty->SetImageSpace(mp_ImageSpace);
    mp_Penalty->SetVerbose(m_verbose);
    mp_Optimizer->SetPenalty(mp_Penalty);
    // Provide configuration file if any
    if (file_options_penalty!="" && mp_Penalty->ReadConfigurationFile(file_options_penalty))
    {
      Cerr("***** oOptimizerManager::ParseOptionsAndInitializeOptimizerAndPenalty() -> A problem occurred while reading and checking penalty's configuration file !" << endl);
      return 1;
    }
    // Provide options if any
    if (list_options_penalty!="" && mp_Penalty->ReadOptionsList(list_options_penalty))
    {
      Cerr("***** oOptimizerManager::ParseOptionsAndInitializeOptimizerAndPenalty() -> A problem occurred while parsing and reading penalty's options !" << endl);
      return 1;
    }
    // Check parameters
    if (mp_Penalty->CheckParameters())
    {
      Cerr("***** oOptimizerManager::ParseOptionsAndInitializeOptimizerAndPenalty() -> A problem occurred while checking penalty parameters !" << endl);
      return 1;
    }
    // Initialize the penalty
    if (mp_Penalty->Initialize())
    {
      Cerr("***** oOptimizerManager::ParseOptionsAndInitializeOptimizerAndPenalty() -> A problem occurred while initializing the penalty !" << endl);
      return 1;
    }
    // Check that the derivatives order of the penalty is compatible with the optimizer
    if (mp_Penalty->GetPenaltyDerivativesOrder() < mp_Optimizer->GetRequiredPenaltyDerivativesOrder())
    {
      Cerr("***** oOptimizerManager::ParseOptionsAndInitializeOptimizerAndPenalty() -> Derivatives order allowed by chosen penalty is not compatible with the order required by the optimizer !" << endl);
      return 1;
    }
  }

  // ---------------------------------------------------------------------------------------------------
  // Finally initialize the optimizer (mother generic function that will call the child specific function at the end)
  // ---------------------------------------------------------------------------------------------------
  if (mp_Optimizer->Initialize())
  {
    Cerr("***** oOptimizerManager::ParseOptionsAndInitializeOptimizerAndPenalty() -> A problem occurred while initializing the optimizer !" << endl);
    return 1;
  }

  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oOptimizerManager::PreDataUpdateStep()
{
  // Simply call the homonymous function from the vOptimizer
  if (mp_Optimizer->PreDataUpdateStep())
  {
    Cerr("***** oOptimizerManager::PreDataUpdateStep() -> A problem occurred while applying the pre-data-update step to the optimizer !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oOptimizerManager::PreImageUpdateStep()
{
  // Simply call the homonymous function from the vOptimizer
  if (mp_Optimizer->PreImageUpdateStep())
  {
    Cerr("***** oOptimizerManager::PreImageUpdateStep() -> A problem occurred while applying the pre-image-update step to the optimizer !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oOptimizerManager::DataUpdateStep( oProjectionLine* ap_Line, vEvent* ap_Event,
                                       int a_bed, int a_timeFrame, int a_respGate, int a_cardGate,
                                       int a_thread)
{
  // ---------------------------------------------------------------------------------
  // Deal with all multiplicative correction factors to be included in the projections
  // ---------------------------------------------------------------------------------

  // Compute the global multiplicative correction factor
  FLTNB multiplicative_correction = ap_Event->GetMultiplicativeCorrections() * mp_ImageDimensionsAndQuantification->GetQuantificationFactor(a_bed,a_timeFrame, a_respGate, a_cardGate);

  // Do nothing if the multiplicative correction factor is null or negative
  if (multiplicative_correction<=0.) return 0;

  // With transmission data, we include a correction factor to convert the default
  // mm-1 unit of CASToR into cm-1 which is the standard attenuation unit in physics
  // (this means that the image is reconstructed in cm-1)
  if (m_dataSpec==SPEC_TRANSMISSION) multiplicative_correction *= 10.;

  // Set the multiplicative correction to the oProjectionLine so that it is part of the system matrix and automatically applied
  ap_Line->SetMultiplicativeCorrection(multiplicative_correction);

  // Set the current attenuation image for SPECT attenuation correction
  if (m_dataType==TYPE_SPECT && mp_ImageSpace->m4p_attenuation!=NULL) mp_Optimizer->SetAttenuationImage(mp_ImageSpace->m4p_attenuation[a_timeFrame][a_respGate][a_cardGate], a_thread);

  // --------------------------------------------------------------------------------
  // Decompose the data update step into 4 steps mandatory steps and 3 optional steps
  // --------------------------------------------------------------------------------

  // Mandatory 1: Compute model (forward-projection + additive background)
  if (mp_Optimizer->DataStep1ForwardProjectModel( ap_Line, ap_Event, a_bed, a_timeFrame, a_respGate, a_cardGate, a_thread ))
  {
    Cerr("***** oOptimizerManager::DataUpdateStep() -> An error occurred while forward projecting !" << endl);
    return 1;
  }

  // Optional 1: Do what is done in the child optimizer
  if (mp_Optimizer->DataStep2Optional( ap_Line, ap_Event, a_bed, a_timeFrame, a_respGate, a_cardGate, a_thread ))
  {
    Cerr("***** oOptimizerManager::DataUpdateStep() -> An error occurred while performing optional step 1 !" << endl);
    return 1;
  }

  // Mandatory 2: Compute sensitivity in histogram mode
  if ( ap_Event->GetDataMode()==MODE_HISTOGRAM && !mp_Optimizer->GetNeedGlobalSensitivity()
    && mp_Optimizer->DataStep3BackwardProjectSensitivity( ap_Line, ap_Event, a_bed, a_timeFrame, a_respGate, a_cardGate, a_thread ))
  {
    Cerr("***** oOptimizerManager::DataUpdateStep() -> An error occurred while backward projecting the sensitivity !" << endl);
    return 1;
  }

  // Optional 2: Do what is done in the child optimizer
  if (mp_Optimizer->DataStep4Optional( ap_Line, ap_Event, a_bed, a_timeFrame, a_respGate, a_cardGate, a_thread ))
  {
    Cerr("***** oOptimizerManager::DataUpdateStep() -> An error occurred while performing optional step 2 !" << endl);
    return 1;
  }

  // Mandatory 3: Compute the correction terms
  if (mp_Optimizer->DataStep5ComputeCorrections( ap_Line, ap_Event, a_bed, a_timeFrame, a_respGate, a_cardGate, a_thread ))
  {
    Cerr("***** oOptimizerManager::DataUpdateStep() -> An error occurred while computing correction terms !" << endl);
    return 1;
  }

  // Optional 3: Do what is done in the child optimizer
  if (mp_Optimizer->DataStep6Optional( ap_Line, ap_Event, a_bed, a_timeFrame, a_respGate, a_cardGate, a_thread ))
  {
    Cerr("***** oOptimizerManager::DataUpdateStep() -> An error occurred while performing optional step 3 !" << endl);
    return 1;
  }

  // Mandatory 4: Backproject correction terms
  if (mp_Optimizer->DataStep7BackwardProjectCorrections( ap_Line, ap_Event, a_bed, a_timeFrame, a_respGate, a_cardGate, a_thread ))
  {
    Cerr("***** oOptimizerManager::DataUpdateStep() -> An error occurred while backward projecting the correction !" << endl);
    return 1;
  }

  // Compute FOM is asked for
  if (m_optimizerFOMFlag && mp_Optimizer->DataStep8ComputeFOM( ap_Line, ap_Event, a_timeFrame, a_respGate, a_cardGate, a_thread ))
  {
    Cerr("***** oOptimizerManager::DataUpdateStep() -> An error occurred while computing FOMs !" << endl);
    return 1;
  }

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oOptimizerManager::ImageUpdateStep()
{
  // Verbose
  if (m_verbose>=2) Cout("oOptimizerManager::ImageUpdateStep() -> Proceed to image update" << endl);
  // Update visited voxels
  if (mp_Optimizer->UpdateVisitedVoxels())
  {
    Cerr("***** oOptimizerManager::ImageUpdateStep() -> Problem while updating visited voxels !" << endl);
    return 1;
  }
  // Image update step
  if (mp_Optimizer->ImageUpdateStep())
  {
    Cerr("***** oOptimizerManager::ImageUpdateStep() -> Problem while updating image space !" << endl);
    return 1;
  }
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

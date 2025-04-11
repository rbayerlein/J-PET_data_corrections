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
  \ingroup datafile

  \brief Implementation of class oDynamicDataManager
*/

#include "oDynamicDataManager.hh"
#include "oImageDimensionsAndQuantification.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

oDynamicDataManager::oDynamicDataManager()
{
  mp_ID = NULL;
  m_verbose = -1;
  m_nbTimeFrames = -1;
  mp_currentFrameIndex = NULL;
  m_respGatingFlag = false;
  m_rMotionCorrFlag = false;
  m_nbRespGates = -1;
  m2p_nbEventsPerRespGate = NULL;
  m2p_indexLastEventRespGate = NULL;
  mp_currentRespGateIndex = NULL;
  m_cardGatingFlag = false;
  m_cMotionCorrFlag = false;
  m_nbCardGates = -1;
  m2p_nbEventsPerCardGate = NULL;
  m2p_indexLastEventCardGate = NULL;
  mp_currentCardGateIndex = NULL;
  m_gateDurationProvidedFlag = false;
  m2p_durationPerGate = NULL;
  m_pMotionCorrFlag = false;
  m_nbPMotionTriggers = -1;
  mp_listPMotionTriggers = NULL;
  mp_MotionTriggersIndexInFrame = NULL;
  m2p_listPMotionWeightInFrame = NULL;
  mp_frameNbPMotionIndexes = NULL;
  mp_framePMotionFirstIndex = NULL;
  mp_framePMotionLastIndex = NULL;
  mp_currentPMotionIndex = NULL;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

oDynamicDataManager::~oDynamicDataManager()
{
  for (int fr=0; fr<m_nbTimeFrames; fr++)
  {
    if (m2p_nbEventsPerRespGate && m2p_nbEventsPerRespGate[fr]) 
      delete m2p_nbEventsPerRespGate[fr];
    if (m2p_nbEventsPerCardGate && m2p_nbEventsPerCardGate[fr]) 
      delete m2p_nbEventsPerCardGate[fr];
    if (m2p_durationPerGate && m2p_durationPerGate[fr]) 
      delete m2p_durationPerGate[fr];
    if (m2p_indexLastEventRespGate && m2p_indexLastEventRespGate[fr]) 
      delete m2p_indexLastEventRespGate[fr];
    if (m2p_indexLastEventCardGate && m2p_indexLastEventCardGate[fr]) 
      delete m2p_indexLastEventCardGate[fr];
    if (m2p_listPMotionWeightInFrame && m2p_listPMotionWeightInFrame[fr])
      delete m2p_listPMotionWeightInFrame[fr];
  }

  if (m2p_nbEventsPerRespGate!=NULL)      delete m2p_nbEventsPerRespGate;
  if (m2p_nbEventsPerCardGate!=NULL)      delete m2p_nbEventsPerCardGate;
  if (m2p_durationPerGate!=NULL)          delete m2p_durationPerGate;
  if (m2p_indexLastEventRespGate!= NULL)  delete m2p_indexLastEventRespGate;
  if (m2p_indexLastEventCardGate!= NULL)  delete m2p_indexLastEventCardGate;
  if (mp_currentFrameIndex!=NULL)         delete mp_currentFrameIndex;
  if (mp_currentRespGateIndex!=NULL)      delete mp_currentRespGateIndex;
  if (mp_currentCardGateIndex!=NULL)      delete mp_currentCardGateIndex;
  if (mp_currentPMotionIndex!=NULL)       delete mp_currentPMotionIndex;
  if (mp_frameNbPMotionIndexes != NULL)    delete mp_frameNbPMotionIndexes;
  if (mp_framePMotionFirstIndex!=NULL)    delete mp_framePMotionFirstIndex;
  if (mp_framePMotionLastIndex!=NULL)     delete mp_framePMotionLastIndex;
  if (mp_listPMotionTriggers!=NULL)       delete mp_listPMotionTriggers;  
  if (m2p_listPMotionWeightInFrame!=NULL) delete[]m2p_listPMotionWeightInFrame;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oDynamicDataManager::InitDynamicData( int a_nbRespGates, int a_nbCardGates, const string& a_pathTo4DDataFile,
                                          int a_rmCorrFlag, int a_cmCorrFlag, int a_pmCorrFlag)
{
  // Verbose
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  if (m_verbose>=VERBOSE_DETAIL) Cout("oDynamicDataManager::InitDynamicData() -> Initialize dynamic data management" << endl);

  // Initialization
  m_nbTimeFrames = mp_ID->GetNbTimeFrames();
  m_nbRespGates = a_nbRespGates;
  m_nbCardGates = a_nbCardGates;
  m_nbPMotionTriggers = 0;
  
  if (m_nbRespGates > 1) m_respGatingFlag = true;
  if (m_nbCardGates > 1) m_cardGatingFlag = true;
  m_rMotionCorrFlag = a_rmCorrFlag;
  m_cMotionCorrFlag = a_cmCorrFlag;
  m_pMotionCorrFlag = a_pmCorrFlag;
  
  // Instanciation & initialization for current indices (used during the loop on events to know which frame/gate the event belongs to)
  mp_currentFrameIndex     = new int[mp_ID->GetNbThreadsForProjection()];
  mp_currentRespGateIndex  = new int[mp_ID->GetNbThreadsForProjection()];
  mp_currentCardGateIndex  = new int[mp_ID->GetNbThreadsForProjection()];
  mp_currentPMotionIndex = new int[mp_ID->GetNbThreadsForProjection()];
  for (int th=0 ; th<mp_ID->GetNbThreadsForProjection() ; th++)
  {
    mp_currentFrameIndex    [0] = 0;
    mp_currentRespGateIndex [0] = 0;
    mp_currentCardGateIndex [0] = 0;
    mp_currentPMotionIndex[0] = 0;
  }
  
  // Instanciation & initialization of number of events per frame for gated reconstructions
  m2p_nbEventsPerRespGate = new int64_t*[m_nbTimeFrames];
  m2p_nbEventsPerCardGate = new int64_t*[m_nbTimeFrames];
  m2p_durationPerGate     = new HPFLTNB*[m_nbTimeFrames];
  m2p_indexLastEventRespGate = new int64_t*[m_nbTimeFrames];
  m2p_indexLastEventCardGate = new int64_t*[m_nbTimeFrames];
  
  for (int fr=0; fr<m_nbTimeFrames; fr++)
  {
    m2p_nbEventsPerRespGate[fr] = new int64_t[m_nbRespGates];
    m2p_nbEventsPerCardGate[fr] = new int64_t[m_nbRespGates*m_nbCardGates];
    m2p_durationPerGate[fr]     = new HPFLTNB[m_nbRespGates*m_nbCardGates];
    m2p_indexLastEventRespGate[fr] = new int64_t[m_nbRespGates];
    m2p_indexLastEventCardGate[fr] = new int64_t[m_nbRespGates*m_nbCardGates];
  
    for (int rg=0; rg<m_nbRespGates; rg++)
    {
      m2p_nbEventsPerRespGate[fr][rg] = -1;
      m2p_indexLastEventRespGate[fr][rg] = -1;
      
      for (int cg=0; cg<m_nbCardGates; cg++)
      {
        m2p_durationPerGate[ fr ][ rg*m_nbCardGates + cg ] = -1;
        m2p_nbEventsPerCardGate[ fr ][ rg*m_nbCardGates + cg ] = -1;
        m2p_indexLastEventCardGate[ fr ][ rg*m_nbCardGates + cg ] = -1;
      }
    }
  }

  // Check coherence of motion initialization
  // Throw error if cardiac gating is enabled alongside respiratory motion correction (and dual-gating motion correction disabled)
  // TODO : warning in the documentation that Respiratory motion correction is not supported for cardiac-gated data in the current implementation
  if (m_rMotionCorrFlag && m_cardGatingFlag && !m_cMotionCorrFlag)
  {
    Cerr("***** oDynamicDataManager::InitDynamicData() -> Respiratory motion correction is enabled for cardiac-gated data. This is not supported in the current implementation  !" << endl);
    return 1;
  }

  // Check iPat motion vs physiological motion
  if (a_pmCorrFlag && (m_rMotionCorrFlag || m_cMotionCorrFlag) )
  {
    Cerr("***** oDynamicDataManager::InitDynamicData() -> Involuntary patient image-based motion correction cannot be used along with respiratory or cardiac motion correction  !" << endl);
    return 1;
  }
  
  // Initialization for respiratory/cardiac gated reconstruction
  // Note : for Analytic Projection, gating flag could be enabled in order to project an image containing several gates, but no description file is required
  //        InitDynamicDataGating() is restricted to reconstruction
  //        Check on the existence of a_pathTo4DDataFile during reconstruction is performed onnly in Grecon
  if ((m_respGatingFlag || m_cardGatingFlag) && !a_pathTo4DDataFile.empty() )
  {
    if(m_verbose>=VERBOSE_DETAIL) Cout("oDynamicDataManager::InitDynamicData() -> Initializing data for gated reconstruction... " << endl);
      
    if (InitDynamicDataGating(a_pathTo4DDataFile))
    {
      Cerr("***** oDynamicDataManager::InitDynamicData() -> A problem occurred during the dynamic gating initialization !" << endl);
      return 1;
    }
  }

  // Initialization with involuntary patient motion corrected reconstruction
  if (m_pMotionCorrFlag && InitDynamicDataPatientMotion(a_pathTo4DDataFile) )
  {
    Cerr("***** oDynamicDataManager::InitDynamicData() -> A problem occurred during the patient involuntary motion correction initialization !" << endl);
    return 1;
  }

  // Some feedback
  if (m_verbose>=VERBOSE_DETAIL)
  {
    if (m_respGatingFlag)  Cout("oDynamicDataManager::InitDynamicData() -> Respiratory gating enabled" << endl);
    if (m_rMotionCorrFlag) Cout("oDynamicDataManager::InitDynamicData() -> Respiratory image-based motion correction enabled" << endl);
    if (m_cardGatingFlag)  Cout("oDynamicDataManager::InitDynamicData() -> Cardiac gating enabled" << endl);
    if (m_cMotionCorrFlag) Cout("oDynamicDataManager::InitDynamicData() -> Cardiac image-based motion correction enabled" << endl);
    if (m_pMotionCorrFlag) Cout("oDynamicDataManager::InitDynamicData() -> Involuntary image-based patient motion correction enabled" << endl);
  }

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oDynamicDataManager::InitDynamicDataGating(const string& a_pathToFile)
{
  // Verbose
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  if (m_verbose>=VERBOSE_DETAIL) Cout("oDynamicDataManager::InitDynamicDataGating() ... " << endl);
    
  // Respiratory gating enabled
  if (m_respGatingFlag)
  {
    // Read information about respiratory gating (number of events per respiratory gate per frame)
    if ( ReadDataASCIIFile(a_pathToFile, 
                          "nb_events_respiratory_gates", 
                           m2p_nbEventsPerRespGate, 
                           m_nbRespGates, 
                           m_nbTimeFrames, 
                           KEYWORD_MANDATORY) )
    {
      Cerr("***** oDynamicDataManager::InitDynamicDataGating() -> Error while trying to read 'nb_events_respiratory_gates' in file '" << a_pathToFile << "' !" << endl);
      return 1;
    }
    
    // Get the index of the last event of each respiratory gate using the number of events inside each gate
    uint64_t event_number_sum = 0;
    if (m_verbose>=VERBOSE_DETAIL) Cout("oDynamicDataManager::InitDynamicDataGating() :" << endl);
            
    for (int fr=0; fr<m_nbTimeFrames; fr++)
      for (int rg=0; rg<m_nbRespGates; rg++)
      { 
        m2p_indexLastEventRespGate[fr][rg] += m2p_nbEventsPerRespGate[fr][rg] + event_number_sum;
        event_number_sum += m2p_nbEventsPerRespGate[fr][rg]; 
        if (m_verbose>=VERBOSE_DETAIL) Cout("  --> Number of events for frame #" << fr << ", respiratory gate #" << rg << " = " << m2p_nbEventsPerRespGate[fr][rg] << endl);
      } 
  }

  // Cardiac gating enabled
  if (m_cardGatingFlag)
  {
    // Read information about cardiac gating (number of events per cardiac gate per respiratory gates times frames)
    if ( ReadDataASCIIFile(a_pathToFile, 
                          "nb_events_cardiac_gates", 
                           m2p_nbEventsPerCardGate,
                           m_nbRespGates*m_nbCardGates, 
                           m_nbTimeFrames, 
                           KEYWORD_MANDATORY) )
    {
      Cerr("***** oDynamicDataManager::InitDynamicDataGating() -> Error while trying to read 'nb_events_cardiac_gates' in file '" << a_pathToFile << "' !" << endl);
      return 1;
    }
    // Get the index of the last event of each cardiac gate using the number of events inside each gate
    uint64_t event_number_sum = 0;
    if (m_verbose>=VERBOSE_DETAIL) Cout("oDynamicDataManager::InitDynamicDataGating() :" << endl);
        
    for (int fr=0; fr<m_nbTimeFrames; fr++)
      for (int g=0; g<m_nbRespGates*m_nbCardGates; g++)
      { 
        m2p_indexLastEventCardGate[fr][g] += m2p_nbEventsPerCardGate[fr][g] + event_number_sum;
        event_number_sum += m2p_nbEventsPerCardGate[fr][g]; 
        if (m_verbose>=VERBOSE_DETAIL) Cout("  --> Number of events for frame #" << fr << ", cardiac gate #" << g << " = " << m2p_nbEventsPerCardGate[fr][g] << endl);
      } 
  }
      
  // Read optional information about respiratory gate duration (if provided)
  if ( ReadDataASCIIFile(a_pathToFile, 
                        "duration_gate", 
                         m2p_durationPerGate, 
                         m_nbRespGates*m_nbCardGates, 
                         m_nbTimeFrames, 
                         KEYWORD_OPTIONAL) == KEYWORD_OPTIONAL_ERROR )
  {
    Cerr("***** oDynamicDataManager::InitDynamicDataGating() -> Error while trying to read 'nb_events_respiratory_gates' in file '" << a_pathToFile << "' !" << endl);
    return 1;
  }
  
  // Flag indicating that gate duration have been provided
  if (m2p_durationPerGate[0][0]>0) 
  { 
    m_gateDurationProvidedFlag = true;
    if (m_verbose>=VERBOSE_DETAIL)
      for (int fr=0; fr<m_nbTimeFrames; fr++)
        for (int g=0; g<m_nbRespGates*m_nbCardGates; g++)
          Cout("  --> Provided gate duration, frame #" << fr << ",  gate #" << g << ": " << m2p_durationPerGate[fr][g] << endl);

  }
      
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oDynamicDataManager::InitDynamicDataPatientMotion(const string& a_pathToFile)
{
  // Verbose
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  if (m_verbose>=VERBOSE_DETAIL) Cout("oDynamicDataManager::InitDynamicDataPatientMotion() ... " << endl);

  // First get the number of time triggers into the dynamic file
  if (ReadDataASCIIFile(a_pathToFile, 
                       "nb_motion_triggers", 
                        &m_nbPMotionTriggers, 
                        1, 
                        KEYWORD_MANDATORY))
    {
      Cerr("***** oDynamicDataManager::InitDynamicDataPatientMotion() -> Involuntary motion correction enabled but no dynamic configuration file has been provided with the -im option," << endl
        << "*****                                                        or the number of triggers could not be found in the dynamic file '" << a_pathToFile << "' !" << endl);
      return 1;
    }
  // Allocate time triggers and read them from the dynamic file
  mp_listPMotionTriggers = (uint32_t*)malloc(m_nbPMotionTriggers*sizeof(uint32_t));
  mp_frameNbPMotionIndexes = (uint16_t*)malloc(m_nbTimeFrames*sizeof(uint16_t));
  mp_framePMotionFirstIndex = (uint16_t*)malloc(m_nbTimeFrames*sizeof(uint16_t));
  mp_framePMotionLastIndex = (uint16_t*)malloc(m_nbTimeFrames*sizeof(uint16_t));
  mp_MotionTriggersIndexInFrame = (int**) malloc(m_nbTimeFrames * sizeof(int*));
  m2p_listPMotionWeightInFrame = (HPFLTNB**)malloc(m_nbTimeFrames*sizeof(HPFLTNB*));
  
  
  for(int fr=0 ; fr<m_nbTimeFrames ; fr++)
  {
    // mp_frameNbPMotionIndexes initialised with 1
    mp_frameNbPMotionIndexes[fr] = 1;
    mp_framePMotionLastIndex[fr] = 0;
    mp_framePMotionFirstIndex[fr] = 0;
  }
  
  // Get timestamp of motion trigger
  if (ReadDataASCIIFile(a_pathToFile, 
                      "timestamp_motion_triggers", 
                        mp_listPMotionTriggers, 
                        m_nbPMotionTriggers, 
                        KEYWORD_MANDATORY))
  {
    Cerr("***** oDynamicDataManager::InitDynamicDataPatientMotion() -> Involuntary motion correction enabled but list of triggers"
      << " not found in the dynamic file '" << a_pathToFile << "' !" << endl);
    return 1;
  }

  // Convert motion trigger in ms as for timestamp
  for (int pm=0 ; pm<m_nbPMotionTriggers ; pm++)
    mp_listPMotionTriggers [ pm ] *= 1000.;
    

  // remove motion triggers taking place after the end of the last frame
  // scan from the end to the first motion trigger and adapt total number of triggers accordingly
  // we assume that the tags are provided in ascending order ( as required )
  int new_m_nbPMotionTriggers=m_nbPMotionTriggers;
  for (int pm=(m_nbPMotionTriggers-1) ; pm>0 ; pm--)
    if (mp_listPMotionTriggers[pm] >= mp_ID->GetFrameTimeStopInMs(0, m_nbTimeFrames-1))
      new_m_nbPMotionTriggers--;
  m_nbPMotionTriggers = new_m_nbPMotionTriggers;


  // --- Recover the number of transformation involved during each frame --- //
  //     The following lines loop on the motion trigger in order to recover,
  //     for each frame, the first and last patient motion trigger 
  
  // First check that the first patient motion trigger is 0s (référence)
  if(mp_listPMotionTriggers[ 0 ] > 0.)
  {
    Cerr("***** oDynamicDataManager::InitDynamicDataPatientMotion() -> "
       <<" Error, the first motion trigger must be 0 ms"
       <<" (the corresponding transformation parameters are usually set to 0,"
       <<" as the acquisition start in the reference position!" << endl);
    return 1;
  }

  // First we need to find the number and index of motion triggers within frames ( for intra-frame motion )
  //
  // Loop over all frames
  for (int fr = 0; fr < m_nbTimeFrames; fr++)
  {
    // Loop on all motion triggers (MT) and get the number motion index within each frame
    for (int t = 0; t < m_nbPMotionTriggers; t++)
    {
      if (mp_ID->GetFrameTimeStartInMs(0, fr) < mp_listPMotionTriggers[t] && mp_listPMotionTriggers[t]  < mp_ID->GetFrameTimeStopInMs(0, fr))
        mp_frameNbPMotionIndexes[fr]++;
    }
    // allocate memory for the motion trigger indexes in each frame for mp_frameNbPMotionIndexes
    mp_MotionTriggersIndexInFrame[fr] = (int*) malloc((mp_frameNbPMotionIndexes[fr] - 1) * sizeof(int));
    // Populate each mp_MotionTriggersIndexInFrame with the motion trigger indexes
    int MT_in_frame = 0;
    for (int t = 0; t < m_nbPMotionTriggers; t++)
    {
      if (mp_ID->GetFrameTimeStartInMs(0, fr) < mp_listPMotionTriggers[t] && mp_listPMotionTriggers[t]  < mp_ID->GetFrameTimeStopInMs(0, fr))
      {
        mp_MotionTriggersIndexInFrame[fr][MT_in_frame] = t;
        MT_in_frame++;
      }
    }
  }

  // Then we set the first and last motion index for each frame
  //  in this loop we also need to consider motion tags that coincide with frame start times.
  //
  // Loop over all frames
  for(int fr=0 ; fr<m_nbTimeFrames ; fr++)
  {
    // Loop on all motion triggers (MT) and get the first and last motion index for each frame
    for (int t = 0; t < m_nbPMotionTriggers; t++)
    {
      // If motion trigger t <= starting time of the current frame, we set it as the first motion index
      if (mp_listPMotionTriggers[t] <= mp_ID->GetFrameTimeStartInMs(0, fr))
      {
        mp_framePMotionFirstIndex[fr] = t;
        mp_framePMotionLastIndex[fr] = t;
      }
      // if higher than frame start time and less that stop time, we update the last motion trigger
      else if (mp_listPMotionTriggers[t] < mp_ID->GetFrameTimeStopInMs(0, fr))
        mp_framePMotionLastIndex[fr] = t;
      // if higher than frame stop time we break and process the next frame
      else if (mp_listPMotionTriggers[t] >= mp_ID->GetFrameTimeStopInMs(0, fr))
        break;
    }
  }

  // Recover number of MT for each frame
  if (m_verbose>=VERBOSE_DETAIL) Cout("oDynamicDataManager::InitDynamicDataPatientMotion() -> First / Last motion trigger index for each dynamic frame: " << endl);
  for(int fr=0 ; fr<m_nbTimeFrames ; fr++)
  {
    mp_frameNbPMotionIndexes[fr] = mp_framePMotionLastIndex[fr] - mp_framePMotionFirstIndex[fr] + 1;
    if (m_verbose>=VERBOSE_DETAIL) Cout("Frm" << fr+1 << " : " << mp_framePMotionFirstIndex[fr]  << " / " << mp_framePMotionLastIndex[fr] << endl);
  }
  
  // ===================================================================
  // Compute the weight of each pm sub-frame within each frame
  // (ratio subset duration/total frame duration)
  // Adjust the decay quantification factors for each taking into account each sub-frame
  // This is only used for list-mode sensitivity image generation where we need to apply the factor to each sub-frame
  //   including application of decay effects

  for(int fr=0 ; fr<m_nbTimeFrames ; fr++)
  {
    m2p_listPMotionWeightInFrame[fr] = (HPFLTNB *) malloc(mp_frameNbPMotionIndexes[fr] * sizeof(HPFLTNB));
    // Initialise sub-frame weights with one
    for (int ipm = 0; ipm < mp_frameNbPMotionIndexes[fr]; ipm++)
      m2p_listPMotionWeightInFrame[fr][ipm]=1;
  }

  // Loop on frames, recover subset start/stop and compute weights
  for(int fr=0 ; fr<m_nbTimeFrames ; fr++)
  {
    if (mp_frameNbPMotionIndexes[fr] > 1)
    {
      // for first sub-frame set start time equal to frame start time
      uint32_t motion_subset_start_time_ms = mp_ID->GetFrameTimeStartInMs(0, fr);
      uint32_t motion_subset_stop_time_ms;
      // loop over motion triggers in frame
      for (int ipm = 0; ipm < mp_frameNbPMotionIndexes[fr]-1; ipm++)
      {
        // set the subset time stop from the stored triggers in Frame
        motion_subset_stop_time_ms = mp_listPMotionTriggers[ mp_MotionTriggersIndexInFrame[fr][ipm] ] ;
        // calculate the weight based on duration
        m2p_listPMotionWeightInFrame[fr][ipm] =
                FLTNB(motion_subset_stop_time_ms - motion_subset_start_time_ms) / mp_ID->GetFrameDurationInMs(0, fr);

        // update motion_subset_start_time_ms for next sub-frame
        motion_subset_start_time_ms = mp_listPMotionTriggers[ mp_MotionTriggersIndexInFrame[fr][ipm] ] ;
      }
      // Treat the last sub-frame in this frame
      motion_subset_stop_time_ms = mp_ID->GetFrameTimeStopInMs(0, fr);
      m2p_listPMotionWeightInFrame[fr][mp_frameNbPMotionIndexes[fr]-1] =
              FLTNB(motion_subset_stop_time_ms - motion_subset_start_time_ms) / mp_ID->GetFrameDurationInMs(0, fr);
    }
  }

  // Check if an isotope has been provided for the reconstruction
  // Loop on frames, and uncorrect for application of frame decay and duration
  if (mp_ID->GetLambda() > 0)
  {
    for (int fr = 0; fr < m_nbTimeFrames; fr++)
      if (mp_frameNbPMotionIndexes[fr] > 1)
      {
        long double lambda = mp_ID->GetLambda();
        long double dstart = lambda * mp_ID->GetFrameTimeStartInSec(0, fr);
        long double dacq = lambda * mp_ID->GetFrameDurationInSec(0, fr);
        for (int ipm = 0; ipm < mp_frameNbPMotionIndexes[fr]; ipm++)
        {
          // Uncorrect for time start decay ( by multiplying with the decay correction factor )
          m2p_listPMotionWeightInFrame[fr][ipm] *= exp(dstart);
          // Uncorrect for frame duration decay ( by multiplying with the decay duration factor )
          m2p_listPMotionWeightInFrame[fr][ipm] *= dacq / (1.0 - exp(-dacq));
        }
      }

    // Loop on frames and sub-frames, recover subset start/stop and apply sub-frame decay factors
    for (int fr = 0; fr < m_nbTimeFrames; fr++)
      if (mp_frameNbPMotionIndexes[fr] > 1)
      {
        uint32_t motion_subset_start_time_ms = mp_ID->GetFrameTimeStartInMs(0, fr);
        uint32_t motion_subset_stop_time_ms;
        long double lambda = mp_ID->GetLambda();
        long double dstart ;
        long double dacq ;
        for (int ipm = 0; ipm < mp_frameNbPMotionIndexes[fr] - 1; ipm++)
        {
          // set the subset time stop from the stored triggers in Frame
          motion_subset_stop_time_ms = mp_listPMotionTriggers[mp_MotionTriggersIndexInFrame[fr][ipm]];
          // calculate the decay factors for the subframe
          dstart = lambda
                   * (long double) motion_subset_start_time_ms
                   * 0.001;
          dacq = lambda
                 * (long double) (motion_subset_stop_time_ms - motion_subset_start_time_ms)
                 * 0.001;
          // Then, apply the correct sub-frame factor the pm subset
          // Time start decay factor
          m2p_listPMotionWeightInFrame[fr][ipm] /= exp(dstart);
          // Frame duration decay factor
          m2p_listPMotionWeightInFrame[fr][ipm] /= dacq / (1.0 - exp(-dacq));
          // update motion_subset_start_time_ms for next sub-frame
          motion_subset_start_time_ms = mp_listPMotionTriggers[mp_MotionTriggersIndexInFrame[fr][ipm]];
        }
        // Treat the last sub-frame in this frame
        motion_subset_stop_time_ms = mp_ID->GetFrameTimeStopInMs(0, fr);
        dstart = lambda
               * (long double)motion_subset_start_time_ms
               * 0.001;
        dacq = lambda
             * (long double)(motion_subset_stop_time_ms - motion_subset_start_time_ms)
             * 0.001 ;
        
        // Time start decay correction
        m2p_listPMotionWeightInFrame[fr][mp_frameNbPMotionIndexes[fr]-1] /= exp(dstart);
        // Frame duration decay correction
        m2p_listPMotionWeightInFrame[fr][mp_frameNbPMotionIndexes[fr]-1] /= dacq / (1.0 - exp(-dacq));

      }
  }

  if (m_verbose>=VERBOSE_DEBUG_NORMAL)
    for(int fr=0 ; fr<m_nbTimeFrames ; fr++)
      for( int ipm=0 ; ipm<mp_frameNbPMotionIndexes[fr] ; ipm++ )
        cout << "m2p_listPMotionWeightInFrame[" <<fr+1 << "] [" <<ipm+1 << "] : " << m2p_listPMotionWeightInFrame[fr][ipm] << endl;



  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oDynamicDataManager::SetDynamicSpecificQuantificationFactors(FLTNB** a2p_quantificationFactors)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  if (m_verbose>=VERBOSE_NORMAL) Cout("oDynamicDataManager::SetDynamicSpecificQuantificationFactors() ... " << endl);
    
  // COMPUTE GATE-SPECIFIC QUANTITATIVE FACTORS
  if (m_nbRespGates>1 || m_nbCardGates>1)
  {
    for(int fr=0 ; fr<m_nbTimeFrames ; fr++)
    {
      // We have cardiac gates and don't correct for motion -> quantification factors required for cardiac gates
      if (m_cardGatingFlag && !m_cMotionCorrFlag)
      {
        // If duration have been provided (value >0, then normalize according to the duration)
        if (m2p_durationPerGate[fr][0]>0)
        {
          for(int g=0 ; g<m_nbRespGates*m_nbCardGates ; g++)
            a2p_quantificationFactors[fr][g] /= m2p_durationPerGate[fr][g];
        }
        else // Default normalization according to the number of events
        {
          uint64_t total_events_in_frame = 0;
          for(int g=0 ; g<m_nbRespGates*m_nbCardGates ; g++) total_events_in_frame += m2p_nbEventsPerCardGate[fr][g];
  
          if (m_verbose>=VERBOSE_DETAIL) 
          
          for(int g=0 ; g<m_nbRespGates*m_nbCardGates ; g++)
            a2p_quantificationFactors[fr][g] *= (FLTNB)total_events_in_frame/(FLTNB)m2p_nbEventsPerCardGate[fr][g];
        }
        
        if (m_verbose>=VERBOSE_DETAIL) 
        {
          Cout("  --> Cardiac gating correction factors :" << endl);
          for(int g=0 ; g<m_nbRespGates*m_nbCardGates ; g++)
            Cout("      Frame #" << fr << ", cardiac gate #" << g << " = " << a2p_quantificationFactors[fr][g] << endl);
        }
      }
      
      // We have resp gates and don't correct for resp motion (and no independent cardiac gate reconstruction) -> quantification factors required for cardiac gates
      else if(m_respGatingFlag && !m_rMotionCorrFlag)
      {
        // If duration have been provided (value >0, then normalize according to the duration)
        if (m2p_durationPerGate[fr][0]>0)
        {
          for(int g=0 ; g<m_nbRespGates*m_nbCardGates ; g++)
          {
            // remove pre-applied frame duration form the quantification factors
            a2p_quantificationFactors[fr][g] *= ((FLTNB)(mp_ID->GetFrameDurationInSec(0,fr)));
            // use the gate duration as the duration of the acquisition for the quantification factor
            a2p_quantificationFactors[fr][g] /= m2p_durationPerGate[fr][g];
          }
        }
        else // Default normalization according to the number of events
        {
          uint64_t total_events_in_frame = 0;
          for(int rg=0 ; rg<m_nbRespGates ; rg++) total_events_in_frame += m2p_nbEventsPerRespGate[fr][rg];
  
          if (m_verbose>=VERBOSE_DETAIL) 
  
          for(int rg=0 ; rg<m_nbRespGates ; rg++)
            a2p_quantificationFactors[fr][rg] *= (FLTNB)total_events_in_frame/(FLTNB)m2p_nbEventsPerRespGate[fr][rg];
        }
        
        if (m_verbose>=VERBOSE_DETAIL)
        {
          Cout("  --> Respiratory gating correction factors :" << endl);
          for(int rg=0 ; rg<m_nbRespGates ; rg++)
            Cout("      Frame #" << fr << ", gate #" << rg << " = " << a2p_quantificationFactors[fr][rg] << endl);
        }
      }
    }
  }
  
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oDynamicDataManager::CheckParameters(int64_t a_nbEvents)
{
  // Verbose
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  if (m_verbose>=VERBOSE_DETAIL) Cout("oDynamicDataManager::CheckParameters() -> Check parameters for dynamic data settings" << endl);
  // Check image dimensions
  if (mp_ID==NULL)
  {
    Cerr("***** oDynamicDataManager::CheckParameters() -> No image dimensions provided !" << endl);
    return 1;
  }
  // Check time frames
  if (m_nbTimeFrames<0)
  {
    Cerr("***** oDynamicDataManager::CheckParameters() -> Wrong number of time frames !" << endl);
    return 1;
  }
  // Check resp gates
  if (m_nbRespGates<0)
  {
    Cerr("***** oDynamicDataManager::CheckParameters() -> Wrong number of respiratory gates !" << endl);
    return 1;
  }
  // Check card gates
  if (m_nbCardGates<0)
  {
    Cerr("***** oDynamicDataManager::CheckParameters() -> Wrong number of respiratory gates !" << endl);
    return 1;
  }
  // Check involuntary motion triggers
  if (m_nbPMotionTriggers<0)
  {
    Cerr("***** oDynamicDataManager::CheckParameters() -> Wrong number of involuntary motion subsets provided !" << endl);
    return 1;
  }
  // Check verbosity
  if (m_verbose<0)
  {
    Cerr("***** oDynamicDataManager::CheckParameters() -> Wrong verbosity level provided !" << endl);
    return 1;
  }  
  
  // Feedback regarding the enabled/disabled options for the reconstruction
  if (m_verbose>=2)
  {
    if (m_respGatingFlag) Cout("oDynamicDataManager::CheckParameters() -> Respiratory gating is enabled" << endl);
    if (m_rMotionCorrFlag) Cout("oDynamicDataManager::CheckParameters() -> Respiratory motion correction enabled" << endl);
    if (m_cardGatingFlag) Cout("oDynamicDataManager::CheckParameters() -> Cardiac gating is enabled" << endl);
    if (m_cMotionCorrFlag) Cout("oDynamicDataManager::CheckParameters() -> Cardiac motion correction is enabled" << endl);
    if (m_pMotionCorrFlag) Cout("oDynamicDataManager::CheckParameters() -> Involuntary motion correction is enabled" << endl);
  }

  // Check consistency between number of events in the datafile and the potential splitting of the dynamic data into respiratory or cardiac gates (if any gating is enabled)
  if (m_respGatingFlag)
  {
    int last_fr = m_nbTimeFrames-1;
    int last_rg = m_nbRespGates-1;
    if (m2p_indexLastEventRespGate[last_fr][last_rg]+1 != a_nbEvents)
    {
      Cerr("***** oDynamicDataManager::CheckParameters() -> Problem while checking consistency of dynamic data !" << endl
        << "                                                The number of events in the datafile (" << a_nbEvents
        << ") is different from the total number of events in respiratory gates (" << m2p_indexLastEventRespGate[last_fr][last_rg]+1 << ") as initialized in the gating file !" << endl);
      return 1;
    }
  }
  if (m_cardGatingFlag)
  {
    int last_fr = m_nbTimeFrames-1;
    int last_cg = (m_nbRespGates*m_nbCardGates) - 1;
    if (m2p_indexLastEventCardGate[last_fr][last_cg]+1 != a_nbEvents)
    {
      Cerr("***** oDynamicDataManager::CheckParameters() -> Problem while checking consistency of dynamic data !" << endl
        << "                                                The number of events in the datafile (" << a_nbEvents
        << ") is different to the total number of events in cardiac gates (" << m2p_indexLastEventCardGate[last_fr][last_cg]+1 << ") as initialized in the gating file !" << endl);
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

void oDynamicDataManager::ResetCurrentDynamicIndices()
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_LIGHT)
  // Reset indices for each thread (this is a thread-safe implementation)
  for (int th=0; th<mp_ID->GetNbThreadsForProjection(); th++)
  {
    // We initialize the frame index to -1 because at the beginning of the events loop, we don't know yet
    // if we are in the first frame (the user can reconstruct only a part of the whole datafile).
    mp_currentFrameIndex[th] = -1;
    mp_currentPMotionIndex[th] = 0;
    mp_currentRespGateIndex[th] = 0;
    mp_currentCardGateIndex[th] = 0;
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int oDynamicDataManager::DynamicSwitch(int64_t a_currentEventIndex, uint32_t a_currentTime, int a_bed, int a_th)
{
  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_EVENT)
  
  // The value that will be returned
  int return_value = DYNAMIC_SWITCH_NOTHING;

  // -------------------------------------------------------------------------------------------------
  // Step 1: involuntary motion management
  // -------------------------------------------------------------------------------------------------

  // Only do this if the involuntary motion correction is enabled (meaning we must proceed to a deformation)
  if (m_pMotionCorrFlag)
  {
    // Search if we passed one of the next motion triggers (starting at current index)
    for (int mt=mp_currentPMotionIndex[a_th]; mt<m_nbPMotionTriggers; mt++)
    {
      // If we passed this trigger, set the return value to DEFORMATION. However, we continue the loop
      // in the case where we passed multiple triggers.
      if (a_currentTime >= mp_listPMotionTriggers[ mt ] // Current timestamp is posterior to the trigger timestamp 
       &&            mt > mp_currentPMotionIndex[a_th] ) // The trigger timestamp changed (posterior to the current trigger) 
      {
        mp_currentPMotionIndex[a_th] = mt;
        
        #ifdef CASTOR_DEBUG
        if (m_verbose>=VERBOSE_DEBUG_EVENT) Cout("oDynamicDataManager::DynamicSwitch() -> Thread " << a_th << ", gate for patient motion correction switched to " << mp_currentPMotionIndex[a_th] << endl
                             << "                                        at event #" << a_currentEventIndex << ", timestamp = " << a_currentTime << endl);
        #endif
        return_value = DYNAMIC_SWITCH_DEFORMATION;
      }
    }
  }

  // -------------------------------------------------------------------------------------------------
  // Step 2: frame management
  // -------------------------------------------------------------------------------------------------

  // Special case if the frame number is still -1
  if (mp_currentFrameIndex[a_th]==-1)
  {
    // We check if we are before the first frame, in this case we directly return a CONTINUE signal to go to the next event
    if (a_currentTime<mp_ID->GetFrameTimeStartInMs(a_bed,0)) return DYNAMIC_SWITCH_CONTINUE;
    // Otherwise, we now are at least in the first frame (which will be updated right after this 'if' section)
    else mp_currentFrameIndex[a_th] = 0;
  }

  // A boolean to know later if the frame index has changed
  bool frame_has_changed = false;

  
  // Now we search for the first frame index in which the event belongs, starting from the current frame index. Note that we do that
  // this way (not simply incrementing the frame index) because we want to deal with the case where the next event managed by this
  // thread could possibly skip multiple frames at once.
  for (int fr=mp_currentFrameIndex[a_th]; fr<m_nbTimeFrames; fr++)
  {
    // If the current time is less than the time stop of the frame 'fr', then the event belongs to this frame
    if (a_currentTime<mp_ID->GetFrameTimeStopInMs(a_bed,fr))
    {
      // Test if the frame has actually changed
      if (mp_currentFrameIndex[a_th]!=fr)
      {
        // Set the boolean to true
        frame_has_changed = true;
        // Update the current frame index
        mp_currentFrameIndex[a_th] = fr;
      }
      break;
    }
    // Otherwise, if in the last frame, then it means that this event is outside of frame definitions.
    // In this case, we return a special value so that this event and the next ones will be ignored.
    else if (fr==m_nbTimeFrames-1)
    {
      return DYNAMIC_SWITCH_END_LOOP;
    }
  }

  #ifdef CASTOR_DEBUG
  if (frame_has_changed && m_verbose >=VERBOSE_DEBUG_EVENT)
    Cout("oDynamicDataManager::DynamicSwitch() -> Thread " << a_th << ", frame switched to " << mp_currentFrameIndex[a_th] << endl);
  #endif


  // -------------------------------------------------------------------------------------------------
  // Step 3: respiratory gating management
  // -------------------------------------------------------------------------------------------------

  // A boolean to know later if the respiratory index has changed
  bool resp_gate_has_changed = false;
  
  // Do this only if respiratory gating is enabled
  if (m_respGatingFlag)
  {
    // Test if the frame index has changed
    if (frame_has_changed)
    {
      // Reset the respiratory gate index
      mp_currentRespGateIndex[a_th] = 0;
    }
    
    // For this frame, search the first gate (from the current gate index) for which the provided event index is below the event-index-stop
    bool resp_gate_found = false;
    for (int rg=mp_currentRespGateIndex[a_th]; rg<m_nbRespGates; rg++)
    {
      // If the current event index is below the last event of this gate, then the event belongs to this gate
      // (We won't enter again in the if statement due to the flag setting to true)
      if (a_currentEventIndex<=m2p_indexLastEventRespGate[mp_currentFrameIndex[a_th]][rg] && resp_gate_found == false)
      {
        // Test if the gate has actually changed
        if (mp_currentRespGateIndex[a_th]!=rg)
        {
          // Update the current gate index
          mp_currentRespGateIndex[a_th] = rg;
          // Verbose
          #ifdef CASTOR_DEBUG
          if (m_verbose>=2) Cout("oDynamicDataManager::DynamicSwitch() -> Thread " << a_th << ", respiratory gate switch to " << mp_currentRespGateIndex[a_th] << endl
                              << "                                                   on event " << a_currentEventIndex << endl
                              << "                                                   current frame : " << mp_currentFrameIndex[a_th] << endl
                              << "                                                   current respiratory gate index " << mp_currentRespGateIndex[a_th] << endl);
          #endif
          // If motion correction is enabled, then we should return a DEFORMATION signal
          if (m_rMotionCorrFlag) return_value = DYNAMIC_SWITCH_DEFORMATION;
          
          // Set the boolean to true
          resp_gate_has_changed = true;
        }
        // Set the boolean to true
        resp_gate_found = true;
      }
    }
    
    
  }
    
  // -------------------------------------------------------------------------------------------------
  // Step 4: cardiac gating management
  // -------------------------------------------------------------------------------------------------

  // A boolean to know later if the respiratory index has changed
  //bool card_gate_has_changed = false;
  
  // Do this only if cardiac gating is enabled
  if (m_cardGatingFlag)
  {
    // Test if the frame or respiratory index have changed
    if (frame_has_changed || resp_gate_has_changed)
    {
      // Reset the cardiac gate index
      mp_currentCardGateIndex[a_th] = 0;
    }
    // For this frame and respiratory gate, search the first gate (from the current gate index) for which the provided event index is below the event-index-stop
    bool card_gate_found = false;
    for (int cg=mp_currentCardGateIndex[a_th]; cg<m_nbCardGates; cg++)
    {
      // If the current event index is below the event-index-stop of this gate, then the event belongs to this gate
      if (a_currentEventIndex<m2p_indexLastEventCardGate[mp_currentFrameIndex[a_th]][mp_currentRespGateIndex[a_th]*m_nbCardGates+cg]  && card_gate_found == false)
      {
        // Test if the gate has actually changed
        if (mp_currentCardGateIndex[a_th]!=cg)
        {
          // Update the current gate index
          mp_currentCardGateIndex[a_th] = cg;
          // Verbose
          #ifdef CASTOR_DEBUG
          if (m_verbose>=3) Cout("oDynamicDataManager::DynamicSwitch() -> Thread " << a_th << ", cardiac gate switch to " <<  mp_currentCardGateIndex[a_th] << endl
                              << "                                                   on event " << a_currentEventIndex << endl
                              << "                                                   current frame : " << mp_currentFrameIndex[a_th] << endl);
          #endif
          // If motion correction is enabled, then we should return a DEFORMATION signal
          if (m_cMotionCorrFlag) return_value = DYNAMIC_SWITCH_DEFORMATION;
          
          // Set the boolean to true
          //card_gate_has_changed = true;
        }
        // Set the boolean to true
        card_gate_found = true;
      }
    }
  }
  
  // -------------------------------------------------------------------------------------------------
  // Step 5: just return the value !
  // -------------------------------------------------------------------------------------------------

  // Return the status of the dynamic switch
  return return_value;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

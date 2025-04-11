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
  \brief    Implementation of class vOptimizer
*/

#include "vOptimizer.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

vOptimizer::vOptimizer()
{
  // Affect default values
  mp_ImageDimensionsAndQuantification = NULL;
  mp_ImageSpace = NULL;
  m_nbBackwardImages = 0;
  m_nbTOFBins = 0;
  m_initialValue = 0.;
  m_listmodeCompatibility = false;
  m_histogramCompatibility = false;
  m_emissionCompatibility = false;
  m_transmissionCompatibility = false;
  m2p_forwardValues = NULL;
  m3p_backwardValues = NULL;
  m_dataMode = MODE_UNKNOWN;
  m_dataType = TYPE_UNKNOWN;
  m_optimizerFOMFlag = false;
  m_optimizerImageStatFlag = false;
  m4p_FOMLogLikelihood = NULL;
  m4p_FOMNbBins = NULL;
  m4p_FOMRMSE = NULL;
  m4p_FOMNbData = NULL;
  m4p_FOMPenalty = NULL;
  m_verbose = 0;
  m2p_attenuationImage = NULL;
  mp_imageStatNbVox = NULL;
  mp_imageStatMin = NULL;
  mp_imageStatMax = NULL;
  mp_imageStatMean = NULL;
  mp_imageStatVariance = NULL;
  mp_correctionStatMean = NULL;
  mp_correctionStatVariance = NULL;
  m_nbIterations = -1;
  mp_nbSubsets = NULL;
  m_currentIteration = -1;
  m_currentSubset = -1;
  mp_Penalty = NULL;
  m_requiredPenaltyDerivativesOrder = -1;
  m_needGlobalSensitivity = false;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

vOptimizer::~vOptimizer()
{
  // Free forward and backward buffers
  if (m2p_forwardValues)
  {
    for (int th=0; th<mp_ImageDimensionsAndQuantification->GetNbThreadsForProjection(); th++)
      if (m2p_forwardValues[th]) free(m2p_forwardValues[th]);
    free(m2p_forwardValues);
  }
  if (m3p_backwardValues)
  {
    for (int th=0; th<mp_ImageDimensionsAndQuantification->GetNbThreadsMax(); th++)
    {
      if (m3p_backwardValues[th])
      {
        for (int tof=0; tof<m_nbTOFBins; tof++) if (m3p_backwardValues[th][tof]) free(m3p_backwardValues[th][tof]);
        free(m3p_backwardValues[th]);
      }
    }
    free(m3p_backwardValues);
  }
  // Free pointers to attenuation images for SPECT
  if (m2p_attenuationImage) free(m2p_attenuationImage);
  // Delete all image update statistics tables
  if (m_optimizerImageStatFlag)
  {
    if (mp_imageStatNbVox) free(mp_imageStatNbVox);
    if (mp_imageStatMin) free(mp_imageStatMin);
    if (mp_imageStatMax) free(mp_imageStatMax);
    if (mp_imageStatMean) free(mp_imageStatMean);
    if (mp_imageStatVariance) free(mp_imageStatVariance);
    if (mp_correctionStatMean) free(mp_correctionStatMean);
    if (mp_correctionStatVariance) free(mp_correctionStatVariance);
  }
  // Delete FOM tables
  if (m_optimizerFOMFlag)
  {
    // Delete the log-likelihood table
    if (m4p_FOMLogLikelihood)
    {
      for (int fr=0; fr<mp_ImageDimensionsAndQuantification->GetNbTimeFrames(); fr++)
      {
        if (m4p_FOMLogLikelihood[fr])
        {
          for (int rg=0; rg<mp_ImageDimensionsAndQuantification->GetNbRespGates(); rg++)
          {
            if (m4p_FOMLogLikelihood[fr][rg])
            {
              for (int cg=0; cg<mp_ImageDimensionsAndQuantification->GetNbCardGates(); cg++)
                if (m4p_FOMLogLikelihood[fr][rg][cg]) free(m4p_FOMLogLikelihood[fr][rg][cg]);
              free(m4p_FOMLogLikelihood[fr][rg]);
            }
          }
          free(m4p_FOMLogLikelihood[fr]);
        }
      }
      free(m4p_FOMLogLikelihood);
    }
    // Delete the RMSE table
    if (m4p_FOMRMSE)
    {
      for (int fr=0; fr<mp_ImageDimensionsAndQuantification->GetNbTimeFrames(); fr++)
      {
        if (m4p_FOMRMSE[fr])
        {
          for (int rg=0; rg<mp_ImageDimensionsAndQuantification->GetNbRespGates(); rg++)
          {
            if (m4p_FOMRMSE[fr][rg])
            {
              for (int cg=0; cg<mp_ImageDimensionsAndQuantification->GetNbCardGates(); cg++)
                if (m4p_FOMRMSE[fr][rg][cg]) free(m4p_FOMRMSE[fr][rg][cg]);
              free(m4p_FOMRMSE[fr][rg]);
            }
          }
          free(m4p_FOMRMSE[fr]);
        }
      }
      free(m4p_FOMRMSE);
    }
    // Delete the nbBins table
    if (m4p_FOMNbBins)
    {
      for (int fr=0; fr<mp_ImageDimensionsAndQuantification->GetNbTimeFrames(); fr++)
      {
        if (m4p_FOMNbBins[fr])
        {
          for (int rg=0; rg<mp_ImageDimensionsAndQuantification->GetNbRespGates(); rg++)
          {
            if (m4p_FOMNbBins[fr][rg])
            {
              for (int cg=0; cg<mp_ImageDimensionsAndQuantification->GetNbCardGates(); cg++)
                if (m4p_FOMNbBins[fr][rg][cg]) free(m4p_FOMNbBins[fr][rg][cg]);
              free(m4p_FOMNbBins[fr][rg]);
            }
          }
          free(m4p_FOMNbBins[fr]);
        }
      }
      free(m4p_FOMNbBins);
    }
    // Delete the nbData table
    if (m4p_FOMNbData)
    {
      for (int fr=0; fr<mp_ImageDimensionsAndQuantification->GetNbTimeFrames(); fr++)
      {
        if (m4p_FOMNbData[fr])
        {
          for (int rg=0; rg<mp_ImageDimensionsAndQuantification->GetNbRespGates(); rg++)
          {
            if (m4p_FOMNbData[fr][rg])
            {
              for (int cg=0; cg<mp_ImageDimensionsAndQuantification->GetNbCardGates(); cg++)
                if (m4p_FOMNbData[fr][rg][cg]) free(m4p_FOMNbData[fr][rg][cg]);
              free(m4p_FOMNbData[fr][rg]);
            }
          }
          free(m4p_FOMNbData[fr]);
        }
      }
      free(m4p_FOMNbData);
    }
    // Delete the penalty table
    if (m4p_FOMPenalty)
    {
      for (int fr=0; fr<mp_ImageDimensionsAndQuantification->GetNbTimeFrames(); fr++)
      {
        if (m4p_FOMPenalty[fr])
        {
          for (int rg=0; rg<mp_ImageDimensionsAndQuantification->GetNbRespGates(); rg++)
          {
            if (m4p_FOMPenalty[fr][rg])
            {
              for (int cg=0; cg<mp_ImageDimensionsAndQuantification->GetNbCardGates(); cg++)
                if (m4p_FOMPenalty[fr][rg][cg]) free(m4p_FOMPenalty[fr][rg][cg]);
              free(m4p_FOMPenalty[fr][rg]);
            }
          }
          free(m4p_FOMPenalty[fr]);
        }
      }
      free(m4p_FOMPenalty);
    }
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void vOptimizer::ShowHelp()
{
  // Call the specific help function from the children
  ShowHelpSpecific();
  // Say if the child optimizer is compatible with histogram and/or list-mode data
  if (m_listmodeCompatibility && m_histogramCompatibility) cout << "This optimizer is compatible with both histogram and list-mode data." << endl;
  else if (m_listmodeCompatibility) cout << "This optimizer is only compatible with list-mode data." << endl;
  else if (m_histogramCompatibility) cout << "This optimizer is only compatible with histogram data." << endl;
  else cout << "This optimizer is a total non-sense, not compatible with histogram nor list-mode data !!!" << endl;
  // Say if the child optimizer is compatible with emission and/or transmission data
  if (m_emissionCompatibility && m_transmissionCompatibility) cout << "This optimizer is compatible with both emission and transmission data." << endl;
  else if (m_emissionCompatibility) cout << "This optimizer is only compatible with emission data." << endl;
  else if (m_transmissionCompatibility) cout << "This optimizer is only compatible with transmission data." << endl;
  else cout << "This optimizer is a total non-sense, not compatible with emission nor transmission data !!!" << endl;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vOptimizer::CheckParameters()
{
  // Check image dimensions
  if (mp_ImageDimensionsAndQuantification==NULL)
  {
    Cerr("***** vOptimizer::CheckParameters() -> oImageDimensionsAndQuantification is null !" << endl);
    return 1;
  }
  // Check image space
  if (mp_ImageSpace==NULL)
  {
    Cerr("***** vOptimizer::CheckParameters() -> oImageSpace is null !" << endl);
    return 1;
  }
  // Check data mode
  if (m_dataMode==MODE_UNKNOWN || (m_dataMode!=MODE_LIST && m_dataMode!=MODE_HISTOGRAM))
  {
    Cerr("***** vOptimizer::CheckParameters() -> No or meaningless data mode provided !" << endl);
    return 1;
  }
  // Check data type
  if (m_dataType==TYPE_UNKNOWN || (m_dataType!=TYPE_PET && m_dataType!=TYPE_SPECT && m_dataType!=TYPE_CT))
  {
    Cerr("***** vOptimizer::CheckParameters() -> No or meaningless data type provided !" << endl);
    return 1;
  }
  // Check data spec
  if (m_dataSpec==SPEC_UNKNOWN || (m_dataSpec!=SPEC_EMISSION && m_dataSpec!=SPEC_TRANSMISSION))
  {
    Cerr("***** vOptimizer::CheckParameters() -> No or meaningless data specificity provided (emission or transmission) !" << endl);
    return 1;
  }
  // Check verbose level
  if (m_verbose<0)
  {
    Cerr("***** vOptimizer::CheckParameters() -> Verbose level is negative !" << endl);
    return 1;
  }
  // Listmode incompatibility
  if (m_dataMode==MODE_LIST && !m_listmodeCompatibility)
  {
    Cerr("***** vOptimizer::CheckParameters() -> The selected optimizer is not compatible with listmode data !" << endl);
    return 1;
  }
  // Histogram incompatibility
  if (m_dataMode==MODE_HISTOGRAM && !m_histogramCompatibility)
  {
    Cerr("***** vOptimizer::CheckParameters() -> The selected optimizer is not compatible with histogram data !" << endl);
    return 1;
  }
  // Emission compatability
  if (m_dataSpec==SPEC_EMISSION && !m_emissionCompatibility)
  {
    Cerr("***** vOptimizer::CheckParameters() -> The selected optimizer is not compatible with emission data !" << endl);
    return 1;
  }
  // Transmission compatibility
  if (m_dataSpec==SPEC_TRANSMISSION && !m_transmissionCompatibility)
  {
    Cerr("***** vOptimizer::CheckParameters() -> The selected optimizer is not compatible with transmission data !" << endl);
    return 1;
  }
  // Send an error if FOM computation is required and verbose is null
  if (m_optimizerFOMFlag && m_verbose==0)
  {
    Cerr("***** vOptimizer::CheckParameters() -> Asked to compute FOMs in the data-space whereas the verbose is not set !" << endl);
    return 1;
  }
  // Send an error if image statistics computation is required and verbose is null
  if (m_optimizerImageStatFlag && m_verbose==0)
  {
    Cerr("***** vOptimizer::CheckParameters() -> Asked to compute image update statistics whereas the verbose is not set !" << endl);
    return 1;
  }
  // Send an error if FOM computation with list-mode data (it has no sense)
  if (m_optimizerFOMFlag && m_dataMode==MODE_LIST)
  {
    Cerr("***** vOptimizer::CheckParameters() -> Computing FOMs in the data-space with list-mode data is not possible !" << endl);
    return 1;
  }
  // Call the CheckSpecificParameters function of the child
  if (CheckSpecificParameters())
  {
    Cerr("***** vOptimizer::CheckParameters() -> A problem occurred while checking parameters specific to the optimizer module !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vOptimizer::Initialize()
{
  // Allocate forward buffers (as many as threads and as many as TOF bins)
  m2p_forwardValues = (FLTNB**)malloc(mp_ImageDimensionsAndQuantification->GetNbThreadsForProjection()*sizeof(FLTNB*));
  for (int th=0; th<mp_ImageDimensionsAndQuantification->GetNbThreadsForProjection(); th++)
    m2p_forwardValues[th] = (FLTNB*)malloc(m_nbTOFBins*sizeof(FLTNB));

  // Allocate backward buffers (as many as threads then as many as TOF bins and as many as backward images)
  m3p_backwardValues = (FLTNB***)malloc(mp_ImageDimensionsAndQuantification->GetNbThreadsMax()*sizeof(FLTNB**));
  for (int th=0; th<mp_ImageDimensionsAndQuantification->GetNbThreadsMax(); th++)
  {
    m3p_backwardValues[th] = (FLTNB**)malloc(m_nbTOFBins*sizeof(FLTNB*));
    for (int tof=0; tof<m_nbTOFBins; tof++)
      m3p_backwardValues[th][tof] = (FLTNB*)malloc(m_nbBackwardImages*sizeof(FLTNB));
  }

  // Allocate pointers to the current attenuation image for SPECT attenuation correction
  m2p_attenuationImage = (FLTNB**)malloc(mp_ImageDimensionsAndQuantification->GetNbThreadsForProjection()*sizeof(FLTNB*));
  for (int th=0; th<mp_ImageDimensionsAndQuantification->GetNbThreadsForProjection(); th++) m2p_attenuationImage[th] = NULL;

  // Allocate the thread-safe tables for image update statistics computation
  if (m_optimizerImageStatFlag)
  {
    mp_imageStatNbVox = (INTNB*)malloc(mp_ImageDimensionsAndQuantification->GetNbThreadsForImageComputation()*sizeof(INTNB));
    mp_imageStatMin = (FLTNB*)malloc(mp_ImageDimensionsAndQuantification->GetNbThreadsForImageComputation()*sizeof(FLTNB));
    mp_imageStatMax = (FLTNB*)malloc(mp_ImageDimensionsAndQuantification->GetNbThreadsForImageComputation()*sizeof(FLTNB));
    mp_imageStatMean = (HPFLTNB*)malloc(mp_ImageDimensionsAndQuantification->GetNbThreadsForImageComputation()*sizeof(HPFLTNB));
    mp_imageStatVariance = (HPFLTNB*)malloc(mp_ImageDimensionsAndQuantification->GetNbThreadsForImageComputation()*sizeof(HPFLTNB));
    mp_correctionStatMean = (HPFLTNB*)malloc(mp_ImageDimensionsAndQuantification->GetNbThreadsForImageComputation()*sizeof(HPFLTNB));
    mp_correctionStatVariance = (HPFLTNB*)malloc(mp_ImageDimensionsAndQuantification->GetNbThreadsForImageComputation()*sizeof(HPFLTNB));
  }

  // We allocate the tables that hold the figures-of-merit if FOM flag is on
  if (m_optimizerFOMFlag)
  {
    m4p_FOMNbBins = (uint64_t****)malloc(mp_ImageDimensionsAndQuantification->GetNbTimeFrames()*sizeof(uint64_t***));
    m4p_FOMNbData = (HPFLTNB****)malloc(mp_ImageDimensionsAndQuantification->GetNbTimeFrames()*sizeof(HPFLTNB***));
    m4p_FOMLogLikelihood = (HPFLTNB****)malloc(mp_ImageDimensionsAndQuantification->GetNbTimeFrames()*sizeof(HPFLTNB***));
    m4p_FOMRMSE = (HPFLTNB****)malloc(mp_ImageDimensionsAndQuantification->GetNbTimeFrames()*sizeof(HPFLTNB***));
    if (mp_Penalty!=NULL) m4p_FOMPenalty = (HPFLTNB****)malloc(mp_ImageDimensionsAndQuantification->GetNbTimeFrames()*sizeof(HPFLTNB***));
    for (int fr=0; fr<mp_ImageDimensionsAndQuantification->GetNbTimeFrames(); fr++)
    {
      m4p_FOMNbBins[fr] = (uint64_t***)malloc(mp_ImageDimensionsAndQuantification->GetNbRespGates()*sizeof(uint64_t**));
      m4p_FOMNbData[fr] = (HPFLTNB***)malloc(mp_ImageDimensionsAndQuantification->GetNbRespGates()*sizeof(HPFLTNB**));
      m4p_FOMLogLikelihood[fr] = (HPFLTNB***)malloc(mp_ImageDimensionsAndQuantification->GetNbRespGates()*sizeof(HPFLTNB**));
      m4p_FOMRMSE[fr] = (HPFLTNB***)malloc(mp_ImageDimensionsAndQuantification->GetNbRespGates()*sizeof(HPFLTNB**));
      if (mp_Penalty!=NULL) m4p_FOMPenalty[fr] = (HPFLTNB***)malloc(mp_ImageDimensionsAndQuantification->GetNbRespGates()*sizeof(HPFLTNB**));
      for (int rg=0; rg<mp_ImageDimensionsAndQuantification->GetNbRespGates(); rg++)
      {
        m4p_FOMNbBins[fr][rg] = (uint64_t**)malloc(mp_ImageDimensionsAndQuantification->GetNbCardGates()*sizeof(uint64_t*));
        m4p_FOMNbData[fr][rg] = (HPFLTNB**)malloc(mp_ImageDimensionsAndQuantification->GetNbCardGates()*sizeof(HPFLTNB*));
        m4p_FOMLogLikelihood[fr][rg] = (HPFLTNB**)malloc(mp_ImageDimensionsAndQuantification->GetNbCardGates()*sizeof(HPFLTNB*));
        m4p_FOMRMSE[fr][rg] = (HPFLTNB**)malloc(mp_ImageDimensionsAndQuantification->GetNbCardGates()*sizeof(HPFLTNB*));
        if (mp_Penalty!=NULL) m4p_FOMPenalty[fr][rg] = (HPFLTNB**)malloc(mp_ImageDimensionsAndQuantification->GetNbCardGates()*sizeof(HPFLTNB*));
        for (int cg=0; cg<mp_ImageDimensionsAndQuantification->GetNbCardGates(); cg++)
        {
          m4p_FOMNbBins[fr][rg][cg] = (uint64_t*)malloc(mp_ImageDimensionsAndQuantification->GetNbThreadsForProjection()*sizeof(uint64_t));
          m4p_FOMNbData[fr][rg][cg] = (HPFLTNB*)malloc(mp_ImageDimensionsAndQuantification->GetNbThreadsForProjection()*sizeof(HPFLTNB));
          m4p_FOMLogLikelihood[fr][rg][cg] = (HPFLTNB*)malloc(mp_ImageDimensionsAndQuantification->GetNbThreadsForProjection()*sizeof(HPFLTNB));
          m4p_FOMRMSE[fr][rg][cg] = (HPFLTNB*)malloc(mp_ImageDimensionsAndQuantification->GetNbThreadsForProjection()*sizeof(HPFLTNB));
          if (mp_Penalty!=NULL) m4p_FOMPenalty[fr][rg][cg] = (HPFLTNB*)malloc(mp_ImageDimensionsAndQuantification->GetNbThreadsForImageComputation()*sizeof(HPFLTNB));
        }
      }
    }
  }

  // Call the specific initialization function of the child
  if (InitializeSpecific())
  {
    Cerr("***** vOptimizer::Initialize() -> A problem occurred while initializing stuff specific to the optimizer module !" << endl);
    return 1;
  }

  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vOptimizer::PreDataUpdateStep()
{
  // Simply reset FOMs if activated
  if (m_optimizerFOMFlag)
  {
    for (int fr=0; fr<mp_ImageDimensionsAndQuantification->GetNbTimeFrames(); fr++)
      for (int rg=0; rg<mp_ImageDimensionsAndQuantification->GetNbRespGates(); rg++)
        for (int cg=0; cg<mp_ImageDimensionsAndQuantification->GetNbCardGates(); cg++)
        {
          for (int th=0; th<mp_ImageDimensionsAndQuantification->GetNbThreadsForProjection(); th++)
          {
            m4p_FOMNbBins[fr][rg][cg][th] = 0;
            m4p_FOMNbData[fr][rg][cg][th] = 0.;
            m4p_FOMLogLikelihood[fr][rg][cg][th] = 0.;
            m4p_FOMRMSE[fr][rg][cg][th] = 0.;
          }
          if (mp_Penalty!=NULL) for (int th=0; th<mp_ImageDimensionsAndQuantification->GetNbThreadsForImageComputation(); th++)
          {
            m4p_FOMPenalty[fr][rg][cg][th] = 0.;
          }
        }
  }
  // Then call the specific pre-update-step function from the child optimizer
  if (PreDataUpdateSpecificStep())
  {
    Cerr("***** vOptimizer::PreDataUpdateStep() -> A problem occurred while applying the specific pre-data-update step of the optimizer module !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vOptimizer::PreDataUpdateSpecificStep()
{
  // By default, just do nothing here. This function is designed to be overloaded by the specific optimizer if needed
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vOptimizer::PreImageUpdateStep()
{
  // Then call the specific post-update-step function from the child optimizer
  if (PreImageUpdateSpecificStep())
  {
    Cerr("***** vOptimizer::PreImageUpdateStep() -> A problem occurred while applying the specific post-data-update step of the optimizer module !" << endl);
    return 1;
  }
  // Compute optimizer figures-of-merit (data fidelity term and penalty)
  if (m_optimizerFOMFlag)
  {
    // Compute data FOM related to the data fildelity term
    for (int fr=0; fr<mp_ImageDimensionsAndQuantification->GetNbTimeFrames(); fr++)
    {
      for (int rg=0; rg<mp_ImageDimensionsAndQuantification->GetNbRespGates(); rg++)
      {
        for (int cg=0; cg<mp_ImageDimensionsAndQuantification->GetNbCardGates(); cg++)
        {
          // Merge thread results
          for (int th=1; th<mp_ImageDimensionsAndQuantification->GetNbThreadsForProjection(); th++)
          {
            m4p_FOMNbBins[fr][rg][cg][0] += m4p_FOMNbBins[fr][rg][cg][th];
            m4p_FOMNbData[fr][rg][cg][0] += m4p_FOMNbData[fr][rg][cg][th];
            m4p_FOMLogLikelihood[fr][rg][cg][0] += m4p_FOMLogLikelihood[fr][rg][cg][th];
            m4p_FOMRMSE[fr][rg][cg][0] += m4p_FOMRMSE[fr][rg][cg][th];
          }
          // Finish mean number of data counts per channel
          m4p_FOMNbData[fr][rg][cg][0] /= ((HPFLTNB)(m4p_FOMNbBins[fr][rg][cg][0]));
          // Finish RMSE computation
          m4p_FOMRMSE[fr][rg][cg][0] /= ((HPFLTNB)(m4p_FOMNbBins[fr][rg][cg][0]));
          m4p_FOMRMSE[fr][rg][cg][0] = sqrt(m4p_FOMRMSE[fr][rg][cg][0]);
        }
      }
    }
    // Compute penalty FOM
    if (mp_Penalty!=NULL)
    {
      // As the penalty will be innocently applied to the current image, it is performed on time basis functions
      // and not frames/gates. As a result, we loop over time basis functions first, compute the penalty, and
      // distribute its value with respect to the contribution of each time basis function to the frames/gates.
      for (int tbf=0; tbf<mp_ImageDimensionsAndQuantification->GetNbTimeBasisFunctions(); tbf++)
      {
        for (int rbf=0; rbf<mp_ImageDimensionsAndQuantification->GetNbRespBasisFunctions(); rbf++)
        {
          for (int cbf=0; cbf<mp_ImageDimensionsAndQuantification->GetNbCardBasisFunctions(); cbf++)
          {
            // In case a problem occurs during the multi-threaded loop
            bool problem = false;
            // The voxel index
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
              // Compute the penalty term
              if (mp_Penalty->LocalPreProcessingStep(tbf,rbf,cbf,v,th))
              {
                Cerr("***** vOptimizer::PreImageUpdateStep() -> A problem occurred while precomputing the local step of the penalty for voxel " << v << " !" << endl);
                problem = true;
              }
              HPFLTNB penalty_value = mp_Penalty->ComputePenaltyValue(tbf,rbf,cbf,v,th);
              // Distribuete the value to frames/gates with respect to basis functions
              for (int fr=0; fr<mp_ImageDimensionsAndQuantification->GetNbTimeFrames(); fr++)
              {
                FLTNB time_basis_coef = mp_ImageDimensionsAndQuantification->GetTimeBasisCoefficient(tbf,fr);
                for (int rg=0; rg<mp_ImageDimensionsAndQuantification->GetNbRespGates(); rg++)
                {
                  FLTNB resp_basis_coef = mp_ImageDimensionsAndQuantification->GetRespBasisCoefficient(rbf,rg);
                  for (int cg=0; cg<mp_ImageDimensionsAndQuantification->GetNbCardGates(); cg++)
                  {
                    FLTNB card_basis_coef = mp_ImageDimensionsAndQuantification->GetCardBasisCoefficient(cbf,cg);
                    m4p_FOMPenalty[fr][rg][cg][th] += time_basis_coef * resp_basis_coef * card_basis_coef * penalty_value;
                  }
                }
              }
            }
            // Potential problems during the loop
            if (problem)
            {
              Cerr("***** vOptimizer::PreImageUpdateStep() -> A problem occurred in the multi-threaded loop for the penalty, stop now !" << endl);
              return 1;
            }
          }
        }
      }
      // Merge thread results
      for (int fr=0; fr<mp_ImageDimensionsAndQuantification->GetNbTimeFrames(); fr++)
        for (int rg=0; rg<mp_ImageDimensionsAndQuantification->GetNbRespGates(); rg++)
          for (int cg=0; cg<mp_ImageDimensionsAndQuantification->GetNbCardGates(); cg++)
          {
            for (int th=1; th<mp_ImageDimensionsAndQuantification->GetNbThreadsForImageComputation(); th++)
              m4p_FOMPenalty[fr][rg][cg][0] += m4p_FOMPenalty[fr][rg][cg][th];
            // Also divide by the number of subsets for balance with respect to the data fidelity term
            m4p_FOMPenalty[fr][rg][cg][0] /= ((HPFLTNB)mp_nbSubsets[m_currentIteration]);
          }
    }
    // Verbose
    Cout("vOptimizer::PreImageUpdateStep() -> Optimizer figures-of-merit related to the objective function" << endl);
    // Number of spaces for printing
    string spaces_fr = ""; if (mp_ImageDimensionsAndQuantification->GetNbTimeFrames()>1) spaces_fr += "  ";
    string spaces_rg = ""; if (mp_ImageDimensionsAndQuantification->GetNbRespGates()>1) spaces_rg += "  ";
    string spaces_cg = ""; if (mp_ImageDimensionsAndQuantification->GetNbCardGates()>1) spaces_cg += "  ";
    // Loop on frames and gates
    for (int fr=0; fr<mp_ImageDimensionsAndQuantification->GetNbTimeFrames(); fr++)
    {
      // Verbose
      if (mp_ImageDimensionsAndQuantification->GetNbTimeFrames()>1)
        Cout(spaces_fr << "Frame " << fr+1 << " on " << mp_ImageDimensionsAndQuantification->GetNbTimeFrames() << endl);
      for (int rg=0; rg<mp_ImageDimensionsAndQuantification->GetNbRespGates(); rg++)
      {
        // Verbose
        if (mp_ImageDimensionsAndQuantification->GetNbRespGates()>1)
          Cout(spaces_fr << spaces_rg << "Respiratory gate " << rg+1 << " on " << mp_ImageDimensionsAndQuantification->GetNbRespGates() << endl);
        for (int cg=0; cg<mp_ImageDimensionsAndQuantification->GetNbCardGates(); cg++)
        {
          // Verbose
          if (mp_ImageDimensionsAndQuantification->GetNbCardGates()>1)
            Cout(spaces_fr << spaces_rg << spaces_cg << "Cardiac gate " << cg+1 << " on " << mp_ImageDimensionsAndQuantification->GetNbCardGates() << endl);
          // Print out
          Cout(spaces_fr << spaces_rg << spaces_cg << "  --> Number of data channels used in optimization: " << m4p_FOMNbBins[fr][rg][cg][0] << endl);
          Cout(spaces_fr << spaces_rg << spaces_cg << "  --> Mean number of data counts per channel: " << m4p_FOMNbData[fr][rg][cg][0] << endl);
          if (mp_Penalty!=NULL)
          {
            Cout(spaces_fr << spaces_rg << spaces_cg << "  --> Penalty: " << m4p_FOMPenalty[fr][rg][cg][0] << endl);
            Cout(spaces_fr << spaces_rg << spaces_cg << "  --> Log-likelihood: " << m4p_FOMLogLikelihood[fr][rg][cg][0]
                           << " (including penalty: " << m4p_FOMLogLikelihood[fr][rg][cg][0] - m4p_FOMPenalty[fr][rg][cg][0] << ")" << endl);
            Cout(spaces_fr << spaces_rg << spaces_cg << "  --> RMSE: " << m4p_FOMRMSE[fr][rg][cg][0]
                           << " (including penalty: " << m4p_FOMRMSE[fr][rg][cg][0] + m4p_FOMPenalty[fr][rg][cg][0] << ")" << endl);
          }
          else
          {
            Cout(spaces_fr << spaces_rg << spaces_cg << "  --> Log-likelihood: " << m4p_FOMLogLikelihood[fr][rg][cg][0] << endl);
            Cout(spaces_fr << spaces_rg << spaces_cg << "  --> RMSE: " << m4p_FOMRMSE[fr][rg][cg][0] << endl);
          }
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

int vOptimizer::PreImageUpdateSpecificStep()
{
  // By default, just do nothing here. This function is designed to be overloaded by the specific optimizer if needed
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vOptimizer::DataStep1ForwardProjectModel( oProjectionLine* ap_Line, vEvent* ap_Event,
                                              int a_bed, int a_timeFrame, int a_respGate, int a_cardGate,
                                              int a_th )
{
  // =======================================================================
  // 4D forward projection including framing, respiratory and cardiac gating
  // =======================================================================

  // Clean buffers
  for (int b=0; b<ap_Line->GetNbTOFBins(); b++)
  {
    m2p_forwardValues[a_th][b] = 0.;
    for (int img=0; img<m_nbBackwardImages; img++) m3p_backwardValues[a_th][b][img] = 0.;
  }

  // Loop over time basis functions
  for (int tbf=0; tbf<mp_ImageDimensionsAndQuantification->GetNbTimeBasisFunctions(); tbf++)
  {
    // Get frame/basis coefficient and continue if null
    FLTNB time_basis_coef = mp_ImageDimensionsAndQuantification->GetTimeBasisCoefficient(tbf,a_timeFrame);
    if (time_basis_coef==0.) continue;
    // Loop over respiratory basis functions
    for (int rbf=0; rbf<mp_ImageDimensionsAndQuantification->GetNbRespBasisFunctions(); rbf++)
    {
      // Get resp_gate/basis coefficient and continue if null
      FLTNB resp_basis_coef = mp_ImageDimensionsAndQuantification->GetRespBasisCoefficient(rbf,a_respGate);
      if (resp_basis_coef==0.) continue;
      // Loop over cardiac basis functions
      for (int cbf=0; cbf<mp_ImageDimensionsAndQuantification->GetNbCardBasisFunctions(); cbf++)
      {
        // Get card_gate_basis coefficient and continue if null
        FLTNB card_basis_coef = mp_ImageDimensionsAndQuantification->GetCardBasisCoefficient(cbf,a_cardGate);
        if (card_basis_coef==0.) continue;
        // Compute global coefficient
        FLTNB global_basis_coef = time_basis_coef * resp_basis_coef * card_basis_coef;
        // Loop over TOF bins
        for (int b=0; b<ap_Line->GetNbTOFBins(); b++)
        {
          // Set the current TOF bin
          ap_Line->SetCurrentTOFBin(b);
          // Project the image and apply dynamic coefficient converions
          m2p_forwardValues[a_th][b] += ForwardProject(ap_Line,mp_ImageSpace->m4p_forwardImage[tbf][rbf][cbf]) * global_basis_coef;
        }
      }
    }
  }

  // ==============================================
  // Differentiate emission and transmission models
  // ==============================================

  // -----------------
  // Case for emission
  // -----------------
  if (m_dataSpec==SPEC_EMISSION)
  {
    // Add additive terms (we multiply by the frame duration because the additive terms are provided as rates)
    if (mp_ImageDimensionsAndQuantification->GateDurationProvided())
    for (int b=0; b<ap_Line->GetNbTOFBins(); b++) m2p_forwardValues[a_th][b] += ap_Event->GetAdditiveCorrections(b) * mp_ImageDimensionsAndQuantification->GetdurationPerGate(a_timeFrame,a_respGate);
    else
    for (int b=0; b<ap_Line->GetNbTOFBins(); b++) m2p_forwardValues[a_th][b] += ap_Event->GetAdditiveCorrections(b) * mp_ImageDimensionsAndQuantification->GetFrameDurationInSec(a_bed, a_timeFrame);
  }
  // ---------------------
  // Case for transmission
  // ---------------------
  else if (m_dataSpec==SPEC_TRANSMISSION)
  {
    int no_tof = 0;
    // Ignore TOF as it has no sense; add scatter rate multiplied by the frame duration; take blank into account
    m2p_forwardValues[a_th][no_tof] = ap_Event->GetBlankValue() * exp(-m2p_forwardValues[a_th][no_tof])
                                    + ap_Event->GetAdditiveCorrections(no_tof) * mp_ImageDimensionsAndQuantification->GetFrameDurationInSec(a_bed, a_timeFrame);
  }

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vOptimizer::DataStep2Optional( oProjectionLine* ap_Line, vEvent* ap_Event,
                                   int a_bed, int a_timeFrame, int a_respGate, int a_cardGate,
                                   int a_th )
{
  // Deliberately do nothing!
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vOptimizer::DataStep3BackwardProjectSensitivity( oProjectionLine* ap_Line, vEvent* ap_Event,
                                                     int a_bed, int a_timeFrame, int a_respGate, int a_cardGate,
                                                     int a_th )
{
  // Loop over TOF bins
  for (int b=0; b<ap_Line->GetNbTOFBins(); b++)
  {
    // Set the current TOF bin of the projection-line
    ap_Line->SetCurrentTOFBin(b);
    // The weight associated to the line to be backprojected
    FLTNB weight = 0.;
    // Call the function to compute specific sensitivity
    if (SensitivitySpecificOperations
    (
      ap_Event->GetEventValue(b),                                                                                         // the measured data
      m2p_forwardValues[a_th][b],                                                                                         // the forward model
      &weight,                                                                                                            // the sensitivity weight
      ap_Event->GetMultiplicativeCorrections(),                                                                           // the multiplicative corrections
      ap_Event->GetAdditiveCorrections(b)*mp_ImageDimensionsAndQuantification->GetFrameDurationInSec(a_bed, a_timeFrame), // the additive corrections
      ap_Event->GetBlankValue(),                                                                                          // the blank value
      mp_ImageDimensionsAndQuantification->GetQuantificationFactor(a_bed,a_timeFrame, a_respGate, a_cardGate),            // the quantification factor
      ap_Line                                                                                                             // the projection line
    ))
    {
      Cerr("***** vOptimizer::DataStep3BackwardProjectSensitivity() -> A problem occurred while performing the sensitivity specific operations !" << endl);
      return 1;
    }
    // Back-project the weight into the sensitivity matrix
    BackwardProject(ap_Line, mp_ImageSpace->m5p_sensitivity[a_th][a_timeFrame][a_respGate][a_cardGate], weight);
  }

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vOptimizer::DataStep4Optional( oProjectionLine* ap_Line, vEvent* ap_Event,
                                   int a_bed, int a_timeFrame, int a_respGate, int a_cardGate,
                                   int a_th )
{
  // Deliberately do nothing!
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vOptimizer::DataStep5ComputeCorrections( oProjectionLine* ap_Line, vEvent* ap_Event,
                                             int a_bed, int a_timeFrame, int a_respGate, int a_cardGate,
                                             int a_th )
{
  // Loop on TOF bins
  for (int b=0; b<ap_Line->GetNbTOFBins(); b++)
  {
    // Set the current TOF bin of the projection-line
    ap_Line->SetCurrentTOFBin(b);
    // Then call the pure virtual function for specific data space operations
    if (DataSpaceSpecificOperations
    (
      ap_Event->GetEventValue(b),                                                                                         // the measured data
      m2p_forwardValues[a_th][b],                                                                                         // the forward model
      m3p_backwardValues[a_th][b],                                                                                        // the backward values (the result)
      ap_Event->GetMultiplicativeCorrections(),                                                                           // the multiplicative corrections
      ap_Event->GetAdditiveCorrections(b)*mp_ImageDimensionsAndQuantification->GetFrameDurationInSec(a_bed, a_timeFrame), // the additive corrections
      ap_Event->GetBlankValue(),                                                                                          // the blank value
      mp_ImageDimensionsAndQuantification->GetQuantificationFactor(a_bed,a_timeFrame, a_respGate, a_cardGate),            // the quantification factor
      ap_Line                                                                                                             // the projection line
    ))
    {
      Cerr("***** vOptimizer::DataStep5ComputeCorrections() -> A problem occurred while performing specific data space operations !" << endl);
      return 1;
    }
  }

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vOptimizer::DataStep6Optional( oProjectionLine* ap_Line, vEvent* ap_Event,
                                   int a_bed, int a_timeFrame, int a_respGate, int a_cardGate,
                                   int a_th )
{
  // Deliberately do nothing!
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vOptimizer::DataStep7BackwardProjectCorrections( oProjectionLine* ap_Line, vEvent* ap_Event,
                                                     int a_bed, int a_timeFrame, int a_respGate, int a_cardGate,
                                                     int a_th )
{
  // ========================================================================
  // 4D backward projection including framing, respiratory and cardiac gating
  // ========================================================================

  // Loop over time basis functions
  for (int tbf=0; tbf<mp_ImageDimensionsAndQuantification->GetNbTimeBasisFunctions(); tbf++)
  {
    // Get frame/basis coefficient and continue if null
    FLTNB time_basis_coef = mp_ImageDimensionsAndQuantification->GetTimeBasisCoefficient(tbf,a_timeFrame);
    if (time_basis_coef==0.) continue;
    // Loop over respiratory basis functions
    for (int rbf=0; rbf<mp_ImageDimensionsAndQuantification->GetNbRespBasisFunctions(); rbf++)
    {
      // Get resp_gate/basis coefficient and continue if null
      FLTNB resp_basis_coef = mp_ImageDimensionsAndQuantification->GetRespBasisCoefficient(rbf,a_respGate);
      if (resp_basis_coef==0.) continue;
      // Loop over cardiac basis functions
      for (int cbf=0; cbf<mp_ImageDimensionsAndQuantification->GetNbCardBasisFunctions(); cbf++)
      {
        // Get card_gate_basis coefficient and continue if null
        FLTNB card_basis_coef = mp_ImageDimensionsAndQuantification->GetCardBasisCoefficient(cbf,a_cardGate);
        if (card_basis_coef==0.) continue;
        // Compute global coefficient
        FLTNB global_basis_coef = time_basis_coef * resp_basis_coef * card_basis_coef;
        // Loop over TOF bins
        for (int b=0; b<ap_Line->GetNbTOFBins(); b++)
        {
          // Set the current TOF bin of the projection-line
          ap_Line->SetCurrentTOFBin(b);
          // Loop over backprojection images
          for (int img=0; img<m_nbBackwardImages; img++)
          {
            // Continue if backproj is null or infinity
            if (m3p_backwardValues[a_th][b][img]==0. || !isfinite(m3p_backwardValues[a_th][b][img])) continue;
            // Backward project into the backward image
            BackwardProject( ap_Line, mp_ImageSpace->m6p_backwardImage[img][a_th][tbf][rbf][cbf] , m3p_backwardValues[a_th][b][img] * global_basis_coef );
          }
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

int vOptimizer::DataStep8ComputeFOM( oProjectionLine* ap_Line, vEvent* ap_Event,
                                     int a_timeFrame, int a_respGate, int a_cardGate,
                                     int a_thread )
{
  // Loop on TOF bins
  for (int b=0; b<ap_Line->GetNbTOFBins(); b++)
  {
    // Get the model and the data
    HPFLTNB model = m2p_forwardValues[a_thread][b];
    HPFLTNB data  = ap_Event->GetEventValue(b);
    // Update log likelihood only if model is strictly positive (to avoid numerical issues, we use a threshold at e-10)
    if (model>1.e-10) m4p_FOMLogLikelihood[a_timeFrame][a_respGate][a_cardGate][a_thread] += data*log(model)-model;
    // Update RMSE
    m4p_FOMRMSE[a_timeFrame][a_respGate][a_cardGate][a_thread] += (data-model)*(data-model);
    // Update number of data channels
    m4p_FOMNbBins[a_timeFrame][a_respGate][a_cardGate][a_thread]++;
    // Update the number of data counts
    m4p_FOMNbData[a_timeFrame][a_respGate][a_cardGate][a_thread] += data;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vOptimizer::UpdateVisitedVoxels()
{
  // Verbose
  if (m_verbose>=3) Cout("vOptimizer::UpdateVisitedVoxels() -> Tick visited voxels based on sensitivity" << endl);

  // The purpose of this mp_visitedVoxelsImage is to keep track within a complete iteration, of which voxels were
  // visited by some LORs. Everything is based on sensitivity. With list-mode data, it means that the visited voxels
  // are the same for all subsets. With histogram data, if the sensitivity is computed on the fly from the data, then
  // the visited voxels may change for each voxel. When updating the voxel values in the ImageUpdateStep, only voxels
  // with non-zero sensitivity are updated, to avoid artificially decreasing the signal in voxels that may not have
  // been visited in this subset but which may be visited in another subset (this is true only for histogram data).
  // But, at the end of each iteration, the mp_visitedVoxelsImage is used to checked for never-visited voxels, in which
  // case they will be set to 0 to avoid that they keep their initial value. In other words, these never-visited voxels
  // are somewhat excluded from the reconstruction.

  // For this matter of visited voxels does not depend on dynamic aspects, we loop on the sensitivity over all frames and only for the first
  // gates and look if voxels were 'visited' by LORs (i.e. the sensitivity is not null);
  // the information is saved in one static 3D image and then applied to all frames/gates;
  int thread0 = 0;
  int resp_gate0 = 0;
  int card_gate0 = 0;
  for (int fr=0; fr<mp_ImageDimensionsAndQuantification->GetNbTimeFrames(); fr++)
    for (int v=0; v<mp_ImageDimensionsAndQuantification->GetNbVoxXYZ(); v++)
    // If the voxel has been visited, we add 1. to the mp_Image->mp_visitedVoxelsImage
    if (mp_ImageSpace->m5p_sensitivity[thread0][fr][resp_gate0][card_gate0][v]!=0.) mp_ImageSpace->mp_visitedVoxelsImage[v] += 1.;

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vOptimizer::ImageUpdateStep()
{
  // Verbose
  if (m_verbose>=2 || m_optimizerImageStatFlag) Cout("vOptimizer::ImageUpdateStep() -> Start image update" << endl);
  // Number of spaces for nice printing
  string spaces_tbf = ""; if (mp_ImageDimensionsAndQuantification->GetNbTimeBasisFunctions()>1) spaces_tbf += "  ";
  string spaces_rbf = ""; if (mp_ImageDimensionsAndQuantification->GetNbRespBasisFunctions()>1) spaces_rbf += "  ";
  string spaces_cbf = ""; if (mp_ImageDimensionsAndQuantification->GetNbCardBasisFunctions()>1) spaces_cbf += "  ";
  // Loop over time basis functions
  for (int tbf=0; tbf<mp_ImageDimensionsAndQuantification->GetNbTimeBasisFunctions(); tbf++)
  {
    // Verbose
    if (mp_ImageDimensionsAndQuantification->GetNbTimeBasisFunctions()>1 && (m_verbose>=3 || m_optimizerImageStatFlag))
    {
      if (mp_ImageDimensionsAndQuantification->GetTimeStaticFlag()) Cout(spaces_tbf << "--> Time frame " << tbf+1 << endl);
      else Cout(spaces_tbf << "--> Time basis function: " << tbf+1 << endl);
    }
    // Loop over respiratory basis functions
    for (int rbf=0; rbf<mp_ImageDimensionsAndQuantification->GetNbRespBasisFunctions(); rbf++)
    {
      // Verbose
      if (mp_ImageDimensionsAndQuantification->GetNbRespBasisFunctions()>1 && (m_verbose>=3 || m_optimizerImageStatFlag))
      {
        if (mp_ImageDimensionsAndQuantification->GetRespStaticFlag()) Cout(spaces_tbf << spaces_rbf << "--> Respiratory gate " << rbf+1 << endl);
        else Cout(spaces_tbf << spaces_rbf << "--> Respiratory basis function: " << rbf+1 << endl);
      }
      // Loop over cardiac basis functions
      for (int cbf=0; cbf<mp_ImageDimensionsAndQuantification->GetNbCardBasisFunctions(); cbf++)
      {
        // Verbose
        if (mp_ImageDimensionsAndQuantification->GetNbCardBasisFunctions()>1 && (m_verbose>=3 || m_optimizerImageStatFlag))
        {
          if (mp_ImageDimensionsAndQuantification->GetCardStaticFlag()) Cout(spaces_tbf << spaces_rbf << spaces_cbf << "--> Cardiac gate " << cbf+1 << endl);
          else Cout(spaces_tbf << spaces_rbf << spaces_cbf << "--> Cardiac basis function: " << cbf+1 << endl);
        }
        // Set the number of threads for image computation
        #ifdef CASTOR_OMP
        omp_set_num_threads(mp_ImageDimensionsAndQuantification->GetNbThreadsForImageComputation());
        #endif
        // Reset counter for image stats
        if (m_optimizerImageStatFlag)
        {
          for (int th=0; th<mp_ImageDimensionsAndQuantification->GetNbThreadsForImageComputation(); th++)
          {
            mp_imageStatNbVox[th] = 0;
            mp_imageStatMin[th] = mp_ImageSpace->m4p_image[tbf][rbf][cbf][0];
            mp_imageStatMax[th] = mp_ImageSpace->m4p_image[tbf][rbf][cbf][0];
            mp_imageStatMean[th] = 0.;
            mp_imageStatVariance[th] = 0.;
            mp_correctionStatMean[th] = 0.;
            mp_correctionStatVariance[th] = 0.;
          }
        }
        // OMP loop over voxels
        INTNB v;
        bool error = false;
        #pragma omp parallel for private(v) schedule(guided)
        for (v=0 ; v<mp_ImageDimensionsAndQuantification->GetNbVoxXYZ() ; v++) 
        {
          // Get the thread index
          int th = 0;
          #ifdef CASTOR_OMP
          th = omp_get_thread_num();
          #endif
          // Compute sensitivity
          FLTNB sensitivity = ComputeSensitivity(mp_ImageSpace->m5p_sensitivity[0],tbf,rbf,cbf,v);
          // If the sensitivity is negative or null, we skip this voxel
          if (sensitivity<=0.) continue;
          // Keep the current image value in a buffer for update statistics
          FLTNB old_image_value = mp_ImageSpace->m4p_image[tbf][rbf][cbf][v];
          // Fill class-local buffer of backward values (thread-safe)
          for (int nb=0; nb<m_nbBackwardImages; nb++)
            m3p_backwardValues[th][0][nb] = mp_ImageSpace->m6p_backwardImage[nb][0][tbf][rbf][cbf][v];
          // Then call the pure virtual function for specific data space operations
          if (ImageSpaceSpecificOperations
          (
            mp_ImageSpace->m4p_image[tbf][rbf][cbf][v],
            // Note: have to be carefull here because we pass rbf and cbf as the actual respgate and cardgate. It may be wrong if one would like
            //       to perform MLAA reconstruction but using PET basis functions different from diracs for example... In this case, we would get
            //       segmentation faults. But actually, this function would be overloaded by MLAA anyway...
            //mp_Image->m4p_attenuation != NULL ? mp_Image->m4p_attenuation[tbf][rbf][cbf][v] : 0. ,
            &mp_ImageSpace->m4p_image[tbf][rbf][cbf][v],
            //mp_Image->m4p_attenuation != NULL ? &mp_Image->m4p_attenuation[tbf][rbf][cbf][v] : NULL ,
            sensitivity,
            m3p_backwardValues[th][0],
            v,
            tbf,
            rbf,
            cbf
          ))
          {
            Cerr("***** vOptimizer::ImageUpdateStep() -> A problem occurred while performing image space specific operations of the optimizer !" << endl);
            error = true;
          }
          // Do some image statistics
          if (m_optimizerImageStatFlag)
          {
            if (mp_imageStatMin[th]>mp_ImageSpace->m4p_image[tbf][rbf][cbf][v]) mp_imageStatMin[th] = mp_ImageSpace->m4p_image[tbf][rbf][cbf][v];
            if (mp_imageStatMax[th]<mp_ImageSpace->m4p_image[tbf][rbf][cbf][v]) mp_imageStatMax[th] = mp_ImageSpace->m4p_image[tbf][rbf][cbf][v];
            // The one pass algorithm is used to compute the variance here for both the image and correction step.
            // Do not change anyhting to the order of the operations !
            mp_imageStatNbVox[th]++;
            HPFLTNB sample = ((HPFLTNB)(mp_ImageSpace->m4p_image[tbf][rbf][cbf][v]));
            HPFLTNB delta = sample - mp_imageStatMean[th];
            mp_imageStatMean[th] += delta / ((HPFLTNB)(mp_imageStatNbVox[th]));
            mp_imageStatVariance[th] += delta * (sample - mp_imageStatMean[th]);
            // Same one-pass variance algorithm for the correction step
            sample = ((HPFLTNB)(mp_ImageSpace->m4p_image[tbf][rbf][cbf][v] - old_image_value));
            delta = sample - mp_correctionStatMean[th];
            mp_correctionStatMean[th] += delta / ((HPFLTNB)(mp_imageStatNbVox[th]));
            mp_correctionStatVariance[th] += delta * (sample - mp_correctionStatMean[th]);
          }
        }
        // Check for errors (we cannot stop the multi-threaded loop so we do it here)
        if (error) return 1;
        // Finish image update statistics computations: these operations are specific to the merging of the
        // multi-threaded one-pass variance estimations used above.
        if (m_optimizerImageStatFlag)
        {
          for (int th=1; th<mp_ImageDimensionsAndQuantification->GetNbThreadsForImageComputation(); th++)
          {
            // Image minimum and maximum
            if (mp_imageStatMin[0]>mp_imageStatMin[th]) mp_imageStatMin[0] = mp_imageStatMin[th];
            if (mp_imageStatMax[0]<mp_imageStatMax[th]) mp_imageStatMax[0] = mp_imageStatMax[th];
            // Cast the counts into HPFLTNB
            HPFLTNB count1 = ((HPFLTNB)(mp_imageStatNbVox[0]));
            HPFLTNB count2 = ((HPFLTNB)(mp_imageStatNbVox[th]));
            // Update the count
            mp_imageStatNbVox[0] += mp_imageStatNbVox[th];
            HPFLTNB count12 = ((HPFLTNB)(mp_imageStatNbVox[0]));
            // Compute the delta mean
            HPFLTNB delta = mp_imageStatMean[0] - mp_imageStatMean[th];
            // Update the variance
            mp_imageStatVariance[0] += mp_imageStatVariance[th] + delta * delta * count1 * count2 / count12;
            // Update the mean
            mp_imageStatMean[0] = (count1*mp_imageStatMean[0] + count2*mp_imageStatMean[th]) / count12;
            // Do the same for the correction step
            delta = mp_correctionStatMean[0] - mp_correctionStatMean[th];
            mp_correctionStatVariance[0] += mp_correctionStatVariance[th] + delta * delta * count1 * count2 / count12;
            mp_correctionStatMean[0] = (count1*mp_correctionStatMean[0] + count2*mp_correctionStatMean[th]) / count12;
          }
          // Finish variance computation
          mp_imageStatVariance[0] /= ((HPFLTNB)(mp_imageStatNbVox[0]));
          mp_correctionStatVariance[0] /= ((HPFLTNB)(mp_imageStatNbVox[0]));
          // Print out
          Cout(spaces_tbf << spaces_rbf << spaces_cbf << "  --> Image update step "
                                                      << " | mean: " << mp_correctionStatMean[0]
                                                      << " | stdv: " << sqrt(mp_correctionStatVariance[0]) << endl);
          Cout(spaces_tbf << spaces_rbf << spaces_cbf << "  --> New image estimate"
                                                      << " | mean: " << mp_imageStatMean[0]
                                                      << " | stdv: " << sqrt(mp_imageStatVariance[0])
                                                      << " | min/max: " << mp_imageStatMin[0] << "/" << mp_imageStatMax[0] << endl);
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

FLTNB vOptimizer::ForwardProject(oProjectionLine* ap_Line, FLTNB* ap_image)
{
  // Particular case for SPECT with attenuation correction
  if (m_dataType==TYPE_SPECT && m2p_attenuationImage[ap_Line->GetThreadNumber()]!=NULL)
    return ap_Line->ForwardProjectWithSPECTAttenuation( m2p_attenuationImage[ap_Line->GetThreadNumber()], ap_image );
  // Otherwise call the function without attenuation
  else
    return ap_Line->ForwardProject( ap_image );
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void vOptimizer::BackwardProject(oProjectionLine* ap_Line, FLTNB* ap_image, FLTNB a_value)
{
  // Particular case for SPECT with attenuation correction
  if (m_dataType==TYPE_SPECT && m2p_attenuationImage[ap_Line->GetThreadNumber()]!=NULL)
    ap_Line->BackwardProjectWithSPECTAttenuation( m2p_attenuationImage[ap_Line->GetThreadNumber()], ap_image, a_value );
  // Otherwise call the function without attenuation
  else
    ap_Line->BackwardProject( ap_image, a_value );
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

FLTNB vOptimizer::ComputeSensitivity( FLTNB**** a4p_sensitivityMatrix,
                                      int a_timeBasisFunction, int a_respBasisFunction, int a_cardBasisFunction,
                                      int a_voxel )
{
  // The sensitivity to be computed
  FLTNB sensitivity = 0.;
  // Loop over time frames
  for (int fr=0; fr<mp_ImageDimensionsAndQuantification->GetNbTimeFrames(); fr++)
  {
    // Get frame/basis coefficient and continue if null
    FLTNB time_basis_coef = mp_ImageDimensionsAndQuantification->GetTimeBasisCoefficient(a_timeBasisFunction,fr);
    if (time_basis_coef==0.) continue;
    // Loop over respiratory gates          
    for (int rg=0; rg<mp_ImageDimensionsAndQuantification->GetNbRespGates(); rg++)
    {
      // Get resp_gate/basis coefficient and continue if null
      FLTNB resp_basis_coef = mp_ImageDimensionsAndQuantification->GetRespBasisCoefficient(a_respBasisFunction,rg);
      if (resp_basis_coef==0.) continue;
      // Loop over cardiac gates
      for (int cg=0; cg<mp_ImageDimensionsAndQuantification->GetNbCardGates(); cg++)
      {
        // Get card_gate_basis coefficient and continue if null
        FLTNB card_basis_coef = mp_ImageDimensionsAndQuantification->GetCardBasisCoefficient(a_cardBasisFunction,cg);
        if (card_basis_coef==0.) continue;
        // Add actual weight contribution
        sensitivity += a4p_sensitivityMatrix[fr][rg][cg][a_voxel] *
                       time_basis_coef * resp_basis_coef * card_basis_coef;
      }
    }
  }
  // If listmode, then divide by the number of subsets
  if (m_dataMode==MODE_LIST || m_needGlobalSensitivity) sensitivity /= ((FLTNB)mp_nbSubsets[m_currentIteration]);
  // Return sensitivity
  return sensitivity;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

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
  \ingroup  projector
  \brief    Implementation of class vProjector
*/

#include "vProjector.hh"
#include "vScanner.hh"
#include "vDataFile.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

vProjector::vProjector()
{
  // Affect default values
  mp_Scanner = NULL;
  mp_ImageDimensionsAndQuantification = NULL;
  mp_sizeVox[0] = -1.;
  mp_sizeVox[1] = -1.;
  mp_sizeVox[2] = -1.;
  mp_nbVox[0] = -1;
  mp_nbVox[1] = -1;
  mp_nbVox[2] = -1;
  m_nbVoxXY = -1;
  mp_halfFOV[0] = -1.;
  mp_halfFOV[1] = -1.;
  mp_halfFOV[2] = -1.;
  m_sensitivityMode = false;
  m_TOFMethod = -1;
  m_TOFNbSigmas = -1.;
  m_TOFMeasurementRangeInMm = -1.;
  m_applyPOI = false;
  m_compatibleWithSPECTAttenuationCorrection = false;
  m_compatibleWithCompression = false;
  m_verbose = 0;
  m_checked = false;
  m_initialized = false;
  mp_mask = NULL;
  m_hasMask = false;
  mp_TOFWeightingFcn = NULL;
  m_TOFWeightingFcnPrecomputedFlag = true;
  m_TOFWeightingFcnNbSamples = -1;
  m_TOFResolutionInMm = -1.;
  m_TOFBinSizeInMm = -1.;
  m_TOFBinProperProcessingFlag = true;
  m_TOFGaussianNormCoef = 0.;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

vProjector::~vProjector()
{
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vProjector::SetImageDimensionsAndQuantification(oImageDimensionsAndQuantification* ap_ImageDimensionsAndQuantification)
{
  // Check that the parameter is not NULL
  if (ap_ImageDimensionsAndQuantification==NULL)
  {
    Cerr("***** vProjector::SetImageDimensionsAndQuantification() -> Input image dimensions object is null !" << endl);
    return 1;
  }
  // Affect image dimensions
  mp_sizeVox[0] = ap_ImageDimensionsAndQuantification->GetVoxSizeX();
  mp_sizeVox[1] = ap_ImageDimensionsAndQuantification->GetVoxSizeY();
  mp_sizeVox[2] = ap_ImageDimensionsAndQuantification->GetVoxSizeZ();
  mp_nbVox[0] = ap_ImageDimensionsAndQuantification->GetNbVoxX();
  mp_nbVox[1] = ap_ImageDimensionsAndQuantification->GetNbVoxY();
  mp_nbVox[2] = ap_ImageDimensionsAndQuantification->GetNbVoxZ();
  m_nbVoxXY = mp_nbVox[0] * mp_nbVox[1];
  mp_halfFOV[0] = mp_sizeVox[0] * ((FLTNB)mp_nbVox[0]) / 2.;
  mp_halfFOV[1] = mp_sizeVox[1] * ((FLTNB)mp_nbVox[1]) / 2.;
  mp_halfFOV[2] = mp_sizeVox[2] * ((FLTNB)mp_nbVox[2]) / 2.;
  mp_ImageDimensionsAndQuantification = ap_ImageDimensionsAndQuantification;
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void vProjector::ShowCommonHelp()
{
  // Return when using MPI and mpi_rank is not 0
  #ifdef CASTOR_MPI
  int mpi_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  if (mpi_rank!=0) return;
  #endif
  // Show help
  cout << "------------------------------------------------------------------" << endl;
  cout << "-----  Common options for all projectors" << endl;
  cout << "------------------------------------------------------------------" << endl;
  cout << "Only options related to TOF implementation are available." << endl;
  cout << "If used, values for all of the following options must be provided as a list of numbers separated by commas." << endl;
  cout << "" << endl;
  cout << "  the number of standard deviations for truncating the nominal TOF Gaussian distribution (-1 for no truncation)." << endl;
  cout << "  whether the TOF weighting function is precomputed (1 for yes, 0 for no)" << endl;
  cout << "  whether the TOF weighting function takes properly into account the TOF bin using convolution or integration (1 for yes, 0 for no)." << endl;
  cout << "" << endl;
  cout << "The default values are -1,1,1" << endl;
  
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void vProjector::ShowHelp()
{
  // Call the specific help function from the children
  ShowHelpSpecific();
  // Then, say if the child projector in use is compatible with SPECT attenuation correction or not
  if (m_compatibleWithSPECTAttenuationCorrection) cout << "This projector is compatible with SPECT attenuation correction." << endl;
  else cout << "This projector is NOT compatible with SPECT attenuation correction." << endl;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vProjector::ReadCommonOptionsList(const string& a_optionsList)
{
  // TODO : make this more explicit, here it is assumed that 3 specific TOF options are provided
  if (a_optionsList!="")
  {
    FLTNB option[3];
    // Read the options
    if (ReadStringOption(a_optionsList, option, 3, ",", "Common options"))
    {
      Cerr("***** vProjector::ReadCommonOptionsList() -> Failed to correctly read the list of options !" << endl);
      return 1;
    }
    m_TOFNbSigmas = option[0];
    m_TOFWeightingFcnPrecomputedFlag = option[1]>0.;
    m_TOFBinProperProcessingFlag = option[2]>0.;
  }
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vProjector::CheckParameters()
{
  // Check scanner
  if (mp_Scanner==NULL)
  {
    Cerr("***** vProjector::CheckParameters() -> Please provide a valid scanner object !" << endl);
    return 1;
  }
  // Check image dimensions
  if (mp_ImageDimensionsAndQuantification==NULL)
  {
    Cerr("***** vProjector::CheckParameters() -> Please provide a valid image dimensions and quantification object !" << endl);
    return 1;
  }
  if (mp_sizeVox[0]<=0. || mp_sizeVox[1]<=0. || mp_sizeVox[2]<=0.)
  {
    Cerr("***** vProjector::CheckParameters() -> One or more voxel sizes is negative or null !" << endl);
    return 1;
  }
  if (mp_nbVox[0]<=0 || mp_nbVox[1]<=0 || mp_nbVox[2]<=0)
  {
    Cerr("***** vProjector::CheckParameters() -> One or more number of voxels is negative or null !" << endl);
    return 1;
  }
  // Check TOF
  if ( m_TOFMethod!=USE_TOFLIST && m_TOFMethod!=USE_TOFHISTO && m_TOFMethod!=USE_NOTOF )
  {
    Cerr("***** vProjector::CheckParameters() -> TOF flag is incorrect or not set !" << endl);
    return 1;
  }
  // Check TOF parameters
  if (m_TOFMethod!=USE_NOTOF)
  {
    if (m_TOFResolutionInMm<=0. || m_TOFMeasurementRangeInMm<=0. || (m_TOFMethod==USE_TOFHISTO && m_TOFBinSizeInMm<=0.))
    {
      Cerr("***** vProjector::CheckParameters() -> Inconsistent TOF related parameters !" << endl);
      return 1;
    }
    if (m_TOFMethod==USE_TOFLIST && m_TOFBinProperProcessingFlag && m_TOFBinSizeInMm<=0.)
    {
      Cout("***** vProjector::CheckParameters() -> Warning: quantization TOF bin size wrong or not provided, so switching to TOF list-mode reconstruction which neglects the quantization TOF bin size!" << endl);
      m_TOFBinProperProcessingFlag = false;
    }
  }
  // Check verbose level
  if (m_verbose<0)
  {
    Cerr("***** vProjector::CheckParameters() -> Verbose level is negative !" << endl);
    return 1;
  }
  // Check parameters of the child class
  if (CheckSpecificParameters())
  {
    Cerr("***** vProjector::CheckParameters() -> An error occurred while checking parameters of the child projector class !" << endl);
    return 1;
  }

  // Normal end
  m_checked = true;
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vProjector::Initialize()
{
  // First check that the parameters have been checked !
  if (!m_checked)
  {
    Cerr("***** vProjector::Initialize() -> Must call the CheckParameters() function before initialization !" << endl);
    return 1;
  }

  // TOF parameters
  if (m_TOFMethod!=USE_NOTOF)
  {
    // Gaussian standard deviation
    HPFLTNB tof_resolution_sigma = m_TOFResolutionInMm / TWO_SQRT_TWO_LN_2;

    // If the provided number of sigmas is negative, set it to a value large enough
    // to approximate no truncation of the Gaussian distribution
    // (the TOF weighting function is in any case limited by the available TOF measurement range)
    if (m_TOFNbSigmas<0.) m_TOFNbSigmas = m_TOFMeasurementRangeInMm / tof_resolution_sigma;

    // Ideal Gaussian normalization coefficient (integral=1)
    m_TOFGaussianNormCoef = INV_SQRT_2_PI / tof_resolution_sigma;

    // Precompute the TOF weighting function if required
    if (m_TOFWeightingFcnPrecomputedFlag)
    {
      // Initialize some variables to default values, TODO make this an optional parameter
      // Sample the TOF weighting function at 1/m_TOFPrecomputedSamplingFactor of the spatial unit (mm)
      m_TOFPrecomputedSamplingFactor = 1000.;

      HPFLTNB tof_half_span = tof_resolution_sigma * m_TOFNbSigmas;

      // number of samples for the Gaussian function truncated at m_TOFNbSigmas standard deviations
      INTNB nb_samples_tof_gauss = (INTNB)ceil(tof_half_span*2.*m_TOFPrecomputedSamplingFactor);
      // make the number of samples odd
      if (nb_samples_tof_gauss % 2 == 0) nb_samples_tof_gauss += 1;

      // normalized Gaussian function for list-mode data (assumption of continuous TOF measurements)
      // normalized Gaussian function for TOF histogram data multiplied by the size of the TOF bin (assumption of constant Gaussian function over a TOF bin)
      if (!m_TOFBinProperProcessingFlag)
      {
        mp_TOFWeightingFcn = new HPFLTNB[nb_samples_tof_gauss];
        // Gaussian mean at the center
        INTNB mu = nb_samples_tof_gauss/2;
        for (INTNB g=0;g<nb_samples_tof_gauss;g++)
        {
          HPFLTNB temp = ((HPFLTNB)(g-mu))/(tof_resolution_sigma*m_TOFPrecomputedSamplingFactor);
          mp_TOFWeightingFcn[g] = m_TOFGaussianNormCoef * exp(-0.5*(temp*temp));
          if (m_TOFMethod==USE_TOFHISTO) mp_TOFWeightingFcn[g] *= m_TOFBinSizeInMm;
        }
        m_TOFWeightingFcnNbSamples = nb_samples_tof_gauss;
      }
      // convolution of the normalized Gaussian with the TOF bin (cumulative or quantization bin)
      else
      {
        // number of samples for the TOF bin door function
        INTNB nb_samples_tof_bin = (INTNB)round(m_TOFBinSizeInMm*m_TOFPrecomputedSamplingFactor);
        // make the number of samples odd
        if (nb_samples_tof_bin % 2 == 0) nb_samples_tof_bin += 1;
    
        // number of samples for the convolution
        INTNB nb_samples_conv = nb_samples_tof_gauss+nb_samples_tof_bin-1;

        // normalized Gaussian function truncated at m_TOFNbSigmas standard deviations
        // padded with zeros to the size of the convolved function
        HPFLTNB* tof_gauss = new HPFLTNB[nb_samples_conv];
        INTNB mu = nb_samples_conv/2;
        for (INTNB sgauss=0;sgauss<nb_samples_conv;sgauss++)
        {
          if (sgauss<mu-nb_samples_tof_gauss/2 || sgauss>mu+nb_samples_tof_gauss/2) tof_gauss[sgauss]=0.;
          else
          {
            HPFLTNB temp = (HPFLTNB)((sgauss-mu))/(tof_resolution_sigma*m_TOFPrecomputedSamplingFactor);
            tof_gauss[sgauss] = m_TOFGaussianNormCoef * exp(-0.5*(temp*temp));
          }
        }

        // TOF bin door function (quantization bin for list-mode, cumulative bin for histogram data)
        HPFLTNB* tof_bin = new HPFLTNB[nb_samples_tof_bin];
        // intensity = 1 for the cumulative bin, 1/binSize for the quantization bin
        HPFLTNB tof_bin_value = (m_TOFMethod==USE_TOFLIST)?(1./(m_TOFBinSizeInMm*m_TOFPrecomputedSamplingFactor)):(1./m_TOFPrecomputedSamplingFactor);
        for (INTNB sbin=0;sbin<nb_samples_tof_bin;sbin++) tof_bin[sbin] = tof_bin_value;

        // the final TOF weighting function
        mp_TOFWeightingFcn = new HPFLTNB[nb_samples_conv];
        for (INTNB s=0;s<nb_samples_conv;s++) mp_TOFWeightingFcn[s] = 0.;

        // convolve the TOF bin with the Gaussian function
        for (INTNB c=0;c<nb_samples_conv;c++)
        {
          for (INTNB ib=0;ib<nb_samples_tof_bin;ib++)
          {
            INTNB temp = c-nb_samples_tof_bin/2+ib;
            if (temp>=0 && temp<nb_samples_conv) mp_TOFWeightingFcn[c] += tof_gauss[temp]*tof_bin[ib];
          }
        }

        m_TOFWeightingFcnNbSamples = nb_samples_conv;

        // clean
        delete []tof_bin;
        delete []tof_gauss;
        tof_bin = NULL;
        tof_gauss = NULL;
      }
    }

    // Verbose
    if (m_verbose>=VERBOSE_NORMAL)
    {
      Cout("  --> TOF weighting function implementation: " << endl);
      Cout("  -->   Gaussian truncation (number of standard deviations) " << m_TOFNbSigmas << endl);
      if (m_TOFWeightingFcnPrecomputedFlag)
      {
        Cout("  -->   Precomputed " << endl);
        if (m_TOFBinProperProcessingFlag )
        {
          Cout("  -->   Convolved with the "<<((m_TOFMethod==USE_TOFHISTO)?"cumulative":"quantization")<<" TOF bin (size " <<m_TOFBinSizeInMm<<"mm)" <<endl);
        }
        else
        {
          if (m_TOFMethod==USE_TOFHISTO) Cout("  -->   Simple multiplication with the cumulative TOF bin (size " <<m_TOFBinSizeInMm<<"mm)" <<endl);
          else if (m_TOFMethod==USE_TOFLIST) Cout("  -->   Simple normalized Gaussian"<<endl);
        }
      }
      else
      {
        Cout("  -->   Computation on the fly " << endl);
        if (m_TOFBinProperProcessingFlag )
        {
          Cout("  -->   Integration of the Gaussian over the "<<((m_TOFMethod==USE_TOFHISTO)?"cumulative":"quantization")<<" TOF bin (size " <<m_TOFBinSizeInMm<<"mm)" <<endl);
        }
        else
        {
          if (m_TOFMethod==USE_TOFHISTO) Cout("  -->   Simple multiplication with the cumulative TOF bin (size " <<m_TOFBinSizeInMm<<"mm)" <<endl);
          else if (m_TOFMethod==USE_TOFLIST) Cout("  -->  Simple normalized Gaussian"<<endl);
        }
      }
    }
  }

  // Call the intialize function specific to the children
  if (InitializeSpecific())
  {
    Cerr("***** vProjector::Initialize() -> A problem occurred while calling the specific initialization of the child projector !" << endl);
    return 1;
  }

  if (m_verbose>=3)
  {
    Cout("vProjector::Initialize() -> Exit function" << endl);
  }

  // Normal end
  m_initialized = true;
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

INTNB vProjector::EstimateMaxNumberOfVoxelsPerLine()
{
  return mp_ImageDimensionsAndQuantification->GetNbVoxXYZ();
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vProjector::Project(int a_direction, oProjectionLine* ap_ProjectionLine, uint32_t* ap_index1, uint32_t* ap_index2, int a_nbIndices)
{
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** vProjector::Project() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif

  DEBUG_VERBOSE(m_verbose,VERBOSE_DEBUG_EVENT)

  // ---------------------------------------------------------------------------------------
  // First: Get cartesian coordinates from the scanner and average positions if compression.
  // Here, we differentiate the case with no compression (i.e. a_nbIndices==1) from the
  // compression case, because it will avoid to perform useless computation without
  // compression. However, it produces some duplication of code parts as a compromise.
  // ---------------------------------------------------------------------------------------

  // ______________________________________________________________________________________
  // Case1: no compression (i.e. a_nbIndices==1)
  if (a_nbIndices==1)
  {
    // Set indices of the line
    ap_ProjectionLine->SetIndex1(((int)(ap_index1[0])));
    ap_ProjectionLine->SetIndex2(((int)(ap_index2[0])));
    // Get positions and orientations from the scanner (mean depth of interaction intrinsicaly taken into account), taking POI into account if any
    if (mp_Scanner->GetPositionsAndOrientations( ((int)(ap_index1[0])), ((int)(ap_index2[0])),
                                                 ap_ProjectionLine->GetPosition1(), ap_ProjectionLine->GetPosition2(),
                                                 ap_ProjectionLine->GetOrientation1(), ap_ProjectionLine->GetOrientation2(),
                                                 ap_ProjectionLine->GetPOI1(), ap_ProjectionLine->GetPOI2() ))
    {
      Cerr("***** vProjector::Project() -> A problem occurred while getting positions and orientations from scanner !" << endl);
      return 1;
    }
  }
  // ______________________________________________________________________________________
  // Case2: compression (i.e. a_nbIndices>1)
  else
  {
    // Set default indices of the line to -1 when compression
    ap_ProjectionLine->SetIndex1(-1);
    ap_ProjectionLine->SetIndex2(-1);
    // Buffer pointers for positions and orientations handled by the projection line
    FLTNB* position1 = ap_ProjectionLine->GetPosition1();
    FLTNB* position2 = ap_ProjectionLine->GetPosition2();
    FLTNB* orientation1 = ap_ProjectionLine->GetOrientation1();
    FLTNB* orientation2 = ap_ProjectionLine->GetOrientation2();
    FLTNB* buffer_position1 = ap_ProjectionLine->GetBufferPosition1();
    FLTNB* buffer_position2 = ap_ProjectionLine->GetBufferPosition2();
    FLTNB* buffer_orientation1 = ap_ProjectionLine->GetBufferOrientation1();
    FLTNB* buffer_orientation2 = ap_ProjectionLine->GetBufferOrientation2();
    // Zero the position and orientation
    for (int i=0; i<3; i++)
    {
      position1[i] = 0.;
      position2[i] = 0.;
      orientation1[i] = 0.;
      orientation2[i] = 0.;
    }
    // Loop on provided indices
    for (int l=0; l<a_nbIndices; l++)
    {
      // Get positions and orientations from scanner (mean depth of interaction intrinsicaly taken into account), taking POI into account if any
      if (mp_Scanner->GetPositionsAndOrientations( ((int)(ap_index1[l])), ((int)(ap_index2[l])),
                                                   buffer_position1, buffer_position2,
                                                   buffer_orientation1, buffer_orientation2,
                                                   ap_ProjectionLine->GetPOI1(), ap_ProjectionLine->GetPOI2() ))
      {
        Cerr("***** vProjector::Project() -> A problem occurred while getting positions and orientations from scanner !" << endl);
        return 1;
      }
      // Add those contributions to the mean position and orientation
      for (int i=0; i<3; i++)
      {
        position1[i] += buffer_position1[i];
        position2[i] += buffer_position2[i];
        orientation1[i] += buffer_orientation1[i];
        orientation2[i] += buffer_orientation2[i];
      }
    }
    // Average positions and orientations
    for (int i=0; i<3; i++)
    {
      position1[i] /= ((FLTNB)a_nbIndices);
      position2[i] /= ((FLTNB)a_nbIndices);
      orientation1[i] /= ((FLTNB)a_nbIndices);
      orientation2[i] /= ((FLTNB)a_nbIndices);
    }
  }

  // --------------------------------------------------------------
  // Second: Modify the end points coordinates from common options,
  // random, offset, and LOR displacements.
  // -----------------------------------------------------------

  // Apply common options TODO

  // Apply LORs displacement TODO

  // Apply global image offset
  ap_ProjectionLine->ApplyOffset();

  // Apply bed position offset
  ap_ProjectionLine->ApplyBedOffset();

  // -----------------------------------------------------------
  // Third: project the line
  // -----------------------------------------------------------

  // Compute LOR length
  ap_ProjectionLine->ComputeLineLength();

  // Switch on different TOF options
  switch (m_TOFMethod)
  {
    case USE_NOTOF:
      if (ProjectWithoutTOF( a_direction, ap_ProjectionLine ))
      {
        Cerr("***** vProjector::Project() -> A problem occurred while projecting a line without time-of-flight !" << endl);
        return 1;
      }
      break;
    case USE_TOFLIST:
      if (ProjectTOFListmode( a_direction, ap_ProjectionLine ))
      {
        Cerr("***** vProjector::Project() -> A problem occurred while projecting a line with time-of-flight position !" << endl);
        return 1;
      }
      break;
    case USE_TOFHISTO:
      if (ProjectTOFHistogram( a_direction, ap_ProjectionLine ))
      {
        Cerr("***** vProjector::Project() -> A problem occurred while projecting a line with binned time-of-flight !" << endl);
        return 1;
      }
      break;
    // No default
  }

  #ifdef CASTOR_VERBOSE
  if (m_verbose>=10)
  {
    Cout("vProjector::Project() -> Exit function" << endl);
  }
  #endif

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

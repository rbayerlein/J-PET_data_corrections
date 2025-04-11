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
  \brief    Implementation of class iProjectorIncrementalSiddonMulti
*/

#include "iProjectorIncrementalSiddonMulti.hh"
#include "sOutputManager.hh"
#ifdef _WIN32
// Avoid compilation errors due to mix up between std::min()/max() and
// min max macros
#undef min
#undef max
#endif

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iProjectorIncrementalSiddonMulti::iProjectorIncrementalSiddonMulti() : vProjector()
{
  m_nbSensLines = 0;
  m_nbRecoLines = 0;
  // This projector is not compatible with SPECT attenuation correction because
  // the voxels contributing to the line are not strictly ordered with respect to
  // their distance to point 2 (due to the use of multiple lines that are
  // stack one after the other)
  m_compatibleWithSPECTAttenuationCorrection = false;
  // This projector is not compatible with compression as it works only 
  // with the detection element indices
  m_compatibleWithCompression = false;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iProjectorIncrementalSiddonMulti::~iProjectorIncrementalSiddonMulti()
{
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorIncrementalSiddonMulti::ReadConfigurationFile(const string& a_configurationFile)
{
  // Read the transaxial FWHM option
  string key_word = "number of lines sensitivity";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_nbSensLines, 1, KEYWORD_OPTIONAL))
  {
    Cerr("***** iProjectorIncrementalSiddonMulti::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  
  key_word = "number of lines reconstruction";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_nbRecoLines, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iProjectorIncrementalSiddonMulti::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  
  // If not initialized, set the number of line used for sens image generation 
  // to the number of lines used for reconstruction
  m_nbSensLines = (m_nbSensLines == 0) ? m_nbRecoLines : m_nbSensLines;
   
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorIncrementalSiddonMulti::ReadOptionsList(const string& a_optionsList)
{
  // Read them
  uint16_t option[2];
  if (ReadStringOption(a_optionsList, option, 2, ",", "Multi-Siddon configuration"))
  {
    Cerr("***** iProjectorIncrementalSiddonMulti::ReadConfigurationFile() -> Failed to correctly read the list of options !" << endl);
    return 1;
  }
  
  m_nbSensLines = option[0];
  m_nbRecoLines = option[1];
  
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iProjectorIncrementalSiddonMulti::ShowHelpSpecific()
{
  cout << "This projector uses multiple ray-tracing for a single event in order to estimate the solid angle contribution." << endl;
  cout << "For each line of an event, the end points of the line are randomly chosen inside the detector element." << endl;
  cout << "For PET systems, if a mean depth of interaction is given, the distribution of the random points will be centered on the provided value" << endl << endl;
  cout << "The ray-tracing is performed with the incremental Siddon algorithm (see incrementalSiddon projector)." << endl;
  cout << "!! IMPORTANT: This implementation has not been assessed for systems using detectors with non-zero axial angle yet." << endl;
  cout << "There are 2 parameters for this projector:" << endl;
  cout << "  number of lines sensitivity: the number of lines used per event for the sensitivity image generation (list-mode only)." << endl;
  cout << "                               this parameter is ignored for histogram format." << endl;
  cout << "  number of lines reconstruction: the number of lines used per event for the reconstruction (and histogram sensitivity image generation)" << endl << endl;
  cout << "This projector can be initialized with a configuration file containing the 'number of lines sensitivity' and 'number of lines reconstruction' keywords and values" << endl;
  cout << "This projector can be initialized with command line options (1st option: number of lines sensitivity, 2nd option: number of lines reconstruction" << endl;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorIncrementalSiddonMulti::CheckSpecificParameters()
{
  // Send an error if less than 1 line
  if (m_nbSensLines<1 || m_nbRecoLines<1)
  {
    Cerr("***** iProjectorIncrementalSiddonMulti::CheckSpecificParameters() -> The provided number of lines is less than 1 !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorIncrementalSiddonMulti::InitializeSpecific()
{
  // Verbose
  if (m_verbose>=2) Cout("iProjectorIncrementalSiddonMulti::Initialize() -> Use incremental Siddon projector with " << m_nbRecoLines << " lines per event in reconstruction" << endl);
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

INTNB iProjectorIncrementalSiddonMulti::EstimateMaxNumberOfVoxelsPerLine()
{
  // Find the maximum number of voxels along a given dimension
  INTNB max_nb_voxels_in_dimension = mp_ImageDimensionsAndQuantification->GetNbVoxX();
  if (mp_ImageDimensionsAndQuantification->GetNbVoxY()>max_nb_voxels_in_dimension) max_nb_voxels_in_dimension = mp_ImageDimensionsAndQuantification->GetNbVoxY();
  if (mp_ImageDimensionsAndQuantification->GetNbVoxZ()>max_nb_voxels_in_dimension) max_nb_voxels_in_dimension = mp_ImageDimensionsAndQuantification->GetNbVoxZ();
  // We should have at most 4 voxels in a given plane, so multiply by 4
  // (note: this is not true however it ensures no overflow and is already quite optimized for RAM usage !)
  max_nb_voxels_in_dimension *= 4;
  // Finally multiply by the number of lines
  max_nb_voxels_in_dimension *= (INTNB)(std::max)( m_nbSensLines, m_nbRecoLines );
  // Return the value
  return max_nb_voxels_in_dimension;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorIncrementalSiddonMulti::ProjectWithoutTOF(int a_direction, oProjectionLine* ap_ProjectionLine )
{
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** iProjectorIncrementalSiddonMulti::ProjectWithoutTOF() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif

  #ifdef CASTOR_VERBOSE
  if (m_verbose>=10)
  {
    string direction = "";
    if (a_direction==FORWARD) direction = "forward";
    else direction = "backward";
    Cout("iProjectorIncrementalSiddonMulti::ProjectWithoutTOF() -> Project line '" << ap_ProjectionLine << "' in " << direction << " direction" << endl);
  }
  #endif

  // For list-mode, the number of lines to use depends on the processing stage (sensitivity generation or reconstruction)
  uint16_t nb_lines = (m_sensitivityMode == true) ?
                       m_nbSensLines : 
                       m_nbRecoLines ;

  // Loop on the number of lines
  for (uint16_t line=0; line<nb_lines; line++)
  {       
    // Get random positions from the scanner (mean depth of interaction intrinsicaly taken into account), taking POI into account if any
    if (mp_Scanner->GetRdmPositionsAndOrientations( ap_ProjectionLine->GetIndex1(), ap_ProjectionLine->GetIndex2(), 
                                                    ap_ProjectionLine->GetPosition1(), ap_ProjectionLine->GetPosition2(),
                                                    ap_ProjectionLine->GetOrientation1(), ap_ProjectionLine->GetOrientation2() ))
    {
      Cerr("***** iProjectorIncrementalSiddonMulti::ProjectWithoutTOF() -> A problem occurred while getting positions and orientations from scanner !" << endl);
      return 1;
    }                                       
    
    // We must called here ApplyOffset() and ComputeLineLength() which 
    // are supposed to be called in vProjector::Project()
    // This is a temporary solution as everyting related to multi-lines 
    // management should be performed in vProjector 
    
    // Apply global image offset
    ap_ProjectionLine->ApplyOffset();

    // Apply bed position offset
    ap_ProjectionLine->ApplyBedOffset();
  
    // Get end points position
    FLTNB* event1Float = ap_ProjectionLine->GetPosition1();
    FLTNB* event2Float = ap_ProjectionLine->GetPosition2();
    HPFLTNB event1[3] = { event1Float[0], event1Float[1], event1Float[2] };
    HPFLTNB event2[3] = { event2Float[0], event2Float[1], event2Float[2] };
  
    // **************************************
    // STEP 1: LOR length calculation
    // **************************************
    // Must recompute LOR length for each new line of response
    ap_ProjectionLine->ComputeLineLength();
    HPFLTNB length_LOR = ap_ProjectionLine->GetLength();
  
    // **************************************
    // STEP 2: Compute entrance voxel indexes
    // **************************************
    HPFLTNB alphaFirst[] = { 0.0, 0.0, 0.0 };
    HPFLTNB alphaLast[] = { 0.0, 0.0, 0.0 };
  
    HPFLTNB alphaMin = 0.0f, alphaMax = 1.0f;
    HPFLTNB delta_pos[] = {
      event2[ 0 ] - event1[ 0 ],
      event2[ 1 ] - event1[ 1 ],
      event2[ 2 ] - event1[ 2 ]
    };
  
    // Computation of alphaMin et alphaMax values (entrance and exit point of the ray)
    for( int i = 0; i < 3; ++i )
    {
      if( delta_pos[ i ] != 0.0 )
      {
        alphaFirst[i] = (-mp_halfFOV[i] - event1[i]) / (delta_pos[ i ]);
        alphaLast[i]  = ( mp_halfFOV[i] - event1[i]) / (delta_pos[ i ]);
        alphaMin = (std::max)(alphaMin,(std::min)(alphaFirst[i],alphaLast[i]));
        alphaMax = (std::min)(alphaMax,(std::max)(alphaFirst[i],alphaLast[i]));
      }
    }
  
    // if alphaMax is less than or equal to alphaMin no intersection
    // and return an empty buffer
    if( alphaMax <= alphaMin ) continue;
  
    // Now we have to find the indices of the particular plane
    // (iMin,iMax), (jMin,jMax), (kMin,kMax)
    int iMin = 0, iMax = 0;
    int jMin = 0, jMax = 0;
    int kMin = 0, kMax = 0;
  
    // For the X-axis
    if( delta_pos[ 0 ] > 0.0f)
    {
      iMin = ::ceil( ( mp_nbVox[ 0 ] + 1 ) - ( mp_halfFOV[ 0 ] - alphaMin * delta_pos[ 0 ] - event1[ 0 ] ) / mp_sizeVox[0] );
      iMax = ::floor( 1 + ( event1[ 0 ] + alphaMax * delta_pos[ 0 ] - (-mp_halfFOV[ 0 ]) ) / mp_sizeVox[0] );
    }
    else if( delta_pos[ 0 ] < 0.0 )
    {
      iMin = ::ceil( ( mp_nbVox[ 0 ] + 1 ) - ( mp_halfFOV[ 0 ] - alphaMax * delta_pos[ 0 ] - event1[ 0 ] ) / mp_sizeVox[0] );
      iMax = ::floor( 1 + ( event1[ 0 ] + alphaMin * delta_pos[ 0 ] - (-mp_halfFOV[ 0 ]) ) / mp_sizeVox[0] );
    }
    if( delta_pos[ 0 ] == 0 )
    {
      iMin = 1, iMax = 0;
    }
  
    // For the Y-axis
    if( delta_pos[ 1 ] > 0 )
    {
      jMin = ::ceil( ( mp_nbVox[ 1 ] + 1 ) - ( mp_halfFOV[ 1 ] - alphaMin * delta_pos[ 1 ] - event1[ 1 ] ) / mp_sizeVox[1] );
      jMax = ::floor( 1 + ( event1[ 1 ] + alphaMax * delta_pos[ 1 ] - (-mp_halfFOV[ 1 ]) ) / mp_sizeVox[1] );
    }
    else if( delta_pos[ 1 ] < 0 )
    {
      jMin = ::ceil( ( mp_nbVox[ 1 ] + 1 ) - ( mp_halfFOV[ 1 ] - alphaMax * delta_pos[ 1 ] - event1[ 1 ] ) / mp_sizeVox[1] );
      jMax = ::floor( 1 + ( event1[ 1 ] + alphaMin * delta_pos[ 1 ] - (-mp_halfFOV[ 1 ]) ) / mp_sizeVox[1] );
    }
    else if( delta_pos[ 1 ] == 0 )
    {
      jMin = 1, jMax = 0;
    }
  
    // For the Z-axis
    if( delta_pos[ 2 ] > 0 )
    {
      kMin = ::ceil( ( mp_nbVox[ 2 ] + 1 ) - ( mp_halfFOV[ 2 ] - alphaMin * delta_pos[ 2 ] - event1[ 2 ] ) / mp_sizeVox[2] );
      kMax = ::floor( 1 + ( event1[ 2 ] + alphaMax * delta_pos[ 2 ] - (-mp_halfFOV[ 2 ]) ) / mp_sizeVox[2] );
    }
    else if( delta_pos[ 2 ] < 0 )
    {
      kMin = ::ceil( ( mp_nbVox[ 2 ] + 1 ) - ( mp_halfFOV[ 2 ] - alphaMax * delta_pos[ 2 ] - event1[ 2 ] ) / mp_sizeVox[2] );
      kMax = ::floor( 1 + ( event1[ 2 ] + alphaMin * delta_pos[ 2 ] - (-mp_halfFOV[ 2 ]) ) / mp_sizeVox[2] );
    }
    else if( delta_pos[ 2 ] == 0 )
    {
      kMin = 1, kMax = 0;
    }
  
    // Computing the last term n number of intersection
    INTNB n = ( iMax - iMin + 1 ) + ( jMax - jMin + 1 )
      + ( kMax - kMin + 1 );
  
    // Array storing the first alpha in X, Y and Z
    HPFLTNB alpha_XYZ[ 3 ] = { 1.0, 1.0, 1.0 };
  
    // Computing the first alpha in X
    if( delta_pos[ 0 ] > 0 )
      alpha_XYZ[ 0 ] = ( ( (-mp_halfFOV[ 0 ]) + ( iMin - 1 ) * mp_sizeVox[0] ) - event1[ 0 ] ) / delta_pos[ 0 ];
    else if( delta_pos[ 0 ] < 0 )
      alpha_XYZ[ 0 ] = ( ( (-mp_halfFOV[ 0 ]) + ( iMax - 1 ) * mp_sizeVox[0] ) - event1[ 0 ] ) / delta_pos[ 0 ];
  
    // Computing the first alpha in Y
    if( delta_pos[ 1 ] > 0 )
      alpha_XYZ[ 1 ] = ( ( (-mp_halfFOV[ 1 ]) + ( jMin - 1 ) * mp_sizeVox[1] ) - event1[ 1 ] ) / delta_pos[ 1 ];
    else if( delta_pos[ 1 ] < 0 )
      alpha_XYZ[ 1 ] = ( ( (-mp_halfFOV[ 1 ]) + ( jMax - 1 ) * mp_sizeVox[1] ) - event1[ 1 ] ) / delta_pos[ 1 ];
  
    // Computing the first alpha in Z
    if( delta_pos[ 2 ] > 0 )
      alpha_XYZ[ 2 ] = ( ( (-mp_halfFOV[ 2 ]) + ( kMin - 1 ) * mp_sizeVox[2] ) - event1[ 2 ] ) / delta_pos[ 2 ];
    else if( delta_pos[ 2 ] < 0 )
      alpha_XYZ[ 2 ] = ( ( (-mp_halfFOV[ 2 ]) + ( kMax - 1 ) * mp_sizeVox[2] ) - event1[ 2 ] ) / delta_pos[ 2 ];
  
    // Computing the alpha updating
    HPFLTNB alpha_u[ 3 ] = {
      mp_sizeVox[0] / std::fabs( delta_pos[ 0 ] ),
      mp_sizeVox[1] / std::fabs( delta_pos[ 1 ] ),
      mp_sizeVox[2] / std::fabs( delta_pos[ 2 ] )
    };
  
    // Computing the index updating 
    INTNB index_u[ 3 ] = {
      (delta_pos[ 0 ] < 0) ? -1 : 1,
      (delta_pos[ 1 ] < 0) ? -1 : 1,
      (delta_pos[ 2 ] < 0) ? -1 : 1
    };
  
    // Check which alpha is the min/max and increment
    if( alpha_XYZ[ 0 ] == alphaMin ) alpha_XYZ[ 0 ] += alpha_u[ 0 ];
    if( alpha_XYZ[ 1 ] == alphaMin ) alpha_XYZ[ 1 ] += alpha_u[ 1 ];
    if( alpha_XYZ[ 2 ] == alphaMin ) alpha_XYZ[ 2 ] += alpha_u[ 2 ];
  
    // Computing the minimum value in the alpha_XYZ buffer
    HPFLTNB const min_alpha_XYZ = (std::min)( alpha_XYZ[ 0 ],
      (std::min)( alpha_XYZ[ 1 ], alpha_XYZ[ 2 ] ) );
  
    // Computing the first index of intersection
    HPFLTNB const alphaMid = ( min_alpha_XYZ + alphaMin ) / 2.0f;
    INTNB index_ijk[ 3 ] = {
      1 + (int)( ( event1[ 0 ] + alphaMid * delta_pos[ 0 ] - (-mp_halfFOV[ 0 ]) ) / mp_sizeVox[0] ),
      1 + (int)( ( event1[ 1 ] + alphaMid * delta_pos[ 1 ] - (-mp_halfFOV[ 1 ]) ) / mp_sizeVox[1] ),
      1 + (int)( ( event1[ 2 ] + alphaMid * delta_pos[ 2 ] - (-mp_halfFOV[ 2 ]) ) / mp_sizeVox[2] )
    };
  
    INTNB const w = mp_nbVox[ 0 ];
    INTNB const wh = w * mp_nbVox[ 1 ];
  
    // Loop over the number of plane to cross
    HPFLTNB alpha_c = alphaMin;
    FLTNB coeff = 0.0f;
    INTNB numVox = 0;
    
    for( int nP = 0; nP < n - 1; ++nP )
    {
      if( ( alpha_XYZ[ 0 ] <= alpha_XYZ[ 1 ] )
       && ( alpha_XYZ[ 0 ] <= alpha_XYZ[ 2 ] ) )
      {
        // Storing values
        if( ( alpha_XYZ[ 0 ] >= alphaMin )
         && ( index_ijk[ 0 ] - 1 >= 0 )
         && ( index_ijk[ 0 ] - 1 <= mp_nbVox[ 0 ] - 1 )
         && ( index_ijk[ 1 ] - 1 >= 0 )
         && ( index_ijk[ 1 ] - 1 <= mp_nbVox[ 1 ] - 1 )
         && ( index_ijk[ 2 ] - 1 >= 0 )
         && ( index_ijk[ 2 ] - 1 <= mp_nbVox[ 2 ] - 1 ) )
        {
          coeff = ( alpha_XYZ[ 0 ] - alpha_c ) * length_LOR;
          numVox = ( index_ijk[ 0 ] - 1 ) + ( ( index_ijk[ 1 ] - 1 ) ) * w + ( ( index_ijk[ 2 ] - 1 ) ) * wh;
          // if this voxel is masked, do nothing
          if (!m_hasMask || mp_mask[numVox]) ap_ProjectionLine->AddVoxel(a_direction, numVox, coeff/((FLTNB)nb_lines) );
        }
        // Increment values
        alpha_c = alpha_XYZ[ 0 ];
        alpha_XYZ[ 0 ] += alpha_u[ 0 ];
        index_ijk[ 0 ] += index_u[ 0 ];
      }
      else if( ( alpha_XYZ[ 1 ] < alpha_XYZ[ 0 ] )
            && ( alpha_XYZ[ 1 ] <= alpha_XYZ[ 2 ] ) )
      {
        // Storing values
        if( ( alpha_XYZ[ 1 ] >= alphaMin )
         && ( index_ijk[ 0 ] - 1 >= 0 )
         && ( index_ijk[ 0 ] - 1 <= mp_nbVox[ 0 ] - 1 )
         && ( index_ijk[ 1 ] - 1 >= 0 )
         && ( index_ijk[ 1 ] - 1 <= mp_nbVox[ 1 ] - 1 )
         && ( index_ijk[ 2 ] - 1 >= 0 )
         && ( index_ijk[ 2 ] - 1 <= mp_nbVox[ 2 ] - 1 ) )
        {
          coeff = ( alpha_XYZ[ 1 ] - alpha_c ) * length_LOR;
          numVox = ( index_ijk[ 0 ] - 1 ) + ( ( index_ijk[ 1 ] - 1 ) ) * w + ( ( index_ijk[ 2 ] - 1 ) ) * wh;
          // if this voxel is masked, do nothing
          if (!m_hasMask || mp_mask[numVox]) ap_ProjectionLine->AddVoxel(a_direction, numVox, coeff/((FLTNB)nb_lines) );
        }
        // Increment values
        alpha_c = alpha_XYZ[ 1 ];
        alpha_XYZ[ 1 ] += alpha_u[ 1 ];
        index_ijk[ 1 ] += index_u[ 1 ];
      }
      else if( ( alpha_XYZ[ 2 ] < alpha_XYZ[ 0 ] )
            && ( alpha_XYZ[ 2 ] < alpha_XYZ[ 1 ] ) )
      {
        // Storing values
        if( ( alpha_XYZ[ 2 ] >= alphaMin )
         && ( index_ijk[ 0 ] - 1 >= 0 )
         && ( index_ijk[ 0 ] - 1 <= mp_nbVox[ 0 ] - 1 )
         && ( index_ijk[ 1 ] - 1 >= 0 )
         && ( index_ijk[ 1 ] - 1 <= mp_nbVox[ 1 ] - 1 )
         && ( index_ijk[ 2 ] - 1 >= 0 )
         && ( index_ijk[ 2 ] - 1 <= mp_nbVox[ 2 ] - 1 ) )
        {
          coeff = ( alpha_XYZ[ 2 ] - alpha_c ) * length_LOR;
          numVox = ( index_ijk[ 0 ] - 1 ) + ( ( index_ijk[ 1 ] - 1 ) ) * w + ( ( index_ijk[ 2 ] - 1 ) ) * wh;
          // if this voxel is masked, do nothing
          if (!m_hasMask || mp_mask[numVox]) ap_ProjectionLine->AddVoxel( a_direction, numVox, coeff/((FLTNB)nb_lines) );
        }
  
        // Increment values
        alpha_c = alpha_XYZ[ 2 ];
        alpha_XYZ[ 2 ] += alpha_u[ 2 ];
        index_ijk[ 2 ] += index_u[ 2 ];
      }
    }                                      
  }

  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorIncrementalSiddonMulti::ProjectTOFListmode(int a_direction, oProjectionLine* ap_ProjectionLine)
{
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** iProjectorIncrementalSiddonMulti::ProjectTOFListmode() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif

  #ifdef CASTOR_VERBOSE
  if (m_verbose>=10)
  {
    string direction = "";
    if (a_direction==FORWARD) direction = "forward";
    else direction = "backward";
    Cout("iProjectorIncrementalSiddonMulti::ProjectTOFListmode() -> Project line '" << ap_ProjectionLine << "' in " << direction << " direction" << endl);
  }
  #endif

  // For list-mode, the number of lines to use depends on the processing stage (sensitivity generation or reconstruction)
  uint16_t nb_lines = (m_sensitivityMode == true) ?
                       m_nbSensLines : 
                       m_nbRecoLines ;

  // Loop on the number of lines
  for (uint16_t line=0; line<nb_lines; line++)
  {       
    // Get random positions from the scanner (mean depth of interaction intrinsicaly taken into account), taking POI into account if any
    if (mp_Scanner->GetRdmPositionsAndOrientations( ap_ProjectionLine->GetIndex1(), ap_ProjectionLine->GetIndex2(), 
                                                    ap_ProjectionLine->GetPosition1(), ap_ProjectionLine->GetPosition2(),
                                                    ap_ProjectionLine->GetOrientation1(), ap_ProjectionLine->GetOrientation2() ))
    {
      Cerr("***** iProjectorIncrementalSiddonMulti::ProjectTOFListmode() -> A problem occurred while getting positions and orientations from scanner !" << endl);
      return 1;
    }                                       
    
    // We must called here ApplyOffset() and ComputeLineLength() which 
    // are supposed to be called in vProjector::Project()
    // This is a temporary solution as everyting related to multi-lines 
    // management should be performed in vProjector 
    
    // Apply global image offset
    ap_ProjectionLine->ApplyOffset();

    // Apply bed position offset
    ap_ProjectionLine->ApplyBedOffset();
  
    // Get end points position
    FLTNB* event1Float = ap_ProjectionLine->GetPosition1();
    FLTNB* event2Float = ap_ProjectionLine->GetPosition2();
    HPFLTNB event1[3] = { event1Float[0], event1Float[1], event1Float[2] };
    HPFLTNB event2[3] = { event2Float[0], event2Float[1], event2Float[2] };
  
    // **************************************
    // STEP 1: LOR length calculation
    // **************************************
    // Must recompute LOR length for each new line of response
    ap_ProjectionLine->ComputeLineLength();
    HPFLTNB length_LOR = ap_ProjectionLine->GetLength();
  
    // **************************************
    // STEP 2: Compute entrance voxel indexes
    // **************************************
    HPFLTNB alphaFirst[] = { 0.0, 0.0, 0.0 };
    HPFLTNB alphaLast[] = { 0.0, 0.0, 0.0 };
  
    HPFLTNB alphaMin = 0.0f, alphaMax = 1.0f;
    HPFLTNB delta_pos[] = {
      event2[ 0 ] - event1[ 0 ],
      event2[ 1 ] - event1[ 1 ],
      event2[ 2 ] - event1[ 2 ]
    };

    // Get TOF info
    HPFLTNB tof_measurement = ap_ProjectionLine->GetTOFMeasurementInPs();

    // convert delta time into delta length
    HPFLTNB tof_delta = tof_measurement * SPEED_OF_LIGHT_IN_MM_PER_PS * 0.5;

    // TOF Gaussian standard deviation and truncation
    HPFLTNB tof_sigma = m_TOFResolutionInMm / TWO_SQRT_TWO_LN_2;
    HPFLTNB tof_sigma_sqrt2 = sqrt(2.)*tof_sigma;
    HPFLTNB tof_half_span = 0.;

    if (m_TOFWeightingFcnPrecomputedFlag) tof_half_span = ((HPFLTNB)m_TOFWeightingFcnNbSamples)/(2.*m_TOFPrecomputedSamplingFactor);
    else if (m_TOFBinProperProcessingFlag) tof_half_span = tof_sigma * m_TOFNbSigmas + 0.5 * m_TOFBinSizeInMm;
    else tof_half_span = tof_sigma * m_TOFNbSigmas;

    // distance between the first event1 and the TOF measurement
    HPFLTNB lor_tof_center = length_LOR * 0.5 + tof_delta;

    // index along each axis of the first voxel falling inside the truncated TOF distribution centered at the TOF measurement
    HPFLTNB tof_edge_low[] = {0,0,0};
    // index along each axis of the last voxel falling inside the truncated TOF distribution centered at the TOF measurement
    HPFLTNB tof_edge_high[] = {0,0,0};
    HPFLTNB tof_center;
    INTNB tof_index_low[] = {0,0,0};
    INTNB tof_index_high[] = {0,0,0};

    // low/high voxel edges (in absolute coordinates) for truncated TOF
    for (int ax=0;ax<3;ax++)
    {
      // absolute coordinate along each axis of the center of the TOF distribution
      tof_center = event1[ax] +  lor_tof_center * delta_pos[ax] / length_LOR;

      // absolute coordinate along each axis of the lowest voxel edge spanned by the TOF distribution, limited by the lowest FOV edge
      tof_edge_low[ax] = tof_center - tof_half_span * fabs(delta_pos[ax]) / length_LOR;
      tof_index_low[ax] = max( (INTNB)::floor( (tof_edge_low[ax] - (-mp_halfFOV[ax])) / mp_sizeVox[ax] ), (INTNB)0);
      // if low TOF edge above the highest FOV edge, return empty line
      if (tof_index_low[ax]>mp_nbVox[ax]-1) return 0;
      tof_edge_low[ax] = (HPFLTNB)tof_index_low[ax] *  mp_sizeVox[ax] - mp_halfFOV[ax];

      // absolute coordinate along each axis of the highest voxel edge spanned by the TOF distribution, limited by the highest FOV edge
      tof_edge_high[ax] = tof_center + tof_half_span * fabs(delta_pos[ax]) / length_LOR;
      tof_index_high[ax] = min( (INTNB)::floor( (tof_edge_high[ax] - (-mp_halfFOV[ax])) / mp_sizeVox[ax] ), mp_nbVox[ax]-1);
      // if high TOF edge below the lowest FOV edge, return empty line
      if (tof_index_high[ax]<0) return 0;
      tof_edge_high[ax] = (HPFLTNB)(tof_index_high[ax]+1) * mp_sizeVox[ax] - mp_halfFOV[ax];
    }

    // Computation of alphaMin et alphaMax values entrance and exit point of the ray at voxel grid (edges),
    // with respect to the truncated TOF distribution centered at the TOF measurement,
    // limited by the FOV, absolute normalized distance from event1
    for( int i = 0; i < 3; ++i )
    {
      if( delta_pos[ i ] != 0.0 )
      {
        alphaFirst[i] = (-mp_halfFOV[i] - event1[i]) / (delta_pos[ i ]);
        alphaLast[i]  = ( mp_halfFOV[i] - event1[i]) / (delta_pos[ i ]);
        alphaMin = (std::max)(alphaMin,(std::min)(alphaFirst[i],alphaLast[i]));
        alphaMax = (std::min)(alphaMax,(std::max)(alphaFirst[i],alphaLast[i]));
      }
    }

    // if alphaMax is less than or equal to alphaMin no intersection
    // and return an empty buffer
    if( alphaMax <= alphaMin ) continue;

    // Now we have to find the indices of the particular plane
    // (iMin,iMax), (jMin,jMax), (kMin,kMax)
    int iMin = 0, iMax = 0;
    int jMin = 0, jMax = 0;
    int kMin = 0, kMax = 0;

    // For the X-axis
    if( delta_pos[ 0 ] > 0.0f )
    {
      iMin = ::ceil( ( mp_nbVox[ 0 ] + 1 ) - ( mp_halfFOV[ 0 ] - alphaMin * delta_pos[ 0 ] - event1[ 0 ] ) / mp_sizeVox[0] );
      iMax = ::floor( 1 + ( event1[ 0 ] + alphaMax * delta_pos[ 0 ] - (-mp_halfFOV[ 0 ]) ) / mp_sizeVox[0] );
    }
    else if( delta_pos[ 0 ] < 0.0f )
    {
      iMin = ::ceil( ( mp_nbVox[ 0 ] + 1 ) - ( mp_halfFOV[ 0 ] - alphaMax * delta_pos[ 0 ] - event1[ 0 ] ) / mp_sizeVox[0] );
      iMax = ::floor( 1 + ( event1[ 0 ] + alphaMin * delta_pos[ 0 ] - (-mp_halfFOV[ 0 ]) ) / mp_sizeVox[0] );
    }
    if( delta_pos[ 0 ] == 0 )
    {
      iMin = 1, iMax = 0;
    }

    // For the Y-axis
    if( delta_pos[ 1 ] > 0 )
    {
      jMin = ::ceil( ( mp_nbVox[ 1 ] + 1 ) - ( mp_halfFOV[ 1 ] - alphaMin * delta_pos[ 1 ] - event1[ 1 ] ) / mp_sizeVox[1] );
      jMax = ::floor( 1 + ( event1[ 1 ] + alphaMax * delta_pos[ 1 ] - (-mp_halfFOV[ 1 ]) ) / mp_sizeVox[1] );
    }
    else if( delta_pos[ 1 ] < 0 )
    {
      jMin = ::ceil( ( mp_nbVox[ 1 ] + 1 ) - ( mp_halfFOV[ 1 ] - alphaMax * delta_pos[ 1 ] - event1[ 1 ] ) / mp_sizeVox[1] );
      jMax = ::floor( 1 + ( event1[ 1 ] + alphaMin * delta_pos[ 1 ] - (-mp_halfFOV[ 1 ]) ) / mp_sizeVox[1] );
    }
    else if( delta_pos[ 1 ] == 0 )
    {
      jMin = 1, jMax = 0;
    }

    // For the Z-axis
    if( delta_pos[ 2 ] > 0 )
    {
      kMin = ::ceil( ( mp_nbVox[ 2 ] + 1 ) - ( mp_halfFOV[ 2 ] - alphaMin * delta_pos[ 2 ] - event1[ 2 ] ) / mp_sizeVox[2] );
      kMax = ::floor( 1 + ( event1[ 2 ] + alphaMax * delta_pos[ 2 ] - (-mp_halfFOV[ 2 ]) ) / mp_sizeVox[2] );
    }
    else if( delta_pos[ 2 ] < 0 )
    {
      kMin = ::ceil( ( mp_nbVox[ 2 ] + 1 ) - ( mp_halfFOV[ 2 ] - alphaMax * delta_pos[ 2 ] - event1[ 2 ] ) / mp_sizeVox[2] );
      kMax = ::floor( 1 + ( event1[ 2 ] + alphaMin * delta_pos[ 2 ] - (-mp_halfFOV[ 2 ]) ) / mp_sizeVox[2] );
    }
    else if( delta_pos[ 2 ] == 0 )
    {
      kMin = 1, kMax = 0;
    }

    // Computing the last term n number of intersection
    INTNB n = ( iMax - iMin + 1 ) + ( jMax - jMin + 1 )
      + ( kMax - kMin + 1 );

    // Array storing the first alpha in X, Y and Z
    HPFLTNB alpha_XYZ[ 3 ] = { 1.0f, 1.0f, 1.0f };

    // Computing the first alpha in X
    if( delta_pos[ 0 ] > 0 )
      alpha_XYZ[ 0 ] = ( ( (-mp_halfFOV[ 0 ]) + ( iMin - 1 ) * mp_sizeVox[0] ) - event1[ 0 ] ) / delta_pos[ 0 ];
    else if( delta_pos[ 0 ] < 0 )
      alpha_XYZ[ 0 ] = ( ( (-mp_halfFOV[ 0 ]) + ( iMax - 1 ) * mp_sizeVox[0] ) - event1[ 0 ] ) / delta_pos[ 0 ];

    // Computing the first alpha in Y
    if( delta_pos[ 1 ] > 0 )
      alpha_XYZ[ 1 ] = ( ( (-mp_halfFOV[ 1 ]) + ( jMin - 1 ) * mp_sizeVox[1] ) - event1[ 1 ] ) / delta_pos[ 1 ];
    else if( delta_pos[ 1 ] < 0 )
      alpha_XYZ[ 1 ] = ( ( (-mp_halfFOV[ 1 ]) + ( jMax - 1 ) * mp_sizeVox[1] ) - event1[ 1 ] ) / delta_pos[ 1 ];

    // Computing the first alpha in Z
    if( delta_pos[ 2 ] > 0 )
      alpha_XYZ[ 2 ] = ( ( (-mp_halfFOV[ 2 ]) + ( kMin - 1 ) * mp_sizeVox[2] ) - event1[ 2 ] ) / delta_pos[ 2 ];
    else if( delta_pos[ 2 ] < 0 )
      alpha_XYZ[ 2 ] = ( ( (-mp_halfFOV[ 2 ]) + ( kMax - 1 ) * mp_sizeVox[2] ) - event1[ 2 ] ) / delta_pos[ 2 ];

    // Computing the alpha updating
    HPFLTNB alpha_u[ 3 ] = {
      mp_sizeVox[0] / std::fabs( delta_pos[ 0 ] ),
      mp_sizeVox[1] / std::fabs( delta_pos[ 1 ] ),
      mp_sizeVox[2] / std::fabs( delta_pos[ 2 ] )
    };

    // Computing the index updating
    INTNB index_u[ 3 ] = {
      (delta_pos[ 0 ] < 0) ? -1 : 1,
      (delta_pos[ 1 ] < 0) ? -1 : 1,
      (delta_pos[ 2 ] < 0) ? -1 : 1
    };
   
    // Check which alpha is the min/max and increment
    if( alpha_XYZ[ 0 ] == alphaMin ) alpha_XYZ[ 0 ] += alpha_u[ 0 ];
    if( alpha_XYZ[ 1 ] == alphaMin ) alpha_XYZ[ 1 ] += alpha_u[ 1 ];
    if( alpha_XYZ[ 2 ] == alphaMin ) alpha_XYZ[ 2 ] += alpha_u[ 2 ];

    // Computing the minimum value in the alpha_XYZ buffer
    HPFLTNB const min_alpha_XYZ = (std::min)( alpha_XYZ[ 0 ], (std::min)( alpha_XYZ[ 1 ], alpha_XYZ[ 2 ] ) );

    // Computing the first index of intersection
    HPFLTNB alphaMid = ( min_alpha_XYZ + alphaMin ) / 2.0f;
    INTNB index_ijk[ 3 ] = {
      1 + (INTNB)( ( event1[ 0 ] + alphaMid * delta_pos[ 0 ] - (-mp_halfFOV[ 0 ]) ) / mp_sizeVox[0] ),
      1 + (INTNB)( ( event1[ 1 ] + alphaMid * delta_pos[ 1 ] - (-mp_halfFOV[ 1 ]) ) / mp_sizeVox[1] ),
      1 + (INTNB)( ( event1[ 2 ] + alphaMid * delta_pos[ 2 ] - (-mp_halfFOV[ 2 ]) ) / mp_sizeVox[2] )
    };

    INTNB const w = mp_nbVox[ 0 ];
    INTNB const wh = w * mp_nbVox[ 1 ];

    //// tof : erf of first voxel edge - Gaussian center
    //HPFLTNB previous_edge_erf = 0.;
    //HPFLTNB next_edge_erf = 0.;
    //bool previous_edge_erf_initialized = false;

    // Loop over the number of plane to cross
    HPFLTNB alpha_c = alphaMin;
    FLTNB coeff = 0.0f;
    INTNB numVox = 0;
    for( int nP = 0; nP < n - 1; ++nP )
    {
      if( ( alpha_XYZ[ 0 ] <= alpha_XYZ[ 1 ] )
       && ( alpha_XYZ[ 0 ] <= alpha_XYZ[ 2 ] ) )
      {
        // Storing values
        if( ( alpha_XYZ[ 0 ] >= alphaMin )
         && ( index_ijk[ 0 ] - 1 >= tof_index_low[0] )
         && ( index_ijk[ 0 ] - 1 <= tof_index_high[0] )
         && ( index_ijk[ 1 ] - 1 >= tof_index_low[1] )
         && ( index_ijk[ 1 ] - 1 <= tof_index_high[1] )
         && ( index_ijk[ 2 ] - 1 >= tof_index_low[2] )
         && ( index_ijk[ 2 ] - 1 <= tof_index_high[2]  ) )
        {
          //// initialize if not initialized yet
          //if (!previous_edge_erf_initialized)
          //{
            //previous_edge_erf = erf( (length_LOR *alpha_c - lor_tof_center) / tof_sigma_sqrt2 );
            //previous_edge_erf_initialized = true;
          //}

          numVox = ( index_ijk[ 0 ] - 1 ) + ( ( index_ijk[ 1 ] - 1 ) ) * w + ( ( index_ijk[ 2 ] - 1 ) ) * wh;

          // voxel center
          alphaMid = ( alpha_XYZ[ 0 ] + alpha_c ) * 0.5;

          // non TOF projection coefficient
          coeff = ( alpha_XYZ[ 0 ] - alpha_c ) * length_LOR;

          if (m_TOFWeightingFcnPrecomputedFlag)
          {
            // fetch the value of the precomputed TOF weighting function (centered at the TOF measurement)
            // nearest to the projection of the voxel center onto the LOR
            INTNB temp = m_TOFWeightingFcnNbSamples/2 + (INTNB)round( (alphaMid * length_LOR - lor_tof_center) * m_TOFPrecomputedSamplingFactor );
            // multiply the non TOF projection coefficient with the TOF weight
            if (temp>=0 && temp<m_TOFWeightingFcnNbSamples) coeff *= mp_TOFWeightingFcn[temp];
            else coeff = 0.;
          }
          else
          {
            if (m_TOFBinProperProcessingFlag)
            {
              // integration between the edges of the current quantization TOF bin (centered at the TOF measurement)
              // of the Gaussian centered at the projection of the voxel center onto the LOR
              // normalized by the size of the bin so that the integral of the final TOF weighting function remains 1
              HPFLTNB temp_erf = std::min(tof_half_span, std::max(-tof_half_span, lor_tof_center - m_TOFBinSizeInMm/2. - alphaMid * length_LOR));
              HPFLTNB prev_erf = erf(temp_erf/tof_sigma_sqrt2);
              temp_erf = std::min(tof_half_span, std::max(-tof_half_span, lor_tof_center + m_TOFBinSizeInMm/2. - alphaMid * length_LOR));
              HPFLTNB new_erf = erf(temp_erf/tof_sigma_sqrt2);
              // multiply the non TOF projection coefficient with the TOF weight
              coeff *= 0.5 * fabs(new_erf - prev_erf) / m_TOFBinSizeInMm;
            }
            else
            {
              //// tof : erf of next voxel edge - Gaussian center
              //next_edge_erf = erf( (length_LOR * alpha_XYZ[ 0 ] - lor_tof_center) / tof_sigma_sqrt2 );
              //// integration of the Gaussian is done on the LOR portion matching the whole current voxel along X axis
              //coeff = 0.5 * std::fabs(previous_edge_erf - next_edge_erf) ;
              //// keeping record of the previous edge, so as to save 1 erf computation
              //previous_edge_erf = next_edge_erf;

              // TOF weight = value of the normalized Gaussian (centered at the TOF measurement)
              // at the projection of the voxel center onto the LOR
              HPFLTNB temp =  (alphaMid * length_LOR - lor_tof_center)  / tof_sigma;
              // multiply the non TOF projection coefficient with the TOF weight
              coeff *=  exp(- 0.5 * temp * temp ) * m_TOFGaussianNormCoef;
            }
          }

          if (!m_hasMask || mp_mask[numVox])  ap_ProjectionLine->AddVoxel(a_direction, numVox, coeff);
        }

        // Increment values
        alpha_c = alpha_XYZ[ 0 ];
        alpha_XYZ[ 0 ] += alpha_u[ 0 ];
        index_ijk[ 0 ] += index_u[ 0 ];
      }
      else if( ( alpha_XYZ[ 1 ] < alpha_XYZ[ 0 ] )
            && ( alpha_XYZ[ 1 ] <= alpha_XYZ[ 2 ] ) )
      {
        // Storing values
        if( ( alpha_XYZ[ 1 ] >= alphaMin )
         && ( index_ijk[ 0 ] - 1 >= tof_index_low[0] )
         && ( index_ijk[ 0 ] - 1 <= tof_index_high[0] )
         && ( index_ijk[ 1 ] - 1 >= tof_index_low[1] )
         && ( index_ijk[ 1 ] - 1 <= tof_index_high[1] )
         && ( index_ijk[ 2 ] - 1 >= tof_index_low[2] )
         && ( index_ijk[ 2 ] - 1 <= tof_index_high[2]  ) )
        {
          //// initialize if not initialized yet
          //if (!previous_edge_erf_initialized)
          //{
            //previous_edge_erf = erf( (length_LOR *alpha_c - lor_tof_center) / tof_sigma_sqrt2 );
            //previous_edge_erf_initialized = true;
          //}

          numVox = ( index_ijk[ 0 ] - 1 ) + ( ( index_ijk[ 1 ] - 1 ) ) * w + ( ( index_ijk[ 2 ] - 1 ) ) * wh;

          alphaMid = ( alpha_XYZ[ 1 ]  + alpha_c ) / 2.;

          // non TOF projection coefficient
          coeff = ( alpha_XYZ[ 1 ] - alpha_c ) * length_LOR;

          if (m_TOFWeightingFcnPrecomputedFlag)
          {
            // fetch the value of the precomputed TOF weighting function (centered at the TOF measurement)
            // nearest to the projection of the voxel center onto the LOR
            INTNB temp = m_TOFWeightingFcnNbSamples/2 + (INTNB)round( (alphaMid * length_LOR - lor_tof_center) * m_TOFPrecomputedSamplingFactor );
            // multiply the non TOF projection coefficient with the TOF weight
            if (temp>=0 && temp<m_TOFWeightingFcnNbSamples) coeff *= mp_TOFWeightingFcn[temp];
            else coeff=0.;
          }
          else
          {
            if (m_TOFBinProperProcessingFlag)
            {
              // integration between the edges of the current quantization TOF bin (centered at the TOF measurement)
              // of the Gaussian centered at the projection of the voxel center onto the LOR
              // normalized by the size of the bin so that the integral of the final TOF weighting function remains 1
              HPFLTNB temp_erf = std::min(tof_half_span, std::max(-tof_half_span, lor_tof_center - m_TOFBinSizeInMm/2. - alphaMid * length_LOR));
              HPFLTNB prev_erf = erf(temp_erf/tof_sigma_sqrt2);
              temp_erf = std::min(tof_half_span, std::max(-tof_half_span, lor_tof_center + m_TOFBinSizeInMm/2. - alphaMid * length_LOR));
              HPFLTNB new_erf = erf(temp_erf/tof_sigma_sqrt2);
              // multiply the non TOF projection coefficient with the TOF weight
              coeff *= 0.5 * fabs(new_erf - prev_erf) / m_TOFBinSizeInMm;
            }
            else
            {
              //// tof : erf of next voxel edge - Gaussian center
              //next_edge_erf = erf( (length_LOR * alpha_XYZ[ 1 ] - lor_tof_center) / tof_sigma_sqrt2 );
              //// integration of the Gaussian is done on the LOR portion matching the whole current voxel along X axis
              //coeff = 0.5 * std::fabs(previous_edge_erf - next_edge_erf) ;
              //// keeping record of the previous edge, so as to save 1 erf computation
              //previous_edge_erf = next_edge_erf;

              // TOF weight = value of the normalized Gaussian (centered at the TOF measurement)
              // at the projection of the voxel center onto the LOR
              HPFLTNB temp = (alphaMid * length_LOR - lor_tof_center) / tof_sigma;
              // multiply the non TOF projection coefficient with the TOF weight
              coeff *=  exp(- 0.5 * temp * temp ) * m_TOFGaussianNormCoef;
            }
          }

          // if the voxel is not masked 
          if (!m_hasMask || mp_mask[numVox]) ap_ProjectionLine->AddVoxel(a_direction, numVox, coeff);
        }

        // Increment values
        alpha_c = alpha_XYZ[ 1 ];
        alpha_XYZ[ 1 ] += alpha_u[ 1 ];
        index_ijk[ 1 ] += index_u[ 1 ];
      }
      else if( ( alpha_XYZ[ 2 ] < alpha_XYZ[ 0 ] )
            && ( alpha_XYZ[ 2 ] < alpha_XYZ[ 1 ] ) )
      {
        // Storing values
        if( ( alpha_XYZ[ 2 ] >= alphaMin )
         && ( index_ijk[ 0 ] - 1 >= tof_index_low[0] )
         && ( index_ijk[ 0 ] - 1 <= tof_index_high[0] )
         && ( index_ijk[ 1 ] - 1 >= tof_index_low[1] )
         && ( index_ijk[ 1 ] - 1 <= tof_index_high[1] )
         && ( index_ijk[ 2 ] - 1 >= tof_index_low[2] )
         && ( index_ijk[ 2 ] - 1 <= tof_index_high[2]  ) )
        {
          //// initialize if not initialized yet
          //if (!previous_edge_erf_initialized)
          //{
            //previous_edge_erf = erf( (length_LOR *alpha_c - lor_tof_center) / tof_sigma_sqrt2 );
            //previous_edge_erf_initialized = true;
          //}

          numVox = ( index_ijk[ 0 ] - 1 ) + ( ( index_ijk[ 1 ] - 1 ) ) * w + ( ( index_ijk[ 2 ] - 1 ) ) * wh;

          alphaMid = ( alpha_XYZ[ 2 ]  + alpha_c ) / 2.;

          // non TOF projection coefficient
          coeff = ( alpha_XYZ[ 2 ] - alpha_c ) * length_LOR;

          if (m_TOFWeightingFcnPrecomputedFlag)
          {
            // fetch the value of the precomputed TOF weighting function (centered at the TOF measurement)
            // nearest to the projection of the voxel center onto the LOR
            INTNB temp = m_TOFWeightingFcnNbSamples/2 + (INTNB)round( (alphaMid * length_LOR - lor_tof_center) * m_TOFPrecomputedSamplingFactor );
            // multiply the non TOF projection coefficient with the TOF weight
            if (temp>=0 && temp<m_TOFWeightingFcnNbSamples) coeff *= mp_TOFWeightingFcn[temp];
            else coeff = 0.;
          }
          else
          {
            if (m_TOFBinProperProcessingFlag)
            {
              // integration between the edges of the current quantization TOF bin (centered at the TOF measurement)
              // of the Gaussian centered at the projection of the voxel center onto the LOR
              // normalized by the size of the bin so that the integral of the final TOF weighting function remains 1
              HPFLTNB temp_erf = std::min(tof_half_span, std::max(-tof_half_span, lor_tof_center - m_TOFBinSizeInMm/2. - alphaMid * length_LOR));
              HPFLTNB prev_erf = erf(temp_erf/tof_sigma_sqrt2);
              temp_erf = std::min(tof_half_span, std::max(-tof_half_span, lor_tof_center + m_TOFBinSizeInMm/2. - alphaMid * length_LOR));
              HPFLTNB new_erf = erf(temp_erf/tof_sigma_sqrt2);
              // multiply the non TOF projection coefficient with the TOF weight
              coeff *= 0.5 * fabs(new_erf - prev_erf) / m_TOFBinSizeInMm;
            }
            else
            {
              //// tof : erf of next voxel edge - Gaussian center
              //next_edge_erf = erf( (length_LOR * alpha_XYZ[ 2 ] - lor_tof_center) / tof_sigma_sqrt2 );
              //// integration of the Gaussian is done on the LOR portion matching the whole current voxel along X axis
              //coeff = 0.5 * std::fabs(previous_edge_erf - next_edge_erf) ;
              //// keeping record of the previous edge, so as to save 1 erf computation
              //previous_edge_erf = next_edge_erf;

              // TOF weight = value of the normalized Gaussian (centered at the TOF measurement)
              // at the projection of the voxel center onto the LOR
              HPFLTNB temp = (alphaMid * length_LOR - lor_tof_center) / tof_sigma;
              // multiply the non TOF projection coefficient with the TOF weight
              coeff *=  exp(- 0.5 * temp * temp ) * m_TOFGaussianNormCoef;
            }
          }

          // if the voxel is not masked 
          if (!m_hasMask || mp_mask[numVox]) ap_ProjectionLine->AddVoxel(a_direction, numVox, coeff);
        }

        // Increment values
        alpha_c = alpha_XYZ[ 2 ];
        alpha_XYZ[ 2 ] += alpha_u[ 2 ];
        index_ijk[ 2 ] += index_u[ 2 ];
      }
    }
  }

  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iProjectorIncrementalSiddonMulti::ProjectTOFHistogram(int a_direction, oProjectionLine* ap_ProjectionLine)
{
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** iProjectorIncrementalSiddonMulti::ProjectTOFHistogram() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif

  #ifdef CASTOR_VERBOSE
  if (m_verbose>=10)
  {
    string direction = "";
    if (a_direction==FORWARD) direction = "forward";
    else direction = "backward";
    Cout("iProjectorIncrementalSiddonMulti::ProjectTOFHistogram() -> Project line '" << ap_ProjectionLine << "' in " << direction << " direction" << endl);
  }
  #endif

  // For list-mode, the number of lines to use depends on the processing stage (sensitivity generation or reconstruction)
  uint16_t nb_lines = (m_sensitivityMode == true) ?
                       m_nbSensLines : 
                       m_nbRecoLines ;

  // Loop on the number of lines
  for (uint16_t line=0; line<nb_lines; line++)
  {       
    // Get random positions from the scanner (mean depth of interaction intrinsicaly taken into account), taking POI into account if any
    if (mp_Scanner->GetRdmPositionsAndOrientations( ap_ProjectionLine->GetIndex1(), ap_ProjectionLine->GetIndex2(), 
                                                    ap_ProjectionLine->GetPosition1(), ap_ProjectionLine->GetPosition2(),
                                                    ap_ProjectionLine->GetOrientation1(), ap_ProjectionLine->GetOrientation2() ))
    {
      Cerr("***** iProjectorIncrementalSiddonMulti::ProjectTOFHistogram() -> A problem occurred while getting positions and orientations from scanner !" << endl);
      return 1;
    }                                       
    
    // We must called here ApplyOffset() and ComputeLineLength() which 
    // are supposed to be called in vProjector::Project()
    // This is a temporary solution as everyting related to multi-lines 
    // management should be performed in vProjector 
    
    // Apply global image offset
    ap_ProjectionLine->ApplyOffset();

    // Apply bed position offset
    ap_ProjectionLine->ApplyBedOffset();
  
    // Get end points position
    FLTNB* event1Float = ap_ProjectionLine->GetPosition1();
    FLTNB* event2Float = ap_ProjectionLine->GetPosition2();
    HPFLTNB event1[3] = { event1Float[0], event1Float[1], event1Float[2] };
    HPFLTNB event2[3] = { event2Float[0], event2Float[1], event2Float[2] };
  
    // **************************************
    // STEP 1: LOR length calculation
    // **************************************
    // Must recompute LOR length for each new line of response
    ap_ProjectionLine->ComputeLineLength();
    HPFLTNB length_LOR = ap_ProjectionLine->GetLength();
    HPFLTNB length_LOR_half = length_LOR * 0.5;
  
    // **************************************
    // STEP 2: Compute entrance voxel indexes
    // **************************************
    HPFLTNB alphaFirst[] = { 0.0, 0.0, 0.0 };
    HPFLTNB alphaLast[] = { 0.0, 0.0, 0.0 };
  
    HPFLTNB alphaMin = 0.0f, alphaMax = 1.0f;
    HPFLTNB delta_pos[] = {
      event2[ 0 ] - event1[ 0 ],
      event2[ 1 ] - event1[ 1 ],
      event2[ 2 ] - event1[ 2 ]
    };

    // Get TOF info
    INTNB tof_nb_bins = ap_ProjectionLine->GetNbTOFBins();
    INTNB tof_half_nb_bins = tof_nb_bins/2;

    // TOF Gaussian standard deviation and truncation
    HPFLTNB tof_sigma = m_TOFResolutionInMm / TWO_SQRT_TWO_LN_2;
    HPFLTNB tof_sigma_sqrt2 = sqrt(2.)*tof_sigma;
    HPFLTNB tof_half_span = 0.;
    if (m_TOFWeightingFcnPrecomputedFlag) tof_half_span = ((HPFLTNB)m_TOFWeightingFcnNbSamples)/(2.*m_TOFPrecomputedSamplingFactor);
    else tof_half_span = tof_sigma * m_TOFNbSigmas;

    // if integration of the Gaussian over the TOF bin, need for help variables to save calls to erf
    HPFLTNB prev_erf = 0., new_erf = 0.;

    // minimum and maximum TOF bins
    INTNB tof_bin_last = tof_half_nb_bins;
    INTNB tof_bin_first = -tof_half_nb_bins;

    // distance between the first event1 and the center of a TOF bin along the LOR
    HPFLTNB lor_tof_center = 0.;
    // the sum of all TOF bin weights for a voxel
    //HPFLTNB tof_norm_coef = 0.;

    // the first and the last relevant TOF bins for a voxel
    INTNB tof_bin_first_for_voxel = 0, tof_bin_last_for_voxel = 0;

    // Computation of alphaMin et alphaMax values (entrance and exit point of the ray)
    for( int i = 0; i < 3; ++i )
    {
      if( delta_pos[ i ] != 0.0 )
      {
        alphaFirst[i] = (-mp_halfFOV[i] - event1[i]) / (delta_pos[ i ]);
        alphaLast[i]  = ( mp_halfFOV[i] - event1[i]) / (delta_pos[ i ]);
        alphaMin = (std::max)(alphaMin,(std::min)(alphaFirst[i],alphaLast[i]));
        alphaMax = (std::min)(alphaMax,(std::max)(alphaFirst[i],alphaLast[i]));
      }
    }

    // if alphaMax is less than or equal to alphaMin no intersection
    // and return an empty buffer
    if( alphaMax <= alphaMin ) continue;

    // temporary storage for TOF bin weights
    // allocation after potential returns
    HPFLTNB* tof_weights_temp = new HPFLTNB[tof_nb_bins];
    for (INTNB tof_bin=0; tof_bin<tof_nb_bins; tof_bin++) tof_weights_temp[tof_bin] = 0.;

    // Now we have to find the indices of the particular plane
    // (iMin,iMax), (jMin,jMax), (kMin,kMax)
    int iMin = 0, iMax = 0;
    int jMin = 0, jMax = 0;
    int kMin = 0, kMax = 0;

    // For the X-axis
    if( delta_pos[ 0 ] > 0.0f )
    {
      iMin = ::ceil( ( mp_nbVox[ 0 ] + 1 ) - ( mp_halfFOV[ 0 ] - alphaMin * delta_pos[ 0 ] - event1[ 0 ] ) / mp_sizeVox[0] );
      iMax = ::floor( 1 + ( event1[ 0 ] + alphaMax * delta_pos[ 0 ] - (-mp_halfFOV[ 0 ]) ) / mp_sizeVox[0] );
    }
    else if( delta_pos[ 0 ] < 0.0f )
    {
      iMin = ::ceil( ( mp_nbVox[ 0 ] + 1 ) - ( mp_halfFOV[ 0 ] - alphaMax * delta_pos[ 0 ] - event1[ 0 ] ) / mp_sizeVox[0] );
      iMax = ::floor( 1 + ( event1[ 0 ] + alphaMin * delta_pos[ 0 ] - (-mp_halfFOV[ 0 ]) ) / mp_sizeVox[0] );
    }
    if( delta_pos[ 0 ] == 0 )
    {
      iMin = 1, iMax = 0;
    }

    // For the Y-axis
    if( delta_pos[ 1 ] > 0 )
    {
      jMin = ::ceil( ( mp_nbVox[ 1 ] + 1 ) - ( mp_halfFOV[ 1 ] - alphaMin * delta_pos[ 1 ] - event1[ 1 ] ) / mp_sizeVox[1] );
      jMax = ::floor( 1 + ( event1[ 1 ] + alphaMax * delta_pos[ 1 ] - (-mp_halfFOV[ 1 ]) ) / mp_sizeVox[1] );
    }
    else if( delta_pos[ 1 ] < 0 )
    {
      jMin = ::ceil( ( mp_nbVox[ 1 ] + 1 ) - ( mp_halfFOV[ 1 ] - alphaMax * delta_pos[ 1 ] - event1[ 1 ] ) / mp_sizeVox[1] );
      jMax = ::floor( 1 + ( event1[ 1 ] + alphaMin * delta_pos[ 1 ] - (-mp_halfFOV[ 1 ]) ) / mp_sizeVox[1] );
    }
    else if( delta_pos[ 1 ] == 0 )
    {
      jMin = 1, jMax = 0;
    }

    // For the Z-axis
    if( delta_pos[ 2 ] > 0 )
    {
      kMin = ::ceil( ( mp_nbVox[ 2 ] + 1 ) - ( mp_halfFOV[ 2 ] - alphaMin * delta_pos[ 2 ] - event1[ 2 ] ) / mp_sizeVox[2] );
      kMax = ::floor( 1 + ( event1[ 2 ] + alphaMax * delta_pos[ 2 ] - (-mp_halfFOV[ 2 ]) ) / mp_sizeVox[2] );
    }
    else if( delta_pos[ 2 ] < 0 )
    {
      kMin = ::ceil( ( mp_nbVox[ 2 ] + 1 ) - ( mp_halfFOV[ 2 ] - alphaMax * delta_pos[ 2 ] - event1[ 2 ] ) / mp_sizeVox[2] );
      kMax = ::floor( 1 + ( event1[ 2 ] + alphaMin * delta_pos[ 2 ] - (-mp_halfFOV[ 2 ]) ) / mp_sizeVox[2] );
    }
    else if( delta_pos[ 2 ] == 0 )
    {
      kMin = 1, kMax = 0;
    }

    // Computing the last term n number of intersection
    INTNB n = ( iMax - iMin + 1 ) + ( jMax - jMin + 1 )
      + ( kMax - kMin + 1 );

    // Array storing the first alpha in X, Y and Z
    HPFLTNB alpha_XYZ[ 3 ] = { 1.0f, 1.0f, 1.0f };

    // Computing the first alpha in X
    if( delta_pos[ 0 ] > 0 )
      alpha_XYZ[ 0 ] = ( ( (-mp_halfFOV[ 0 ]) + ( iMin - 1 ) * mp_sizeVox[0] ) - event1[ 0 ] ) / delta_pos[ 0 ];
    else if( delta_pos[ 0 ] < 0 )
      alpha_XYZ[ 0 ] = ( ( (-mp_halfFOV[ 0 ]) + ( iMax - 1 ) * mp_sizeVox[0] ) - event1[ 0 ] ) / delta_pos[ 0 ];

    // Computing the first alpha in Y
    if( delta_pos[ 1 ] > 0 )
      alpha_XYZ[ 1 ] = ( ( (-mp_halfFOV[ 1 ]) + ( jMin - 1 ) * mp_sizeVox[1] ) - event1[ 1 ] ) / delta_pos[ 1 ];
    else if( delta_pos[ 1 ] < 0 )
      alpha_XYZ[ 1 ] = ( ( (-mp_halfFOV[ 1 ]) + ( jMax - 1 ) * mp_sizeVox[1] ) - event1[ 1 ] ) / delta_pos[ 1 ];

    // Computing the first alpha in Z
    if( delta_pos[ 2 ] > 0 )
      alpha_XYZ[ 2 ] = ( ( (-mp_halfFOV[ 2 ]) + ( kMin - 1 ) * mp_sizeVox[2] ) - event1[ 2 ] ) / delta_pos[ 2 ];
    else if( delta_pos[ 2 ] < 0 )
      alpha_XYZ[ 2 ] = ( ( (-mp_halfFOV[ 2 ]) + ( kMax - 1 ) * mp_sizeVox[2] ) - event1[ 2 ] ) / delta_pos[ 2 ];

    // Computing the alpha updating
    HPFLTNB alpha_u[ 3 ] = {
      mp_sizeVox[0] / std::fabs( delta_pos[ 0 ] ),
      mp_sizeVox[1] / std::fabs( delta_pos[ 1 ] ),
      mp_sizeVox[2] / std::fabs( delta_pos[ 2 ] )
    };

    // Computing the index updating
    INTNB index_u[ 3 ] = {
      (delta_pos[ 0 ] < 0) ? -1 : 1,
      (delta_pos[ 1 ] < 0) ? -1 : 1,
      (delta_pos[ 2 ] < 0) ? -1 : 1
    };

    // Check which alpha is the min/max and increment
    if( alpha_XYZ[ 0 ] == alphaMin ) alpha_XYZ[ 0 ] += alpha_u[ 0 ];
    if( alpha_XYZ[ 1 ] == alphaMin ) alpha_XYZ[ 1 ] += alpha_u[ 1 ];
    if( alpha_XYZ[ 2 ] == alphaMin ) alpha_XYZ[ 2 ] += alpha_u[ 2 ];

    // Computing the minimum value in the alpha_XYZ buffer
    HPFLTNB const min_alpha_XYZ = (std::min)( alpha_XYZ[ 0 ], (std::min)( alpha_XYZ[ 1 ], alpha_XYZ[ 2 ] ) );

    // Computing the first index of intersection
    HPFLTNB alphaMid = ( min_alpha_XYZ + alphaMin ) / 2.0f;
    INTNB index_ijk[ 3 ] = {
      1 + (INTNB)( ( event1[ 0 ] + alphaMid * delta_pos[ 0 ] - (-mp_halfFOV[ 0 ]) ) / mp_sizeVox[0] ),
      1 + (INTNB)( ( event1[ 1 ] + alphaMid * delta_pos[ 1 ] - (-mp_halfFOV[ 1 ]) ) / mp_sizeVox[1] ),
      1 + (INTNB)( ( event1[ 2 ] + alphaMid * delta_pos[ 2 ] - (-mp_halfFOV[ 2 ]) ) / mp_sizeVox[2] )
    };

    INTNB const w = mp_nbVox[ 0 ];
    INTNB const wh = w * mp_nbVox[ 1 ];

    // Loop over the number of plane to cross
    HPFLTNB alpha_c = alphaMin;
    FLTNB coeff = 0.0f;
    INTNB numVox = 0;

    for( int nP = 0; nP < n - 1; ++nP )
    {
      if( ( alpha_XYZ[ 0 ] <= alpha_XYZ[ 1 ] )
       && ( alpha_XYZ[ 0 ] <= alpha_XYZ[ 2 ] ) )
      {
        // Storing values
        if( ( alpha_XYZ[ 0 ] >= alphaMin )
         && ( index_ijk[ 0 ] - 1 >= 0 )
         && ( index_ijk[ 0 ] - 1 <= mp_nbVox[ 0 ] - 1 )
         && ( index_ijk[ 1 ] - 1 >= 0 )
         && ( index_ijk[ 1 ] - 1 <= mp_nbVox[ 1 ] - 1 )
         && ( index_ijk[ 2 ] - 1 >= 0 )
         && ( index_ijk[ 2 ] - 1 <= mp_nbVox[ 2 ] - 1 ) )
        {
          numVox = ( index_ijk[ 0 ] - 1 ) + ( ( index_ijk[ 1 ] - 1 ) ) * w + ( ( index_ijk[ 2 ] - 1 ) ) * wh;
          
          // if this voxel is masked, do nothing
          if (!m_hasMask || mp_mask[numVox])
          {
            coeff = ( alpha_XYZ[ 0 ] - alpha_c ) * length_LOR;

            alphaMid = ( alpha_XYZ[ 0 ] + alpha_c ) * 0.5;

            // compute the first and the last relevant TOF bin for this voxel
            if (!m_TOFWeightingFcnPrecomputedFlag && m_TOFBinProperProcessingFlag)
            {
              // taking the TOF bins reached by the truncated Gaussian centered at the voxel center projection
              tof_bin_first_for_voxel = max(tof_bin_first , (INTNB)(floor( (alphaMid * length_LOR - tof_half_span - length_LOR_half - m_TOFBinSizeInMm/2.) / m_TOFBinSizeInMm )));
              tof_bin_last_for_voxel = min(tof_bin_last , (INTNB)(ceil( (alphaMid * length_LOR + tof_half_span - length_LOR_half + m_TOFBinSizeInMm/2.) / m_TOFBinSizeInMm )));
            }
            else
            {
              // taking the TOF bins whose TOF weighting function, centered at the bin center, reaches the voxel center projection 
              tof_bin_first_for_voxel = max(tof_bin_first , (INTNB)(ceil( (alphaMid * length_LOR - tof_half_span - length_LOR_half ) / m_TOFBinSizeInMm )));
              tof_bin_last_for_voxel = min(tof_bin_last , (INTNB)(floor( (alphaMid * length_LOR + tof_half_span - length_LOR_half) / m_TOFBinSizeInMm )));
            }

            lor_tof_center = length_LOR_half + tof_bin_first_for_voxel * m_TOFBinSizeInMm;

            // shift tof bin indices from -/+ to 0:nbBins range
            tof_bin_first_for_voxel += tof_half_nb_bins;
            tof_bin_last_for_voxel += tof_half_nb_bins;

            // initialization of help variables for reducing calls to erf
            if (!m_TOFWeightingFcnPrecomputedFlag && m_TOFBinProperProcessingFlag)
            {
              // bound the integration to the Gaussian truncation
              HPFLTNB temp = std::min(tof_half_span,std::max(-tof_half_span, lor_tof_center - m_TOFBinSizeInMm/2. - alphaMid * length_LOR));
              prev_erf = erf(temp/tof_sigma_sqrt2);
            }

            // compute TOF bin weights for the current voxel for all relevant TOF bins
            //tof_norm_coef = 0.;
            for (INTNB tof_bin=tof_bin_first_for_voxel; tof_bin<=tof_bin_last_for_voxel; tof_bin++)
            {
              if (m_TOFWeightingFcnPrecomputedFlag)
              {
                // fetch the value of the precomputed TOF weighting function (centered at the TOF bin center)
                // nearest to the projection of the voxel center onto the LOR
                INTNB temp = m_TOFWeightingFcnNbSamples/2 + (INTNB)round( ( alphaMid * length_LOR - lor_tof_center) * m_TOFPrecomputedSamplingFactor);
                if (temp>=0 && temp<m_TOFWeightingFcnNbSamples) tof_weights_temp[tof_bin] = mp_TOFWeightingFcn[temp];
                // add the weight to the sum for normalization
                //tof_norm_coef += tof_weights_temp[tof_bin];
                // update TOF bin center along the LOR for the next TOF bin
                lor_tof_center += m_TOFBinSizeInMm;
              }
              else
              {
                if (m_TOFBinProperProcessingFlag)
                {
                  // integration between the edges of the current TOF bin of the Gaussian centered at the projection of the voxel center onto the LOR
                  // reuse of integration from the previous bin to save calls to erf
                  // bound the integration to the Gaussian truncation
                  HPFLTNB temp = std::min(tof_half_span, std::max(-tof_half_span, lor_tof_center + m_TOFBinSizeInMm/2. - alphaMid * length_LOR));
                  new_erf = erf(temp/tof_sigma_sqrt2);
                  tof_weights_temp[tof_bin] = 0.5 * fabs(new_erf - prev_erf);
                  // add the weight to the sum
                 // tof_norm_coef += tof_weights_temp[tof_bin];
                  prev_erf = new_erf;
                  // update TOF bin center along the LOR for the next TOF bin
                  lor_tof_center += m_TOFBinSizeInMm;
                }
                else
                {
                  // TOF weight = TOF bin size * value of the normalized Gaussian (centered at the TOF bin center) at the projection of the voxel center onto the LOR
                  HPFLTNB temp = ( alphaMid * length_LOR - lor_tof_center) / tof_sigma;
                  // save the weight temporarily
                  tof_weights_temp[tof_bin] = exp(- 0.5 * temp * temp ) * m_TOFGaussianNormCoef * m_TOFBinSizeInMm;
                  // add the weight to the sum
                  //tof_norm_coef += tof_weights_temp[tof_bin];
                  // update TOF center along the LOR for the next TOF bin
                  lor_tof_center += m_TOFBinSizeInMm;
                }
              }
            }
  /*
            if (tof_norm_coef>0.)
            {
              // first normalize TOF bin weights so that they sum to 1
              for (INTNB tof_bin=tof_bin_first_for_voxel; tof_bin<=tof_bin_last_for_voxel; tof_bin++) tof_weights_temp[tof_bin] /= tof_norm_coef;
            }
  */
            // write the final TOF bin weighted projection coefficients for the current voxel
            ap_ProjectionLine->AddVoxelAllTOFBins(a_direction, numVox, coeff, tof_weights_temp, tof_bin_first_for_voxel, tof_bin_last_for_voxel);

          }
        }
        // Increment values
        alpha_c = alpha_XYZ[ 0 ];
        alpha_XYZ[ 0 ] += alpha_u[ 0 ];
        index_ijk[ 0 ] += index_u[ 0 ];
      }
      else if( ( alpha_XYZ[ 1 ] < alpha_XYZ[ 0 ] )
            && ( alpha_XYZ[ 1 ] <= alpha_XYZ[ 2 ] ) )
      {
        // Storing values
        if( ( alpha_XYZ[ 1 ] >= alphaMin )
         && ( index_ijk[ 0 ] - 1 >= 0 )
         && ( index_ijk[ 0 ] - 1 <= mp_nbVox[ 0 ] - 1 )
         && ( index_ijk[ 1 ] - 1 >= 0 )
         && ( index_ijk[ 1 ] - 1 <= mp_nbVox[ 1 ] - 1 )
         && ( index_ijk[ 2 ] - 1 >= 0 )
         && ( index_ijk[ 2 ] - 1 <= mp_nbVox[ 2 ] - 1 ) )
        {
          numVox = ( index_ijk[ 0 ] - 1 ) + ( ( index_ijk[ 1 ] - 1 ) ) * w + ( ( index_ijk[ 2 ] - 1 ) ) * wh;

          // if this voxel is masked, do nothing
          if (!m_hasMask || mp_mask[numVox])
          {
            coeff = ( alpha_XYZ[ 1 ] - alpha_c ) * length_LOR;

            alphaMid = ( alpha_XYZ[ 1 ] + alpha_c ) * 0.5;

            // compute the first and the last relevant TOF bin for this voxel
            if (!m_TOFWeightingFcnPrecomputedFlag && m_TOFBinProperProcessingFlag)
            {
              // taking the TOF bins reached by the truncated Gaussian centered at the voxel center projection
              tof_bin_first_for_voxel = max(tof_bin_first , (INTNB)(floor( (alphaMid * length_LOR - tof_half_span - length_LOR_half - m_TOFBinSizeInMm/2.) / m_TOFBinSizeInMm )));
              tof_bin_last_for_voxel = min(tof_bin_last , (INTNB)(ceil( (alphaMid * length_LOR + tof_half_span - length_LOR_half + m_TOFBinSizeInMm/2.) / m_TOFBinSizeInMm )));
            }
            else
            {
              // taking the TOF bins whose TOF weighting function, centered at the bin center, reaches the voxel center projection 
              tof_bin_first_for_voxel = max(tof_bin_first , (INTNB)(ceil( (alphaMid * length_LOR - tof_half_span - length_LOR_half ) / m_TOFBinSizeInMm )));
              tof_bin_last_for_voxel = min(tof_bin_last , (INTNB)(floor( (alphaMid * length_LOR + tof_half_span - length_LOR_half) / m_TOFBinSizeInMm )));
            }
            
            lor_tof_center = length_LOR_half + tof_bin_first_for_voxel * m_TOFBinSizeInMm;

            // shift tof bin indices from -/+ to 0:nbBins range
            tof_bin_first_for_voxel += tof_half_nb_bins;
            tof_bin_last_for_voxel += tof_half_nb_bins;

            // initialization of help variables for reducing calls to erf
            if (!m_TOFWeightingFcnPrecomputedFlag && m_TOFBinProperProcessingFlag)
            {
              // bound the integration to the Gaussian truncation
              HPFLTNB temp = std::min(tof_half_span,std::max(-tof_half_span, lor_tof_center - m_TOFBinSizeInMm/2. - alphaMid * length_LOR));
              prev_erf = erf(temp/tof_sigma_sqrt2);
            }

            // compute TOF bin weights for the current voxel for all relevant TOF bins
            //tof_norm_coef = 0.;
            for (INTNB tof_bin=tof_bin_first_for_voxel; tof_bin<=tof_bin_last_for_voxel; tof_bin++)
            {
              if (m_TOFWeightingFcnPrecomputedFlag)
              {
                // fetch the value of the precomputed TOF weighting function (centered at the TOF bin center)
                // nearest to the projection of the voxel center onto the LOR
                INTNB temp = m_TOFWeightingFcnNbSamples/2 + (INTNB)round( ( alphaMid * length_LOR - lor_tof_center) * m_TOFPrecomputedSamplingFactor );
                if (temp>=0 && temp<m_TOFWeightingFcnNbSamples) tof_weights_temp[tof_bin] = mp_TOFWeightingFcn[temp];
                // add the weight to the sum
                //tof_norm_coef += tof_weights_temp[tof_bin];
                // update TOF center along the LOR for the next TOF bin
                lor_tof_center += m_TOFBinSizeInMm;
              }
              else
              {
                if (m_TOFBinProperProcessingFlag)
                {
                  // integration between the edges of the current TOF bin of the Gaussian centered at the projection of the voxel center onto the LOR
                  // reuse of integration from the previous bin to save calls to erf
                  // bound the integration to the Gaussian truncation
                  HPFLTNB temp = std::min(tof_half_span, std::max(-tof_half_span, lor_tof_center + m_TOFBinSizeInMm/2. - alphaMid * length_LOR));
                  new_erf = erf(temp/tof_sigma_sqrt2);
                  tof_weights_temp[tof_bin] = 0.5 * fabs(new_erf - prev_erf);
                  // add the weight to the sum for normalization
                  //tof_norm_coef += tof_weights_temp[tof_bin];
                  prev_erf = new_erf;
                  // update TOF center along the LOR for the next TOF bin
                  lor_tof_center += m_TOFBinSizeInMm;
                }
                else
                {
                  // TOF weight = TOF bin size * value of the normalized Gaussian (centered at the TOF bin center) at the projection of the voxel center onto the LOR
                  HPFLTNB temp = ( alphaMid * length_LOR - lor_tof_center) / tof_sigma;
                  // save the weight temporarily
                  tof_weights_temp[tof_bin] = exp(- 0.5 * temp * temp ) * m_TOFGaussianNormCoef * m_TOFBinSizeInMm;
                  // add the weight to the sum for normalization
                  //tof_norm_coef += tof_weights_temp[tof_bin];
                  // update TOF center along the LOR for the next TOF bin
                  lor_tof_center += m_TOFBinSizeInMm;
                } 
              }
            }
  /*
            if (tof_norm_coef>0.)
            {
              // first normalize TOF bin weights so that they sum to 1
              for (INTNB tof_bin=tof_bin_first_for_voxel; tof_bin<=tof_bin_last_for_voxel; tof_bin++) tof_weights_temp[tof_bin] /= tof_norm_coef;
            }
  */
            // write the final TOF bin weighted projection coefficients for the current voxel
            ap_ProjectionLine->AddVoxelAllTOFBins(a_direction, numVox, coeff, tof_weights_temp, tof_bin_first_for_voxel, tof_bin_last_for_voxel);
          }
        }

        // Increment values
        alpha_c = alpha_XYZ[ 1 ];
        alpha_XYZ[ 1 ] += alpha_u[ 1 ];
        index_ijk[ 1 ] += index_u[ 1 ];
      }
      else if( ( alpha_XYZ[ 2 ] < alpha_XYZ[ 0 ] )
            && ( alpha_XYZ[ 2 ] < alpha_XYZ[ 1 ] ) )
      {
        // Storing values
        if( ( alpha_XYZ[ 2 ] >= alphaMin )
         && ( index_ijk[ 0 ] - 1 >= 0 )
         && ( index_ijk[ 0 ] - 1 <= mp_nbVox[ 0 ] - 1 )
         && ( index_ijk[ 1 ] - 1 >= 0 )
         && ( index_ijk[ 1 ] - 1 <= mp_nbVox[ 1 ] - 1 )
         && ( index_ijk[ 2 ] - 1 >= 0 )
         && ( index_ijk[ 2 ] - 1 <= mp_nbVox[ 2 ] - 1 ) )
        {
          numVox = ( index_ijk[ 0 ] - 1 ) + ( ( index_ijk[ 1 ] - 1 ) ) * w + ( ( index_ijk[ 2 ] - 1 ) ) * wh;

          // if this voxel is masked, do nothing
          if (!m_hasMask || mp_mask[numVox])
          {
            coeff = ( alpha_XYZ[ 2 ] - alpha_c ) * length_LOR;

            alphaMid = ( alpha_XYZ[ 2 ] + alpha_c ) * 0.5;

            // compute the first and the last relevant TOF bin for this voxel
            if (!m_TOFWeightingFcnPrecomputedFlag && m_TOFBinProperProcessingFlag)
            {
              // taking the TOF bins reached by the truncated Gaussian centered at the voxel center projection
              tof_bin_first_for_voxel = max(tof_bin_first , (INTNB)(floor( (alphaMid * length_LOR - tof_half_span - length_LOR_half - m_TOFBinSizeInMm/2.) / m_TOFBinSizeInMm )));
              tof_bin_last_for_voxel = min(tof_bin_last , (INTNB)(ceil( (alphaMid * length_LOR + tof_half_span - length_LOR_half + m_TOFBinSizeInMm/2.) / m_TOFBinSizeInMm )));
            }
            else
            {
              // taking the TOF bins whose TOF weighting function, centered at the bin center, reaches the voxel center projection 
              tof_bin_first_for_voxel = max(tof_bin_first , (INTNB)(ceil( (alphaMid * length_LOR - tof_half_span - length_LOR_half ) / m_TOFBinSizeInMm )));
              tof_bin_last_for_voxel = min(tof_bin_last , (INTNB)(floor( (alphaMid * length_LOR + tof_half_span - length_LOR_half) / m_TOFBinSizeInMm )));
            }

            lor_tof_center = length_LOR_half + tof_bin_first_for_voxel * m_TOFBinSizeInMm;

            // shift tof bin indices from -/+ to 0:nbBins range
            tof_bin_first_for_voxel += tof_half_nb_bins;
            tof_bin_last_for_voxel += tof_half_nb_bins;


            // initialization of help variables for reducing calls to erf
            if (!m_TOFWeightingFcnPrecomputedFlag && m_TOFBinProperProcessingFlag)
            {
              // bound the integration to the Gaussian truncation
              HPFLTNB temp = std::min(tof_half_span,std::max(-tof_half_span, lor_tof_center - m_TOFBinSizeInMm/2. - alphaMid * length_LOR));
              prev_erf = erf(temp/tof_sigma_sqrt2);
            }

            // compute TOF bin weights for the current voxel for all relevant TOF bins
            //tof_norm_coef = 0.;
            for (INTNB tof_bin=tof_bin_first_for_voxel; tof_bin<=tof_bin_last_for_voxel; tof_bin++)
            {
              if (m_TOFWeightingFcnPrecomputedFlag)
              {
                // fetch the value of the precomputed TOF weighting function (centered at the TOF bin center)
                // nearest to the projection of the voxel center onto the LOR
                INTNB temp = m_TOFWeightingFcnNbSamples/2 + (INTNB)round( ( alphaMid * length_LOR - lor_tof_center) * m_TOFPrecomputedSamplingFactor );
                if (temp>=0 && temp<m_TOFWeightingFcnNbSamples) tof_weights_temp[tof_bin] = mp_TOFWeightingFcn[temp];
                // add the weight to the sum
                //tof_norm_coef += tof_weights_temp[tof_bin];
                // update TOF bin center along the LOR for the next TOF bin
                lor_tof_center += m_TOFBinSizeInMm;
              }
              else
              {
                if (m_TOFBinProperProcessingFlag)
                {
                  // integration between the edges of the current TOF bin of the Gaussian centered at the projection of the voxel center onto the LOR
                  // reuse of integration from the previous bin to save calls to erf
                  // bound the integration to the Gaussian truncation
                  HPFLTNB temp = std::min(tof_half_span, std::max(-tof_half_span, lor_tof_center + m_TOFBinSizeInMm/2. - alphaMid * length_LOR));
                  new_erf = erf(temp/tof_sigma_sqrt2);
                  tof_weights_temp[tof_bin] = 0.5 * fabs(new_erf - prev_erf);
                  // add the weight to the sum
                  //tof_norm_coef += tof_weights_temp[tof_bin];
                  prev_erf = new_erf;
                  // update TOF center along the LOR for the next TOF bin
                  lor_tof_center += m_TOFBinSizeInMm;
                }
                else
                {
                  // TOF weight = TOF bin size * value of the normalized Gaussian (centered at the TOF bin center) at the projection of the voxel center onto the LOR
                  HPFLTNB temp = ( alphaMid * length_LOR - lor_tof_center) / tof_sigma;
                  tof_weights_temp[tof_bin] = exp(- 0.5 * temp * temp ) * m_TOFGaussianNormCoef * m_TOFBinSizeInMm;
                  // add the weight to the sum
                  //tof_norm_coef += tof_weights_temp[tof_bin];
                  // update TOF center along the LOR for the next TOF bin
                  lor_tof_center += m_TOFBinSizeInMm;
                }
              }
            }
  /*
            if (tof_norm_coef>0.)
            {
              // first normalize TOF bin weights so that they sum to 1
              for (INTNB tof_bin=tof_bin_first_for_voxel; tof_bin<=tof_bin_last_for_voxel; tof_bin++) tof_weights_temp[tof_bin] /= tof_norm_coef;
            }
  */
            // write the final TOF bin weighted projection coefficients for the current voxel
            ap_ProjectionLine->AddVoxelAllTOFBins(a_direction, numVox, coeff, tof_weights_temp, tof_bin_first_for_voxel, tof_bin_last_for_voxel);
          }
        }

        // Increment values
        alpha_c = alpha_XYZ[ 2 ];
        alpha_XYZ[ 2 ] += alpha_u[ 2 ];
        index_ijk[ 2 ] += index_u[ 2 ];
      }
    }
    delete[] tof_weights_temp;
  }

  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

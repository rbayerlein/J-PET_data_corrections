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
  \brief    Implementation of class iPenaltyMedianRootPrior
*/

#include "iPenaltyMedianRootPrior.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iPenaltyMedianRootPrior::iPenaltyMedianRootPrior() : vPenalty()
{
  // --------------------------
  // Specific member parameters
  // --------------------------
  m_neighborhoodMaxNbVoxels = -1;
  m_neighborhoodShape = MRF_NOT_DEFINED;
  m_neighborhoodSphereRadius = -1.;
  m_neighborhoodBoxOrder = -1;
  m_neighborhoodBoxExcludeCorners = 0; // Keep them by default
  mp_neighborhoodNbVoxels = NULL;
  m2p_neighborhoodIndices = NULL;
  m2p_neighborhoodKernel = NULL;
  mp_medianValue = NULL;
  // Derivative order is infinite
  m_penaltyDerivativesOrder = INT_MAX;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iPenaltyMedianRootPrior::~iPenaltyMedianRootPrior()
{
  // Destroy mp_neighborhoodNbVoxels
  if (mp_neighborhoodNbVoxels)
  {
    free(mp_neighborhoodNbVoxels);
  }
  // Destroy m2p_neighborhoodIndices
  if (m2p_neighborhoodIndices)
  {
    for (int th=0; th<mp_ImageDimensionsAndQuantification->GetNbThreadsForImageComputation(); th++)
    {
      if (m2p_neighborhoodIndices[th]) free(m2p_neighborhoodIndices[th]);
    }
    free(m2p_neighborhoodIndices);
  }
  // Destroy m2p_neighborhoodKernel
  if (m2p_neighborhoodKernel)
  {
    for (INTNB v=0; v<m_neighborhoodMaxNbVoxels; v++)
    {
      if (m2p_neighborhoodKernel[v]) free(m2p_neighborhoodKernel[v]);
    }
    free(m2p_neighborhoodKernel);
  }
  // Destroy mp_medianValue
  if (mp_medianValue) free(mp_medianValue);
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iPenaltyMedianRootPrior::ShowHelpSpecific()
{
  cout << "The Median Root Prior (MRP) is implemented for several types of neighborhood. The parameters are described below, and a" << endl;
  cout << "template configuration file for parameters selection can be found in the configuration folder (config/optimizer/MRP.conf)." << endl;
  cout << "The configuration of this penalty can only be done through a configuration file. It is not possible to parameterize it" << endl;
  cout << "directly from the command line. The current implementation is the simplest one, based on the following reference:" << endl;
  cout << "S. Alenius and U Ruotsalainen, \"Bayesian image reconstruction for emission tomography based on the median root prior\"," << endl;
  cout << "European Journal of Nuclear Medicine, vol. 24, pp. 258-265, 1997." << endl;
  cout << "------------" << endl;
  cout << "Neighborhood" << endl;
  cout << "------------" << endl;
  cout << "The neighborhood is set by the 'neighborhood shape' keyword and can be described using one of the following setting:" << endl;
  cout << "[6-nearest]: Consider the 6 nearest neighbors, i.e. 2 along each dimension. In 2D, only 4 in the plane are used." << endl;
  cout << "[box]:       Consider all voxels included in a box centered on the voxel of interest. The size of the box is parameterized" << endl;
  cout << "             using the 'box neighborhood order' keyword. The side of the box is equal to '2*order+1'. The 8 corner voxels" << endl;
  cout << "             can also be excluded from the neighborhood by setting the 'exclude box neighborhood corners' keyword to 1." << endl;
  cout << "[sphere]:    Consider all voxels whose center is included in a sphere centered on the voxel of interest and of radius" << endl;
  cout << "             provided by the 'sphere neighborhood radius (mm)' keyword." << endl;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iPenaltyMedianRootPrior::ReadConfigurationFile(const string& a_configurationFile)
{
  string buffer = "";
  string key_word = "";

  // -------------------------------------------------------------------
  // Read the type of neighborhood
  // -------------------------------------------------------------------
  key_word = "neighborhood shape";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &buffer, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iPenaltyMedianRootPrior::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Sphere neighborhood
  if (buffer == "sphere")
  {
    m_neighborhoodShape = MRF_NEIGHBORHOOD_SPHERE;
    // Radius of the sphere
    key_word = "sphere neighborhood radius (mm)";
    if (ReadDataASCIIFile(a_configurationFile, key_word, &m_neighborhoodSphereRadius, 1, KEYWORD_MANDATORY))
    {
      Cerr("***** iPenaltyMedianRootPrior::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
      return 1;
    }
  }
  // Box neighborhood
  else if (buffer == "box")
  {
    m_neighborhoodShape = MRF_NEIGHBORHOOD_BOX;
    // Order of the box (number of voxels which will be included in each direction)
    key_word = "box neighborhood order";
    if (ReadDataASCIIFile(a_configurationFile, key_word, &m_neighborhoodBoxOrder, 1, KEYWORD_MANDATORY))
    {
      Cerr("***** iPenaltyMedianRootPrior::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
      return 1;
    }
    // Include corners of the box or not (some people do that) (not mandatory as we keep them by default)
    key_word = "exclude box neighborhood corners";
    if (ReadDataASCIIFile(a_configurationFile, key_word, &m_neighborhoodBoxExcludeCorners, 1, KEYWORD_OPTIONAL) == KEYWORD_OPTIONAL_ERROR)
    {
      Cerr("***** iPenaltyMedianRootPrior::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
      return 1;
    }
  }
  // 6-nearest
  else if (buffer == "6-nearest")
  {
    m_neighborhoodShape = MRF_NEIGHBORHOOD_6_NEAREST;
  }
  // Unknown
  else
  {
    Cerr("***** iPenaltyMedianRootPrior::ReadConfigurationFile() -> Unknown neighborhood type '" << buffer << "' !" << endl);
    return 1;
  }

  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
 
int iPenaltyMedianRootPrior::ReadOptionsList(const string& a_optionsList)
{ 
  // To complicated to do it that way
  Cerr("***** iPenaltyMedianRootPrior::ReadOptionsList() -> Options can be specified only using a configuration file !" << endl);
  return 1;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iPenaltyMedianRootPrior::CheckSpecificParameters()
{
  // Check neighborhood type
  if (m_neighborhoodShape==MRF_NOT_DEFINED)
  {
    Cerr("***** iPenaltyMarkovRandomField::CheckSpecificParameters() -> Should provide a neighborhood type !" << endl);
    return 1;
  }
  // Check some parameters of the neighborhood
  if (m_neighborhoodShape==MRF_NEIGHBORHOOD_SPHERE && m_neighborhoodSphereRadius<0.)
  {
    Cerr("***** iPenaltyMarkovRandomField::CheckSpecificParameters() -> Provided radius of the sphere neighborhood (" << m_neighborhoodSphereRadius << ") is negative" << endl);
    return 1;
  }
  else if (m_neighborhoodShape==MRF_NEIGHBORHOOD_BOX && m_neighborhoodBoxOrder<0)
  {
    Cerr("***** iPenaltyMarkovRandomField::CheckSpecificParameters() -> Provided order of the box neighborhood (" << m_neighborhoodBoxOrder << ") is negative" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iPenaltyMedianRootPrior::InitializeSpecific()
{
  // Build the neighborhood kernel that will be used to facilitate the computation of the neighborhood of each voxel later on
  if (BuildNeighborhoodKernel())
  {
    Cerr("***** iPenaltyMedianRootPrior::InitializeSpecific() -> Failed to build the neighborhood !" << endl);
    return 1;  
  }
  // Some allocations for the specific neighborhood of each voxel used during computations
  mp_neighborhoodNbVoxels = (INTNB*)malloc(mp_ImageDimensionsAndQuantification->GetNbThreadsForImageComputation()*sizeof(INTNB));
  m2p_neighborhoodIndices = (INTNB**)malloc(mp_ImageDimensionsAndQuantification->GetNbThreadsForImageComputation()*sizeof(INTNB*));
  for (int th=0; th<mp_ImageDimensionsAndQuantification->GetNbThreadsForImageComputation(); th++)
  {
    m2p_neighborhoodIndices[th]=(INTNB*)malloc((m_neighborhoodMaxNbVoxels)*sizeof(INTNB));
  }
  // Allocate the table to store median values for each thread
  mp_medianValue = (FLTNB*)malloc(mp_ImageDimensionsAndQuantification->GetNbThreadsForImageComputation()*sizeof(FLTNB));
  // Verbose
  if (m_verbose>=VERBOSE_NORMAL)
  {
    Cout("iPenaltyMedianRootPrior::InitializeSpecific() -> Use the Median Root Prior" << endl);
    if (m_verbose>=VERBOSE_DETAIL)
    {
      // Type of neighborhood
      if (m_neighborhoodShape == MRF_NEIGHBORHOOD_SPHERE)
      {
        Cout("  --> shape of neighborhood: sphere" << endl);
        Cout("  --> radius of the sphere neighborhood: " << m_neighborhoodSphereRadius << " mm " << endl);
      }
      else if (m_neighborhoodShape == MRF_NEIGHBORHOOD_BOX)
      {
        Cout("  --> shape of neighborhood: box" << endl);
        Cout("  --> order of the box (number of voxels in each direction): " << m_neighborhoodBoxOrder << endl);
        if (m_neighborhoodBoxExcludeCorners) Cout("  --> exclude corner voxels" << endl);
      }
      else if (m_neighborhoodShape == MRF_NEIGHBORHOOD_6_NEAREST)
      {
        Cout("  --> shape of neighborhood: 6-nearest" << endl);
      }
      Cout("  --> Max number of voxels in the neighborhood: " << m_neighborhoodMaxNbVoxels << endl);
      Cout("  --> penalty energy function derivatives order: " << m_penaltyDerivativesOrder << endl);
    }
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iPenaltyMedianRootPrior::BuildNeighborhoodKernel()
{
  // Check if the neighborhood has already been created
  if (m2p_neighborhoodKernel)
  {
    Cerr("***** iPenaltyMedianRootPrior::BuildNeighborhoodKernel() -> Neighborhood kernel has already been created !" << endl);
    return 1;
  }
  // Get voxel sizes locally
  FLTNB vox_size_x = mp_ImageDimensionsAndQuantification->GetVoxSizeX();
  FLTNB vox_size_y = mp_ImageDimensionsAndQuantification->GetVoxSizeY();
  FLTNB vox_size_z = mp_ImageDimensionsAndQuantification->GetVoxSizeZ();
  // ========================================================================================================
  // Note
  // ========================================================================================================
  // As opposed to the markov random fields in which we discard the voxel itself from its neighborhood,
  // when computing the median in MRP, the voxel itself is included in the neighborhood.
  // ========================================================================================================
  // Case where the neighborhood is defined by a sphere with radius provided in mm
  // ========================================================================================================
  if (m_neighborhoodShape == MRF_NEIGHBORHOOD_SPHERE)
  {
    // Compute the size of the box including the sphere (minimum of 1 voxel radius)
    INTNB radius_vox_x = max( (INTNB)1 , ((INTNB)(m_neighborhoodSphereRadius/vox_size_x)) );
    INTNB radius_vox_y = max( (INTNB)1 , ((INTNB)(m_neighborhoodSphereRadius/vox_size_y)) );
    INTNB radius_vox_z = max( (INTNB)1 , ((INTNB)(m_neighborhoodSphereRadius/vox_size_z)) );
    // Particular case for 2D reconstruction
    if (mp_ImageDimensionsAndQuantification->GetNbVoxZ()==1) radius_vox_z = 0;
    // Precompute the maximum number of voxels in the neighborhood as the box including the sphere
    m_neighborhoodMaxNbVoxels = (radius_vox_x*2+1) * (radius_vox_y*2+1) * (radius_vox_z*2+1);
    // Allocate the neighborhood kernel
    m2p_neighborhoodKernel = (INTNB**)malloc(m_neighborhoodMaxNbVoxels*sizeof(INTNB*));
    for (INTNB v = 0; v<m_neighborhoodMaxNbVoxels; v++) m2p_neighborhoodKernel[v] = (INTNB*)malloc(3*sizeof(INTNB));
    // Counter for the number of voxels in the neighborhood
    INTNB nb_voxels_in_neighborhood = 0;
    // Loop on increments of voxels indices in box neighborhood
    for (INTNB vox_z = -radius_vox_z; vox_z<=radius_vox_z; vox_z++)
      for (INTNB vox_y = -radius_vox_y; vox_y<=radius_vox_y; vox_y++)
        for (INTNB vox_x = -radius_vox_x; vox_x<=radius_vox_x; vox_x++)
        {
          // Compute the square distance between the two voxels
          FLTNB squareDistance = ((FLTNB)(vox_x*vox_x))*vox_size_x*vox_size_x
                               + ((FLTNB)(vox_y*vox_y))*vox_size_y*vox_size_y
                               + ((FLTNB)(vox_z*vox_z))*vox_size_z*vox_size_z;
          // Include voxel if the distance from the central voxel is less than the radius
          if (squareDistance<m_neighborhoodSphereRadius*m_neighborhoodSphereRadius)
          {
            m2p_neighborhoodKernel[nb_voxels_in_neighborhood][MRF_NEIGHBOR_X] = vox_x; 
            m2p_neighborhoodKernel[nb_voxels_in_neighborhood][MRF_NEIGHBOR_Y] = vox_y;
            m2p_neighborhoodKernel[nb_voxels_in_neighborhood][MRF_NEIGHBOR_Z] = vox_z;
            nb_voxels_in_neighborhood++;
          }
        }
    // Update the maximum number of voxels in the neighborhood by the exact value
    m_neighborhoodMaxNbVoxels = nb_voxels_in_neighborhood; 
    m2p_neighborhoodKernel = (INTNB**)realloc(m2p_neighborhoodKernel,m_neighborhoodMaxNbVoxels*sizeof(INTNB*));
  }
  // ========================================================================================================
  // Case where the neighborhood is defined by a box as described in many papers
  // ========================================================================================================
  else if (m_neighborhoodShape == MRF_NEIGHBORHOOD_BOX)
  {
    // Particular case for 2D reconstruction
    INTNB neighborhood_box_order_z = m_neighborhoodBoxOrder;
    if (mp_ImageDimensionsAndQuantification->GetNbVoxZ()==1) neighborhood_box_order_z = 0;
    // Maximum number of voxels in the neighborhood of a voxel with no image boundary constraints
    m_neighborhoodMaxNbVoxels = (m_neighborhoodBoxOrder*2+1)*(m_neighborhoodBoxOrder*2+1)*(m_neighborhoodBoxOrder*2+1);
    // Allocate the neighborhood kernel
    m2p_neighborhoodKernel = (INTNB**)malloc(m_neighborhoodMaxNbVoxels*sizeof(INTNB*));
    for (INTNB v = 0; v<m_neighborhoodMaxNbVoxels; v++) m2p_neighborhoodKernel[v] = (INTNB*)malloc(3*sizeof(INTNB));
    // Counter for the number of voxels in the neighborhood
    INTNB nb_voxels_in_neighborhood = 0;
    // Loop on increments of voxels indices in the box neighborhood
    for (INTNB vox_z = -neighborhood_box_order_z; vox_z<=neighborhood_box_order_z; vox_z++)
      for (INTNB vox_y = -m_neighborhoodBoxOrder; vox_y<=m_neighborhoodBoxOrder; vox_y++)
        for (INTNB vox_x = -m_neighborhoodBoxOrder; vox_x<=m_neighborhoodBoxOrder; vox_x++)
        {
          // Exclude corner voxels if requested
          if (!m_neighborhoodBoxExcludeCorners || (vox_x!=-m_neighborhoodBoxOrder && vox_x!=m_neighborhoodBoxOrder)
                                               || (vox_y!=-m_neighborhoodBoxOrder && vox_y!=m_neighborhoodBoxOrder)
                                               || (vox_z!=-m_neighborhoodBoxOrder && vox_z!=m_neighborhoodBoxOrder))
          {
            m2p_neighborhoodKernel[nb_voxels_in_neighborhood][MRF_NEIGHBOR_X] = vox_x;
            m2p_neighborhoodKernel[nb_voxels_in_neighborhood][MRF_NEIGHBOR_Y] = vox_y;
            m2p_neighborhoodKernel[nb_voxels_in_neighborhood][MRF_NEIGHBOR_Z] = vox_z;
            nb_voxels_in_neighborhood++;
          }
        }
    // Update the maximum number of voxels in the neighborhood by the exact value
    m_neighborhoodMaxNbVoxels = nb_voxels_in_neighborhood; 
    m2p_neighborhoodKernel = (INTNB**)realloc(m2p_neighborhoodKernel,m_neighborhoodMaxNbVoxels*sizeof(INTNB*));
  }
  // ========================================================================================================
  // Case where the neighborhood is defined simply by the 6-nearest neighbors
  // ========================================================================================================
  else if (m_neighborhoodShape == MRF_NEIGHBORHOOD_6_NEAREST)
  {
    // Maximum number of voxels in the neighborhood is 6 + 1 (the voxel itself)
    m_neighborhoodMaxNbVoxels = 7;
    // Particular case for 2D reconstruction
    if (mp_ImageDimensionsAndQuantification->GetNbVoxZ()==1) m_neighborhoodMaxNbVoxels = 5;
    // Allocate the neighborhood kernel
    m2p_neighborhoodKernel = (INTNB**)malloc(m_neighborhoodMaxNbVoxels*sizeof(INTNB*));
    for (INTNB v = 0; v<m_neighborhoodMaxNbVoxels; v++) m2p_neighborhoodKernel[v] = (INTNB*)malloc(3*sizeof(INTNB));
    // The 7 voxels
    for (INTNB v=0; v<m_neighborhoodMaxNbVoxels; v++) for (INTNB dim=0; dim<3; dim++) m2p_neighborhoodKernel[v][dim] = 0;
    m2p_neighborhoodKernel[1][MRF_NEIGHBOR_X] = -1;
    m2p_neighborhoodKernel[2][MRF_NEIGHBOR_X] = 1;
    m2p_neighborhoodKernel[3][MRF_NEIGHBOR_Y] = -1;
    m2p_neighborhoodKernel[4][MRF_NEIGHBOR_Y] = 1;
    if (mp_ImageDimensionsAndQuantification->GetNbVoxZ()!=1)
    {
      m2p_neighborhoodKernel[5][MRF_NEIGHBOR_Z] = -1;
      m2p_neighborhoodKernel[6][MRF_NEIGHBOR_Z] = 1;
    }
  }
  // ========================================================================================================
  // Unknown neighborhood type
  // ========================================================================================================
  else
  {
    Cerr("***** iPenaltyMedianRootPrior::BuildNeighborhoodKernel() -> The provided shape of neighborhood (" << m_neighborhoodShape << ") is unknown !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iPenaltyMedianRootPrior::BuildSpecificNeighborhood(INTNB a_voxel, int a_th)
{
  // Get number of voxels locally for code readability
  INTNB nb_vox_x  = mp_ImageDimensionsAndQuantification->GetNbVoxX();
  INTNB nb_vox_y  = mp_ImageDimensionsAndQuantification->GetNbVoxY();
  INTNB nb_vox_z  = mp_ImageDimensionsAndQuantification->GetNbVoxZ();
  INTNB nb_vox_xy = mp_ImageDimensionsAndQuantification->GetNbVoxXY();
  // Compute the X, Y and Z components of the current voxel
  INTNB index_x = a_voxel % nb_vox_x;
  INTNB index_y = (a_voxel/nb_vox_x) % nb_vox_y;
  INTNB index_z = a_voxel / nb_vox_xy;
  // Count the number of valid neighbors
  mp_neighborhoodNbVoxels[a_th] = 0;
  // Loop on all possible neighbors in the kernel
  for (INTNB neigh = 0; neigh<m_neighborhoodMaxNbVoxels; neigh++)
  {
    // Compute X, Y and Z indices of this potential neighbor
    INTNB neighbor_x = index_x + m2p_neighborhoodKernel[neigh][MRF_NEIGHBOR_X];
    INTNB neighbor_y = index_y + m2p_neighborhoodKernel[neigh][MRF_NEIGHBOR_Y];
    INTNB neighbor_z = index_z + m2p_neighborhoodKernel[neigh][MRF_NEIGHBOR_Z];
    // Check if the neighbour is inside the image
    if ( neighbor_x>=0 && neighbor_x<nb_vox_x && neighbor_y>=0 && neighbor_y<nb_vox_y && neighbor_z>=0 && neighbor_z<nb_vox_z )
    {
      // Add the global index of this neighbor to the neighborhood list
      m2p_neighborhoodIndices[a_th][neigh] = neighbor_x + neighbor_y*nb_vox_x + neighbor_z*nb_vox_xy;
      // One more neighbor
      mp_neighborhoodNbVoxels[a_th]++;
    }
    // Otherwise, discard it by setting -1
    else
    {
      m2p_neighborhoodIndices[a_th][neigh] = -1;
    }
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iPenaltyMedianRootPrior::LocalPreProcessingStep(int a_tbf, int a_rbf, int a_cbf, INTNB a_voxel, int a_th)
{
  // First build the list of neighbours of this voxel
  if (BuildSpecificNeighborhood(a_voxel,a_th))
  {
    Cerr("***** iPenaltyMedianRootPrior::LocalPreProcessingStep() -> A problem occurred while building specific neighborhood of voxel " << a_voxel << " !" << endl);
    return 1;
  }
  // Compute the median value in this neighborhood
  if (ComputeMedianInNeighborhood(a_tbf,a_rbf,a_cbf,a_voxel,a_th))
  {
    Cerr("***** iPenaltyMedianRootPrior::LocalPreProcessingStep() -> A problem occurred while computing the median in neighborhood of voxel " << a_voxel << " !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iPenaltyMedianRootPrior::ComputeMedianInNeighborhood(int a_tbf, int a_rbf, int a_cbf, INTNB a_voxel, int a_th)
{
  // Sort the neighborhood indices according to increasing intensity values
  std::sort( &(m2p_neighborhoodIndices[a_th][0]), &(m2p_neighborhoodIndices[a_th][m_neighborhoodMaxNbVoxels]),
  [a_voxel,a_tbf,a_rbf,a_cbf,a_th,this](INTNB i1, INTNB i2) 
  {
    // This two lines pushes the -1 indices (voxels out of image) at the end
    if (i1 == -1) return false;
    if (i2 == -1) return true;
    return mp_ImageSpace->m4p_image[a_tbf][a_rbf][a_cbf][i1] < mp_ImageSpace->m4p_image[a_tbf][a_rbf][a_cbf][i2];
  });
  // There can be an even number of voxels if we are near the border of the image (restrained neighborhood)
  if (mp_neighborhoodNbVoxels[a_th]%2==0)
  {
    // In this case, we take the mean of the tow median voxels
    mp_medianValue[a_th] = (mp_ImageSpace->m4p_image[a_tbf][a_rbf][a_cbf][m2p_neighborhoodIndices[a_th][mp_neighborhoodNbVoxels[a_th]/2]]
                         +  mp_ImageSpace->m4p_image[a_tbf][a_rbf][a_cbf][m2p_neighborhoodIndices[a_th][mp_neighborhoodNbVoxels[a_th]/2-1]])
                         * 0.5;
  }
  // Even number of voxels
  else
  {
    // The median is in the middle
    mp_medianValue[a_th] = mp_ImageSpace->m4p_image[a_tbf][a_rbf][a_cbf][m2p_neighborhoodIndices[a_th][mp_neighborhoodNbVoxels[a_th]/2]];
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

FLTNB iPenaltyMedianRootPrior::ComputePenaltyValue(int a_tbf, int a_rbf, int a_cbf, INTNB a_voxel, int a_th)
{
  // Compute the value of the penalty
  FLTNB penalty = 0.;
  if (mp_medianValue[a_th]!=0.) penalty = m_penaltyStrength * 0.5
                                        * (mp_ImageSpace->m4p_image[a_tbf][a_rbf][a_cbf][a_voxel]-mp_medianValue[a_th])
                                        * (mp_ImageSpace->m4p_image[a_tbf][a_rbf][a_cbf][a_voxel]-mp_medianValue[a_th])
                                        / mp_medianValue[a_th];
  // Return result
  return penalty;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

FLTNB iPenaltyMedianRootPrior::ComputeFirstDerivative(int a_tbf, int a_rbf, int a_cbf, INTNB a_voxel, int a_th)
{
  // Compute the first derivative
  FLTNB first_derivative = 0.;
  if (mp_medianValue[a_th]!=0.) first_derivative = m_penaltyStrength * (mp_ImageSpace->m4p_image[a_tbf][a_rbf][a_cbf][a_voxel]-mp_medianValue[a_th])
                                                                     / mp_medianValue[a_th];
  return first_derivative;
}


// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

FLTNB iPenaltyMedianRootPrior::ComputeSecondDerivative(int a_tbf, int a_rbf, int a_cbf, INTNB a_voxel, int a_th)
{
  // Compute the second derivative of the penalty
  FLTNB second_derivative = 0.;
  if (mp_medianValue[a_th]!=0.) second_derivative = m_penaltyStrength / mp_medianValue[a_th];
  // Return result
  return second_derivative;
}

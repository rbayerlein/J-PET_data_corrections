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
  \brief    Implementation of class iPenaltyMarkovRandomField
*/

#include "iPenaltyMarkovRandomField.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iPenaltyMarkovRandomField::iPenaltyMarkovRandomField() : vPenalty()
{
  // --------------------------
  // Specific member parameters
  // --------------------------
  m_penaltyID = "MRF";
  m_potentialType = MRF_NOT_DEFINED;
  m_potentialRelativeDifferenceGamma = -1.;
  m_potentialGemanMcClureDelta = -1.;
  m_potentialGreenLogCoshDelta = -1.;
  m_potentialHebertLeahyMu = -1.;
  m_potentialHuberDelta = -1.;
  m_proximityType = MRF_NOT_DEFINED;
  m_proximityCharacteristicDistance = -1.;
  mp_proximityKernel = NULL;
  m_similarityType = MRF_NOT_DEFINED;
  m_similarityBowsherThreshold = -1;
  m_similarityBowsherNbVoxels = -1;
  m2p_similarityFactors = NULL;
  m_neighborhoodMaxNbVoxels = -1;
  m_neighborhoodShape = MRF_NOT_DEFINED;
  m_neighborhoodSphereRadius = -1.;
  m_neighborhoodBoxOrder = -1;
  m_neighborhoodBoxExcludeCorners = 0; // Keep them by default
  mp_neighborhoodNbVoxels = NULL;
  m2p_neighborhoodIndices = NULL;
  m2p_neighborhoodKernel = NULL;
  // The derivative order depends on the potential function that we do not know yet, so we let it by default to 0
  m_penaltyDerivativesOrder = 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iPenaltyMarkovRandomField::~iPenaltyMarkovRandomField()
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
  // Destroy m2p_similarityFactors
  if (m2p_similarityFactors)
  {
    for (int th=0; th<mp_ImageDimensionsAndQuantification->GetNbThreadsForImageComputation(); th++)
    {
      if (m2p_similarityFactors[th]) free(m2p_similarityFactors[th]);
    }
    free(m2p_similarityFactors);
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
  // Destroy mp_proximityKernel
  if (mp_proximityKernel)
  {
    free(mp_proximityKernel);
  }
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iPenaltyMarkovRandomField::ShowHelpSpecific()
{
  cout << "The Markov Random Field (MRF) penalty is implemented for several types of neighborhood, potential functions," << endl;
  cout << "similarity and proximity factors. All these parameters are described below, and a template configuration file" << endl;
  cout << "for parameters selection can be found in the configuration folder (config/optimizer/MRF.conf). The configuration" << endl;
  cout << "of this penalty can only be done through a configuration file. It is not possible to parameterize it directly" << endl;
  cout << "from the command line. The MRF penalty has the following expression:" << endl;
  cout << "penalty = beta * sum_on_voxels * sum_on_neighborhood * proximity_factor * similarity_factor * potential_function" << endl;
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
  cout << "-----------------" << endl;
  cout << "Proximity factors" << endl;
  cout << "-----------------" << endl;
  cout << "The proximity factors are used to weight the contribution of each neighbor to the global penalty, based on proximity to the" << endl;
  cout << "voxel of interest. These factors are always normalized so that their sum is 1. They can be set using the 'proximity factor'" << endl;
  cout << "keyword based on one of the following setting:" << endl;
  cout << "[none]:      Consider uniform proximity factors." << endl;
  cout << "[voxel]:     Consider factors inversely proportional to the distance in voxels from the voxel of interest (i.e. the unit" << endl;
  cout << "             is voxels, whatever the axis and the voxel dimensions)." << endl;
  cout << "[euclidian]: Consider factors inversely proportional to the euclidian distance from the voxel of interest (i.e. the unit" << endl;
  cout << "             is mm)." << endl;
  cout << endl;
  cout << "------------------" << endl;
  cout << "Similarity factors" << endl;
  cout << "------------------" << endl;
  cout << "The similarity factors are used to weight the contribution of each neighbor to the global penalty, based on similarity to the" << endl;
  cout << "voxel of interest. These factors are not normalized by default. They can be set using the 'similarity factor' keyword based on" << endl;
  cout << "one of the following setting:" << endl;
  cout << "[none]:      Consider uniform similarity factors, i.e. the value of 1 for each voxel in the neighborhood." << endl;
  cout << "[aBowsher]:  Consider similarity factors based on the asymmetrical Bowsher's method and an additional image. The additional" << endl;
  cout << "             image must be provided using the '-multimodal' option in the command line. Based on additional image values," << endl;
  cout << "             voxels of the neighborhood most similar to the voxel of interest will have a similarity factor of 1, and the" << endl;
  cout << "             other voxels will have a similarity factor of 0. The number of most similar voxels is parameterized by a percentage" << endl;
  cout << "             of the number of voxels in the neighborhood, provided with the 'similarity threshold Bowsher (%)' keyword." << endl;
  cout << "             For an explanation of asymmetrical Bowsher, see Schramm et al, IEEE Trans. Med. Imaging, vol. 37, pp. 590, 2018." << endl;
  cout << "-------------------" << endl;
  cout << "Potential functions" << endl;
  cout << "-------------------" << endl;
  cout << "The potential function actually penalizes the difference between the voxel of interest and a neighbor. It can be set using the" << endl;
  cout << "'potential function' keyword based on one of the following setting:" << endl;
  cout << "[quadratic]:       p(u,v) = 0.5*(u-v)^2" << endl;
  cout << "                   Reference: Geman and Geman, IEEE Trans. Pattern Anal. Machine Intell., vol. PAMI-6, pp. 721-741, 1984." << endl;
  cout << "[geman mcclure]:   p(u,v,d) = (u-v)^2 / (d^2+(u-v)^2)" << endl;
  cout << "                   The parameter 'd' can be set using the 'deltaGMC' keyword." << endl;
  cout << "                   Reference: Geman and McClure, Proc. Amer. Statist. Assoc., 1985." << endl;
  cout << "[hebert leahy]:    p(u,v,m) = log(1+(u-v)^2/m^2)" << endl;
  cout << "                   The parameter 'm' can be set using the 'muHL' keyword." << endl;
  cout << "                   Reference: Hebert and Leahy, IEEE Trans. Med. Imaging, vol. 8, pp. 194-202, 1989." << endl;
  cout << "[green logcosh]:   p(u,v,d) = log(cosh((u-v)/d))" << endl;
  cout << "                   The parameter 'd' can be set using the 'deltaLogCosh' keyword." << endl;
  cout << "                   Reference: Green, IEEE Trans. Med. Imaging, vol. 9, pp. 84-93, 1990." << endl;
  cout << "[huber piecewise]: " << endl;
  cout << "                   p(u,v,d) = d*abs(u-v)-0.5*d^2  if  abs(u-v) >  d" << endl;
  cout << "                            = 0.5*(u-v)^2         if  abs(u-v) <= d" << endl;
  cout << "                   The parameter 'd' can be set using the 'deltaHuber' keyword." << endl;
  cout << "                   Reference: e.g. Mumcuoglu et al, Phys. Med. Biol., vol. 41, pp. 1777-1807, 1996." << endl;
  cout << "[nuyts relative]:  p(u,v,g) = (u-v)^2 / (u+v+g*abs(u-v))" << endl;
  cout << "                   The parameter 'g' can be set using the 'gammaRD' keyword." << endl;
  cout << "                   Reference: Nuyts et al, IEEE Trans. Nucl. Sci., vol. 49, pp. 56-60, 2002." << endl;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iPenaltyMarkovRandomField::ReadConfigurationFile(const string& a_configurationFile)
{
  string buffer = "";
  string key_word = "";

  // -------------------------------------------------------------------
  // Read the type of neighborhood
  // -------------------------------------------------------------------
  key_word = "neighborhood shape";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &buffer, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iPenaltyMarkovRandomField::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
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
      Cerr("***** iPenaltyMarkovRandomField::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
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
      Cerr("***** iPenaltyMarkovRandomField::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
      return 1;
    }
    // Include corners of the box or not (some people do that) (not mandatory as we keep them by default)
    key_word = "exclude box neighborhood corners";
    if (ReadDataASCIIFile(a_configurationFile, key_word, &m_neighborhoodBoxExcludeCorners, 1, KEYWORD_OPTIONAL) == KEYWORD_OPTIONAL_ERROR)
    {
      Cerr("***** iPenaltyMarkovRandomField::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
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
    Cerr("***** iPenaltyMarkovRandomField::ReadConfigurationFile() -> Unknown neighborhood type '" << buffer << "' !" << endl);
    return 1;
  }

  // -------------------------------------------------------------------
  // Read and check the type of the proximity factor
  // -------------------------------------------------------------------
  key_word = "proximity factor";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &buffer, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iPenaltyMarkovRandomField::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Inverse of euclidian distance
  if (buffer == "euclidian")
  {
    m_proximityType = MRF_PROXIMITY_EUCLIDIAN;
  }
  // Inverse of voxel distance (in numbers of voxels)
  else if (buffer == "voxel")
  {
    m_proximityType = MRF_PROXIMITY_VOXEL;
  }
  // Exponential distance squared
  /*
  else if (buffer == "exponential distance squared")
  {
    m_proximityType = MRF_PROXIMITY_EXP_DIST_SQUARED;
    key_word = "proximity characteristic distance (mm)";
    if (ReadDataASCIIFile(a_configurationFile, key_word, &m_proximityCharacteristicDistance, 1, KEYWORD_MANDATORY))
    {
      Cerr("***** iPenaltyMarkovRandomField::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
      return 1;
    }
  }
  */
  // No proximity factor
  else if (buffer == "none")
  {
    m_proximityType = MRF_PROXIMITY_NONE;
  }
  // Unknown
  else
  {
    Cerr("***** iPenaltyMarkovRandomField::ReadConfigurationFile() -> Unknown proximity factor type '" << buffer << "' !" << endl);
    return 1;
  }

  // -------------------------------------------------------------------
  // Read and check the type of the similarity factor
  // -------------------------------------------------------------------
  key_word = "similarity factor";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &buffer, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iPenaltyMarkovRandomField::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // No similarity factors
  if (buffer == "none")
  {
    m_similarityType = MRF_SIMILARITY_NONE;
  }
  // Bowsher's similarity factors
  else if (buffer == "aBowsher")
  {
    m_similarityType = MRF_SIMILARITY_BOWSHER;
    // Read the Bowsher threshold
    key_word = "similarity threshold Bowsher (%)";
    if (ReadDataASCIIFile(a_configurationFile, key_word, &m_similarityBowsherThreshold, 1, KEYWORD_MANDATORY))
    {
      Cerr("***** iPenaltyMarkovRandomField::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
      return 1;
    }
  }
  // Unknown
  else
  {
    Cerr("***** iPenaltyMarkovRandomField::ReadConfigurationFile() -> Unknown similarity factor type '" << buffer << "' !" << endl);
    return 1;
  }

  // -------------------------------------------------------------------
  // Read and check the type of the potential energy function
  // -------------------------------------------------------------------
  key_word = "potential function";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &buffer, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iPenaltyMarkovRandomField::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Quadratic potential function
  if (buffer == "quadratic")
  {
    m_potentialType = MRF_POTENTIAL_QUADRATIC;
    m_penaltyDerivativesOrder = INT_MAX;
  }
  // Relative differences potential function
  else if (buffer == "nuyts relative")
  {
    m_potentialType = MRF_POTENTIAL_RELATIVE_DIFFERENCE;
    m_penaltyDerivativesOrder = INT_MAX;
    key_word = "gammaRD";
    if (ReadDataASCIIFile(a_configurationFile, key_word, &m_potentialRelativeDifferenceGamma, 1, KEYWORD_MANDATORY))
    {
      Cerr("***** iPenaltyMarkovRandomField::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
      return 1;
    }
  }
  // Geman and McClure function
  else if (buffer == "geman mcclure")
  {
    m_potentialType = MRF_POTENTIAL_GEMAN_MCCLURE;
    m_penaltyDerivativesOrder = INT_MAX;
    key_word = "deltaGMC";
    if (ReadDataASCIIFile(a_configurationFile, key_word, &m_potentialGemanMcClureDelta, 1, KEYWORD_MANDATORY))
    {
      Cerr("***** iPenaltyMarkovRandomField::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
      return 1;
    }
  }
  // Green's logcosh
  else if (buffer == "green logcosh")
  {
    m_potentialType = MRF_POTENTIAL_GREEN;
    m_penaltyDerivativesOrder = INT_MAX;
    key_word = "deltaLogCosh";
    if (ReadDataASCIIFile(a_configurationFile, key_word, &m_potentialGreenLogCoshDelta, 1, KEYWORD_MANDATORY))
    {
      Cerr("***** iPenaltyMarkovRandomField::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
      return 1;
    }
  }
  // Hebert and Leahy function
  else if (buffer == "hebert leahy")
  {
    m_potentialType = MRF_POTENTIAL_HEBERT_LEAHY;
    m_penaltyDerivativesOrder = INT_MAX;
    key_word = "muHL";
    if (ReadDataASCIIFile(a_configurationFile, key_word, &m_potentialHebertLeahyMu, 1, KEYWORD_MANDATORY))
    {
      Cerr("***** iPenaltyMarkovRandomField::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
      return 1;
    }
  }
  // Huber function
  else if (buffer == "huber piecewise")
  {
    m_potentialType = MRF_POTENTIAL_HUBER;
    m_penaltyDerivativesOrder = INT_MAX;
    key_word = "deltaHuber";
    if (ReadDataASCIIFile(a_configurationFile, key_word, &m_potentialHuberDelta, 1, KEYWORD_MANDATORY))
    {
      Cerr("***** iPenaltyMarkovRandomField::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
      return 1;
    }
  }
  // Unknown
  else
  {
    Cerr("***** iPenaltyMarkovRandomField::ReadConfigurationFile() -> Unknown potential function '" << buffer << "' !" << endl);
    return 1;
  }

  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
 
int iPenaltyMarkovRandomField::ReadOptionsList(const string& a_optionsList)
{ 
  // Too complicated to do it that way
  Cerr("***** iPenaltyMarkovRandomField::ReadOptionsList() -> Options can be specified only using a configuration file !" << endl);
  return 1;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iPenaltyMarkovRandomField::CheckSpecificParameters()
{
  // Check potential function type
  if (m_potentialType==MRF_NOT_DEFINED)
  {
    Cerr("***** iPenaltyMarkovRandomField::CheckSpecificParameters() -> Should provide a potential function type !" << endl);
    return 1;
  }
  // Check proximity factor type
  if (m_proximityType==MRF_NOT_DEFINED)
  {
    Cerr("***** iPenaltyMarkovRandomField::CheckSpecificParameters() -> Should provide a proximity factor type !" << endl);
    return 1;
  }
  // Check similarity factor type
  if (m_similarityType==MRF_NOT_DEFINED)
  {
    Cerr("***** iPenaltyMarkovRandomField::CheckSpecificParameters() -> Should provide a similarity factor type !" << endl);
    return 1;
  }
  // Check neighborhood type
  if (m_neighborhoodShape==MRF_NOT_DEFINED)
  {
    Cerr("***** iPenaltyMarkovRandomField::CheckSpecificParameters() -> Should provide a neighborhood type !" << endl);
    return 1;
  }
  /*
  // check parameters of the proximity factors
  if (m_proximityType == MRF_PROXIMITY_EXP_DIST_SQUARED && m_proximityCharacteristicDistance <=0.)
  {
    Cerr("***** iPenaltyMarkovRandomField::CheckSpecificParameters() -> Provided proximity characteristic distance (" << m_proximityCharacteristicDistance << ") is non-positive" << endl);
    return 1;    
  }
  */
  // Check parameters of the similarity factors
  if (m_similarityType == MRF_SIMILARITY_BOWSHER)
  {
    if (m_similarityBowsherThreshold<0. || m_similarityBowsherThreshold>100.)
    {
      Cerr("***** iPenaltyMarkovRandomField::CheckSpecificParameters() -> Provided threshold parameter for the Bowsher similarity " << m_similarityBowsherThreshold << " does not fall into [0,100]% !" << endl);
      return 1;
    }
    else if (mp_ImageDimensionsAndQuantification->GetNbMultiModalImages()==0)
    {
      Cerr("***** iPenaltyMarkovRandomField::CheckSpecificParameters() -> Bowsher similarity requires a multimodal image !"<< endl);
      return 1;
    }
    else if (mp_ImageDimensionsAndQuantification->GetNbMultiModalImages()>1)
    {
      Cout("***** iPenaltyMarkovRandomField::CheckSpecificParameters() -> Warning : More than one multimodal image, Bowsher similarity will use only the first one !"<< endl);
    }
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
  // Check for gamma value of the relative distance potential function
  if (m_potentialType==MRF_POTENTIAL_RELATIVE_DIFFERENCE && m_potentialRelativeDifferenceGamma<0.)
  {
    Cerr("***** iPenaltyMarkovRandomField::CheckSpecificParameters() -> Provided gamma parameter of relative differences potential function must be positive or null !" << endl);
    return 1;
  }
  // Check for delta value of the Geman and McClure potential function
  if (m_potentialType==MRF_POTENTIAL_GEMAN_MCCLURE && m_potentialGemanMcClureDelta<=0.)
  {
    Cerr("***** iPenaltyMarkovRandomField::CheckSpecificParameters() -> Provided delta parameter of Geman and McClure's potential function must be strictly positive !" << endl);
    return 1;
  }
  // Check for delta value of the Green's log-cosh potential function
  if (m_potentialType==MRF_POTENTIAL_GREEN && m_potentialGreenLogCoshDelta<=0.)
  {
    Cerr("***** iPenaltyMarkovRandomField::CheckSpecificParameters() -> Provided delta parameter of Green's log-cosh potential function must be strictly positive !" << endl);
    return 1;
  }
  // Check for mu value of the Hebert and Leahy potential function
  if (m_potentialType==MRF_POTENTIAL_HEBERT_LEAHY && m_potentialHebertLeahyMu<=0.)
  {
    Cerr("***** iPenaltyMarkovRandomField::CheckSpecificParameters() -> Provided mu parameter of Hebert and Leahy's potential function must be strictly positive !" << endl);
    return 1;
  }
  // Check for delta value of the Huber potential function
  if (m_potentialType==MRF_POTENTIAL_HUBER && m_potentialHuberDelta<=0.)
  {
    Cerr("***** iPenaltyMarkovRandomField::CheckSpecificParameters() -> Provided delta parameter of Huber's piecewise potential function must be strictly positive !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iPenaltyMarkovRandomField::InitializeSpecific()
{
  // Build the neighborhood kernel that will be used to facilitate the computation of the neighborhood of each voxel later on
  if (BuildNeighborhoodKernel())
  {
    Cerr("***** iPenaltyMarkovRandomField::InitializeSpecific() -> Failed to build the neighborhood !" << endl);
    return 1;  
  }
  // Build the proximity factors of the neighborhood that remain constant for any voxel
  if (BuildProximityFactors())
  {
    Cerr("***** iPenaltyMarkovRandomField::InitializeSpecific() -> Failed to build the proximity factors !" << endl);
    return 1;  
  }
  // Some allocations for the specific neighborhood of each voxel used during computations
  mp_neighborhoodNbVoxels = (INTNB*)malloc(mp_ImageDimensionsAndQuantification->GetNbThreadsForImageComputation()*sizeof(INTNB));
  m2p_neighborhoodIndices = (INTNB**)malloc(mp_ImageDimensionsAndQuantification->GetNbThreadsForImageComputation()*sizeof(INTNB*));
  m2p_similarityFactors = (FLTNB**)malloc(mp_ImageDimensionsAndQuantification->GetNbThreadsForImageComputation()*sizeof(FLTNB*));
  for (int th=0; th<mp_ImageDimensionsAndQuantification->GetNbThreadsForImageComputation(); th++)
  {
    m2p_neighborhoodIndices[th]=(INTNB*)malloc((m_neighborhoodMaxNbVoxels)*sizeof(INTNB));
    m2p_similarityFactors[th]=(FLTNB*)malloc((m_neighborhoodMaxNbVoxels)*sizeof(FLTNB));
    // Initialize the similarity factors to 1, in case they stay uniform for the reconstruction
    for (INTNB v=0; v<m_neighborhoodMaxNbVoxels; v++) m2p_similarityFactors[th][v] = 1.;
  }
  // Verbose
  if (m_verbose>=VERBOSE_NORMAL)
  {
    Cout("iPenaltyMarkovRandomField::InitializeSpecific() -> Use the MRF penalty" << endl);
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
    }
    Cout("  --> Max number of voxels in the neighborhood: " << m_neighborhoodMaxNbVoxels << endl);
    if (m_verbose>=VERBOSE_DETAIL)
    {
      // Proximity factors
      if (m_proximityType == MRF_PROXIMITY_NONE)
      {
        Cout("  --> proximity factor: uniform" << endl);
      }
      else if (m_proximityType == MRF_PROXIMITY_EUCLIDIAN)
      {
        Cout("  --> proximity factor: inverse of euclidian distance" << endl);
      }
      else if (m_proximityType == MRF_PROXIMITY_VOXEL)
      {
        Cout("  --> proximity factor: inverse of voxel distance (i.e. distance in numbers of voxels)" << endl);
      }
/*
      else if (m_proximityType == MRF_PROXIMITY_EXP_DIST_SQUARED)
      {
        Cout("  --> proximity factor: exponential square distances"<< endl);
        Cout("  --> characteristic distance for the proximity factors: "<< m_proximityCharacteristicDistance << endl);
      }
*/
    }
    // Similarity factors
    if (m_similarityType == MRF_SIMILARITY_NONE )
    {
      Cout("  --> similarity factor: uniform" << endl);
    }    
    else if (m_similarityType == MRF_SIMILARITY_BOWSHER)
    {
      Cout("  --> similarity factor: asymmetric Bowsher" << endl);
      if (m_verbose>=VERBOSE_DETAIL)
      {
        Cout("  --> Bowsher's threshold: " << m_similarityBowsherThreshold << " %" << endl);
        Cout("  --> Bowsher's most similar voxels kept: " << m_similarityBowsherNbVoxels << endl);
      }
    }
    // Potential function
    if (m_potentialType == MRF_POTENTIAL_QUADRATIC)
    {
      Cout("  --> potential function: quadratic" << endl);
    }
    else if (m_potentialType == MRF_POTENTIAL_RELATIVE_DIFFERENCE)
    {
      Cout("  --> potential function: Nuyts relative differences 2002" << endl);
      if (m_verbose>=VERBOSE_DETAIL) Cout("  --> gamma: " << m_potentialRelativeDifferenceGamma << endl);
    }
    else if (m_potentialType == MRF_POTENTIAL_GEMAN_MCCLURE)
    {
      Cout("  --> potential function: Geman and McClure 1985" << endl);
      if (m_verbose>=VERBOSE_DETAIL) Cout("  --> delta: " << m_potentialGemanMcClureDelta << endl);
    }
    else if (m_potentialType == MRF_POTENTIAL_GREEN)
    {
      Cout("  --> potential function: Green log-cosh 1990" << endl);
      if (m_verbose>=VERBOSE_DETAIL) Cout("  --> delta: " << m_potentialGreenLogCoshDelta << endl);
    }
    else if (m_potentialType == MRF_POTENTIAL_HEBERT_LEAHY)
    {
      Cout("  --> potential function: Hebert and Leahy 1989" << endl);
      if (m_verbose>=VERBOSE_DETAIL) Cout("  --> mu: " << m_potentialHebertLeahyMu << endl);
    }
    else if (m_potentialType == MRF_POTENTIAL_HUBER)
    {
      Cout("  --> potential function: Huber piecewise linear-quadratic" << endl);
      if (m_verbose>=VERBOSE_DETAIL) Cout("  --> delta: " << m_potentialHuberDelta << endl);
    }
    if (m_verbose>=VERBOSE_DETAIL) Cout("  --> penalty energy function derivatives order: " << m_penaltyDerivativesOrder << endl);
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iPenaltyMarkovRandomField::BuildNeighborhoodKernel()
{
  // Check if the neighborhood has already been created
  if (m2p_neighborhoodKernel)
  {
    Cerr("***** iPenaltyMarkovRandomField::BuildNeighborhoodKernel() -> Neighborhood kernel has already been created !" << endl);
    return 1;
  }
  // Get voxel sizes locally
  FLTNB vox_size_x = mp_ImageDimensionsAndQuantification->GetVoxSizeX();
  FLTNB vox_size_y = mp_ImageDimensionsAndQuantification->GetVoxSizeY();
  FLTNB vox_size_z = mp_ImageDimensionsAndQuantification->GetVoxSizeZ();
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
          // Exclude the voxel itself
          if (vox_x!=0 || vox_y!=0 || vox_z!=0)
          {
            // Compute the square distance between the two voxels
            FLTNB squareDistance = ((FLTNB)(vox_x*vox_x))*vox_size_x*vox_size_x
                                 + ((FLTNB)(vox_y*vox_y))*vox_size_y*vox_size_y
                                 + ((FLTNB)(vox_z*vox_z))*vox_size_z*vox_size_z;
            // Include voxel if the distance from the central voxel is less than the radius
            if (squareDistance<=m_neighborhoodSphereRadius*m_neighborhoodSphereRadius)
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
    m_neighborhoodMaxNbVoxels = (m_neighborhoodBoxOrder*2+1)*(m_neighborhoodBoxOrder*2+1)*(m_neighborhoodBoxOrder*2+1)-1;
    // Allocate the neighborhood kernel
    m2p_neighborhoodKernel = (INTNB**)malloc(m_neighborhoodMaxNbVoxels*sizeof(INTNB*));
    for (INTNB v = 0; v<m_neighborhoodMaxNbVoxels; v++) m2p_neighborhoodKernel[v] = (INTNB*)malloc(3*sizeof(INTNB));
    // Counter for the number of voxels in the neighborhood
    INTNB nb_voxels_in_neighborhood = 0;
    // Loop on increments of voxels indices in the box neighborhood
    for (INTNB vox_z = -neighborhood_box_order_z; vox_z<=neighborhood_box_order_z; vox_z++)
      for (INTNB vox_y = -m_neighborhoodBoxOrder; vox_y<=m_neighborhoodBoxOrder; vox_y++)
        for (INTNB vox_x = -m_neighborhoodBoxOrder; vox_x<=m_neighborhoodBoxOrder; vox_x++)
          // Exclude the voxel itself
          if (vox_x!=0 || vox_y!=0 || vox_z!=0)
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
    // Maximum number of voxels in the neighborhood is 6
    m_neighborhoodMaxNbVoxels = 6;
    // Particular case for 2D reconstruction
    if (mp_ImageDimensionsAndQuantification->GetNbVoxZ()==1) m_neighborhoodMaxNbVoxels = 4;
    // Allocate the neighborhood kernel
    m2p_neighborhoodKernel = (INTNB**)malloc(m_neighborhoodMaxNbVoxels*sizeof(INTNB*));
    for (INTNB v = 0; v<m_neighborhoodMaxNbVoxels; v++) m2p_neighborhoodKernel[v] = (INTNB*)malloc(3*sizeof(INTNB));
    // The 6 voxels
    for (INTNB v=0; v<m_neighborhoodMaxNbVoxels; v++) for (INTNB dim=0; dim<3; dim++) m2p_neighborhoodKernel[v][dim] = 0;
    m2p_neighborhoodKernel[0][MRF_NEIGHBOR_X] = -1;
    m2p_neighborhoodKernel[1][MRF_NEIGHBOR_X] = 1;
    m2p_neighborhoodKernel[2][MRF_NEIGHBOR_Y] = -1;
    m2p_neighborhoodKernel[3][MRF_NEIGHBOR_Y] = 1;
    if (mp_ImageDimensionsAndQuantification->GetNbVoxZ()!=1)
    {
      m2p_neighborhoodKernel[4][MRF_NEIGHBOR_Z] = -1;
      m2p_neighborhoodKernel[5][MRF_NEIGHBOR_Z] = 1;
    }
  }
  // ========================================================================================================
  // Unknown neighborhood type
  // ========================================================================================================
  else
  {
    Cerr("***** iPenaltyMarkovRandomField::BuildNeighborhoodKernel() -> The provided shape of neighborhood (" << m_neighborhoodShape << ") is unknown !" << endl);
    return 1;
  }
  // Check that there are some voxels in the neighborhood
  if (m_neighborhoodMaxNbVoxels<1)
  {
    Cerr("***** iPenaltyMarkovRandomField::BuildNeighborhoodKernel() -> There is no voxel in the neighborhood ! Check your neighborhood definition." << endl);
    return 1;
  }  
  // ========================================================================================================
  // With Bowsher's similarity factors, we compute the number of voxels that will be kept in the penalty
  // ========================================================================================================
  if (m_similarityType == MRF_SIMILARITY_BOWSHER)
  {
    m_similarityBowsherNbVoxels = (INTNB)round(m_similarityBowsherThreshold*((FLTNB)m_neighborhoodMaxNbVoxels)/100.);
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iPenaltyMarkovRandomField::BuildProximityFactors()
{
  // Check for allocations
  if (mp_proximityKernel)
  {
    Cerr("***** iPenaltyMarkovRandomField::BuildProximityFactors() -> Proximity kernel has already been created somewhere else !" << endl);
    return 1;
  }
  if (!m2p_neighborhoodKernel)
  {
    Cerr("***** iPenaltyMarkovRandomField::BuildProximityFactors() -> Neighborhood kernel must have been created and initialized before !" << endl);
    return 1;
  }
  // Compute the sum of proximity factors for normalization
  FLTNB proximity_factor_sum = 0.;
  // Allocate the proximity kernel
  mp_proximityKernel = (FLTNB*)malloc(m_neighborhoodMaxNbVoxels * sizeof(FLTNB));
  // Loop on the maximal neighborhood
  for (INTNB neigh=0; neigh<m_neighborhoodMaxNbVoxels; neigh++)
  {
    // Proximity factors based on inverse of euclidian distance
    if (m_proximityType == MRF_PROXIMITY_EUCLIDIAN)
    {
      mp_proximityKernel[neigh]  = fabs(m2p_neighborhoodKernel[neigh][MRF_NEIGHBOR_X])
                                 * mp_ImageDimensionsAndQuantification->GetVoxSizeX()
                                 * fabs(m2p_neighborhoodKernel[neigh][MRF_NEIGHBOR_X])
                                 * mp_ImageDimensionsAndQuantification->GetVoxSizeX();
      mp_proximityKernel[neigh] += fabs(m2p_neighborhoodKernel[neigh][MRF_NEIGHBOR_Y])
                                 * mp_ImageDimensionsAndQuantification->GetVoxSizeY()
                                 * fabs(m2p_neighborhoodKernel[neigh][MRF_NEIGHBOR_Y])
                                 * mp_ImageDimensionsAndQuantification->GetVoxSizeY();
      mp_proximityKernel[neigh] += fabs(m2p_neighborhoodKernel[neigh][MRF_NEIGHBOR_Z])
                                 * mp_ImageDimensionsAndQuantification->GetVoxSizeZ()
                                 * fabs(m2p_neighborhoodKernel[neigh][MRF_NEIGHBOR_Z])
                                 * mp_ImageDimensionsAndQuantification->GetVoxSizeZ();
      mp_proximityKernel[neigh] = 1./sqrt(mp_proximityKernel[neigh]);
    }
    // Proximity factors based on inverse of voxel distance (as described in many papers)
    else if (m_proximityType == MRF_PROXIMITY_VOXEL)
    {
      mp_proximityKernel[neigh]  = fabs(m2p_neighborhoodKernel[neigh][MRF_NEIGHBOR_X])
                                 * fabs(m2p_neighborhoodKernel[neigh][MRF_NEIGHBOR_X]);
      mp_proximityKernel[neigh] += fabs(m2p_neighborhoodKernel[neigh][MRF_NEIGHBOR_Y])
                                 * fabs(m2p_neighborhoodKernel[neigh][MRF_NEIGHBOR_Y]);
      mp_proximityKernel[neigh] += fabs(m2p_neighborhoodKernel[neigh][MRF_NEIGHBOR_Z])
                                 * fabs(m2p_neighborhoodKernel[neigh][MRF_NEIGHBOR_Z]);
      mp_proximityKernel[neigh] = 1./sqrt(mp_proximityKernel[neigh]);
    }
    /*
    // Proximity factors based on decreasing exponential of squared distance
    else if (m_proximityType == MRF_PROXIMITY_EXP_DIST_SQUARED)
    {
      FLTNB squareDistance = abs(m2p_neighborhoodKernel[neigh][MRF_NEIGHBOR_X])
                           * mp_ImageDimensionsAndQuantification->GetVoxSizeX()
                           * abs(m2p_neighborhoodKernel[neigh][MRF_NEIGHBOR_X])
                           * mp_ImageDimensionsAndQuantification->GetVoxSizeX();
      squareDistance += abs(m2p_neighborhoodKernel[neigh][MRF_NEIGHBOR_Y])
                      * mp_ImageDimensionsAndQuantification->GetVoxSizeY()
                      * abs(m2p_neighborhoodKernel[neigh][MRF_NEIGHBOR_Y])
                      * mp_ImageDimensionsAndQuantification->GetVoxSizeY();
      squareDistance += abs(m2p_neighborhoodKernel[neigh][MRF_NEIGHBOR_Z])
                      * mp_ImageDimensionsAndQuantification->GetVoxSizeZ()
                      * abs(m2p_neighborhoodKernel[neigh][MRF_NEIGHBOR_Z])
                      * mp_ImageDimensionsAndQuantification->GetVoxSizeZ();
      mp_proximityKernel[neigh] = exp(-squareDistance/(m_proximityCharacteristicDistance*m_proximityCharacteristicDistance))/squareDistance;
    }
    */
    // Uniform proximity factors
    else if (m_proximityType == MRF_PROXIMITY_NONE)
    {
      mp_proximityKernel[neigh] = 1.;
    }
    else
    {
      Cerr("***** iPenaltyMarkovRandomField::BuildProximityFactors() -> The provided type of proximity factor (" << m_proximityType << ") is unknown !" << endl);
      return 1;
    }
    // Add to the sum
    proximity_factor_sum += mp_proximityKernel[neigh];
  }
  // Normalize the proximity factors
  for (INTNB neigh=0; neigh<m_neighborhoodMaxNbVoxels; neigh++) mp_proximityKernel[neigh] /= proximity_factor_sum;
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iPenaltyMarkovRandomField::BuildSpecificNeighborhood(INTNB a_voxel, int a_th)
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

int iPenaltyMarkovRandomField::ComputeSimilarityFactors(int a_tbf, int a_rbf, int a_cbf, INTNB a_voxel, int a_th)
{
  // No similarity factor (uniform)
  if (m_similarityType == MRF_SIMILARITY_NONE)
  {
    // We just return as the factors were set to one after allocation
    return 0;
  }
  // Bowsher's similarity factors
  else if (m_similarityType == MRF_SIMILARITY_BOWSHER)
  {
    // Apply Bowsher's similarity only if there are more voxels in the neighborhood of this voxel than Bowsher's voxels
    if (m_similarityBowsherNbVoxels<mp_neighborhoodNbVoxels[a_th])
    {
      // Sort the neighborhood indices according to ascending absolute intensity difference from the voxel v
      std::sort( &(m2p_neighborhoodIndices[a_th][0]), &(m2p_neighborhoodIndices[a_th][m_neighborhoodMaxNbVoxels]),
      [a_voxel,a_tbf,a_rbf,a_cbf,a_th,this](INTNB i1, INTNB i2) 
      {
        // This two lines pushes the -1 indices (voxels out of image) at the end
        if (i1 == -1) return false;
        if (i2 == -1) return true;
        return (abs(mp_ImageSpace-> m2p_multiModalImage[0][i1] - mp_ImageSpace-> m2p_multiModalImage[0][a_voxel]))
             < (abs(mp_ImageSpace-> m2p_multiModalImage[0][i2] - mp_ImageSpace-> m2p_multiModalImage[0][a_voxel]));
      });
      // Consider the most similar neighbors as defined by the number of Bowsher's voxels
      for (INTNB n=0; n<m_similarityBowsherNbVoxels; n++) m2p_similarityFactors[a_th][n] = 1.;
      // And discard all the others
      for (INTNB n=m_similarityBowsherNbVoxels; n<mp_neighborhoodNbVoxels[a_th]; n++) m2p_similarityFactors[a_th][n] = 0.;
    }
  }
  // Unknown
  else
  {
    Cerr("***** iPenaltyMarkovRandomField::ComputeSimilarityFactors() -> Unknown similarity type provided !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iPenaltyMarkovRandomField::LocalPreProcessingStep(int a_tbf, int a_rbf, int a_cbf, INTNB a_voxel, int a_th)
{
  // First build the list of neighbours of this voxel
  if (BuildSpecificNeighborhood(a_voxel,a_th))
  {
    Cerr("***** iPenaltyMarkovRandomField::LocalPreProcessingStep() -> A problem occurred while building specific neighborhood of voxel " << a_voxel << " !" << endl);
    return 1;
  }
  // Then compute the similarity factors
  if (ComputeSimilarityFactors(a_tbf,a_rbf,a_cbf,a_voxel,a_th))
  {
    Cerr("***** iPenaltyMarkovRandomField::LocalPreProcessingStep() -> A problem occurred while computing the similirity factors of the neighborhood of voxel " << a_voxel << " !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

FLTNB iPenaltyMarkovRandomField::ComputePenaltyValue(int a_tbf, int a_rbf, int a_cbf, INTNB a_voxel, int a_th)
{
  // Compute the value of the penalty
  FLTNB result = 0.;
  // Sum on the contribution of the different neighbors
  for (INTNB n=0; n<m_neighborhoodMaxNbVoxels; n++)
  {
    // Get the neighbor index
    INTNB neighbor = m2p_neighborhoodIndices[a_th][n];
    // If the voxel is not in the neighborhood, skip it
    if (neighbor==-1) continue;
    else
    {
      // Pointer to the image
      FLTNB* p_image = mp_ImageSpace->m4p_image[a_tbf][a_rbf][a_cbf];
      // Proximity factor
      FLTNB proximity_factor = mp_proximityKernel[n];
      // Similarity factor
      FLTNB similarity_factor = m2p_similarityFactors[a_th][n];
      // Precompute the difference
      FLTNB difference = p_image[a_voxel] - p_image[neighbor];
      // Compute value of the potential function (consider the potential function to be symmetric)
      FLTNB value = 0.;
      switch(m_potentialType)
      {
        // Quadratic
        case MRF_POTENTIAL_QUADRATIC:
        {
          value = 0.5 * difference * difference;
          break;
        }
        // Relative differences
        case MRF_POTENTIAL_RELATIVE_DIFFERENCE:
        {
          FLTNB denominator = p_image[a_voxel] + p_image[neighbor] + m_potentialRelativeDifferenceGamma * abs(difference);
          if (denominator==0.) value = 0.;
          else value = difference * difference / denominator;
          break;
        }
        // Geman and McClure + 1 (we add 1 to the original definition to make it positive)
        case MRF_POTENTIAL_GEMAN_MCCLURE:
        {
          FLTNB squared_difference = difference * difference;
          value = squared_difference / (m_potentialGemanMcClureDelta*m_potentialGemanMcClureDelta + squared_difference);
          break;
        }
        // Green logcosh
        case MRF_POTENTIAL_GREEN:
        {
          value = log(cosh( difference / m_potentialGreenLogCoshDelta));
          break;
        }
        // Hebert and Leahy
        case MRF_POTENTIAL_HEBERT_LEAHY:
        {
          value = log( 1. + difference * difference / (m_potentialHebertLeahyMu * m_potentialHebertLeahyMu));
          break;
        }
        // Huber piecewise
        case MRF_POTENTIAL_HUBER:
        {
          FLTNB abs_difference = fabs(difference);
          if (abs_difference>m_potentialHuberDelta) value = m_potentialHuberDelta*abs_difference - 0.5*m_potentialHuberDelta*m_potentialHuberDelta;
          else value = 0.5 * abs_difference * abs_difference;
          break;
        }
      }
      // Check that the value is a normal number
      if (fpclassify(value) != FP_NORMAL) value = 0.;
      // Add the final contribution to the penalty value
      else result += proximity_factor * similarity_factor * value;
    }
  }
  // Apply the penalty strength
  result *= m_penaltyStrength;
  // Check for Inf, Nan, etc
  if (fpclassify(result) != FP_NORMAL) result = 0.;
  // Return result
  return result;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

FLTNB iPenaltyMarkovRandomField::ComputeFirstDerivative(int a_tbf, int a_rbf, int a_cbf, INTNB a_voxel, int a_th)
{
  // Compute the first derivative of the penalty
  FLTNB result = 0.;
  // Sum on the contribution of the different neighbors
  for (INTNB n=0; n<m_neighborhoodMaxNbVoxels; n++)
  {
    // Get the neighbor index
    INTNB neighbor = m2p_neighborhoodIndices[a_th][n];
    // If the voxel is not in the neighborhood, skip it
    if (neighbor==-1) continue;
    else
    {
      // Pointer to the image
      FLTNB* p_image = mp_ImageSpace->m4p_image[a_tbf][a_rbf][a_cbf];
      // Proximity factor
      FLTNB proximity_factor = mp_proximityKernel[n];
      // Similarity factor
      FLTNB similarity_factor = m2p_similarityFactors[a_th][n];
      // Precompute the difference
      FLTNB difference = p_image[a_voxel] - p_image[neighbor];
      // Compute first derivative of the potential function (consider the potential function to be symmetric)
      FLTNB first_derivative = 0.;
      switch(m_potentialType)
      {
        // Quadratic
        case MRF_POTENTIAL_QUADRATIC:
        {
          first_derivative = difference;
          break;
        }
        // Relative differences
        case MRF_POTENTIAL_RELATIVE_DIFFERENCE:
        {
          FLTNB term = p_image[a_voxel] + p_image[neighbor] + m_potentialRelativeDifferenceGamma * fabs(difference);
          FLTNB term_square = term * term;
          if (term_square==0.) first_derivative = 0.;
          else first_derivative = difference * (term + 2.*p_image[neighbor]) / (term_square);
          break;
        }
        // Geman and McClure
        case MRF_POTENTIAL_GEMAN_MCCLURE:
        {
          FLTNB difference = p_image[a_voxel] - p_image[neighbor];
          FLTNB delta_square = m_potentialGemanMcClureDelta * m_potentialGemanMcClureDelta;
          first_derivative = 2. * delta_square * difference
                           / pow( difference*difference + delta_square , 2. );
          break;
        }
        // Green logcosh
        case MRF_POTENTIAL_GREEN:
        {
          FLTNB rel_diff = difference / m_potentialGreenLogCoshDelta;
          first_derivative = tanh(rel_diff) / m_potentialGreenLogCoshDelta;
          break;
        }
        // Hebert and Leahy
        case MRF_POTENTIAL_HEBERT_LEAHY:
        {
          first_derivative = 2. * difference / (difference*difference + m_potentialHebertLeahyMu*m_potentialHebertLeahyMu);
          break;
        }
        // Huber piecewise
        case MRF_POTENTIAL_HUBER:
        {
          if (difference==0.) first_derivative = 0.;
          else
          {
            FLTNB abs_difference = fabs(difference);
            if (abs_difference>m_potentialHuberDelta) first_derivative = m_potentialHuberDelta * difference / abs_difference;
            else first_derivative = difference;
          }
          break;
        }
      }
      // Check that the derivative is a normal number
      if (fpclassify(first_derivative) != FP_NORMAL) first_derivative = 0.;
      // Add the final contribution to the derivative value
      else result += proximity_factor * similarity_factor * first_derivative;
    }
  }
  // Apply the penalty strength
  result *= m_penaltyStrength;
  // Check for Inf, Nan, etc
  if (fpclassify(result) != FP_NORMAL) result = 0.;
  // Return result
  return result;
}


// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

FLTNB iPenaltyMarkovRandomField::ComputeSecondDerivative(int a_tbf, int a_rbf, int a_cbf, INTNB a_voxel, int a_th)
{
  // Compute the second derivative of the penalty
  FLTNB result = 0.;
  // Sum on the contribution of the different neighbors
  for (INTNB n=0; n<m_neighborhoodMaxNbVoxels; n++)
  {
    // Get the neighbor index
    INTNB neighbor = m2p_neighborhoodIndices[a_th][n];
    // If the voxel is not in the neighborhood, skip it
    if (neighbor==-1) continue;
    else
    {
      // Pointer to the image
      FLTNB* p_image = mp_ImageSpace->m4p_image[a_tbf][a_rbf][a_cbf];
      // Proximity factor
      FLTNB proximity_factor = mp_proximityKernel[n];
      // Similarity factor
      FLTNB similarity_factor = m2p_similarityFactors[a_th][n];
      // Compute second derivative of the potential function (consider the potential function to be symmetric)
      FLTNB second_derivative = 0.;
      switch(m_potentialType)
      {
        // Quadratic
        case MRF_POTENTIAL_QUADRATIC:
        {
          second_derivative = 1.;
          break;
        }
        // Relative differences
        case MRF_POTENTIAL_RELATIVE_DIFFERENCE:
        {
          FLTNB term = p_image[a_voxel] + p_image[neighbor] + m_potentialRelativeDifferenceGamma * fabs(p_image[a_voxel] - p_image[neighbor]);
          FLTNB term_cube = term * term * term;
          if (term_cube==0.) second_derivative = 0.;
          else second_derivative = (8.*p_image[neighbor]*p_image[neighbor])/(term_cube);
          break;
        }
        // Geman and McClure
        case MRF_POTENTIAL_GEMAN_MCCLURE:
        {
          FLTNB difference = p_image[a_voxel] - p_image[neighbor];
          FLTNB difference_square = difference * difference;
          FLTNB delta_square = m_potentialGemanMcClureDelta * m_potentialGemanMcClureDelta;
          second_derivative = -2. * delta_square * (3.*difference_square - delta_square)
                            / pow( difference_square + delta_square , 3. );
          break;
        }
        // Green logcosh
        case MRF_POTENTIAL_GREEN:
        {
          second_derivative = 1. / (m_potentialGreenLogCoshDelta * cosh((p_image[a_voxel]-p_image[neighbor])/m_potentialGreenLogCoshDelta));
          second_derivative *= second_derivative;
          break;
        }
        // Hebert and Leahy
        case MRF_POTENTIAL_HEBERT_LEAHY:
        {
          FLTNB difference = p_image[a_voxel] - p_image[neighbor];
          second_derivative = - 2. * (difference+m_potentialHebertLeahyMu) * (difference-m_potentialHebertLeahyMu)
                            / pow( difference*difference + m_potentialHebertLeahyMu*m_potentialHebertLeahyMu , 2. );
          break;
        }
        // Huber piecewise
        case MRF_POTENTIAL_HUBER:
        {
          // Assume null derivative in all domain
          second_derivative = 0.;
          break;
        }
      }
      // Check that the derivative is a normal number
      if (fpclassify(second_derivative) != FP_NORMAL) second_derivative = 0.;
      // Add the final contribution to the derivative value
      else result += proximity_factor * similarity_factor * second_derivative;
    }
  }
  // Apply the penalty strength
  result *= m_penaltyStrength;
  // Check for Inf, Nan, etc
  if (fpclassify(result) != FP_NORMAL) result = 0.;
  // Return result
  return result;
}

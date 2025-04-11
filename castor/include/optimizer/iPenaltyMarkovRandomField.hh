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
  \brief    Declaration of class iPenaltyMarkovRandomField
*/

#ifndef IPENALTYMARKOVRANDOMFIELD_HH
#define IPENALTYMARKOVRANDOMFIELD_HH 1

#include "vPenalty.hh"
#include "sAddonManager.hh"
#include "sOutputManager.hh"

#define MRF_POTENTIAL_QUADRATIC 0
#define MRF_POTENTIAL_RELATIVE_DIFFERENCE 1
#define MRF_POTENTIAL_GEMAN_MCCLURE 2
#define MRF_POTENTIAL_GREEN 3
#define MRF_POTENTIAL_HEBERT_LEAHY 4
#define MRF_POTENTIAL_HUBER 5

#define MRF_PROXIMITY_NONE 0
#define MRF_PROXIMITY_EUCLIDIAN 1
#define MRF_PROXIMITY_VOXEL 2
//#define MRF_PROXIMITY_EXP_DIST_SQUARED 3

#define MRF_SIMILARITY_NONE 0
#define MRF_SIMILARITY_BOWSHER 1

#define MRF_NEIGHBORHOOD_SPHERE 0
#define MRF_NEIGHBORHOOD_BOX 1
#define MRF_NEIGHBORHOOD_6_NEAREST 2

#define MRF_NOT_DEFINED -1

#define MRF_NEIGHBOR_X 0
#define MRF_NEIGHBOR_Y 1
#define MRF_NEIGHBOR_Z 2

/*!
  \class   iPenaltyMarkovRandomField
  \brief   This class implements the general Markov Random Field kind of penalties
  \details This class inherits from vPenalty and implements the general MRF penalties.
*/
class iPenaltyMarkovRandomField : public vPenalty
{
  // -------------------------------------------------------------------
  // Constructor & Destructor
  public:

    /*!
      \fn      public iPenaltyMarkovRandomField::iPenaltyMarkovRandomField()
      \brief   The constructor of iPenaltyMarkovRandomField
      \details This is the default and unique constructor. It does not take any parameter and
               its role is only to affect default values to each member of the class.
    */
    iPenaltyMarkovRandomField();
    /*!
      \fn      public iPenaltyMarkovRandomField::~iPenaltyMarkovRandomField()
      \brief   The destructor of iPenaltyMarkovRandomField
      \details This is the default and unique destructor. It does not take any parameter and
               its role is only to free or delete all structures that were built by this class.
    */
    ~iPenaltyMarkovRandomField();

  // -------------------------------------------------------------------
  // Public member functions
  public:
    // Function for automatic insertion (put the class name as the parameter and do not add semi-column at the end of the line)
    FUNCTION_PENALTY(iPenaltyMarkovRandomField)
    /*!
      \fn      public int iPenaltyMarkovRandomField::ReadConfigurationFile()
      \param   const string& a_configurationFile
      \brief   A function used to read options from a configuration file
      \details This function implements the reading of all options associated to the child penalty, from
               a configuration file. It is the implementation of the pure virtual function inherited
               from the abstract class vPenalty. It checks the reading status but not
               the options values that will be checked by the CheckSpecificParameters() function.
      \return  An integer reflecting the reading success; 0 if success, another value otherwise.
    */
    int ReadConfigurationFile(const string& a_configurationFile);
    /*!
      \fn      public int iPenaltyMarkovRandomField::ReadOptionsList()
      \param   const string& a_optionsList
      \brief   A function used to read options from a list of options
      \details This function implements the reading of all options associated to the child penalty, from
               a list of options. It is the implementation of the pure virtual function inherited
               from the abstract class vPenalty. It checks the reading status but not
               the options values that will be checked by the CheckSpecificParameters() function.
      \return  An integer reflecting the reading success; 0 if success, another value otherwise.
    */
    int ReadOptionsList(const string& a_optionsList);
    /*!
      \fn      public int iPenaltyMarkovRandomField::LocalPreProcessingStep()
      \param   int a_tbf
      \param   int a_rbf
      \param   int a_cbf
      \param   INTNB a_voxel
      \param   int a_th
      \brief   A public function computing a local pre-processing step for the penalty
      \details This function overloads the mother class function that does nothing. The idea here is
               to build the specific neighborhood of the given voxel and to compute all related weights.
      \return  An integer reflecting the initialization status; 0 if no problem, another value otherwise.
    */
    int LocalPreProcessingStep(int a_tbf, int a_rbf, int a_cbf, INTNB a_voxel, int a_th);
    /*!
      \fn      public int iPenaltyMarkovRandomField::ComputePenaltyValue()
      \param   int a_tbf
      \param   int a_rbf
      \param   int a_cbf
      \param   INTNB a_voxel
      \param   int a_th
      \brief   Implementation of the pure virtual vPenalty::ComputePenaltyValue()
      \details This function computes the value of the penalty function for the provided indices.
               It is the implementation of the pure virtual vPenalty::ComputePenaltyValue().
      \return  The penalty value
    */
    FLTNB ComputePenaltyValue(int a_tbf, int a_rbf, int a_cbf, INTNB a_voxel, int a_th);
    /*!
      \fn      public int iPenaltyMarkovRandomField::ComputeFirstDerivative()
      \param   int a_tbf
      \param   int a_rbf
      \param   int a_cbf
      \param   INTNB a_voxel
      \param   int a_th
      \brief   Implementation of the pure virtual vPenalty::ComputeFirstDerivative()
      \details This function computes the first derivative of the penalty.
               It is the implementation of the pure virtual vPenalty::ComputeFirstDerivative().
      \return  The derivative
    */
    FLTNB ComputeFirstDerivative(int a_tbf, int a_rbf, int a_cbf, INTNB a_voxel, int a_th);
    /*!
      \fn      public int iPenaltyMarkovRandomField::ComputeSecondDerivative()
      \param   int a_tbf
      \param   int a_rbf
      \param   int a_cbf
      \param   INTNB a_voxel
      \param   int a_th
      \brief   Implementation of the pure virtual vPenalty::ComputeSecondDerivative()
      \details This function computes the second derivative of the penalty.
               It is the implementation of the pure virtual vPenalty::ComputeSecondDerivative().
      \return  The derivative
    */ 
    FLTNB ComputeSecondDerivative(int a_tbf, int a_rbf, int a_cbf, INTNB a_voxel, int a_th);
    /*!
      \fn      public int iPenaltyMarkovRandomField::BuildSpecificNeighborhood()
      \param   INTNB a_voxel
      \param   int a_th
      \brief   Computes the specific neighborhood of a voxel, m2p_neighborhoodIndices, as well as mp_neighborhoodNbVoxels
      \details This function computes the specific neighborhood of a voxel, from the neighborhood kernel and the provided index.
               From the neighborhood kernel that contains indices shifts, it computes the actual voxel indices of each
               neighbor and check if it is outside of boundaries. It records the index voxel if valid and -1 if not.
               This function use a thread index to be thread safe during computations.
      \return  An integer reflecting the check status; 0 if no problem, another value otherwise.
    */   
    int BuildSpecificNeighborhood(INTNB a_voxel, int a_th);
    
  // -------------------------------------------------------------------
  // Public Get & Set functions
  public:
    /*!
      \fn      public inline int iPenaltyMarkovRandomField::GetPotentialType()
      \return  m_potentialType
    */
    inline int GetPotentialType() {return m_potentialType;}
    /*!
      \fn      public inline FLTNB* iPenaltyMarkovRandomField::GetProximityKernel()
      \return  mp_proximityKernel
    */
    inline FLTNB* GetProximityKernel() {return mp_proximityKernel;}
    /*!
      \fn      public inline FLTNB* iPenaltyMarkovRandomField::GetSimilarityFactors()
      \param   int a_th
      \return  m2p_similarityFactors[a_th]
    */
    inline FLTNB* GetSimilarityFactors(int a_th) {return m2p_similarityFactors[a_th];}
    /*!
      \fn      public inline INTNB iPenaltyMarkovRandomField::GetNeighborhoodMaxNbVoxels()
      \return  m_neighborhoodMaxNbVoxels
    */
    inline INTNB GetNeighborhoodMaxNbVoxels() {return m_neighborhoodMaxNbVoxels;}
    /*!
      \fn      public inline INTNB** iPenaltyMarkovRandomField::GetNeighborhoodIndices()
      \return  m2p_neighborhoodIndices
    */
    inline INTNB** GetNeighborhoodIndices() {return m2p_neighborhoodIndices;}

  // -------------------------------------------------------------------
  // Private member functions
  private:
    /*!
      \fn      private void iPenaltyMarkovRandomField::ShowHelpSpecific()
      \brief   A function used to show help about the child penalty
      \details This function must describe what the module does and how to use it. It describes in
               details the different parameters of the penalty, and how to set them through the use
               of a configuration file or a list of options. It is the implementation of the pure
               virtual function inherited from the abstract class vPenalty. It is called by the
               public ShowHelp() function.
    */
    void ShowHelpSpecific();
    /*!
      \fn      private int iPenaltyMarkovRandomField::CheckSpecificParameters()
      \brief   A private function used to check the parameters settings specific to the child penalty
      \details This function is used to check that all parameters specific to the penalty are correctly set
               within allowed values. It is called by the CheckParameters() function of the mother class.
               It is the implementation of the pure virtual function inherited from the abstract mother
               class vPenalty.
      \return  An integer reflecting the check status; 0 if no problem, another value otherwise.
    */
    int CheckSpecificParameters();
    /*!
      \fn      private int iPenaltyMarkovRandomField::InitializeSpecific()
      \brief   This function is used to initialize specific stuff to the child penalty.
      \details It is called by the public Initialize() function from the mother.
      \return  An integer reflecting the initialization status; 0 if no problem, another value otherwise.
    */
    int InitializeSpecific();
    /*!
      \fn      private int iPenaltyMarkovRandomField::BuildNeighborhoodKernel()
      \brief   A function used to build the neighborhood kernel
      \details This function build the table m2p_neighborhoodKernel, that prepare neighborhood indices operations.
               It is called only one time at initialization. The kernel will be used when computing the specific
               neighborhood of any voxel later on.
      \return  An integer reflecting the reading success; 0 if success, another value otherwise.
    */
    int BuildNeighborhoodKernel();
    /*!
      \fn      private int iPenaltyMarkovRandomField::BuildProximityFactors()
      \brief   A function used to build the kernel of the proximity factors
      \details This function build the table FLTNB* mp_proximityKernel
      \return  An integer reflecting the reading success; 0 if success, another value otherwise.
    */
    int BuildProximityFactors();
    /*!
      \fn      private int iPenaltyMarkovRandomField::ComputeSimilarityFactors()
      \param   int a_tbf
      \param   int a_rbf
      \param   int a_cbf
      \param   INTNB a_voxel
      \param   int a_th
      \brief   Computes the similarity factors of this specific voxel with respect to the similarity type used.
      \details Store the results in m2p_similarityFactors
      \return  An integer reflecting the check status; 0 if no problem, another value otherwise.
    */
    int ComputeSimilarityFactors(int a_tbf, int a_rbf, int a_cbf, INTNB a_voxel, int a_th);

  // -------------------------------------------------------------------
  // Data members
  private:
    // Potential functions
    int m_potentialType;                      /*!< Type of the potential function used (e.g. quadratic) */
    FLTNB m_potentialRelativeDifferenceGamma; /*!< Parameter of the relative difference potential function */
    FLTNB m_potentialGemanMcClureDelta;       /*!< Parameter of the Geman and McClure's potential function */
    FLTNB m_potentialGreenLogCoshDelta;       /*!< Parameter of the Green's log-cosh potential function */
    FLTNB m_potentialHebertLeahyMu;           /*!< Parameter of the Hebert and Leahy's potential function */
    FLTNB m_potentialHuberDelta;              /*!< Parameter of the piecewise Huber potential function */
    // Proximity factors
    int m_proximityType;                      /*!< Type of weighting based on proximity (e.g. the inverse of the euclidian distance between voxels) */
    FLTNB m_proximityCharacteristicDistance;  /*!< When the proximity factors are the exponential distance squared, this is the characteristic distance */
    FLTNB* mp_proximityKernel;                /*!< Precomputation of the proximity factors. Index is on neighboring voxels, as for m2p_neighborhoodKernel */
    // Similarity factors
    int m_similarityType;                     /*!< Type of weighting based on similarity (e.g. joint entropy) */
    FLTNB m_similarityBowsherThreshold;       /*!< Percentage of the neighborhood voxels for defining the most similar Bowsher neighbor voxels */
    FLTNB** m2p_similarityFactors;            /*!< Similarity factors of the neighborhood. First index on threads, second on neighboring voxels, as for m2p_neighborhoodIndices */
    INTNB m_similarityBowsherNbVoxels;        /*!< The number of voxels kept from the Bowsher's threashold */
    // Neighborhood
    int       m_neighborhoodShape;             /*!< Type of neighborhood (sphere, box) */
    FLTNB     m_neighborhoodSphereRadius;      /*!< In case of a spherical neighborhood, radius of the sphere */
    int       m_neighborhoodBoxOrder;          /*!< In case of a box neighborhood, order of the box */
    int       m_neighborhoodBoxExcludeCorners; /*!< In case of a box neighborhood, remove the corners as in some papers, or not (keep them by default) */
    INTNB     m_neighborhoodMaxNbVoxels;       /*!< Maximum number of voxels in the neighborhood */
    INTNB** m2p_neighborhoodKernel;            /*!< Neighborhood of a virtual voxel of coordinates 0, 0, 0 without boundary. First index on neighboring voxels, second one on indexes */
    INTNB*   mp_neighborhoodNbVoxels;          /*!< Number of voxels in a specific neighborhood. The index is on threads */
    INTNB** m2p_neighborhoodIndices;           /*!< Neighborhood indices of a specific voxel. First index is on the threads, second one contains the list of the voxels in the neighborhood (as many as m_neighborhoodMaxNbVoxels). It follows strictly the same order than m2p_neighborhoodKernel, with the convention -1: out of image boundaries */
};

// Class for automatic insertion (set here the visible optimizer's name as the first parameter,
// put the class name as the second parameter and do NOT add semi-colon at the end of the line)
CLASS_PENALTY(MRF,iPenaltyMarkovRandomField)

#endif


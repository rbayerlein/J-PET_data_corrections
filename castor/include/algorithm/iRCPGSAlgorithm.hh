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
  \file iRCPGSAlgorithm.hh
  \class iRCPGSAlgorithm
  \brief RCP-GS : Random Clustering Prior - Gibbs Sampler
  \details Sampling of the posterior probability distribution (Bayesian inference)
  Image prior : ddCRP, with intensity ~ Gamma (conjugate prior for Poisson likelihood)
  Multinomial distribution for the backprojection from acquired data into (latent) emissions from voxel j into projection i
*/

#ifndef IRCPGSALGORITHM_HH
#define IRCPGSALGORITHM_HH 1


#include "vAlgorithm.hh"
#include <functional>


class iRCPGSAlgorithm : public vAlgorithm
{
  // Constructor & Destructor
  public:
    iRCPGSAlgorithm();
    ~iRCPGSAlgorithm();
    
    int InitSpecificOptions(string a_specificOptions);
    void ShowHelpSpecific();

  // Protected overriden virtual member functions
  protected:
    int StepBeforeIterationLoop();
    int StepAfterIterationLoop();
    int StepBeforeSubsetLoop(int a_iteration);
    int StepAfterSubsetLoop(int a_iteration);
    int StepPreProcessInsideSubsetLoop(int a_iteration, int a_subset);
    int StepInnerLoopInsideSubsetLoop(int a_iteration, int a_subset, int a_bed);
    int StepPostProcessInsideSubsetLoop(int a_iteration, int a_subset);
    int StepImageOutput(int a_iteration, int a_subset = -1);
    int ReadConfigurationFile(const string& a_configurationFile);
    
  // Private member functions
  private:
    /*!
      \fn private int iRCPGSAlgorithm::SampleConditionalCompleteData()
      \brief Gibbs sampler : first conditional probability (backprojection)
      \return 0 if success, positive value otherwise.
    */
    int SampleConditionalCompleteData(int a_iteration, int a_subset, int a_bed);
    /*!
      \fn private int iRCPGSAlgorithm::SampleConditionalClustering()
      \brief Gibbs sampler : second conditional probability (ddCRP links)
      \return 0 if success, positive value otherwise.
    */
    int SampleConditionalClustering(int a_iteration);
    /*!
      \fn private int iRCPGSAlgorithm::SampleConditionalClusterIntensity()
      \brief Gibbs sampler : third conditional probability (cluster intensity)
      \return 0 if success, positive value otherwise.
    */
    int SampleConditionalClusterIntensity();
    /*!
      \fn private int iRCPGSAlgorithm::UpdateVisitedVoxels()
      \brief Check for voxels which do not contribute to recorded observed data
      \return 0 if success, positive value otherwise.
    */
    int UpdateVisitedVoxels();
    /*!
      \fn private int iRCPGSAlgorithm::InitializeHelperVar()
      \brief Allocate and initialize temporary variables, after the main variables have been initialized
      Assumption : the initial image and the sensitivity image have already been initialized 
      \return 0 if success, positive value otherwise.
    */
    int InitializeHelperVar();
    /*!
      \fn private int iRCPGSAlgorithm::ProcessMultiModalInfo()
      \brief Check input multimodal images
      \return 0 if success, positive value otherwise.
    */
    int ProcessMultiModalInfo();
    /*!
      \fn private int iRCPGSAlgorithm::GenerateCurrentImage()
      \brief Generate the current image estimation from the current sample of partition/clustering and cluster intensities
      \return 0 if success, positive value otherwise.
    */
    int GenerateCurrentImage();
    /*!
      \fn private int iRCPGSAlgorithm::ComputeFellowVoxelsList()
      \brief Fill the input list of fellow voxels for the current type of neighbourhood and the input voxel, for ddCRP
      \return 0 if success, positive value otherwise.
    */
    int ComputeFellowVoxelsList(vector<INTNB>& a_fellow_voxels, int a_current_voxel);
    /*!
      \fn private int iRCPGSAlgorithm::ComputeSumsPerClusters()
      \brief Compute sums of voxel features for each cluster
      \return 0 if success, positive value otherwise.
    */
    int ComputeSumsPerClusters(int a_iteration);
    /*!
      \fn private int iRCPGSAlgorithm::SaveIntermediaryData()
      \brief Save potentially useful by-product variables specific to this algorithm
      \return 0 if success, positive value otherwise.
    */
    int SaveIntermediaryData(int a_iteration);
    
    
  // Private data members
  private:
    // algorithm parameters
    int m_ddCRP; /*!< flag describing the ddCRP prior : 0 = no ddCRP, 1 = original ddCRP, 2 = modified ddCRP */
    int m_backprojection; /*!< flag describing the multinomial backprojection of the current iteration/subset data (sampling of conditional complete data):
                              0 = the previous backprojection state is cleared, as for ML-EM
                              1 = update of the previous backprojection state
                              2 = backprojection marginalized over cluster intensity, implies update of the previous backprojection state */
    int m_neighbourhood; /*!< Number of voxels in the neighbourhood */
    HPFLTNB m_gammaShape; /*!< Gamma prior distribution shape parameter */
    HPFLTNB m_gammaRate; /*!< Gamma prior distribution rate parameter */
    HPFLTNB m_ddcrpAlpha; /*!< ddCRP : the unnormalized probability of drawing a self link */
    INTNB m_multiModalLag; /*!< Number of iterations after which the multimodal images start affecting voxels clustering */
    HPFLTNB m_ddcrpAlphaIncrement; /*!< Multiplicative increment for optimizing the ddCRP alpha parameter through iterations */

    HPFLTNB* mp_multiModalNoiseSigma; /*!< Standard deviation of Gaussian noise in multimodal images */
    HPFLTNB* mp_multiModalParam; /*!< Parameter for the standard deviation of the prior on multimodal images */

    INTNB* mp_voxelClusterMapping;/*!< Mapping : voxel index -> cluster index */
    HPFLTNB* mp_clusterValues;/*!< Values (here voxel intensities) associated with each cluster */
    HPFLTNB** mp_clusterN;/*!< Ns : sum over cluster voxels and all LORs of latent backprojected variable (complete data), requires HPFLTNB precision, can be threaded  */
    HPFLTNB* mp_clusterSensitivity;/*!< Sensitivity : sum over cluster voxels and all LORs of system matrix components, requires HPFLTNB precision  */
    HPFLTNB** mpp_clusterMultiModal;/*!< Multimodal image(s) : sum over cluster voxels, requires HPFLTNB precision */
    INTNB* mp_clusterCount;/*!< Number of voxels per cluster */
    INTNB* mp_nextLink;/*!< Mapping : voxel index -> next voxel index (link from this voxel to another voxel)  */
    vector<INTNB>* mpv_parentLinks;/*!< Mapping : voxel index -> list of previous voxels (links from other voxels to this voxel)   */
    vector<INTNB> mv_newClusters;/*!< Indices of free available (empty) clusters (the number of clusters <= the number of voxels + 1)  */
    INTNB* mp_listVoxelIndices;/*!< List of voxel indices, used for random shuffling  */
    INTNB** mp_listEventsIndices;/*!< List of events indices, per bed, used for random shuffling     */
    INTNB**** m4p_EventsBackprojectionTrace;/*!< Trace of the backprojection destination for each observed count at previous iteration, 
                                              Number of beds x Number of events x Number of TOF bins x Number of counts */
    bool* mp_listRelevantVoxelIndices;/*!< Indices of voxels which will be taken into account in the algorithm; all the other voxels are regarded as background and set to 0 */
    HPFLTNB* temp_count_multimodal; /*!< helper variable for the second sampling step, just to avoid allocating repeatedly >*/
    HPFLTNB m_ddcrpLogAlpha; /*!< helper variable for precomputed log of ddCRP : the probability of drawing a self link (not normalized) */
    HPFLTNB m_currentMeanClusterVolume; /*!< helper variable for the average cluster volume after each iteration/subset */
    HPFLTNB m_meanClusterVolumeMin; /*!< Minimum threshold for average cluster volume (in mm3), used for adapting ddCRP alpha automatically through iterations */
    HPFLTNB m_meanClusterVolumeMax; /*!< Maximum threshold for average cluster volume (in mm3), used for adapting ddCRP alpha automatically through iterations */
    FLTNB* mp_permanentBackwardImage; /*!< A buffer image used for holding the complete up-to-date backward image at any time.
                                             It is updated at each iteration/subset, and used in backprojection modes 1 and 2 */

};

#endif














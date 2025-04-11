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
  \ingroup  algorithm
  \brief    Declaration of class vAlgorithm
*/

#ifndef VALGORITHM_HH
#define VALGORITHM_HH 1

#include "gVariables.hh"
#include "oImageDimensionsAndQuantification.hh"
#include "oProjectorManager.hh"
#include "oImageConvolverManager.hh"
#include "oImageProcessingManager.hh"
#include "oOptimizerManager.hh"
#include "oDeformationManager.hh"
#include "oDynamicModelManager.hh"
#include "oImageSpace.hh"
#include "vDataFile.hh"
#include "sChronoManager.hh"


/*!
  \class  vAlgorithm
  \brief  This is the base class for reconstructions, containing a framework with iteration and data subset loops. \n
          It contains all the managers and the images.
*/
class vAlgorithm
{
  // Constructor & Destructor
  public:
    /*!
      \fn    public vAlgorithm::vAlgorithm
      \brief vAlgorithm constructor. 
             Initialize the member variables to their default values.
    */
    vAlgorithm();
    /*!
      \fn    public vAlgorithm::~vAlgorithm
      \brief vAlgorithm destructor. 
    */
    virtual ~vAlgorithm();



  // -------------------------------------------------------------------
  // Public member functions
  public:
    /*!
      \fn public vAlgorithm::Run
      \brief Just call either the RunCPU or the RunGPU function as asked for
      \return 0 if success, positive value otherwise
    */
    int Run();
    /*!
      \fn public vAlgorithm::RunCPU
      \brief Perform the iterative loop of the algorithm.
             Function designed to be executed on the CPU only.
      \details Loops over the iterations, data subsets, bed positions \n
                \n
               StepBeforeIterationLoop() \n
               / Loop on iterations \n
               | StepBeforeSubsetLoop(iteration) \n
               |  / Loop on data subsets \n
               |  | StepPreProcessInsideSubsetLoop(iteration,subset) \n
               |  | / Loop on bed positions \n
               |  | | StepInnerLoopInsideSubsetLoop(iteration,subset,bed) \n
               |  | StepPostProcessInsideSubsetLoop(iteration,subset) \n
               |  StepAfterSubsetLoop(iteration) \n
               StepAfterIterationLoop() \n
          
      \return 0 if success, positive value otherwise
    */
    int RunCPU();
    #ifdef CASTOR_GPU
    /*!
      \fn public vAlgorithm::RunGPU
      \brief Perform the iterative loop of the algorithm.
             Function designed to be executed on the GPU only. \n
             This function is NOT yet implemented.
      \return 0 if success, positive value otherwise
    */
    int RunGPU();
    #endif
    /*!
      \fn     inline public void vAlgorithm::SetSaveSensitivityHistoFlag()
      \param  bool a_saveSensitivityHistoFlag
      \brief  Set the flag that specifies if the sensitivity image in histogram mode has to be saved for each subset/iteration
    */
    inline void SetSaveSensitivityHistoFlag(bool a_saveSensitivityHistoFlag)
           {m_saveSensitivityHistoFlag = a_saveSensitivityHistoFlag;}
    /*!
      \fn     inline public void vAlgorithm::SetSaveSubsetImageFlag()
      \param  bool a_saveImageAfterSubsets
      \brief  Set the flag that specifies if the image has to be saved for each subset
    */
    inline void SetSaveSubsetImageFlag(bool a_saveImageAfterSubsets)
           {m_saveImageAfterSubsets = a_saveImageAfterSubsets;}
    /*!
      \fn     inline public void vAlgorithm::SetSaveDynamicBasisCoefficientImages()
      \param  bool a_saveDynamicBasisCoefficients
      \brief  Set the flag that specifies if the dynamic basis functions coefficients images have to be saved
    */
    inline void SetSaveDynamicBasisCoefficientImages(bool a_saveDynamicBasisCoefficients)
           {m_saveDynamicBasisCoefficients = a_saveDynamicBasisCoefficients;}
    /*!
      \fn inline public vAlgorithm::SetOptimizerManager
      \param ap_OptimizerManager
      \brief Set the Optimizer Manager Object 
    */
    inline void SetOptimizerManager(oOptimizerManager* ap_OptimizerManager)
           {mp_OptimizerManager = ap_OptimizerManager;}
    /*!
      \fn inline public vAlgorithm::SetImageDimensionsAndQuantification
      \param ap_ImageDimensionsAndQuantification
      \brief Set the Image Dimensions and Quantification Object 
    */
    inline void SetImageDimensionsAndQuantification(oImageDimensionsAndQuantification* ap_ImageDimensionsAndQuantification)
           {mp_ID = ap_ImageDimensionsAndQuantification;}
    /*!
      \fn inline public vAlgorithm::SetImageSpace
      \param ap_ImageSpace
      \brief Set the Image Space Object 
    */
    inline void SetImageSpace(oImageSpace* ap_ImageSpace)
           {mp_ImageSpace = ap_ImageSpace;}
    /*!
      \fn inline public vAlgorithm::SetProjectorManager
      \param ap_ProjectorManager
      \brief Set the Projector Manager Object 
    */
    inline void SetProjectorManager(oProjectorManager* ap_ProjectorManager)
           {mp_ProjectorManager = ap_ProjectorManager;}
    /*!
      \fn inline public vAlgorithm::SetImageConvolverManager
      \param ap_ImageConvolverManager
      \brief Set the Image Convolver Manager Object 
    */
    inline void SetImageConvolverManager(oImageConvolverManager* ap_ImageConvolverManager)
           {mp_ImageConvolverManager = ap_ImageConvolverManager;}
    /*!
      \fn inline public vAlgorithm::SetImageProcessingManager
      \param ap_ImageProcessingManager
      \brief Set the Image Processing Manager Object 
    */
    inline void SetImageProcessingManager(oImageProcessingManager* ap_ImageProcessingManager)
           {mp_ImageProcessingManager = ap_ImageProcessingManager;}
    /*!
      \fn inline public vAlgorithm::SetDynamicModelManager
      \param ap_DynamicModelManager
      \brief Set the Dynamic Model Manager Object 
    */
    inline void SetDynamicModelManager(oDynamicModelManager* ap_DynamicModelManager)
           {mp_DynamicModelManager = ap_DynamicModelManager;}
    /*!
      \fn inline public vAlgorithm::SetDeformationManager
      \param ap_DeformationManager
      \brief Set the Deformation Manager Object 
    */
    inline void SetDeformationManager(oDeformationManager* ap_DeformationManager)
           {mp_DeformationManager = ap_DeformationManager;}
    /*!
      \fn inline public vAlgorithm::SetDataFile
      \param a2p_DataFile
      \brief Set the list of DataFile
    */
    inline void SetDataFile(vDataFile** a2p_DataFile)
           {m2p_DataFile = a2p_DataFile;}
    /*!
      \fn inline public vAlgorithm::SetGPUflag
      \param a_flagGPU
      \brief Set the GPU flag
    */
    inline void SetGPUflag(bool a_flagGPU)
           {m_flagGPU = a_flagGPU;}
    /*!
      \fn inline public vAlgorithm::SetVerbose
      \param a_verboseLevel
      \brief Set Verbosity
    */
    inline void SetVerbose(int a_verboseLevel)
           {m_verbose = a_verboseLevel;}
    /*!
      \fn inline public vAlgorithm::SetNbBeds
      \param a_nbBeds
      \brief Set number of beds (bed positions)
    */
    inline void SetNbBeds(int a_nbBeds)
           {m_nbBeds = a_nbBeds;}
    /*!
      \fn inline public vAlgorithm::SetPathInitImage
      \param a_pathToInitialImage
      \brief Set path to an initial image
    */
    inline void SetPathInitImage(string a_pathToInitialImage)
           {m_pathToInitialImg = a_pathToInitialImage;}
    /*!
      \fn      inline public vAlgorithm::SetPathToAttenuationImage
      \param   string a_pathToAttenuationImage
      \brief   This function is used to set the path to the attenuation image.
    */
    inline void SetPathToAttenuationImage(string a_pathToAttenuationImage)
           {m_pathToAtnImg = a_pathToAttenuationImage;}
    /*!
      \fn inline public vAlgorithm::SetPathToSensitivityImage
      \param a_pathToSensitivityImage
      \brief Set path to the sensitivity image
    */
    inline void SetPathToSensitivityImage(string a_pathToSensitivityImage)
           {m_pathToSensitivityImg = a_pathToSensitivityImage;}
    /*!
      \fn inline public vAlgorithm::SetPathToMultiModalImage
      \param a_pathToMultiModalImage
      \brief Set path to multimodal images
    */
    inline void SetPathToMultiModalImage(vector<string> a_pathToMultiModalImage)
           {m_pathToMultiModalImg = a_pathToMultiModalImage;}
    /*!
      \fn inline public vAlgorithm::SetPathToMaskImage
      \param a_pathToMaskImage
      \brief Set path to a mask image
    */
    inline void SetPathToMaskImage(string a_pathToMaskImage)
           {m_pathToMaskImg = a_pathToMaskImage;}
    /*!
      \fn public vAlgorithm::SetNbIterationsAndSubsets
      \param a_nbIterationsSubsets
      \brief Set the number of iterations and subsets
      \details The provided string is a list of iteration:subset pairs 
               separated by commas.
      \return 0 if success, positive value otherwise
    */
    int SetNbIterationsAndSubsets(const string& a_nbIterationsSubsets);
    /*!
      \fn public vAlgorithm::SetOutputIterations
      \param a_outputIterations
      \brief Set the selected output iterations
      \details The provided string is a list of couple a:b separated by commas. \n
              It means that we save one iteration over a until b is reached. \n
               "b" must be incrementing for each successive couples.
              If the list is empty, we save all iterations by default. \n
              If the list is equal to "-1", then we save only the last iteration.
      \return 0 if success, positive value otherwise
    */
    int SetOutputIterations(const string& a_outputIterations);
    /*!
      \fn public virtual vAlgorithm::InitSpecificOptions
      \param a_specificOptions
      \return 0
    */
    virtual int InitSpecificOptions(string a_specificOptions);
    /*!
      \fn public virtual vAlgorithm::ShowHelpSpecific
      \brief Show help for the child algorithm
    */
    virtual void ShowHelpSpecific(){}


  // -------------------------------------------------------------------
  // Protected member functions
  protected:
    /*!
      \fn vAlgorithm::StepBeforeIterationLoop
      \brief This function is called at the beginning of the RunCPU function.
      \return 0 if success, positive value otherwise.
    */
    virtual int StepBeforeIterationLoop();
    /*!
      \fn vAlgorithm::StepAfterIterationLoop
      \brief This function is called at the end of the RunCPU function.
      \return 0 if success, positive value otherwise.
    */
    virtual int StepAfterIterationLoop();
    /*!
      \fn vAlgorithm::StepBeforeSubsetLoop
      \param a_iteration : iteration index
      \brief This function is called before starting the data subset loop.
      \return 0 if success, positive value otherwise.
    */
    virtual int StepBeforeSubsetLoop(int a_iteration);
    /*!
      \fn vAlgorithm::StepAfterSubsetLoop
      \param a_iteration : iteration index
      \brief This function is called after finishing the data subset loop.
      \return 0 if success, positive value otherwise.
    */
    virtual int StepAfterSubsetLoop(int a_iteration);
    /*!
      \fn vAlgorithm::StepPreProcessInsideSubsetLoop
      \param a_iteration : iteration index
      \param a_subset : subset index
      \brief This function is called right after starting the data subset loop.
      \return 0 if success, positive value otherwise.
    */
    virtual int StepPreProcessInsideSubsetLoop(int a_iteration, int a_subset);
    /*!
      \fn vAlgorithm::StepPreProcessInsideSubsetLoop
      \param a_iteration : iteration index
      \param a_subset : subset index
      \brief This function is called right after starting the subset loop. \n
      \return 0 if success, positive value otherwise.
    */
    virtual int StepPostProcessInsideSubsetLoop(int a_iteration, int a_subset);
    /*!
      \fn vAlgorithm::StepInnerLoopInsideSubsetLoop
      \param a_iteration : iteration index
      \param a_subset : subset index
      \param a_bed : bed position
      \brief This function is called inside the subset loop. \n
             It contains the core operations of the algorithm and must be implemented by child classes.
      \return 0 if success, positive value otherwise.
    */
    virtual int StepInnerLoopInsideSubsetLoop(int a_iteration, int a_subset, int a_bed) = 0;
    /*!
      \fn vAlgorithm::StepImageOutput
      \brief This function deals with everything about saving output images from the reconstruction
    */
    virtual int StepImageOutput(int a_iteration, int a_subset = -1) = 0;


  // -------------------------------------------------------------------
  // Data members
  protected:
    int m_nbIterations;                                 /*!< Number of iterations (default=1)*/
    int* mp_nbSubsets;                                  /*!< Number of subsets (default=1)*/
    bool* mp_outputIterations;                          /*!< A boolean for each iteration saying if we save it or not */
    int m_verbose;                                      /*!< Verbosity (default=-1)*/
    bool m_flagGPU;                                     /*!< Do we use GPU or not (default=false) */
    oImageDimensionsAndQuantification* mp_ID;           /*!< Pointer to the oImageDimensionsAndQuantification object */
    vDataFile** m2p_DataFile;                           /*!< Pointer to the array of vDataFile object */
    oProjectorManager* mp_ProjectorManager;             /*!< Pointer to the Projector Manager object */
    oOptimizerManager* mp_OptimizerManager;             /*!< Pointer to the Optimizer Manager object */
    oDeformationManager* mp_DeformationManager;         /*!< Pointer to the Deformation Manager object */
    oDynamicModelManager* mp_DynamicModelManager;       /*!< Pointer to the Dynamic Model Manager object */
    oImageSpace* mp_ImageSpace;                         /*!< Pointer to the Image Space object */
    oImageConvolverManager* mp_ImageConvolverManager;   /*!< Pointer to the Image Convolver Manager object */
    oImageProcessingManager* mp_ImageProcessingManager; /*!< Pointer to the Image Processing Manager object */
    int m_nbBeds;                                       /*!< number of bed FOVs (1 datafile by bed) (default=-1) */
    string m_pathToInitialImg;                          /*!< String containing the path to an initialization image */
    string m_pathToAtnImg;                              /*!< String of the path to the attenuation images */
    string m_pathToSensitivityImg;                      /*!< String containing the path to a sensitivity image */
    vector<string> m_pathToMultiModalImg;               /*!< String vector containing paths to multimodal images */
    string m_pathToMaskImg;                             /*!< String containing the path to a mask image */
    bool m_saveSensitivityHistoFlag;                    /*!< Flag specifying that the sensitivity image has to be saved for each subset/iteration in histogram mode */
    bool m_saveImageAfterSubsets;                       /*!< Flag specifying that the image has to be saved after each subset */
    bool m_saveDynamicBasisCoefficients;                /*!< Flag specifying that the dynamic basis coefficient images will be saved */


};

#endif














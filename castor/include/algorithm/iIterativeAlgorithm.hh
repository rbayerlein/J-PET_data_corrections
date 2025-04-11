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
  \brief    Declaration of class oIterativeAlgorithm
*/

#ifndef IITERATIVEALGORITHM_HH
#define IITERATIVEALGORITHM_HH 1

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
#include "vAlgorithm.hh"



/*!
  \class  oIterativeAlgorithm
  \brief  This is the main class for iterative optimization reconstruction algorithms. \n
*/
class iIterativeAlgorithm : public vAlgorithm
{
  // Constructor & Destructor
  public:
    iIterativeAlgorithm();
    ~iIterativeAlgorithm();

  // -------------------------------------------------------------------
  // Public member functions
  public:
    int InitSpecificOptions(string a_specificOptions);

  // -------------------------------------------------------------------
  // Protected member functions
  protected:
    /*!
      \fn iIterativeAlgorithm::StepBeforeIterationLoop
      \brief This function is called at the beginning of the RunCPU function.
      \details Initialization and memory allocation for the imageSpace and
               some managers
      \return 0 if success, positive value otherwise.
    */
    int StepBeforeIterationLoop();
    /*!
      \fn iIterativeAlgorithm::StepAfterIterationLoop
      \brief This function is called at the end of the RunCPU function.
      \details Free Memory for the imageSpace and
               some managers
      \return 0 if success, positive value otherwise.
    */
    int StepAfterIterationLoop();
    /*!
      \fn iIterativeAlgorithm::StepBeforeSubsetLoop
      \param a_iteration : iteration index
      \brief This function is called before starting the subset loop.
      \return 0 if success, positive value otherwise.
    */
    int StepBeforeSubsetLoop(int a_iteration);
    /*!
      \fn iIterativeAlgorithm::StepAfterSubsetLoop
      \param a_iteration : iteration index
      \brief This function is called after finishing the subset loop.
      \details Clean the main images from never visited voxels \n
               Apply post-convolution/post-processing if needed \n
               Write output images on disk as requested by the user
      \todo manage output of coeff  (save parametric images using intrinsic basis functions)
      \return 0 if success, positive value otherwise.
    */
    int StepAfterSubsetLoop(int a_iteration);
    /*!
      \fn iIterativeAlgorithm::StepPreProcessInsideSubsetLoop
      \param a_iteration : iteration index
      \param a_subset : subset index
      \brief This function is called right after starting the subset loop.
             Apply any kind of processing on the forward image before projections
      \details Copy current main image into forward image matrix \n
               Reinitialize backward image and 4D gating indices \n
               Apply image processing/convolution on the forward image matrix
               (image to be projeted)
      \todo    mp_DeformationManager->ApplyDeformationsToForwardImage() : Do we keep this function ?
      \return 0 if success, positive value otherwise.
    */
    int StepPreProcessInsideSubsetLoop(int a_iteration, int a_subset);
    /*!
      \fn iIterativeAlgorithm::StepPreProcessInsideSubsetLoop
      \param a_iteration : iteration index
      \param a_subset : subset index
      \brief This function is called right after starting the subset loop. \n
             Apply any kind of image processing on the backward image and
             main image after backprojections
      \details Synchronize parallel results \n
               Apply image deformation/processing/convolution on the backward image \n
               Update the image using the optimizer functions \n
               Apply dynamic model/image processing/convolution on the main image
      \return 0 if success, positive value otherwise.
    */
    int StepPostProcessInsideSubsetLoop(int a_iteration, int a_subset);
    /*!
      \fn iIterativeAlgorithm::StepInnerLoopInsideSubsetLoop
      \param a_iteration : iteration index
      \param a_subset : subset index
      \param a_bed : bed position
      \brief This function is called inside the subset loop and manages the main loop over the events \n
             The loop over the events is multithreaded, and involves a thread lock in the case of image-based deformation
      \details Step 1: Get the current event for that thread index \n
               Step 2: Call the dynamic switch function that updates the current frame and gate numbers, and also detects involuntary patient motion \n
                       Perform image-based deformation if needed (thread lock is required) \n
               Step 3: Compute the projection line \n
               Step 4: Optimize in the data space (forward-proj, update, backward-proj)
      \todo  Check implementation of thread-locked operation on different architectures
      \return 0 if success, positive value otherwise.
    */
    int StepInnerLoopInsideSubsetLoop(int a_iteration, int a_subset, int a_bed);
    /*!
      \fn iIterativeAlgorithm::StepImageOutput
      \brief This function deals with everything about saving output images from the reconstruction
      \details It applies post-smoothing, post-processing, output masking, and save dynamic basis coefficients if asked for.
    */
    int StepImageOutput(int a_iteration, int a_subset = -1);

  // -------------------------------------------------------------------
  // Data members
  protected:
};

#endif














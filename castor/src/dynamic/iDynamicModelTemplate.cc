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
  \ingroup  dynamic
  \brief    Implementation of class iDynamicModelTemplate
*/

#include "iDynamicModelTemplate.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn iDynamicModelTemplate
  \brief Constructor of iDynamicModelTemplate. Simply set all data members to default values.
*/
iDynamicModelTemplate::iDynamicModelTemplate() : vDynamicModel() 
{
  // --- Parameters inherited from vDynamicModel class --- //
  
  m_nbTimeBF = 1; // Number of basis functions in the model
  m_nbModelParam = 1; // Number of model parameters
  
  // Image matrix containing the parametric images.
  // 2 pointers:
  // 1: Parametric image related to the dynamic model basis functions.
  // 2: 3D voxels
  // Should be allocated in InitializeSpecific() function
  m2p_parametricImages = NULL;

  // Vector containing the Model temporal basis functions
  // 2 pointers: 
  // 1: index of the temporal function 
  // 2: coefficient of the functions for each time points (frame) of the dynamic acquisition
  // Dimensions : [m_nbModelParam][voxels] 
  // Should be allocated in InitializeSpecific() function
  m2p_nestedModelTimeBasisFunctions = NULL;

  // Image matrix to gather the parametric image before writing on disk
  // By default it will point directly to the parametric m2p_parametricImages.
  // They are allocated if post-processing are enabled before writing the image (i.e FOV masking)
  // 2 pointers:  \n
  // 1: Parametric image related to the dynamic model basis functions.
  // 2: 3D voxels 
  m2p_outputParImages = NULL; 

  // Image matrix containing the voxels which the model cannot fit
  // 1 pointer: 3D voxels 
  mp_blackListedvoxelsImage = NULL; 
  
  // Flag indicating if the blacklisted voxels mask image should be written on disk
  m_saveBlacklistedImageMaskFlag = false;  
  
  // Input image containing a mask defining in which voxels the model must be applied (1) or not (0). Default: all voxels to 1
  // 1 pointer: 3D voxels 
  mp_maskModel = NULL;
  
  // Number of voxels in mask
  m_nbVoxelsMask = -1;
    
  // Path to a configuration file. To be used for initialization in ReadAndCheckConfigurationFileSpecific() and InitializeSpecific() functions
  m_fileOptions = ""; 
  
  // String containing a list of options. To be used for initialization in ReadAndCheckOptionsList() and InitializeSpecific() functions
  m_listOptions = ""; 
  
  // Boolean indicating whether the parameters were checked or not
  m_checked = false;                
  
  // Boolean indicating whether the manager was initialized or not
  m_initialized = false;            
  
  // Flag indicating if parametric images should be written on disk (default=true)
  // Could be disabled with the keyword 'Save_parametric_images: 0' in a configuration file (managed by vDynamicModel)
  m_saveParImageFlag = true;
  
  // If true, the EstimateImageWithModel() functions is not called, so the reconstructed images are not estimated from the parametric images, 
  // The class will only estimate parametric images from each current estimation of the images (default=false)
  // Could be disabled with the keyword 'No_image_update: 1' in a configuration file (managed by vDynamicModel)
  m_noImageUpdateFlag = false; 
                                           
  // If true, the EstimateModelParameters() functions is not called, so the parameters are not estimated from the serie of dynamic images (default=false)
  // Could be disabled with the keyword 'No_parameters_update: 1' in a configuration file (managed by vDynamicModel)
  m_noParametersUpdateFlag = false; /*!<  */


  // Number of iterations after which the reconstructed images are estimated from the parametric images. (default or negative value =0)
  // Could be modified with the keyword 'Number of iterations before image update: xxx' in a configuration file (managed by vDynamicModel)
  m_startIteUpdateFlag = 0;
  
  //Verbosity
  m_verbose = -1;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ~iDynamicModelTemplate
  \brief Destructor of iDynamicModelTemplate
*/
iDynamicModelTemplate::~iDynamicModelTemplate() 
{
  if(m_initialized)
  {
    // Free variables
  }
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ShowHelpModelSpecific
  \brief Print out specific help about the implementation of the model 
         and its initialization
*/
void iDynamicModelTemplate::ShowHelpModelSpecific()
{
  // ===================================================================
  // Here, display some help and guidance to how to use this dynamic model and what it does
  // ===================================================================
  cout << "This class is a template class dedicated to add your own dynamic model." << endl;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ReadAndCheckConfigurationFileSpecific
  \brief This function is used to read options from a configuration file.
  \return 0 if success, other value otherwise.
*/
int iDynamicModelTemplate::ReadAndCheckConfigurationFileSpecific()
{
  if(m_verbose >=3) Cout("iDynamicModelTemplate::ReadAndCheckConfigurationFileSpecific ..."<< endl); 
  
  // ===================================================================
  // Implement here the reading of any options specific to this dynamic model 
  // (i.e : parameters or path to deformation files), through a configuration file
  // The ReadDataASCIIFile() functions could be helpful to recover data from a file
  // (check other dynamicModels for examples)
  // ===================================================================
    
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ReadAndCheckOptionsList
  \param a_optionsList : a list of parameters separated by commas
  \brief This function is used to read parameters from a string.
  \return 0 if success, other value otherwise.
*/
int iDynamicModelTemplate::ReadAndCheckOptionsList(string a_listOptions)
{
  // ===================================================================
  // Implement here the reading of any options specific to this dynamic model,
  // through a list of options separated by commas
  // The ReadStringOption() function could be helpful to parse the list of parameters in an array
  // ===================================================================
  
  // Normal end
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn CheckSpecificParameters
  \brief This function is used to check whether all member variables
         have been correctly initialized or not.
  \return 0 if success, positive value otherwise.
*/
int iDynamicModelTemplate::CheckSpecificParameters()
{
  // ===================================================================
  // Implement here checks over parameters which should be read using either
  // ReadAndCheckConfigurationFile() and ReadAndCheckOptionsList() functions
  // ===================================================================
    
  if(m_verbose >=3) Cout("iDynamicModelTemplate::CheckSpecificParameters ..."<< endl); 
  
  // Normal end
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn InitializeSpecific
  \brief This function is used to initialize the model parametric images and basis functions
  \return 0 if success, other value otherwise.
*/
int iDynamicModelTemplate::InitializeSpecific()
{
  if(m_verbose >=3) Cout("iDynamicModelTemplate::InitializeSpecific ..."<< endl); 


  // ===================================================================
  // Implement here the allocation/initialization of whatever member
  // variables specifically used by this dynamic model
  // ===================================================================
  
  
  // Forbid initialization without check
  if (!m_checked)
  {
    Cerr("***** oDynamicModelManager::InitializeSpecific() -> Must call CheckParameters functions before Initialize() !" << endl);
    return 1;
  }
  
  // Normal end
  m_initialized = true;
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn EstimateModelParameters
  \param ap_ImageS : pointer to the ImageSpace
  \param a_ite : index of the actual iteration (not used)
  \param a_sset : index of the actual subset (not used)
  \brief Estimate parametric images
  \return 0 if success, other value otherwise.
*/
int iDynamicModelTemplate::EstimateModelParameters(oImageSpace* ap_ImageS, int a_ite, int a_sset) 
{
  if(m_verbose >=3) Cout("iDynamicModelTemplate::EstimateModelParameters ..." <<endl);

  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** iDynamicModelTemplate::EstimateModelParameters() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif
  
  // ===================================================================
  // This function can be used to implement the estimation any parameters of the model
  // and parametric images 
  // This function is called dring the iterative image reconstruction 
  // during the image update which occurs at the end of the iteration/subset
  //
  // 
  // The main image matrices are stored in the ImageSpace object, passed
  // in argument.
  // The main image contains 4 dimensions :
  // ap_ImageS->m4p_image[fr][rg][cg][v]
  // fr = time frames
  // rg = respiratory gates
  // cg = cardiac gates
  //  v = actual voxel of the 3D volume
  
  /* IMAGE DIMENSIONS :
  *  For code efficiency and readability, the spatial index of a voxel is a cumulative 1D index. That is to say, given a voxel [indexX,indexY,indexZ],
  *  its cumulative 1D index is computed by 'index = indexZ*nbVoxXY + indexY*nbVoxX + indexX'.
  *
  *  The image dimensions can be recovered from the mp_ID class
  *  Total number of voxels         : mp_ID->GetNbVoxXYZ()
  *  Number of voxels in a slice    : mp_ID->GetNbVoxXY()
  *  Number of voxels on the X-axis : mp_ID->GetNbVoxX()
  */
  
  // Any error should return a value >0.
    
  // ===================================================================
  
  return 0;
}
  
  
  
// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn EstimateImageWithModel
  \param ap_ImageS : pointer to the ImageSpace
  \param a_ite : index of the actual iteration (not used)
  \param a_sset : index of the actual subset (not used)
  \brief Estimate image using the model parametric images and basis functions
  \return 0 if success, other value otherwise.
*/
int iDynamicModelTemplate::EstimateImageWithModel(oImageSpace* ap_ImageS, int a_ite, int a_sset) 
{
  if(m_verbose >= 3) Cout("iDynamicModelTemplate::EstimateImageWithModel ... " <<endl);
  
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** iDynamicModelTemplate::EstimateImageWithModel() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif


  // ===================================================================
  // This function can be used to generate the serie of dynamic images 
  // from the model, after estimation of the model parameters in EstimateModelParameters()
  // It is called right after EstimateModelParameters()
  //
  // The main image matrices are stored in the ImageSpace object, passed
  // in argument.
  // The main image contains 4 dimensions :
  // ap_ImageS->m4p_image[fr][rg][cg][v]
  // fr = time frames
  // rg = respiratory gates
  // cg = cardiac gates
  //  v = actual voxel of the 3D volume
  
  /* IMAGE DIMENSIONS :
  * For code efficiency and readability, the spatial index of a voxel is a cumulative 1D index. That is to say, given a voxel [indexX,indexY,indexZ],
  * its cumulative 1D index is computed by 'index = indexZ*nbVoxXY + indexY*nbVoxX + indexX'.
  *
  * The image dimensions can be recovered from the mp_ID class
  * Total number of voxels         : mp_ID->GetNbVoxXYZ()
  * Number of voxels in a slice    : mp_ID->GetNbVoxXY()
  * Number of voxels on the X-axis : mp_ID->GetNbVoxX()
  */
  
  // Any error should return a value >0.
    
  // ===================================================================
        
  return 0;
}

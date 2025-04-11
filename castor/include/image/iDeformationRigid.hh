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
  \ingroup  image
  \brief    Declaration of class iDeformationRigid
*/

#ifndef IDEFORMATIONRIGID_HH
#define IDEFORMATIONRIGID_HH 1

#include "oImageSpace.hh"
#include "sAddonManager.hh"
#include "vDeformation.hh"
#include "oMatrix.hh"

/*!
  \class    iDeformationRigid
  \brief    This class performs rigid transformation based on a trilinear interpolation
            It requires ASCII parameter files containing x,y,z transformation vectors
  \details  TODO
*/
class iDeformationRigid : public vDeformation
{
  // -----------------------------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \fn      iDeformationRigid::iDeformationRigid
      \brief   Constructor of iDeformationRigid. Simply set all data members to default values.
    */
    iDeformationRigid();
    /*!
      \fn      iDeformationRigid::~iDeformationRigid
      \brief   Destructor of iDeformationRigid. Free memory from all allocated tabs.
    */
    ~iDeformationRigid();


  // -----------------------------------------------------------------------------------------
  // Public member functions
  public:
    // Function for automatic insertion (put the class name as the parameters and do not add semi-colon at the end of the line)
    FUNCTION_DEFORMATION(iDeformationRigid)
    /*!
      \fn      iDeformationRigid::ReadAndCheckConfigurationFile
      \param   a_configurationFile
      \brief   This function is used to read options from a configuration file.
      \return  0 if success, other value otherwise.
    */
    int ReadAndCheckConfigurationFile(const string& a_fileOptions);
    /*!
      \fn      iDeformationRigid::ReadAndCheckOptionsList
      \param   a_optionsList
      \brief   This function is used to read options from a list of options.
               Throw error by defaut for this method, as a file has to be used for initialization
      \return  0 if success, other value otherwise.
    */
    int ReadAndCheckOptionsList(const string& a_listOptions);
    /*!
      \fn      iDeformationRigid::CheckSpecificParameters
      \brief   This function is used to check parameters after the latter
               have been all set using Set functions.
      \return  0 if success, positive value otherwise.
    */
    int CheckSpecificParameters();
    /*!
      \fn      iDeformationRigid::Initialize
      \brief   This function is used to initialize specific stuff in the deformation model.
      \details Allocate memory for image matrices, and initialize images containing transformation parameters
      \return  0 if success, other value otherwise.
    */
    int Initialize();
    /*!
      \fn      iDeformationRigid::ShowHelp
      \brief   This function is used to print out specific help about the deformation and its options.
    */
    void ShowHelp();
    /*!
      \fn      iDeformationRigid::ApplyDeformations
      \param   ap_inputImage : input image to deform
      \param   ap_outputImage : image in which the output of the deformation should be recovered
      \param   a_direction : a direction for the deformation to perform (forward or backward)
      \param   a_defIdx : index of the deformation
      \brief   This function prepares the deformation to perform 
      \details 1. Selects the right deformation parameters file according to the direction and deformation index argument \n
               2. Copy the image to deform in the buffer image of this class \n
               3. Call the deformation function (TransformImage) \n
               4. Copy back the output of the deformation image to the output image matrice passed in argument
      \return  0 if success, other value otherwise.
    */
    int ApplyDeformations(FLTNB* ap_inputImage, FLTNB* ap_outputImage, int a_direction, int a_defIdx);


  // -----------------------------------------------------------------------------------------
  // Private member functions
  private:
    /*!
      \fn      iDeformationRigid::TransformImage
      \param   a_direction : direction of deformation
      \param   a_defIdx : index of deformation
      \param   ap_inputImage : image on which the deformation should be performed
      \param   ap_outputImage : matrice which recovers the image after deformation
      \brief   This function performs a rigid transformation using x,y,z displacement vectors and trilinear interpolation \n
      \details 1. Load the X,Y,Z transformation vectors from the file passed in argument
               2. Loop on the voxels and compute the new voxels value using trilinear interpolation
      \todo    perhaps add options to load all the vector parameters in RAM during initialization, instead on reading files during reconstruction
      \todo    perhaps use padded image in order to avoid if statements
      \return  0 if success, other value otherwise.
    */
    int TransformImage(int a_direction, int a_defIdx, HPFLTNB *ap_inputImage, HPFLTNB *ap_outputImage);

    /*!
      \fn ComputeTransformationMatrices()
      \brief Initialize transformation matrices from parameters
      \return 0 if success, positive value otherwise.
    */
    int ComputeTransformationMatrices();

    /*!
      \fn SetRotationMatrix(oMatrix* apRotMtx, FLTNB a_ang1, FLTNB a_ang2, FLTNB a_ang3, string a_cvt)
      \param apRotMtx : Rotation oMatrix to initialize
      \param a_ang1 : 1st rotation angle (rad)
      \param a_ang2 : 2nd rotation angle (rad)
      \param a_ang3 : 3rd rotation angle (rad)
      \param a_cvt : rotation convention (must be any combinations of 'x', 'y', and 'z'. (Default: xyz) )
      \brief This function set the rotation matrix passed in parameter with the provided angles in radian
             and rotation convention
      \todo Check this function, add missing conventions
      \return 0 if success, other value otherwise.
    */
    int SetRotationMatrix(oMatrix* apRotMtx, FLTNB a_ang1, FLTNB a_ang2, FLTNB a_ang3, string a_cvt);

    //int Tlerp(HPFLTNB *ap_inputImage, HPFLTNB *ap_outputImage, uint32_t iov, uint32_t iiv, FLTNB dX, FLTNB dY, FLTNB dZ);

  // -----------------------------------------------------------------------------------------
  // Data members
  private:
    vector<string> 
      m_pathToFwdDeformationFiles; /*!< Containers of the location of forward deformation file.
                                        Nb of elements should be equal to the number of transformations 'm_nbTransformations' */
    vector<string> 
      m_pathToBwdDeformationFiles; /*!< Containers of the location of backward deformation file.
                                        Nb of elements should be equal to the number of transformations 'm_nbTransformations' */
    HPFLTNB* mp_OriginalArray;       /*!< Image matrice to get the input image for deformation */
    HPFLTNB* mp_OutputArray;         /*!< Image matrice to recover the output image after deformation */
    HPFLTNB* mp_tX;                /*!< Translation on the X-axis */
    HPFLTNB* mp_tY;                /*!< Translation on the Y-axis */
    HPFLTNB* mp_tZ;                /*!< Translation on the Z-axis */
    HPFLTNB* mp_rA;                /*!< First rotation */
    HPFLTNB* mp_rB;                /*!< Second rotation */
    HPFLTNB* mp_rC;                /*!< Third rotation */
    oMatrix** m2p_FTmtx;          /*!< Forward transformation matrix */
    oMatrix** m2p_BTmtx;          /*!< Backward transformation matrix */
    string m_rotConvention;       /*!< Convention for rotation. Default : XYZ */
    bool m_cmpTransfoFlag;        /*!< Flag defining if transformations from reference to each position must be computed (provided
                                       transformation parameters from one position to the next one, default) or not (already computed by the user)*/
};

// Class for automatic insertion (set here the visible image deformation's name, put the class name as the parameters and do not add semi-colon at the end of the line)
CLASS_DEFORMATION(deformationRigid,iDeformationRigid)

#endif

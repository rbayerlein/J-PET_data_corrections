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
  \brief    Implementation of class iDeformationRigid
*/

#include "iDeformationRigid.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn iDeformationRigid
  \brief Constructor of iDeformationRigid. Simply set all data members to default values.
*/
iDeformationRigid::iDeformationRigid() : vDeformation()
{
  mp_OriginalArray = NULL;
  mp_OutputArray = NULL;
  mp_tX = NULL;
  mp_tY = NULL;
  mp_tZ = NULL;
  mp_rA = NULL;
  mp_rB = NULL;
  mp_rC = NULL;
  m2p_FTmtx = NULL;
  m2p_BTmtx = NULL;
  m_rotConvention = "XYZ";
  m_cmpTransfoFlag = true;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ~iDeformationRigid
  \brief Destructor of iDeformationRigid. Free memory from all allocated tabs.
*/
iDeformationRigid::~iDeformationRigid() 
{
  if(m_initialized)
  {
    if(mp_OriginalArray) delete[] mp_OriginalArray;
    if(mp_OutputArray) delete[] mp_OutputArray;  
  }

  for(int t=0; t<m_nbTransformations; t++)
  {
    if (m2p_FTmtx[t]) delete m2p_FTmtx[t];
    if (m2p_BTmtx[t]) delete m2p_BTmtx[t];
  }
  
  if(m2p_FTmtx) delete[] m2p_FTmtx;
  if(m2p_BTmtx) delete[] m2p_BTmtx;
  
  
  if(mp_tX) delete[] mp_tX;
  if(mp_tY) delete[] mp_tY;
  if(mp_tZ) delete[] mp_tZ;
  if(mp_rA) delete[] mp_rA;
  if(mp_rB) delete[] mp_rB;
  if(mp_rC) delete[] mp_rC;

}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ShowHelp
  \brief This function is used to print out specific help about the deformation and its options.
*/
void iDeformationRigid::ShowHelp()
{
  // todo : really we should rather use param files in binary format here
  cout << "-- Rigid transformation algorithm based on trilinear interpolation --" << endl;
  cout << "Each transformation requires 3 translation (x,y,z) and 3 rotation parameters, provided by command-line option or in a text file" << endl;
  cout << "The transformation parameters can be either defined between: " << endl;
  cout << "- A position n and the following position n+1 " << endl;
  cout << "- A reference position and the actual position n (default) " << endl;
  cout << "  (Check the 'Transformation mode' keyword below for more information) " << endl;
  cout << endl;
  cout << "--------------------------------------------------------" << endl;
  cout << "Command-line options:" <<endl;
  cout << "Parameters for each transformation must be provided sequentially and separated by commas" << endl;
  cout << "Example for a dataset containing 3 transformations:" << endl;
  cout << "-rm deformationRigid:trX1,trY1,trZ1,rotX1,rotY1,rotZ1,trX2,trY2,trZ2,rotX2,rotY2,rotZ2,trX3,trY3,trZ3,rotX3,rotY3,rotZ3" << endl;
  cout << "--------------------------------------------------------" << endl;
  cout << "Text file:" <<endl;
  cout << "Parameters for each transformation must be provided after a 'Transformation parameters' key, on a separate line, and separated by commas" << endl;
  cout << "Translations must be provided in mm, for each axis." << endl;
  cout << "Rotations must be provided in degree. Use the Rotation convention optional parameter in order to set up the rotation orientation." << endl;
  cout << "Example for a dataset containing 3 transformations:" << endl;
  cout << "Transformation_parameters:" << endl;
  cout << "trX1,trY1,trZ1,rotA1,rotB1,rotC1" << endl;
  cout << "trX2,trY2,trZ2,rotA2,rotB2,rotC2" << endl;
  cout << "trX3,trY3,trZ3,rotA3,rotB3,rotC3" << endl;
  cout << "etc..." << endl;
  cout << endl;
  cout << "--------------------------------------------------------" << endl;
  cout << "Optional parameters (configuration file only):" << endl;
  cout << endl;
  cout << "Rotation_convention: Any combinations of 'x', 'y', and 'z'" << endl;
  cout << "Define the axe on which the rotation will be performed using rotA, rotB and rotC angle in 1 (default): XYZ" << endl;
  cout << endl;
  cout << "Transformation_mode: (1 or 0)" << endl;
  cout << "0 (default): Transformation between the reference position to the n position." << endl;
  cout << "1          : Transformation between the n-1 position to the n position." << endl;
  cout << "             They will be recomputed to represent transformations from a reference position to i  as required by CASToR" << endl;

}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ReadAndCheckConfigurationFile
  \param a_configurationFile
  \brief This function is used to read options from a configuration file.
  \return 0 if success, other value otherwise.
*/
int iDeformationRigid::ReadAndCheckConfigurationFile(const string& a_fileOptions)
{
  if (m_verbose >=VERBOSE_DETAIL) Cout("iDeformationRigid::ReadAndCheckConfigurationFile ..."<< endl); 
  
  // First, check the nb of transformations has correctly been initializeds
  if (m_nbTransformations<0)
  {
    Cerr("***** iDeformationRigid::ReadAndCheckConfigurationFile() -> Number of transformations in the deformation has not been initialized  !" << endl);
    return 1;
  }
  
  
  ifstream in_file(a_fileOptions.c_str(), ios::in);
  
  if(in_file)
  {
    // Local array to recover rigid transformation parameters
    HPFLTNB* p_coeffs = new HPFLTNB[m_nbTransformations*6];
    
    // Read 6 parameters for each transformations, on each line
    if (ReadDataASCIIFile(a_fileOptions, "Transformation_parameters", p_coeffs, 6, m_nbTransformations, KEYWORD_MANDATORY))
    {
      Cerr("***** iDeformationRigid::ReadAndCheckConfigurationFile -> A problem occurred while trying to recover transformation parameters from file !" << endl);
      return 1;
    }

    // Read 6 parameters for each transformations, on each line
    if (ReadDataASCIIFile(a_fileOptions, "Transformation_mode", &m_cmpTransfoFlag, 1, KEYWORD_OPTIONAL) == KEYWORD_OPTIONAL_ERROR)
    {
      Cerr("***** iDeformationRigid::ReadAndCheckConfigurationFile -> A problem occurred while trying to read 'compute_transformations' flag (must be 1 (true) or 0 (false))!" << endl);
      return 1;
    }

    // Read 6 parameters for each transformations, on each line
    if (ReadDataASCIIFile(a_fileOptions, "Rotation_convention", &m_rotConvention, 1, KEYWORD_OPTIONAL) == KEYWORD_OPTIONAL_ERROR)
    {
      Cerr("***** iDeformationRigid::ReadAndCheckConfigurationFile -> A problem occurred while trying to read rotation convention (string) !" << endl);
      return 1;
    }
    
    // Allocate memory for each parameter
    mp_tX = new HPFLTNB[m_nbTransformations];
    mp_tY = new HPFLTNB[m_nbTransformations];
    mp_tZ = new HPFLTNB[m_nbTransformations];
    mp_rA = new HPFLTNB[m_nbTransformations];
    mp_rB = new HPFLTNB[m_nbTransformations];
    mp_rC = new HPFLTNB[m_nbTransformations];
    
    // Recover parameters
    for(int t=0 ; t<m_nbTransformations ; t++)
    {
      mp_tX[t] = p_coeffs[t*6];
      mp_tY[t] = p_coeffs[t*6+1];
      mp_tZ[t] = p_coeffs[t*6+2];
      mp_rA[t] = p_coeffs[t*6+3];
      mp_rB[t] = p_coeffs[t*6+4];
      mp_rC[t] = p_coeffs[t*6+5];
    }

    if(ComputeTransformationMatrices())
    {
      Cerr("***** iDeformationRigid::ReadAndCheckOptionsList() -> An error occurred while building Transformation matrices !" << endl);
      return 1;
    }
  
    // Delete local objects
    delete[] p_coeffs;
  }
  else
  {
    Cerr("***** iDeformationRigid::ReadAndCheckConfigurationFile -> Error while trying to read configuration file : " << a_fileOptions << endl);
    return 1;
  }
  
  return 0;  
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ReadAndCheckOptionsList
  \param a_optionsList
  \brief This function is used to read options from a list of options.
         Throw error by defaut for this method, as a file has to be used for initialization
  \return 0 if success, other value otherwise.
*/
int iDeformationRigid::ReadAndCheckOptionsList(const string& a_listOptions)
{
  if(m_verbose >=VERBOSE_DETAIL) Cout("iDeformationRigid::ReadAndCheckOptionsList ..."<< endl); 

  // First, check the nb of transformations has correctly been initializeds
  if (m_nbTransformations<0)
  {
    Cerr("***** iDeformationRigid::ReadAndCheckOptionsList() -> Number of transformations in the deformation has not been initialized  !" << endl);
    return 1;
  }
  
  // Local array to recover rigid transformation parameters
  HPFLTNB* p_coeffs = new HPFLTNB[m_nbTransformations*6];
  
  if (ReadStringOption(a_listOptions,
                       p_coeffs,
                       m_nbTransformations*6,
                       ",",
                       "Rigid deformation configuration"))
  {
    Cerr("***** iDeformationRigid::ReadAndCheckOptionsList() -> Failed to correctly read the list of parameters in command-line options !" << endl);
    Cerr("*****                                                 "<<m_nbTransformations*6<<" parameters were expected (6 parameters for each transformation)" << endl);
    return 1;
  }
  
  // Allocate memory for each parameter
  
  mp_tX = new HPFLTNB[m_nbTransformations];
  mp_tY = new HPFLTNB[m_nbTransformations];
  mp_tZ = new HPFLTNB[m_nbTransformations];
  mp_rA = new HPFLTNB[m_nbTransformations];
  mp_rB = new HPFLTNB[m_nbTransformations];
  mp_rC = new HPFLTNB[m_nbTransformations];
  
  // Recover parameters
  for(int t=0 ; t<m_nbTransformations ; t++)
  {
    mp_tX[t] = p_coeffs[t*6];
    mp_tY[t] = p_coeffs[t*6+1];
    mp_tZ[t] = p_coeffs[t*6+2];
    mp_rA[t] = p_coeffs[t*6+3];
    mp_rB[t] = p_coeffs[t*6+4];
    mp_rC[t] = p_coeffs[t*6+5];
  }
  
  if(ComputeTransformationMatrices())
  {
    Cerr("***** iDeformationRigid::ReadAndCheckOptionsList() -> An error occurred while building Transformation matrices !" << endl);
    return 1;
  }
  
  delete[] p_coeffs;
  
  // Normal end
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ComputeTransformationMatrices
  \brief Initialize transformation matrices from parameters
  \return 0 if success, positive value otherwise.
*/
int iDeformationRigid::ComputeTransformationMatrices()
{
  // Check initialization has been done
  if(    mp_tX==NULL || mp_tY==NULL || mp_tZ==NULL ||
         mp_rA==NULL || mp_rB==NULL || mp_rC==NULL )
  {
    Cerr("***** iDeformationRigid::ComputeTransformationMatrices() -> Parameters not initialized  !" << endl);
    return 1;
  }
  
  // Allocate tmp matrices (translation, rotation, and backup)
  oMatrix *p_trs_mtx = new oMatrix(4,4); 
  oMatrix *p_rot_mtx = new oMatrix(4,4); 
  oMatrix *p_tmp_mtx = new oMatrix(4,4);
  
  // Allocate transformation matrices
  m2p_FTmtx = new oMatrix*[m_nbTransformations];
  m2p_BTmtx = new oMatrix*[m_nbTransformations];
  
  
  for(int t=0; t<m_nbTransformations; t++)
  {
    // Allocate matrices for this transformation
    m2p_FTmtx[t] = new oMatrix(4,4);
    m2p_BTmtx[t] = new oMatrix(4,4);
   
    FLTNB tx = mp_tX[t];
    FLTNB ty = mp_tY[t];
    FLTNB tz = mp_tZ[t];
    FLTNB x_rad = mp_rA[t] * M_PI/180.;
    FLTNB y_rad = mp_rB[t] * M_PI/180.;
    FLTNB z_rad = mp_rC[t] * M_PI/180.;
    
    // For both direction
    for(int d=0 ; d<2 ; d++)
    {
      if (d == 1)
      {
        tx = -tx;
        ty = -ty;
        tz = -tz;
        x_rad = -mp_rA[t] * M_PI/180.;
        y_rad = -mp_rB[t] * M_PI/180.;
        z_rad = -mp_rC[t] * M_PI/180.;
      }
      
      // Set translation matrix
      p_trs_mtx->SetMatriceElt(0,0, 1 );
      p_trs_mtx->SetMatriceElt(0,1, 0 );
      p_trs_mtx->SetMatriceElt(0,2, 0);
      p_trs_mtx->SetMatriceElt(0,3, -tx);
      p_trs_mtx->SetMatriceElt(1,0, 0 );
      p_trs_mtx->SetMatriceElt(1,1, 1 );
      p_trs_mtx->SetMatriceElt(1,2, 0);
      p_trs_mtx->SetMatriceElt(1,3, -ty);
      p_trs_mtx->SetMatriceElt(2,0, 0 );
      p_trs_mtx->SetMatriceElt(2,1, 0 );
      p_trs_mtx->SetMatriceElt(2,2, 1);
      p_trs_mtx->SetMatriceElt(2,3, -tz);
      p_trs_mtx->SetMatriceElt(3,0, 0 );
      p_trs_mtx->SetMatriceElt(3,1, 0 );
      p_trs_mtx->SetMatriceElt(3,2, 0);
      p_trs_mtx->SetMatriceElt(3,3, 1);
    
      // Set rotation matrix (ZYX)
      if(SetRotationMatrix(p_rot_mtx, x_rad, y_rad, z_rad, m_rotConvention))
      {
        Cerr("***** iDeformationRigid::ComputeTransformationMatrices() -> Error occurred when initializing rotation matrices  !" << endl);
        return 1;
      }

      // Compute forward or backward transformation matrix from t to t+1
      (d == 0) ? 
        p_trs_mtx->Multiplication(p_rot_mtx, m2p_FTmtx[t]):
        p_trs_mtx->Multiplication(p_rot_mtx, m2p_BTmtx[t]);
    }
    
    // Compute transformation matrices from and toward reference position
    
    if(m_cmpTransfoFlag && t>0)
    {
      // forward transformation matrix
      m2p_FTmtx[t]->Multiplication( m2p_FTmtx[t-1], p_tmp_mtx );
      
      for(int l=0 ; l<4 ; l++)
        for(int c=0 ; c<4 ; c++)
          m2p_FTmtx[t]->SetMatriceElt(l,c, p_tmp_mtx->GetMatriceElt(l,c) );
      
      // backward transformation matrix
      m2p_BTmtx[t]->Multiplication( m2p_BTmtx[t-1], p_tmp_mtx );
            
      for(int l=0 ; l<4 ; l++)
        for(int c=0 ; c<4 ; c++)
          m2p_BTmtx[t]->SetMatriceElt(l,c, p_tmp_mtx->GetMatriceElt(l,c) );
    }
  }

  delete p_trs_mtx; 
  delete p_rot_mtx; 
  delete p_tmp_mtx;
  
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn CheckSpecificParameters
  \brief This function is used to check parameters after the latter
         have been all set using Set functions.
  \return 0 if success, positive value otherwise.
*/
int iDeformationRigid::CheckSpecificParameters()
{
  if(m_verbose >= VERBOSE_DETAIL) Cout("iDeformationRigid::CheckSpecificParameters ..."<< endl); 
  
  // Check vector files
  if (!mp_tX || 
      !mp_tY || 
      !mp_tZ || 
      !mp_rA || 
      !mp_rB || 
      !mp_rC  )
  {
    Cerr("***** iDeformationRigid::CheckSpecificParameters() -> Transformation parameters not initialized !" << endl);
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
/*
  \fn Initialize
  \brief This function is used to initialize specific stuff in the deformation model.
  \details Allocate memory for image matrices, and initialize images containing transformation parameters
  \return 0 if success, other value otherwise.
*/
int iDeformationRigid::Initialize()
{
  if(m_verbose >=VERBOSE_NORMAL) Cout("iDeformationRigid::Initialize ..."<< endl); 
  
  // Forbid initialization without check
  if (!m_checked)
  {
    Cerr("***** iDeformationRigid::Initialize() -> Parameters should be checked before Initialize() !" << endl);
    return 1;
  }
  
  // Memory allocation for image matrices
  mp_OriginalArray = new HPFLTNB [mp_ID->GetNbVoxXYZ()];
  mp_OutputArray   = new HPFLTNB [mp_ID->GetNbVoxXYZ()];

  // Normal end
  m_initialized = true;
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ApplyDeformations
  \param ap_inputImage : input image to deform
  \param ap_outputImage : image in which the output of the deformation should be recovered
  \param a_direction : a direction for the deformation to perform (forward or backward)
  \param a_defIdx : index of the deformation
  \brief This function prepares the deformation to perform 
  \details 1. Selects the right deformation parameters file according to the direction and deformation index argument
           2. Copy the image to deform in the buffer image of this class
           3. Call the deformation function (TransformImage)
           4. Copy back the output of the deformation image to the output image matrice passed in argument
  \return 0 if success, other value otherwise.
*/
int iDeformationRigid::ApplyDeformations(FLTNB* ap_inputImage, FLTNB* ap_outputImage, int a_direction, int a_defIdx)
{
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** iDeformationRigid::ApplyDeformations() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif

  if(m_verbose >= VERBOSE_DETAIL) Cout("iDeformationRigid::ApplyDeformations #" << a_defIdx+1 << " ("
                                 << (string)((a_direction==FORWARD_DEFORMATION)?"forward":"backward") << "), with parameters "
                                 << mp_tX[a_defIdx] << ";" << mp_tY[a_defIdx] << ";" << mp_tZ[a_defIdx] << ";"
                                 << mp_rA[a_defIdx] << ";" << mp_rB[a_defIdx] << ";" << mp_rC[a_defIdx] << ";" << endl);

  for(int v=0; v<mp_ID->GetNbVoxXYZ(); v++)
  {
    mp_OriginalArray[v] = (HPFLTNB) ap_inputImage[v];
    mp_OutputArray[v] = 0.;
  }
  
  // Apply all required deformations
  TransformImage(a_direction, a_defIdx, mp_OriginalArray, mp_OutputArray); 

  for(int v=0; v<mp_ID->GetNbVoxXYZ(); v++)
  {
    if (mp_OutputArray[v] < 0) mp_OutputArray[v] = 0.;
    ap_outputImage[v] = (FLTNB) mp_OutputArray[v];
  }
  
  return 0;
}




// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn TransformImage
  \param a_vectorFileIdx : deformation vector file
  \param ap_inputImage : image on which the deformation should be performed
  \param ap_outputImage : matrice which recovers the image after deformation
  \brief This function performs a rigid transformation using x,y,z displacement vectors and trilinear interpolation
  \details 1. Load the X,Y,Z transformation vectors from the file passed in argument
           2. Loop on the voxels and compute the new voxels value using trilinear interpolation
  \todo : perhaps add options to load all the vector parameters in RAM during initialization, instead on reading files during reconstruction
  \todo : perhaps use padded image in order to avoid if statements
  \return 0 if success, other value otherwise.
*/
int iDeformationRigid::TransformImage(int a_direction, int a_defIdx, HPFLTNB *ap_inputImage, HPFLTNB *ap_outputImage)
{
  #ifdef CASTOR_VERBOSE
  if(m_verbose >=VERBOSE_DETAIL) Cout("iDeformationRigid::TransformImage() ... : "<< endl); 
  #endif

  // Set vectors to recover voxel locations before/after transformation
  oMatrix *iVect = new oMatrix(4,1);
  oMatrix *oVect = new oMatrix(4,1);
  
  // Recover image dimensions variables
  FLTNB sX = mp_ID->GetVoxSizeX();
  FLTNB sY = mp_ID->GetVoxSizeY();
  FLTNB sZ = mp_ID->GetVoxSizeZ();
  FLTNB dX = mp_ID->GetNbVoxX();
  FLTNB dY = mp_ID->GetNbVoxY();
  FLTNB dZ = mp_ID->GetNbVoxZ();
  FLTNB dXY = mp_ID->GetNbVoxXY();

    
  // Loop on voxels
  for(int v=0 ; v<mp_ID->GetNbVoxXYZ() ; v++)
  {
    // Get voxel index by axis
    uint32_t iZ = v/mp_ID->GetNbVoxXY();
    uint32_t iY = (v - iZ*mp_ID->GetNbVoxXY()) / mp_ID->GetNbVoxX();
    uint32_t iX =  v - iZ*mp_ID->GetNbVoxXY() - iY*mp_ID->GetNbVoxX();
    
    // Set input vector
    iVect->SetMatriceElt(0,0, iX*sX - dX*sX/2 + sX*0.5);
    iVect->SetMatriceElt(1,0, iY*sY - dY*sY/2 + sY*0.5);
    iVect->SetMatriceElt(2,0, iZ*sZ - dZ*sZ/2 + sZ*0.5);
    iVect->SetMatriceElt(3,0, 1);
    
    // Compute output vector with the correct related matrix
    if (a_direction==FORWARD_DEFORMATION)
      m2p_FTmtx[a_defIdx]->Multiplication(iVect, oVect);
    else if (a_direction==BACKWARD_DEFORMATION)
      m2p_BTmtx[a_defIdx]->Multiplication(iVect, oVect);
    else
    {
      Cerr("*****iDeformationRigid::ApplyDeformationToForwardImage()->Error, unknown direction type( must be 'forward' or 'backward' # " << endl);
      return 1;
    }
    
    // Get index of the voxel to be interpolated in the input image
    uint32_t ivZ = (oVect->GetMatriceElt(2,0) + dZ*sZ/2) / sZ;
    uint32_t ivY = (oVect->GetMatriceElt(1,0) + dY*sY/2) / sY;
    uint32_t ivX = (oVect->GetMatriceElt(0,0) + dX*sX/2) / sX;
    uint32_t iv =  ivZ*dXY + ivY*dX + ivX;
    
    // Compute difference between cartesian position after transformation and center of estimated vox position
    FLTNB dvX, dvY, dvZ;
    dvZ = (oVect->GetMatriceElt(2,0) - (ivZ*sZ - dZ*sZ/2 + sZ*0.5) ) / sZ;
    dvY = (oVect->GetMatriceElt(1,0) - (ivY*sY - dY*sY/2 + sY*0.5) ) / sY;
    dvX = (oVect->GetMatriceElt(0,0) - (ivX*sX - dX*sX/2 + sX*0.5) ) / sX;

    // Set constraints
    if (dvX>1) dvX =1;
    if (dvY>1) dvY =1;
    if (dvZ>1) dvZ =1;
    if (dvX<-1) dvX =-1;
    if (dvY<-1) dvY =-1;
    if (dvZ<-1) dvZ =-1;
    
    // Get Tlerp of destination voxel
    Tlerp(ap_inputImage, ap_outputImage, v, iv, dvX, dvY, dvZ);
  }  

  delete iVect;
  delete oVect;
  
  return 0;
}




/*
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
int iDeformationRigid::SetRotationMatrix(oMatrix* apRotMtx, FLTNB a_ang1, FLTNB a_ang2, FLTNB a_ang3, string a_cvt)
{
  oMatrix *r1 = new oMatrix(3,3);
  oMatrix *r2 = new oMatrix(3,3);
  oMatrix *r3 = new oMatrix(3,3);
  oMatrix *rtmp = new oMatrix(3,3);
  

  // Set up rotation matrix according to the convention
  string a = a_cvt.substr(0, 1);
  string b = a_cvt.substr(1, 1);
  string c = a_cvt.substr(2, 1);

  
  if     (a == "X" || a == "x")
    r1->SetXRotMtx(a_ang1);
  else if(a == "Y" || a == "y")
    r1->SetYRotMtx(a_ang1);
  else if(a == "Z" || a == "z")
    r1->SetZRotMtx(a_ang1);
  else
  {
    Cerr("*****iDeformationRigid::SetRotationMatrix()-> Symbol'"<< a << "' unknown for rotation convention! " << endl);
    return 1;
  }
  
  // Set up rotation matrix according to the convention
  if     (b == "X" || b == "x")
    r2->SetXRotMtx(a_ang2);
  else if(b == "Y" || b == "y")
    r2->SetYRotMtx(a_ang2);
  else if(b == "Z" || b == "z")
    r2->SetZRotMtx(a_ang2);
  else
  {
    Cerr("*****iDeformationRigid::SetRotationMatrix()-> Symbol'"<< b << "' unknown for rotation convention! " << endl);
    return 1;
  }

  // Set up rotation matrix according to the convention
  if     (c == "X" || c == "x")
    r3->SetXRotMtx(a_ang3);
  else if(c == "Y" || c == "y")
    r3->SetYRotMtx(a_ang3);
  else if(c == "Z" || c == "z")
    r3->SetZRotMtx(a_ang3);
  else
  {
    Cerr("*****iDeformationRigid::SetRotationMatrix()-> Symbol'"<< c << "' unknown for rotation convention! " << endl);
    return 1;
  }
  
  // r1*r2. Recover result in rtmp
  r1->Multiplication(r2, rtmp);

  // rtmp*r3. Recover result in r1
  rtmp->Multiplication(r3, r1);
  
  
  // Set up rotation matrix with r1
  apRotMtx->SetMatriceElt(0,0, r1->GetMatriceElt(0,0) );
  apRotMtx->SetMatriceElt(0,1, r1->GetMatriceElt(0,1) );
  apRotMtx->SetMatriceElt(0,2, r1->GetMatriceElt(0,2));
  apRotMtx->SetMatriceElt(0,3, 0);
  apRotMtx->SetMatriceElt(1,0, r1->GetMatriceElt(1,0) );
  apRotMtx->SetMatriceElt(1,1, r1->GetMatriceElt(1,1) );
  apRotMtx->SetMatriceElt(1,2, r1->GetMatriceElt(1,2));
  apRotMtx->SetMatriceElt(1,3, 0);
  apRotMtx->SetMatriceElt(2,0, r1->GetMatriceElt(2,0));
  apRotMtx->SetMatriceElt(2,1, r1->GetMatriceElt(2,1) );
  apRotMtx->SetMatriceElt(2,2, r1->GetMatriceElt(2,2));
  apRotMtx->SetMatriceElt(2,3, 0);
  apRotMtx->SetMatriceElt(3,0, 0 );
  apRotMtx->SetMatriceElt(3,1, 0 );
  apRotMtx->SetMatriceElt(3,2, 0);
  apRotMtx->SetMatriceElt(3,3, 1);
  
  /*
  HPFLTNB cx = cos(a_ang1);
  HPFLTNB cy = cos(a_ang2);
  HPFLTNB cz = cos(a_ang3);
  HPFLTNB sx = sin(a_ang1);
  HPFLTNB sy = sin(a_ang2);
  HPFLTNB sz = sin(a_ang3);
  
  if(a_cvt == "XYZ" || a_cvt == "xyz")
  {
    apRotMtx->SetMatriceElt(0,0, cy*cz );
    apRotMtx->SetMatriceElt(0,1, -cy*sz );
    apRotMtx->SetMatriceElt(0,2, sy);
    apRotMtx->SetMatriceElt(0,3, 0);
    apRotMtx->SetMatriceElt(1,0, sx*sy*cz + cx*sz );
    apRotMtx->SetMatriceElt(1,1, -sx*sy*sz + cx*cz );
    apRotMtx->SetMatriceElt(1,2, -sx*cy);
    apRotMtx->SetMatriceElt(1,3, 0);
    apRotMtx->SetMatriceElt(2,0, -cx*sy*cz + sx*sz );
    apRotMtx->SetMatriceElt(2,1, cx*sy*sz + sx*cz );
    apRotMtx->SetMatriceElt(2,2, cx*cy);
    apRotMtx->SetMatriceElt(2,3, 0);
    apRotMtx->SetMatriceElt(3,0, 0 );
    apRotMtx->SetMatriceElt(3,1, 0 );
    apRotMtx->SetMatriceElt(3,2, 0);
    apRotMtx->SetMatriceElt(3,3, 1);
  }
  */
  
  delete r1;
  delete r2;
  delete r3;
  delete rtmp;
  
  return 0;
}

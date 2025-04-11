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
  \brief    Implementation of class iImageConvolverStationaryIsotropicGaussian
*/

#include "vImageConvolver.hh"
#include "oImageDimensionsAndQuantification.hh"
#include "oImageSpace.hh"
#include "iImageConvolverStationaryIsotropicGaussian.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iImageConvolverStationaryIsotropicGaussian::iImageConvolverStationaryIsotropicGaussian() : vImageConvolver()
{
  m_nbSigmas = -1.;
  m_FWHM = -1.;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iImageConvolverStationaryIsotropicGaussian::~iImageConvolverStationaryIsotropicGaussian()
{
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iImageConvolverStationaryIsotropicGaussian::ShowHelp()
{
  cout << "This convolver is based on a classic stationary isotropic Gaussian kernel." << endl;
  cout << "To speed up computation time, three consecutive 1D convolutions are performed." << endl;
  cout << "The following options are used (in this particular order when provided as a list, separated by commas):" << endl;
  cout << "  FWHM: to set the isotropic FWHM (in mm)" << endl;
  cout << "  number of sigmas: to set the number of sigmas included in each dimension of the kernel (recommendations: at least 3. and maximum 5.)" << endl;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iImageConvolverStationaryIsotropicGaussian::ReadConfigurationFile(const string& a_configurationFile)
{
  string key_word = "";
  // Read the FWHM option
  key_word = "FWHM";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_FWHM, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iImageConvolverStationaryIsotropicGaussian::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Read the number of sigma option
  key_word = "number of sigmas";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_nbSigmas, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iImageConvolverStationaryIsotropicGaussian::ReadConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iImageConvolverStationaryIsotropicGaussian::ReadOptionsList(const string& a_optionsList)
{
  // There are 2 floating point variables as options
  FLTNB options[2];
  // Read them
  if (ReadStringOption(a_optionsList, options, 2, ",", "Stationary Isotropic Gaussian convolver configuration"))
  {
    Cerr("***** iImageConvolverStationaryIsotropicGaussian::ReadConfigurationFile() -> Failed to correctly read the list of options !" << endl);
    return 1;
  }
  // Affect options
  m_FWHM = options[0];
  m_nbSigmas = options[1];
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iImageConvolverStationaryIsotropicGaussian::CheckSpecificParameters()
{
  // Check number of sigmas
  if (m_nbSigmas<=0.)
  {
    Cerr("***** iImageConvolverStationaryIsotropicGaussian::CheckSpecificParameters() -> Number of sigmas included in the kernel must be strictly positive !" << endl);
    return 1;
  }
  // Check the FWHM
  if (m_FWHM<0.)
  {
    Cerr("***** iImageConvolverStationaryIsotropicGaussian::CheckSpecificParameters() -> FWHM is negative !" << endl);
    return 1;
  }
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iImageConvolverStationaryIsotropicGaussian::BuildConvolutionKernel()
{
  // Verbose
  if (m_verbose>=2)
  {
    Cout("iImageConvolverStationaryIsotropicGaussian::BuildConvolutionKernel() -> Compute convolution kernel" << endl);
    Cout("  --> Isotropic FWHM: " << m_FWHM << " mm" << endl);
    Cout("  --> Included sigmas: " << m_nbSigmas << endl);
  }

  // --------------------------------
  // Stationary kernel
  // --------------------------------

  // Set the number of kernels to 3
  m_nbKernels = 3;
  // This convolver is stationary but we overload the Convolve() function (however, the ConvolveTranspose() function calls the Convolve() function when stationary)
  m_stationary = true;
  // Compute kernel half size based on the number of sigmas included in the kernel
  FLTNB two_sqrt_two_ln_two = 2.*sqrt(2.*log(2.));
  FLTNB half_kern_floatX = m_FWHM * m_nbSigmas / (two_sqrt_two_ln_two*mp_ImageDimensionsAndQuantification->GetVoxSizeX());
  FLTNB half_kern_floatY = m_FWHM * m_nbSigmas / (two_sqrt_two_ln_two*mp_ImageDimensionsAndQuantification->GetVoxSizeY());
  FLTNB half_kern_floatZ = m_FWHM * m_nbSigmas / (two_sqrt_two_ln_two*mp_ImageDimensionsAndQuantification->GetVoxSizeZ());
  // For each half size, if the decimal part is below 0.5, then it will be taken into account into the central kernel's voxel,
  // otherwise, we had a voxel to the half size
  INTNB half_kern_dimX = ((INTNB)(half_kern_floatX));
  if (half_kern_floatX-((FLTNB)half_kern_dimX)>0.5) half_kern_dimX++;
  INTNB half_kern_dimY = ((INTNB)(half_kern_floatY));
  if (half_kern_floatY-((FLTNB)half_kern_dimY)>0.5) half_kern_dimY++;
  INTNB half_kern_dimZ = ((INTNB)(half_kern_floatZ));
  if (half_kern_floatZ-((FLTNB)half_kern_dimZ)>0.5) half_kern_dimZ++;
  // Allocate kernel dimensions
  mp_dimKernelX = (INTNB*)malloc(m_nbKernels*sizeof(INTNB));
  mp_dimKernelY = (INTNB*)malloc(m_nbKernels*sizeof(INTNB));
  mp_dimKernelZ = (INTNB*)malloc(m_nbKernels*sizeof(INTNB));
  // Compute kernel dimensions
  for (INTNB k=0; k<m_nbKernels; k++)
  {
    mp_dimKernelX[k] = half_kern_dimX * 2 + 1;
    mp_dimKernelY[k] = half_kern_dimY * 2 + 1;
    mp_dimKernelZ[k] = half_kern_dimZ * 2 + 1;
  }
  // Verbose
  if (m_verbose>=2) Cout ("  --> Kernel dimensions: [" << mp_dimKernelX[0] << ";" << mp_dimKernelY[0] << ";" << mp_dimKernelZ[0] << "]" << endl);
  // Allocate kernel
  m2p_kernel = (FLTNB**)malloc(m_nbKernels*sizeof(FLTNB*));
  m2p_kernel[KERNEL_1D_X] = (FLTNB*)malloc(mp_dimKernelX[0]*sizeof(FLTNB));
  m2p_kernel[KERNEL_1D_Y] = (FLTNB*)malloc(mp_dimKernelY[0]*sizeof(FLTNB));
  m2p_kernel[KERNEL_1D_Z] = (FLTNB*)malloc(mp_dimKernelZ[0]*sizeof(FLTNB));

  // --------------------------------
  // Compute the kernel
  // --------------------------------

  // Variables
  FLTNB sigma = m_FWHM / two_sqrt_two_ln_two;
  FLTNB mu_x = (FLTNB)(half_kern_dimX);
  FLTNB mu_y = (FLTNB)(half_kern_dimY);
  FLTNB mu_z = (FLTNB)(half_kern_dimZ);
  // Compute kernel by integration of the Gaussian distribution (this method is not exactly exact... but is closer to it compared to single Gaussian points)
  FLTNB sum_kernel_X = 0.;
  for (int x=0; x<mp_dimKernelX[0]; x++)
  {
    m2p_kernel[KERNEL_1D_X][x] = 0.5 * ( erf(  ( (((FLTNB)x)-mu_x+0.5)*mp_ImageDimensionsAndQuantification->GetVoxSizeX()) / (sqrt(2.)*sigma)  )
                                       - erf(  ( (((FLTNB)x)-mu_x-0.5)*mp_ImageDimensionsAndQuantification->GetVoxSizeX()) / (sqrt(2.)*sigma)  ) );
    sum_kernel_X += m2p_kernel[KERNEL_1D_X][x];
  }
  FLTNB sum_kernel_Y = 0.;
  for (int y=0; y<mp_dimKernelY[0]; y++)
  {
    m2p_kernel[KERNEL_1D_Y][y] = 0.5 * ( erf(  ( (((FLTNB)y)-mu_y+0.5)*mp_ImageDimensionsAndQuantification->GetVoxSizeY()) / (sqrt(2.)*sigma)  )
                                       - erf(  ( (((FLTNB)y)-mu_y-0.5)*mp_ImageDimensionsAndQuantification->GetVoxSizeY()) / (sqrt(2.)*sigma)  ) );
    sum_kernel_Y += m2p_kernel[KERNEL_1D_Y][y];
  }
  FLTNB sum_kernel_Z = 0.;
  for (int z=0; z<mp_dimKernelZ[0]; z++)
  {
    m2p_kernel[KERNEL_1D_Z][z] = 0.5 * ( erf(  ( (((FLTNB)z)-mu_z+0.5)*mp_ImageDimensionsAndQuantification->GetVoxSizeZ()) / (sqrt(2.)*sigma)  )
                                       - erf(  ( (((FLTNB)z)-mu_z-0.5)*mp_ImageDimensionsAndQuantification->GetVoxSizeZ()) / (sqrt(2.)*sigma)  ) );
    sum_kernel_Z += m2p_kernel[KERNEL_1D_Z][z];
  }
  // Normalize kernels
  for (int x=0; x<mp_dimKernelX[0]; x++) m2p_kernel[KERNEL_1D_X][x] /= sum_kernel_X;
  for (int y=0; y<mp_dimKernelY[0]; y++) m2p_kernel[KERNEL_1D_Y][y] /= sum_kernel_Y;
  for (int z=0; z<mp_dimKernelZ[0]; z++) m2p_kernel[KERNEL_1D_Z][z] /= sum_kernel_Z;
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iImageConvolverStationaryIsotropicGaussian::Convolve(FLTNB* ap_outputImage)
{
  #ifdef CASTOR_DEBUG
  if (!m_initialized)
  {
    Cerr("***** iImageConvolverStationaryIsotropicGaussian::Convolve() -> Called while not initialized !" << endl);
    Exit(EXIT_DEBUG);
  }
  #endif

  INTNB stationary_kernel = 0;
  
  // Convolve along Z
  INTNB z;
  #pragma omp parallel for private(z) schedule(guided)
  for (z=0; z<mp_ImageDimensionsAndQuantification->GetNbVoxZ(); z++)
  {
    // As a first step, we virtually go through the kernel along Z to compute the actual kernel integral that will
    // be really used in the convolution. For each virtual slice contributing to the convolved slice, we test if it
    // is inside or outside the image; in the former case we count it into the integral and in the latter we don't.
    FLTNB kernel_integral = 0.;
    for (INTNB zk=0; zk<mp_dimKernelZ[stationary_kernel]; zk++)
    {
      // Index of the virtual slice which will (or not) contribute to the convolved slice
      INTNB virtualZ = z + zk - m_offsetZ;
      // Test if it is not in the image, then we move onto the next virtual slice
      if (virtualZ<0 || virtualZ>=mp_ImageDimensionsAndQuantification->GetNbVoxZ()) continue;
      // Otherwise, we update the kernel integral
      kernel_integral += m2p_kernel[KERNEL_1D_Z][zk];
    }
    // Some precomputation for maximum efficiency
    INTNB indexZ_img_base = z*mp_ImageDimensionsAndQuantification->GetNbVoxXY();
    INTNB indexZ_pad_base = z*m_dimPadXY;
    for (INTNB y=0; y<mp_ImageDimensionsAndQuantification->GetNbVoxY(); y++)
    {
      // Some precomputation for maximum efficiency
      INTNB indexY_img_base = y*mp_ImageDimensionsAndQuantification->GetNbVoxX();
      INTNB indexY_pad_base = (y+mp_dimKernelY[stationary_kernel]/2)*m_dimPadX;
      for (INTNB x=0; x<mp_ImageDimensionsAndQuantification->GetNbVoxX(); x++)
      {
        // Output image index
        INTNB index_img = indexZ_img_base + indexY_img_base + x;
        // Base of padded index
        INTNB index_pad_base = indexZ_pad_base + indexY_pad_base + x + mp_dimKernelX[stationary_kernel]/2;
        // Inner loop on convolution kernel dimension
        for (INTNB zk=0; zk<mp_dimKernelZ[stationary_kernel]; zk++)
        {
          // Padded image index
          INTNB index_pad = index_pad_base + zk * m_dimPadXY;
          // Apply contribution
          ap_outputImage[index_img] += mp_paddedImage[index_pad] * m2p_kernel[KERNEL_1D_Z][zk] / kernel_integral;
        }
      }
    }
  }
  // Copy back to padded image
  CopyToPaddedImage(ap_outputImage);
  // Convolve along Y
  #pragma omp parallel for private(z) schedule(guided)
  for (z=0; z<mp_ImageDimensionsAndQuantification->GetNbVoxZ(); z++)
  {
    // Some precomputation for maximum efficiency
    INTNB indexZ_img_base = z*mp_ImageDimensionsAndQuantification->GetNbVoxXY();
    INTNB indexZ_pad_base = (z+m_offsetZ)*m_dimPadXY;
    for (INTNB y=0; y<mp_ImageDimensionsAndQuantification->GetNbVoxY(); y++)
    {
      // Some precomputation for maximum efficiency
      INTNB indexY_img_base = y*mp_ImageDimensionsAndQuantification->GetNbVoxX();
      INTNB indexY_pad_base = y*m_dimPadX;
      for (INTNB x=0; x<mp_ImageDimensionsAndQuantification->GetNbVoxX(); x++)
      {
        // Output image index
        INTNB index_img = indexZ_img_base + indexY_img_base + x;
        // Base of padded index
        INTNB index_pad_base = indexZ_pad_base + indexY_pad_base + x + m_offsetX;
        // Inner loop on convolution kernel dimension
        for (INTNB yk=0; yk<mp_dimKernelY[stationary_kernel]; yk++)
        {
          // Padded image index
          INTNB index_pad = index_pad_base + yk * m_dimPadX;
          // Apply contribution
          ap_outputImage[index_img] += mp_paddedImage[index_pad] * m2p_kernel[KERNEL_1D_Y][yk];
        }
      }
    }
  }
  // Copy back to padded image
  CopyToPaddedImage(ap_outputImage);
  // Convolve along X
  #pragma omp parallel for private(z) schedule(guided)
  for (z=0; z<mp_ImageDimensionsAndQuantification->GetNbVoxZ(); z++)
  {
    // Some precomputation for maximum efficiency
    INTNB indexZ_img_base = z*mp_ImageDimensionsAndQuantification->GetNbVoxXY();
    INTNB indexZ_pad_base = (z+mp_dimKernelZ[stationary_kernel]/2)*m_dimPadXY;
    for (INTNB y=0; y<mp_ImageDimensionsAndQuantification->GetNbVoxY(); y++)
    {
      // Some precomputation for maximum efficiency
      INTNB indexY_img_base = y*mp_ImageDimensionsAndQuantification->GetNbVoxX();
      INTNB indexY_pad_base = (y+mp_dimKernelY[stationary_kernel]/2)*m_dimPadX;
      for (INTNB x=0; x<mp_ImageDimensionsAndQuantification->GetNbVoxX(); x++)
      {
        // Output image index
        INTNB index_img = indexZ_img_base + indexY_img_base + x;
        // Base of padded index
        INTNB index_pad_base = indexZ_pad_base + indexY_pad_base + x;
        // Inner loop on convolution kernel dimension
        for (INTNB xk=0; xk<mp_dimKernelX[stationary_kernel]; xk++)
        {
          // Padded image index
          INTNB index_pad = index_pad_base + xk;
          // Apply contribution
          ap_outputImage[index_img] += mp_paddedImage[index_pad] * m2p_kernel[KERNEL_1D_X][xk];
        }
      }
    }
  }
  // Normal end
  return 0;
}

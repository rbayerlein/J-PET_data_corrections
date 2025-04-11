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
  \brief    Declaration of class iImageConvolverStationaryIsotropicGaussian
*/

#ifndef IIMAGECONVOLVERSTATIONARYISOTROPICGAUSSIAN_HH
#define IIMAGECONVOLVERSTATIONARYISOTROPICGAUSSIAN_HH 1

#include "gVariables.hh"
#include "gOptions.hh"
#include "oImageSpace.hh"
#include "vImageConvolver.hh"
#include "sAddonManager.hh"

#define KERNEL_1D_X 0
#define KERNEL_1D_Y 1
#define KERNEL_1D_Z 2

/*!
  \class   iImageConvolverStationaryIsotropicGaussian
  \brief   This class is an image convolver module implementing stationary and isotropic gaussian filtering
  \details This class is a child of vImageConvolver. It implements a stationary isotropic gaussian
           3D convolution filter as three separated 1D filters, parameterized by an FWHM, as well as the
           number of sigmas to be included in the convolution kernel.
*/
class iImageConvolverStationaryIsotropicGaussian : public vImageConvolver
{
  // -----------------------------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \fn      public iImageConvolverStationaryIsotropicGaussian::iImageConvolverStationaryIsotropicGaussian()
      \brief   The constructor of iImageConvolverStationaryIsotropicGaussian
      \details This is the default and unique constructor. It does not take any parameter and
               its role is to affect default values to parameters specific to this module.
    */
    iImageConvolverStationaryIsotropicGaussian();
    /*!
      \fn      public iImageConvolverStationaryIsotropicGaussian::~iImageConvolverStationaryIsotropicGaussian()
      \brief   The destructor of iImageConvolverStationaryIsotropicGaussian
      \details This is the default and unique destructor. It does not take any parameter and
               its role is only to free or delete all structures that were built specifically
               by this module. Here it does nothing.
    */
    ~iImageConvolverStationaryIsotropicGaussian();


  // -----------------------------------------------------------------------------------------
  // Public member functions
  public:
    // Function for automatic insertion (put the class name as the parameter and do not add semi-colon at the end of the line)
    FUNCTION_IMAGE_CONVOLVER(iImageConvolverStationaryIsotropicGaussian)
    /*!
      \fn      public int iImageConvolverStationaryIsotropicGaussian::ReadConfigurationFile()
      \param   const string& a_configurationFile
      \brief   A function used to read options from a configuration file
      \details This function implements the reading of all options associated to a child module, from
               a configuration file. It is the implementation of the pure virtual function inherited
               from the abstract class vImageConvolver. It checks the reading status but not
               the options values that will be checked by the CheckSpecificParameters() function.
      \return  An integer reflecting the reading success; 0 if success, another value otherwise.
    */
    int ReadConfigurationFile(const string& a_fileOptions);
    /*!
      \fn      public int iImageConvolverStationaryIsotropicGaussian::ReadOptionsList()
      \param   const string& a_configurationFile
      \brief   A function used to read options from a list of options
      \details This function implements the reading of all options associated to a child module, from
               a list of options. It is the implementation of the pure virtual function inherited
               from the abstract class vImageConvolver. It checks the reading status but not
               the options values that will be checked by the CheckSpecificParameters() function.
      \return  An integer reflecting the reading success; 0 if success, another value otherwise.
    */
    int ReadOptionsList(const string& a_listOptions);
    /*!
      \fn      public int iImageConvolverStationaryIsotropicGaussian::ShowHelp()
      \brief   A function used to show help about the child module
      \details This function must describe what the module does and how to use it. It describes in
               details the different parameters of the module, and how to set them through the use
               of a configuration file or a list of options. It is the implementation of the pure
               virtual function inherited from the abstract class vImageConvolverModule.
    */
    void ShowHelp();


  // -----------------------------------------------------------------------------------------
  // Private member functions
  private:
    /*!
      \fn      private int iImageConvolverStationaryIsotropicGaussian::CheckSpecificParameters()
      \brief   A private function used to check the parameters settings specific to the child module
      \details This function is used to check that all parameters specific to the module are correctly set
               within allowed values. It is called by the CheckParameters() function of the mother class.
               It is the implementation of the pure virtual function inherited from the abstract mother
               class vImageConvolver.
      \return  An integer reflecting the check status; 0 if no problem, another value otherwise.
    */
    int CheckSpecificParameters();
    /*!
      \fn      private int iImageConvolverStationaryIsotropicGaussian::BuildConvolutionKernel()
      \brief   A private function used to build the convolution kernel specific to the child convolver
      \details This function is used to build the convolution kernels associated to the child convolver.
               It is called by the Initialize() function. It is the implementation of the pure virtual
               function inherited from the abstract mother class vImageConvolver. In this case, the
               kernel is stationary so m_nbKernels = 1. We keep the implementation of the convolution
               from the mother class which already implements stationary convolutions.
      \return  An integer reflecting the building status; 0 if no problem, another value otherwise.
    */
    int BuildConvolutionKernel();
    /*!
      \fn      private int iImageConvolverStationaryIsotropicGaussian::Convolve()
      \param   FLTNB* ap_outputImage
      \brief   A private function used to apply the convolution on the padded image to the provided output image
      \details This is an overload of the mother function which applies 3D stationary convolutions. Here, we assume
               isotropic gaussian convolutions which allows us to perform three separated by consecutive 1D
               convolution to enhance computational efficiency.
      \return  An integer reflecting the convolution status; 0 if no problem, another value
               otherwise.
    */
    int Convolve(FLTNB* ap_outputImage);


  // -----------------------------------------------------------------------------------------
  // Data members
  protected:
    FLTNB m_FWHM;    /*!< The isotropic FWHM in mm */
    FLTNB m_nbSigmas;     /*!< The number of sigmas of the Gaussian distribution included in the kernel */
};


// Class for automatic insertion (set here the visible image convolver's name as the first parameter,
// put the class name as the second parameter and do NOT add semi-colon at the end of the line)
CLASS_IMAGE_CONVOLVER(isotropicGaussian,iImageConvolverStationaryIsotropicGaussian)

#endif

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
  \brief    Declaration of class iLinearSpectralModel
*/

#ifndef ILINEARSPECTRALMODEL_HH
#define ILINEARSPECTRALMODEL_HH 1

// Mother class of dynamic model
#include "vDynamicModel.hh"
// Required to automatically add the class in the CASToR code
#include "sAddonManager.hh"

// Inherit from Generic Linear Model
#include "iLinearModel.hh"

/*!
  \class   iLinearSpectralModel
  \brief   This class is a child of the iLinearModel class implementing the Nested Spectral reconsutction
  \details
*/
class iLinearSpectralModel : public iLinearModel
{
  // -----------------------------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \fn      iLinearSpectralModel::iLinearSpectralModel
      \brief   Constructor of iLinearSpectralModel. Simply set all data members to default values.
    */
    iLinearSpectralModel();
    /*!
      \fn      iLinearSpectralModel::~iLinearSpectralModel
      \brief   Destructor of iLinearSpectralModel
    */
    ~iLinearSpectralModel();


  // -----------------------------------------------------------------------------------------
  // Public member functions related to the initialization of the model
  public:
    // Function for automatic insertion (put the class name as the parameters and do not add semi-colon at the end of the line)
    FUNCTION_DYNAMICMODEL(iLinearSpectralModel)
    /*!
      \fn      iLinearSpectralModel::CheckSpecificParameters
      \brief   This function is used to check whether all member variables
               have been correctly initialized or not.
      \return  0 if success, positive value otherwise.
    */
    int CheckSpecificParameters();
    /*!
      \fn      iLinearSpectralModel::ReadAndCheckConfigurationFileSpecific
      \brief   This function is used to read options from a configuration file, specific to this model.
      \return  0 if success, other value otherwise.
    */
    int ReadAndCheckConfigurationFileSpecific();
    /*!
      \fn      iLinearSpectralModel::ReadAndCheckOptionsList
      \param   a_optionsList : a list of parameters separated by commas
      \brief   This function is used to read parameters from a string.
      \return  0 if success, other value otherwise.
    */
    int ReadAndCheckOptionsList(string a_listOptions);
    /*!
      \fn      iLinearSpectralModel::InitializeSpecific
      \brief   This function is used to initialize the model parametric images and basis functions
      \todo    Read Interfile for parametric images initialization
      \return  0 if success, other value otherwise.
    */
    int InitializeSpecific();
    /*!
      \fn      iLinearSpectralModel::ShowHelp
      \brief   This function is used to print out specific help about the model and its options.
    */
    void ShowHelpModelSpecific();


  // -----------------------------------------------------------------------------------------
  // Public member functions called by the main iterative algorithm class


  // -----------------------------------------------------------------------------------------
  // Data members
  protected:

    //int  m_OptimisationMethod;         /*!<Flag indicating the method to estimate Patlak parameters. */
    HPFLTNB m_fast_exp;                  /*!<Fastest spectral function. */
    HPFLTNB m_slow_exp;                  /*!<Fastest spectral function. */
    bool m_full_trapping_basis_flag;     /*!<Value for constant basis function to be added. */
    bool m_blood_fraction_fasis_flag;    /*!<Value for constant basis function to be added. */
    int m_spectral_bank_size;            /*!<requested number for size Bank of spectral function decay rates. */
    int m_additionalBF_size;             /*!<Size of additional basis functions required for the model (eg. constant etc) . */
    HPFLTNB* mp_spectral_bank;           /*!<Bank of spectral function decay rates. */

};

// Class for automatic insertion (set here the visible dynamic model's name, put the class name as the parameters and do not add semi-colon at the end of the line)
CLASS_DYNAMICMODEL(Spectral,iLinearSpectralModel)

#endif

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
  \brief    Declaration of class iDynamicModelTemplate
*/

#ifndef IDYNAMICMODELTEMPLATE_HH
#define IDYNAMICMODELTEMPLATE_HH 1

// Mother class of dynamic model
#include "vDynamicModel.hh"
// Required to automatically add the class in the CASToR code
#include "sAddonManager.hh"

// Inherit from Generic Linear Model
#include "iLinearModel.hh"

// This class inherits from iLinearModel,
// which  itself inherits from the virtual vDynamicModel
//  vDynamicModel
//  |
//  |-->> iLinearModel
//        |
//        |-->> iLinearModelTemplate

/*!
  \class   iLinearModelTemplate
  \brief   This class is a child of the vDynamicModel class implementing a template squeleton
  \details Use this class to implement your own custom deformation model.
*/

class iLinearModelTemplate : public iLinearModel
{
  // -----------------------------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \fn      iLinearModelTemplate::iLinearModelTemplate
      \brief   Constructor of iLinearModelTemplate. Simply set all data members to default values.
    */
    iLinearModelTemplate();
    /*!
      \fn      iLinearModelTemplate::~iLinearModelTemplate
      \brief   Destructor of iLinearModelTemplate
    */
    ~iLinearModelTemplate();


  // -----------------------------------------------------------------------------------------
  // Public member functions related to the initialization of the model
  public:
    // Function for automatic insertion (put the class name as the parameters and do not add semi-colon at the end of the line)
    FUNCTION_DYNAMICMODEL(iLinearModelTemplate)
    /*!
      \fn      iLinearModelTemplate::CheckSpecificParameters
      \brief   This function is used to check whether all member variables
               have been correctly initialized or not.
      \return  0 if success, positive value otherwise.
    */
    int CheckSpecificParameters();
    /*!
      \fn      iLinearModelTemplate::ReadAndCheckConfigurationFileSpecific
      \brief   This function is used to read options from a configuration file, specific to this model.
      \return  0 if success, other value otherwise.
    */
    int ReadAndCheckConfigurationFileSpecific();
    /*!
      \fn      iLinearModelTemplate::ReadAndCheckOptionsList
      \param   a_optionsList : a list of parameters separated by commas
      \brief   This function is used to read parameters from a string.
      \return  0 if success, other value otherwise.
    */
    int ReadAndCheckOptionsList(string a_listOptions);
    /*!
      \fn      iLinearModelTemplate::InitializeSpecific
      \brief   This function is used to initialize the model parametric images and basis functions
      \todo    Read Interfile for parametric images initialization
      \return  0 if success, other value otherwise.
    */
    int InitializeSpecific();
    /*!
      \fn      iLinearModelTemplate::ShowHelp
      \brief   This function is used to print out specific help about the model and its options.
    */
    void ShowHelpModelSpecific();


  // -----------------------------------------------------------------------------------------
  // Public member functions called by the main iterative algorithm class


  // -----------------------------------------------------------------------------------------
  // Data members
  protected:

};

// Class for automatic insertion (set here the visible dynamic model's name, put the class name as the parameters and do not add semi-colon at the end of the line)
CLASS_DYNAMICMODEL(template,iLinearModelTemplate)

#endif
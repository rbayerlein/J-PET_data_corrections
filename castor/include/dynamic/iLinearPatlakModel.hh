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
  \brief    Declaration of class iLinearPatlakModel
*/

#ifndef ILINEARPATLAKMODEL_HH
#define ILINEARPATLAKMODEL_HH 1


// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

#include "vDynamicModel.hh"
#include "iLinearModel.hh"
#include "sAddonManager.hh"

/*!
  \class   iLinearPatlakModel
  \brief   This class implements the Patlak model, to model kinetics of irreversible radiotracers
*/
class iLinearPatlakModel : public iLinearModel
{
  // -----------------------------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \fn      iLinearPatlakModel::iLinearPatlakModel
      \brief   Constructor of iLinearPatlakModel. Simply set all data members to default values.
    */
    iLinearPatlakModel();
    /*!
      \fn      iLinearPatlakModel::~iLinearPatlakModel
      \brief   Destructor of iLinearPatlakModel
    */
    ~iLinearPatlakModel();


  // -----------------------------------------------------------------------------------------
  // Public member functions related to the initialization of the model
  public:
    // Function for automatic insertion (put the class name as the parameters and do not add semi-colon at the end of the line)
    FUNCTION_DYNAMICMODEL(iLinearPatlakModel)
    /*!
      \fn      iLinearPatlakModel::CheckSpecificParameters
      \brief   This function is used to check whether all member variables
               have been correctly initialized or not.
      \return  0 if success, positive value otherwise.
    */
    int CheckSpecificParameters();
    /*!
      \fn      iLinearPatlakModel::ReadAndCheckConfigurationFileSpecific
      \brief   This function is used to read options from a configuration file.
      \return  0 if success, other value otherwise.
    */
    int ReadAndCheckConfigurationFileSpecific();
    /*!
      \fn      iLinearPatlakModel::ReadAndCheckOptionsList
      \param   const string& a_optionsList : a list of parameters separated by commas
      \brief   This function is used to read parameters from a string.
      \return  0 if success, other value otherwise.
    */
   int ReadAndCheckOptionsList(string a_listOptions);
    /*!
      \fn      iLinearPatlakModel::InitializeSpecific
      \brief   This function is used to initialize Patlak parametric images and basis functions
      \todo    Read Interfile for parametric images initialization
      \return  0 if success, other value otherwise.
    */
    int InitializeSpecific();
    /*!
      \fn      iLinearPatlakModel::ShowHelpModelSpecific
      \brief   Print out specific help about the implementation of the Patlak
               model and its initialization
    */
    void ShowHelpModelSpecific();


  // -----------------------------------------------------------------------------------------


  // -----------------------------------------------------------------------------------------
  // Data members
  protected:


};

// Class for automatic insertion (set here the visible dynamic model's name, put the class name as the parameters and do not add semi-colon at the end of the line)
CLASS_DYNAMICMODEL(Patlak,iLinearPatlakModel)

#endif

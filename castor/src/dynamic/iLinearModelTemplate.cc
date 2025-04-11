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
  \brief    Implementation of class iLinearModelTemplate
*/

#include "iLinearModelTemplate.hh"



// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn iLinearModelTemplate
  \brief Constructor of iLinearModelTemplate. Simply set all data members to default values.
*/
iLinearModelTemplate::iLinearModelTemplate() : iLinearModel()
{
  // --- Parameters inherited from vDynamicModel class --- //

  m_nbTimeBF = 1; // Number of basis functions in the model
  m_nbModelParam = 1; // Number of model parameters
  m_nnlsN = 2; // Number of parameters for NNLS optimisation
  // Set default Card and Resp gating parameters to 1 as this is only a dynamic model
  m_nbRgateBF = 1; m_nbRGModelParam=1;
  m_nbCgateBF = 1; m_nbCGModelParam=1;

}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ~iLinearPatlakModel
  \brief Destructor of iLinearPatlakModel
*/
iLinearModelTemplate::~iLinearModelTemplate()
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
  \brief Print out specific help about the implementation of this linear dynamic
         model and its initialization
*/

void iLinearModelTemplate::ShowHelpModelSpecific()
{
  // ===================================================================
  // Here, display some help and guidance to how to use this linear dynamic model and what it does
  // ===================================================================
  cout << "This class is a template class dedicated to add your own linear dynamic model." << endl;
}





// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ReadAndCheckConfigurationFileSpecific
  \param const string& a_configurationFile : ASCII file containing information about the linear dynamic model
  \brief This function is used to read options from a configuration file.
  \return 0 if success, other value otherwise.
*/

int iLinearModelTemplate::ReadAndCheckConfigurationFileSpecific()
{

  if(m_verbose >=3) Cout("iLinearModelTemplate::ReadAndCheckConfigurationFileSpecific ..."<< endl);

  // Apply the generic linear parameter for all Linear Models
  if(ReadAndCheckConfigurationFileSpecificToAllLinearModels())
  {
    Cerr("***** iLinearModelTemplate::ReadAndCheckConfigurationFileSpecific -> Error while trying to read configuration file for generic options of all linear models" << endl);
    return 1;
  }

  // ===================================================================
  // Implement here the reading of any options specific to this linear dynamic model
  // Generic parameters that apply to genera linear models will be read by the
  // function ReadAndCheckConfigurationFileSpecificToAllLinearModels()
  // The ReadDataASCIIFile() functions could be helpful to recover data from a file
  // (check other linear dynamicModels for examples)
  // ===================================================================


  // Normal End
  return 0;
}


// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn ReadAndCheckOptionsList
  \brief This function is used to read parameters from a string.
  \return 0 if success, other value otherwise.
*/

int iLinearModelTemplate::ReadAndCheckOptionsList(string a_listOptions)
{
  if(m_verbose >=3) Cout("iLinearModelTemplate::ReadAndCheckOptionsList ..."<< endl);

  // ===================================================================
  // Implement here the reading of any options specific to this linear dynamic model,
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

int iLinearModelTemplate::CheckSpecificParameters()
{
  if(m_verbose >=3) Cout("iLinearModelTemplate::CheckSpecificParameters ..."<< endl);

  // Perform generic checks that apply for the Linear Models
  if(CheckSpecificParametersForAllLinearModels())
  {
    Cerr("***** iLinearModelTemplate::CheckSpecificParameters -> A problem occurred while checking specific parameters ! " << endl);
    return 1;
  }

  // ===================================================================
  // Implement here checks over parameters for this specific linear
  // dynamic model which should be read using either ReadAndCheckConfigurationFile()
  // or ReadAndCheckOptionsList() functions.
  // Generic parameters applying to all linear models will be checked
  // using the function CheckSpecificParametersForAllLinearModels()
  // ===================================================================


  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
/*
  \fn InitializeSpecific
  \brief This function is used to initialize the linear model parametric images and basis functions
  \return 0 if success, other value otherwise.
*/
int iLinearModelTemplate::InitializeSpecific()
{
  if(m_verbose >=3) Cout("iLinearModelTemplate::InitializeSpecific ..."<< endl);

  // Forbid initialization without check
  if (!m_checked)
  {
    Cerr("***** iLinearModelTemplate::InitializeSpecific() -> Must call CheckParameters functions before Initialize() !" << endl);
    return 1;
  }

  // Run generic Initialization for all Linear Models
  if (InitializeSpecificToAllLinearModels())
  {
    Cerr("***** iLinearPatlakModel::InitializeSpecific() -> Error while performing generic initialisations for linear models !" << endl);
    return 1;
  }

  // ===================================================================
  // Implement here the allocation/initialization of whatever member
  // variables specifically used by this linear dynamic model.
  // The generic variables will be initialised by the
  // InitializeSpecificToAllLinearModels() function.
  // ===================================================================

  // ===================================================================
  // if your model makes use of an Input function curve, implement here
  // specific calculations to be applied on the interpolated input curve;
  // Input and interpolation of the curve is handled by oArterialInputCurve
  // See other implemented linear models for specific examples.
  // ===================================================================

 // --- Default Initialization of time basis functions --- //
 // The basis function imaging ( parametric images ) can be initialised here
 // By default the InitializeSpecificToAllLinearModels() function
 // will initialise them with 1

  // Normal end
  m_initialized = true;
  return 0;
}

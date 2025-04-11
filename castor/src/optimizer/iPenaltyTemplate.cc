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
  \ingroup  optimizer
  \brief    Implementation of class iPenaltyTemplate
*/

#include "iPenaltyTemplate.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iPenaltyTemplate::iPenaltyTemplate() : vPenalty()
{
  // ---------------------------
  // Mandatory member parameters
  // (inherited from vPenalty)
  // ---------------------------

  // Specify here the derivative order of the penalty.
  // Most algorithms able to handle penalties require 1 derivative,
  // and some 2 derivatives. If infinite, then set INT_MAX.
  m_penaltyDerivativesOrder = INT_MAX;

  // --------------------------
  // Specific member parameters
  // --------------------------

  // Set here the default values of the parameters specific to this optimizer

}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

iPenaltyTemplate::~iPenaltyTemplate()
{
  // Delete or free ONLY things that were built into this class
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iPenaltyTemplate::ShowHelpSpecific()
{
  cout << "This penalty is only a squeleton template to explain how to add a penalty into CASToR. If you" << endl;
  cout << "want to implement your own penalty, start from here and look at the specific documentation." << endl;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iPenaltyTemplate::ReadConfigurationFile(const string& a_configurationFile)
{
  // This function is designed to read parameters specific to the optimizer through a configuration file.
  // To do that, use the ReadDataASCIIFile() function, and take a look at other optimizers to see how it is done.
  // Do not check the parameters' values, the CheckSpecificParameters() function is designed to do that.
  // Return 1 if any problem. See other penalties' implementation to get guidance.

  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================
 
int iPenaltyTemplate::ReadOptionsList(const string& a_optionsList)
{ 
  // This function is designed to read parameters specific to the optimizer through a list of options in a string.
  // To do that, use the ReadStringOption() function, and take a look at other optimizers to see how it is done.
  // Do not check the parameters' values, the CheckSpecificParameters() function is designed to do that.
  // Return 1 if any problem. See other penalties' implementation to get guidance.

  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iPenaltyTemplate::CheckSpecificParameters()
{
  // This function is designed to check that all parameters specific to the optimizer have been set properly.
  // Return 1 if any problem.

  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iPenaltyTemplate::InitializeSpecific()
{
  // Verbose
  if (m_verbose>=2) Cout("iOptimizerTemplate::InitializeSpecific() -> Use the template optimizer" << endl);

  // This function is designed to initialize everything that should be initialized before being able to launch the iteration.
  // Return 1 if any problem.

  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

FLTNB iPenaltyTemplate::ComputePenaltyValue(int a_tbf, int a_rbf, int a_cbf, INTNB a_voxel, int a_th)
{
  // This is where you must implement the computation of the penalty value, for the provided parameters.
  // It is used to compute the overall cost function.
  // This function is called inside a parallel loop. The index of the thread is provided.
  // See other penalties' implementation to get guidance.
  FLTNB penalty = 0.;
  return penalty;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

FLTNB iPenaltyTemplate::ComputeFirstDerivative(int a_tbf, int a_rbf, int a_cbf, INTNB a_voxel, int a_th)
{
  // This is where you must implement the computation of the penalty's first derivative value, for the
  // provided parameters. It is used by the optimizer to compute the next update.
  // This function is called inside a parallel loop. The index of the thread is provided.
  // See other penalties' implementation to get guidance.
  FLTNB first_derivative = 0.;
  return first_derivative;
}


// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

FLTNB iPenaltyTemplate::ComputeSecondDerivative(int a_tbf, int a_rbf, int a_cbf, INTNB a_voxel, int a_th)
{
  // This is where you must implement the computation of the penalty's second derivative value, for the
  // provided parameters. It may be used by the optimizer to compute the next update.
  // If the penalty does not admit a second derivative, then just let the function as is.
  // This function is called inside a parallel loop. The index of the thread is provided.
  // See other penalties' implementation to get guidance.
  FLTNB second_derivative = 0.;
  return second_derivative;
}

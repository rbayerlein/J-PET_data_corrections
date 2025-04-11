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
  \brief    Declaration of class iPenaltyTemplate
*/

#ifndef IPENALTYTEMPLATE_HH
#define IPENALTYTEMPLATE_HH 1

#include "vPenalty.hh"
#include "sAddonManager.hh"
#include "sOutputManager.hh"

/*!
  \class   iPenaltyTemplate
  \brief   This class is a template for penalties
  \details This class inherits from vPenalty and provides details on how to implement a penalty.
*/
class iPenaltyTemplate : public vPenalty
{
  // -------------------------------------------------------------------
  // Constructor & Destructor
  public:

    /*!
      \fn      public iPenaltyTemplate::iPenaltyTemplate()
      \brief   The constructor of iPenaltyTemplate
      \details This is the default and unique constructor. It does not take any parameter and
               its role is only to affect default values to each member of the class.
    */
    iPenaltyTemplate();
    /*!
      \fn      public iPenaltyTemplate::~iPenaltyTemplate()
      \brief   The destructor of iPenaltyTemplate
      \details This is the default and unique destructor. It does not take any parameter and
               its role is only to free or delete all structures that were built by this class.
    */
    ~iPenaltyTemplate();

  // -------------------------------------------------------------------
  // Public member functions
  public:
    // Function for automatic insertion (put the class name as the parameter and do not add semi-column at the end of the line)
    FUNCTION_PENALTY(iPenaltyTemplate)
    /*!
      \fn      public int iPenaltyTemplate::ReadConfigurationFile()
      \param   const string& a_configurationFile
      \brief   A function used to read options from a configuration file
      \details This function implements the reading of all options associated to the child penalty, from
               a configuration file. It is the implementation of the pure virtual function inherited
               from the abstract class vPenalty. It checks the reading status but not
               the options values that will be checked by the CheckSpecificParameters() function.
      \return  An integer reflecting the reading success; 0 if success, another value otherwise.
    */
    int ReadConfigurationFile(const string& a_configurationFile);
    /*!
      \fn      public int iPenaltyTemplate::ReadOptionsList()
      \param   const string& a_optionsList
      \brief   A function used to read options from a list of options
      \details This function implements the reading of all options associated to the child penalty, from
               a list of options. It is the implementation of the pure virtual function inherited
               from the abstract class vPenalty. It checks the reading status but not
               the options values that will be checked by the CheckSpecificParameters() function.
      \return  An integer reflecting the reading success; 0 if success, another value otherwise.
    */
    int ReadOptionsList(const string& a_optionsList);
    /*!
      \fn      public int iPenaltyTemplate::ComputePenaltyValue()
      \param   int a_tbf
      \param   int a_rbf
      \param   int a_cbf
      \param   INTNB a_voxel
      \param   int a_th
      \brief   Implementation of the pure virtual vPenalty::ComputePenaltyValue()
      \details This function computes the value of the penalty function for the provided indices.
               It is the implementation of the pure virtual vPenalty::ComputePenaltyValue().
      \return  The penalty value
    */
    FLTNB ComputePenaltyValue(int a_tbf, int a_rbf, int a_cbf, INTNB a_voxel, int a_th);
    /*!
      \fn      public int iPenaltyTemplate::ComputeFirstDerivative()
      \param   int a_tbf
      \param   int a_rbf
      \param   int a_cbf
      \param   INTNB a_voxel
      \param   int a_th
      \brief   Implementation of the pure virtual vPenalty::ComputeFirstDerivative()
      \details This function computes the first derivative of the penalty.
               It is the implementation of the pure virtual vPenalty::ComputeFirstDerivative().
      \return  The derivative
    */
    FLTNB ComputeFirstDerivative(int a_tbf, int a_rbf, int a_cbf, INTNB a_voxel, int a_th);
    /*!
      \fn      public int iPenaltyTemplate::ComputeSecondDerivative()
      \param   int a_tbf
      \param   int a_rbf
      \param   int a_cbf
      \param   INTNB a_voxel
      \param   int a_th
      \brief   Implementation of the pure virtual vPenalty::ComputeSecondDerivative()
      \details This function computes the second derivative of the penalty.
               It is the implementation of the pure virtual vPenalty::ComputeSecondDerivative().
      \return  The derivative
    */
    FLTNB ComputeSecondDerivative(int a_tbf, int a_rbf, int a_cbf, INTNB a_voxel, int a_th);

  // -------------------------------------------------------------------
  // Private member functions
  private:
    /*!
      \fn      private void iPenaltyTemplate::ShowHelpSpecific()
      \brief   A function used to show help about the child penalty
      \details This function must describe what the module does and how to use it. It describes in
               details the different parameters of the penalty, and how to set them through the use
               of a configuration file or a list of options. It is the implementation of the pure
               virtual function inherited from the abstract class vPenalty. It is called by the
               public ShowHelp() function.
    */
    void ShowHelpSpecific();
    /*!
      \fn      private int iPenaltyTemplate::CheckSpecificParameters()
      \brief   A private function used to check the parameters settings specific to the child penalty
      \details This function is used to check that all parameters specific to the penalty are correctly set
               within allowed values. It is called by the CheckParameters() function of the mother class.
               It is the implementation of the pure virtual function inherited from the abstract mother
               class vPenalty.
      \return  An integer reflecting the check status; 0 if no problem, another value otherwise.
    */
    int CheckSpecificParameters();
    /*!
      \fn      private int iPenaltyTemplate::InitializeSpecific()
      \brief   This function is used to initialize specific stuff to the child penalty.
      \details It is called by the public Initialize() function from the mother.
      \return  An integer reflecting the initialization status; 0 if no problem, another value otherwise.
    */
    int InitializeSpecific();

  // -------------------------------------------------------------------
  // Data members
  private:

};

// Class for automatic insertion (set here the visible optimizer's name as the first parameter,
// put the class name as the second parameter and do NOT add semi-colon at the end of the line)
CLASS_PENALTY(template,iPenaltyTemplate)

#endif


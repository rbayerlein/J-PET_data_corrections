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
  \brief    Declaration of class vPenalty
*/

#ifndef VPENALTY_HH
#define VPENALTY_HH 1

#include "gVariables.hh"
#include "oImageDimensionsAndQuantification.hh"
#include "oImageSpace.hh"
#include "vDataFile.hh"
#include "oProjectionLine.hh"

/*!
  \class   vPenalty
  \brief   This class is designed to generically described any penalty applied to MAP algorithms
  \details This class is an abstract one, in the sense that it cannot be used on its own
           because several pure virtual functions belong to it. Its children are
           implementations of actual penalties. Everywhere in the code, this parent class
           should be used instead of any of its children.
           Nothing is yet implemented. To be designed.
*/
class vPenalty
{
  // -------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \fn      public vPenalty::vPenalty()
      \brief   The constructor of vPenalty
      \details This is the default and unique constructor. It does not take any parameter and
               its role is only to affect default values to each member of the class.
    */
    vPenalty();
    /*!
      \fn      virtual public vPenalty::~vPenalty()
      \brief   The destructor of vPenalty
      \details This is the default and unique destructor. It does not take any parameter and
               its role is only to free or delete all structures that were built by this class.
               It is virtual, so that it is automatically called when a child object is deleted.
    */
    virtual ~vPenalty();


  // -------------------------------------------------------------------
  // Public member functions
  public:
  
    /*!
      \fn      public int vPenalty::ShowHelp()
      \brief   A function used to show help about the penalty
      \details This function simply calls the ShowHelpSpecific() function implemented by children.
    */
    void ShowHelp();
    /*!
      \fn      public int vPenalty::CheckParameters()
      \brief   A public function used to check the parameters settings
      \details This function does not take any parameter and is used to check that all mandatory
               members were correctly parameterized. At the end, it calls the pure virtual
               CheckSpecificParameters() function implemented by children.
      \return  An integer reflecting the check status; 0 if no problem, another value otherwise.
    */
    int CheckParameters();
    /*!
      \fn      public int vPenalty::Initialize()
      \brief   A public function used to initialize the penalty
      \details This function does not take any parameter and is used to initialize everything that
               should be initialized. At the end, it calls the pure virtual InitializeSpecific()
               function implemented by children.
      \return  An integer reflecting the initialization status; 0 if no problem, another value otherwise.
    */
    int Initialize();
    /*!
      \fn      public virtual int vPenalty::GlobalPreProcessingStep()
      \brief   A public function computing a global pre-processing step for the penalty
      \details This function of this mother class does nothing but can be overloaded by children in order
               to compute some required pre-processing step before going through the loops on dimensions.
      \return  An integer reflecting the initialization status; 0 if no problem, another value otherwise.
    */  
    virtual int GlobalPreProcessingStep();
    /*!
      \fn      public virtual int vPenalty::LocalPreProcessingStep()
      \param   int a_timeBasisFunction
      \param   int a_respBasisFunction
      \param   int a_cardBasisFunction
      \param   INTNB a_voxel
      \param   int a_th
      \brief   A public function computing a local pre-processing step for the penalty
      \details This function of this mother class does nothing but can be overloaded by children in order
               to compute some required pre-processing step before being able to compute the first and
               second derivatives for a given voxel, or the actual penalty value.
      \return  An integer reflecting the initialization status; 0 if no problem, another value otherwise.
    */  
    virtual int LocalPreProcessingStep(int a_tbf, int a_rbf, int a_cbf, INTNB a_voxel, int a_th);
    /*!  
      \fn      public virtual int vPenalty::ComputePenaltyValue() = 0
      \param   int a_timeBasisFunction
      \param   int a_respBasisFunction
      \param   int a_cardBasisFunction
      \param   INTNB a_voxel
      \param   int a_th
      \brief   A public function computing the value of the penalty function
      \details This function computes the value of the penalty function for the provided indices.
               This function is supposed to be called within a multi-threaded loop, so the thread index is provided.
      \return  The penalty value
    */  
    virtual FLTNB ComputePenaltyValue(int a_tbf, int a_rbf, int a_cbf, INTNB a_voxel, int a_th) = 0;
    /*!  
      \fn      public virtual int vPenalty::ComputeFirstDerivative() = 0
      \param   int a_timeBasisFunction
      \param   int a_respBasisFunction
      \param   int a_cardBasisFunction
      \param   INTNB a_voxel
      \param   int a_th
      \brief   A public function computing the derivative of the penalty
      \details This function computes the first derivative of the penalty for the provided indices.
               This function is supposed to be called within a multi-threaded loop, so the thread index is provided.
      \return  The derivative
    */  
    virtual FLTNB ComputeFirstDerivative(int a_tbf, int a_rbf, int a_cbf, INTNB a_voxel, int a_th) = 0;
    /*!  
      \fn      public virtual int vPenalty::ComputeSecondDerivative() = 0
      \param   int a_timeBasisFunction
      \param   int a_respBasisFunction
      \param   int a_cardBasisFunction
      \param   INTNB a_voxel
      \param   int a_th
      \brief   A public function computing the second derivative of the penalty (the two derivatives are according to the same variable)
      \details This function compute the second derivative of the penalty for the provided indices. (the two derivatives are according to the same variable)
               This function is supposed to be called within a multi-threaded loop, so the thread index is provided.
      \return  The second derivative
    */  
    virtual FLTNB ComputeSecondDerivative(int a_tbf, int a_rbf, int a_cbf, INTNB a_voxel, int a_th) = 0;
  

  // -------------------------------------------------------------------
  // Public Get & Set functions
  public:
    /*!
      \fn      public inline void vPenalty::SetVerbose()
      \param   int a_verboseLevel
      \brief   Set the verbose level.
    */
    inline void SetVerbose(int a_verbose)
           {m_verbose = a_verbose;}
    /*!
      \fn      public inline void vPenalty::SetImageDimensionsAndQuantification()
      \param   oImageDimensionsAndQuantification* ap_ImageDimensionsAndQuantification
      \brief   Set the pointer to the image dimensions in use.
    */
    inline void SetImageDimensionsAndQuantification(oImageDimensionsAndQuantification* ap_ImageDimensionsAndQuantification)
           {mp_ImageDimensionsAndQuantification = ap_ImageDimensionsAndQuantification;}
    /*!
      \fn      public inline void vPenalty::SetImageSpace()
      \param   oImageSpace* ap_ImageSpace
      \brief   Set the pointer to the image space in use.
    */
    inline void SetImageSpace(oImageSpace* ap_ImageSpace)
           {mp_ImageSpace = ap_ImageSpace;}
    /*!
      \fn      public inline void vPenalty::SetPenaltyStrength()
      \param   FLTNB a_penaltyStrength
      \brief   Set the penalty strength.
    */
    inline void SetPenaltyStrength(FLTNB a_penaltyStrength)
           {m_penaltyStrength = a_penaltyStrength;}
    /*!
      \fn      public inline int vPenalty::GetPenaltyStrength()
      \brief   Get the penalty strength.
      \return  m_penaltyStrength
    */
    inline FLTNB GetPenaltyStrength()
           {return m_penaltyStrength;}
    /*!
      \fn      public inline int vPenalty::GetPenaltyDerivativesOrder()
      \brief   Get the penalty deratives order.
      \return  m_penaltyDerivativesOrder
    */
    inline int GetPenaltyDerivativesOrder()
           {return m_penaltyDerivativesOrder;}
    /*!
      \fn      public inline void vPenalty::SetPenaltyID()
      \param   const string& a_penaltyID
      \brief   Set the penalty ID
    */
    inline void SetPenaltyID(const string& a_penaltyID)
           {m_penaltyID = a_penaltyID;}
    /*!
      \fn      public inline const string& vPenalty::GetPenaltyID()
      \return  m_penaltyID
    */
    inline const string& GetPenaltyID()
           {return m_penaltyID;}


  // -----------------------------------------------------------------------------------------
  // Pure virtual private member functions that need to be implemented by child optimizers
  private:
    /*!
      \fn      private virtual int vPenalty::CheckSpecificParameters() = 0
      \brief   A private function used to check the parameters settings specific to the child penalty
      \details This function is used to check that all parameters specific to the penalty are correctly set
               within allowed values. It is called by the CheckParameters() function. It is pure virtual
               so is implemented by children.
      \return  An integer reflecting the check status; 0 if no problem, another value otherwise.
    */
    virtual int CheckSpecificParameters() = 0;
    /*!
      \fn      private virtual void vPenalty::ShowHelpSpecific() = 0
      \brief   A function used to show help about the child module
      \details This function must describe what the penalty does and how to use it. It describes in
               details the different parameters of the penalty, and how to set them through the use
               of a configuration file or a list of options. It is pure virtual so is implemented by
               children. It is private because called by the public ShowHelp() function.
    */
    virtual void ShowHelpSpecific() = 0;
    /*!
      \fn      private virtual int vPenalty::InitializeSpecific() = 0
      \brief   A private function used to initialize everything specific to the child penalty
      \details This function is used to initialize everything specific to the penalty that should be
               initialized. It is called by the Initialize() function. It is pure virtual so is
               implemented by children.
      \return  An integer reflecting the initialization status; 0 if no problem, another value
               otherwise.
    */
    virtual int InitializeSpecific() = 0;


  // -------------------------------------------------------------------
  // Virtual public member functions that may be overloaded by specific child penalties
  public:
  
    /*!
      \fn      public virtual int vPenalty::ReadConfigurationFile() = 0
      \param   const string& a_configurationFile
      \brief   A function used to read options from a configuration file and check that the corresponding values are correct
      \details This function implements the reading of all options associated to a child penalty, from
               a configuration file. It is pure virtual so is implemented by children. It checks the reading
               status but not the options values that will be checked by the CheckSpecificParameters()
               function.
      \return  An integer reflecting the reading success; 0 if success, another value otherwise.
    */
    virtual int ReadConfigurationFile( const string& a_configurationFile ) = 0;
    /*!
      \fn      public virtual int vPenalty::ReadOptionsList() = 0
      \param   const string& a_optionsList
      \brief   A function used to read options from a list of options
      \details This function implements the reading of all options associated to a child penalty, from
               a list of options. It is pure virtual so is implemented by children. It checks the reading
               status but not the options values that will be checked by the CheckSpecificParameters()
               function.
      \return  An integer reflecting the reading success; 0 if success, another value otherwise.
    */
    virtual int ReadOptionsList( const string& a_optionsList ) = 0;
    
  // -------------------------------------------------------------------
  // Data members
  protected:
    // Penalty ID
    string m_penaltyID;                    /*!< String containing the name provided as the class identifer in the children classes */
    // Image dimensions
    oImageDimensionsAndQuantification* 
      mp_ImageDimensionsAndQuantification; /*!< The pointer to the image dimensions and quantification object */
    oImageSpace* mp_ImageSpace;            /*!< The pointer to the image object */
    // Verbosity
    int m_verbose;                         /*!< The verbose level */
    // Order of penalty derivatives
    int m_penaltyDerivativesOrder;         /*!< The derivative order of the penalty */
    FLTNB m_penaltyStrength;               /*!< Regularization parameter (penalty strength) */
};


// ---------------------------------------------------------------------
// Part of code that manages the auto declaration of children classes
// ---------------------------------------------------------------------

// Macro for the function that creates the object
#define FUNCTION_PENALTY(CLASS) \
  static vPenalty *make_penalty() { return new CLASS(); };

// Macro for the class that links the appropriate function to the map of objects
#define CLASS_PENALTY(NAME,CLASS)                                                          \
  class NAME##PenaltyCreator                                                               \
  {                                                                                        \
    public:                                                                                \
      NAME##PenaltyCreator()                                                               \
        { sAddonManager::GetInstance()->mp_listOfPenalties[#NAME] = CLASS::make_penalty; } \
  };                                                                                       \
  static NAME##PenaltyCreator PenaltyCreator##NAME;

#endif

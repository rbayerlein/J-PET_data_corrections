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
  \ingroup  projector
  \brief    Declaration of class vProjector
*/

#ifndef VPROJECTOR_HH
#define VPROJECTOR_HH 1

#include "gVariables.hh"
#include "vEvent.hh"
#include "oProjectionLine.hh"
#include "vScanner.hh"
#include "oImageDimensionsAndQuantification.hh"

/**
 * @defgroup TOF_USE TOF use
 *
 *    \brief Keywords defining the type of Time of Flight enabled or not \n
 *           Defined in vProjector.hh
 * @{
 */
/** Constant integer saying TOF is currently being used through histogrammed data */
#define USE_TOFHISTO 1
/** Constant integer saying TOF is currently being used through list-mode data */
#define USE_TOFLIST 2
/** Constant integer saying TOF is not currently being used */
#define USE_NOTOF  3
/** @} */

/*!
  \file    vProjector.hh
  \class   vProjector
  \brief   This class is designed to generically described any on-the-fly projector
  \details This class is an abstract one, in the sense that it cannot be used on its own
           because several pure virtual functions belong to it. Its children are
           implementations of actual on-the-fly projectors. Everywhere in the code,
           this parent class should be used instead of any of its children. It can be used
           during the projection/reconstruction process by the oProjectorManager through the
           use of the Project() function that cannot be overloaded. This function is called
           from the Project() function of the oProjectionManager to get a oProjectionLine
           associated to an vEvent. \n
           All children must implement the following pure virtual functions: \n
            - ReadAndCheckConfigurationFile(): read specific options from a configuration file \n
            - ReadAndCheckOptionsList(): read specific options from a string \n
            - ShowHelp(): print helps about the projector specifications \n
            - Initialize(): initialize specific stuff of the projector \n
            - ProjectWithoutTOF(): project without the TOF info \n
            - ProjectTOFListmode(): project with TOF info from list-mode data \n
            - ProjectTOFHistogram(): project with TOF info from histogram data
*/
class vProjector
{
  // -------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \fn      public vProjector::vProjector()
      \brief   The constructor of vProjector
      \details This is the default and unique constructor. It does not take any parameter and
               its role is only to affect default values to each member of the class.
    */
    vProjector();
    /*!
      \fn      virtual public vProjector::~vProjector()
      \brief   The destructor of vProjector
      \details This is the default and unique destructor. It does not take any parameter and
               its role is only to free or delete all structures that were built by this class.
               It is virtual, so that it is automatically called when a child object is deleted.
    */
    virtual ~vProjector();


  // -------------------------------------------------------------------
  // Public member functions
  public:
    /*!
      \fn      public int vProjector::ShowHelp()
      \brief   A function used to show help about the projector
      \details This function simply calls the ShowHelpSpecific() function implemented by children.
    */
    void ShowHelp();
    /*!
      \fn      public static int vProjector::ShowCommonHelp()
      \brief   This function is used to print out some help about the use of options common to all projectors.
               It is static because it is called in main without instantiating an object.
    */
    static void ShowCommonHelp();
    /*!
      \fn      public int vProjector::ReadCommonOptionsList()
      \brief   This function is used to read options common to all projectors given as a string.
      \return  An integer reflecting the reading success; 0 if success, another value otherwise.
    */
    int ReadCommonOptionsList(const string& a_optionsList);
    /*!
      \fn      public int vProjector::CheckParameters()
      \brief   A public function used to check the parameters settings
      \details This function does not take any parameter and is used to check that all mandatory
               members were correctly parameterized. At the end, it calls the pure virtual
               CheckSpecificParameters() function implemented by children.
      \return  An integer reflecting the check status; 0 if no problem, another value otherwise.
    */
    int CheckParameters();
    /*!
      \fn      public int vProjector::Initialize()
      \brief   A public function used to initialize the projector
      \details This function does not take any parameter and is used to initialize everything that
               should be initialized. At the end, it calls the pure virtual InitializeSpecific()
               function implemented by children.
      \return  An integer reflecting the initialization status; 0 if no problem, another value
               otherwise.
    */
    int Initialize();
    /*!
      \fn      public int vProjector::Project()
      \param   int a_direction
      \param   oProjectionLine* ap_ProjectionLine
      \param   int* ap_index1
      \param   int* ap_index2
      \param   int a_nbIndices
      \brief   A function use to computed the projection elements with respect to the provided parameters
      \details This function is used to fill the provided oProjectionLine with the system matrix elements associated
               to the provided indices. The filling is done by calling one of the three projection functions (w/o TOF)
               implemented by its children. This function cannot be overloaded.
               From a collection of scanner elements indices associated to a oProjectionLine and to a projection
               direction (FORWARD or BACKWARD), it first uses the scanner to get associated cartesian coordinates
               taking compression, mean depth of interaction or POI into account, it then applies offsets (general
               one from oImageDimensionsAndQuantification), LOR displacement (associated to each event).
               The computed coordinates are embedded into the oProjectionLine which are passed
               to the different pure virtual projection functions (w/o TOF). If the number of provided indices is
               one (no compression), then the provided indices are also passed to the oProjectionLine if needed by
               a highly specified projector. In the case of compression, the indices are set to -1 in the oProjectionLine,
               while the positions and orientations are averaged.
      \return  An integer reflecting the projection status; 0 if no problem, another value otherwise.
    */
    int Project(int a_direction, oProjectionLine* ap_ProjectionLine, uint32_t* ap_index1, uint32_t* ap_index2, int a_nbIndices);


  // -------------------------------------------------------------------
  // Virtual but not pure for children
  public:
    /*!
      \fn      public INTNB vProjector::EstimateMaxNumberOfVoxelsPerLine()
      \brief   This function is used to compute and provide an estimate of the maximum number of voxels that could
               contribute to a projected line.
      \details The vProjector implementation simply returns the total image's number of voxels, but it can be
               overloaded by children to provide a better estimate in order to optimize and reduce memory requirements
               of the oProjectionLine buffers when using the FIXED_LIST_STRATEGY.
      \return  The estimate of the maximum number of voxels contributing to a line.
    */
    virtual INTNB EstimateMaxNumberOfVoxelsPerLine();


  // -------------------------------------------------------------------
  // Pure virtual public member functions that need to be implemented by children
  public:
    /*!
      \fn      public virtual int vProjector::ReadConfigurationFile() = 0
      \param   const string& a_configurationFile
      \brief   A function used to read options from a configuration file
      \details This function implements the reading of all options associated to a child projector, from
               a configuration file. It is pure virtual so is implemented by children. It checks the reading
               status but not the options values that will be checked by the CheckSpecificParameters()
               function.
      \return  An integer reflecting the reading success; 0 if success, another value otherwise.
    */
    virtual int ReadConfigurationFile(const string& a_configurationFile) = 0;
    /*!
      \fn      public virtual int vProjector::ReadOptionsList() = 0
      \param   const string& a_optionsList
      \brief   A function used to read options from a list of options
      \details This function implements the reading of all options associated to a child projector, from
               a list of options. It is pure virtual so is implemented by children. It checks the reading
               status but not the options values that will be checked by the CheckSpecificParameters()
               function.
      \return  An integer reflecting the reading success; 0 if success, another value otherwise.
    */
    virtual int ReadOptionsList(const string& a_optionsList) = 0;


  // -------------------------------------------------------------------
  // Pure virtual private member functions that need to be implemented by children
  private:
    /*!
      \fn      private virtual void vProjector::ShowHelpSpecific() = 0
      \brief   A function used to show help about the child module
      \details This function must describe what the projector does and how to use it. It describes in
               details the different parameters of the projector, and how to set them through the use
               of a configuration file or a list of options. It is pure virtual so is implemented by
               children. It is private because called by the public ShowHelp() function.
    */
    virtual void ShowHelpSpecific() = 0;
    /*!
      \fn      private virtual int vProjector::CheckSpecificParameters() = 0
      \brief   A private function used to check the parameters settings specific to the child projector
      \details This function is used to check that all parameters specific to the projector are correctly set
               within allowed values. It is called by the CheckParameters() function. It is pure virtual
               so is implemented by children.
      \return  An integer reflecting the check status; 0 if no problem, another value otherwise.
    */
    virtual int CheckSpecificParameters() = 0;
    /*!
      \fn      private virtual int vProjector::InitializeSpecific() = 0
      \brief   A private function used to initialize everything specific to the child projector
      \details This function is used to initialize everything specific to the projector that should be
               initialized. It is called by the Initialize() function. It is pure virtual so is
               implemented by children.
      \return  An integer reflecting the initialization status; 0 if no problem, another value
               otherwise.
    */
    virtual int InitializeSpecific() = 0;
    /*!
      \fn      private virtual int vProjector::ProjectWithoutTOF() = 0
      \param   int a_direction
      \param   oProjectionLine* ap_ProjectionLine
      \brief   A function to project without TOF.
      \details Projects the provided line following the provided direction, without TOF. It fills the provided
               oProjectionLine. This is a pure virtual function that must be implemented by children.
      \return  An integer reflecting the projection status; 0 if no problem, another value otherwise.
    */
    virtual int ProjectWithoutTOF( int a_direction, oProjectionLine* ap_ProjectionLine ) = 0;
    /*!
      \fn      private virtual int vProjector::ProjectTOFListmode() = 0
      \param   int a_direction
      \param   oProjectionLine* ap_ProjectionLine
      \brief   A function to project with TOF continuous information.
      \details Projects the provided line following the provided direction, with TOF described as a continuous
               measurement. It fills the provided oProjectionLine. This is a pure virtual function that must
               be implemented by children.
      \return  An integer reflecting the projection status; 0 if no problem, another value otherwise.
    */
    virtual int ProjectTOFListmode( int a_direction, oProjectionLine* ap_ProjectionLine ) = 0;
    /*!
      \fn      private virtual int vProjector::ProjectTOFHistogram() = 0
      \param   int a_direction
      \param   oProjectionLine* ap_ProjectionLine
      \brief   A function to project with TOF binned information.
      \details Projects the provided line following the provided direction, with TOF information describe as a
               histogram bin. It fills the provided oProjectionLine. This is a pure virtual function that must
               be implemented by children.
      \return  An integer reflecting the projection status; 0 if no problem, another value otherwise.
    */
    virtual int ProjectTOFHistogram( int a_direction, oProjectionLine* ap_ProjectionLine ) = 0;


  // -------------------------------------------------------------------
  // Public Get & Set functions
  public:
    /*!
      \fn      public inline void vProjector::SetVerbose()
      \param   int a_verboseLevel
      \brief   Set the verbose level.
    */
    inline void SetVerbose(int a_verbose)
           {m_verbose = a_verbose;}
    /*!
      \fn      public inline void vProjector::SetScanner()
      \param   vScanner* ap_Scanner
      \brief   Set the pointer to the scanner in use.
    */
    inline void SetScanner(vScanner* ap_Scanner)
           {mp_Scanner = ap_Scanner;}
    /*!
      \fn      public inline void vProjector::SetSensitivityMode()
      \param   bool a_sensitivityMode
      \brief   Set the sensitivity mode on or off
    */
    inline void SetSensitivityMode(bool a_sensitivityMode)
           {m_sensitivityMode = a_sensitivityMode;}
    /*!
      \fn      public inline void vProjector::SetApplyTOF()
      \param   int a_applyTOF
      \brief   Set the TOF mode
    */
    inline void SetApplyTOF(int a_applyTOF)
           {m_TOFMethod = a_applyTOF;}
    /*!
      \fn      public inline void vProjector::SetApplyPOI()
      \param   bool a_applyPOI
      \brief   Set the POI mode
    */
    inline void SetApplyPOI(bool a_applyPOI)
           {m_applyPOI = a_applyPOI;}
    /*!
      \fn      public void vProjector::SetImageDimensionsAndQuantification()
      \param   oImageDimensionsAndQuantification* ap_ImageDimensionsAndQuantification
      \brief   Set the pointer to the image dimensions in use and copy locally some often use variables.
      \return  0 if success, other value otherwise.
    */
    int SetImageDimensionsAndQuantification(oImageDimensionsAndQuantification* ap_ImageDimensionsAndQuantification);
    /*!
      \fn      inline public bool vProjector::GetCompatibilityWithSPECTAttenuationCorrection()
      \return  m_compatibleWithSPECTAttenuationCorrection
    */
    inline bool GetCompatibilityWithSPECTAttenuationCorrection()
           {return m_compatibleWithSPECTAttenuationCorrection;}
    /*!
      \fn      inline public bool vProjector::GetCompatibilityWithCompression()
      \return  m_compatibleWithCompression
    */
    inline bool GetCompatibilityWithCompression()
           {return m_compatibleWithCompression;}
    /*!
      \fn      inline public bool vProjector::SetMask()
      \param   bool* ap_mask
      \brief   Set a mask for voxels
    */
    inline void SetMask(bool* ap_mask){ mp_mask = ap_mask; m_hasMask = true;}
    /*!
      \fn      public inline int vProjector::GetTOFResolutionInMm()
      \brief   This function is used to get the TOF resolution in mm.
      \return  TOF resolution.
    */
    inline FLTNB GetTOFResolutionInMm()
           {return m_TOFResolutionInMm;}
    /*!
      \fn      public inline void vProjector::SetTOFResolutionInMm()
      \param   FLTNB a_TOFResolutionInMm
      \brief   This function is used to set the TOF resolution in use.
    */
    inline void SetTOFResolutionInMm(FLTNB a_TOFResolutionInMm)
           {m_TOFResolutionInMm = a_TOFResolutionInMm;}
    /*!
      \fn      public inline int vProjector::GetTOFMeasurementRangeInMm()
      \brief   This function is used to get the TOF measurement range in mm.
      \return  TOF measurement range.
    */
    inline FLTNB GetTOFMeasurementRangeInMm()
           {return m_TOFMeasurementRangeInMm;}
    /*!
      \fn      public inline void vProjector::SetTOFMeasurementRangeInMm()
      \param   FLTNB a_TOFMeasurementRangeInMm
      \brief   This function is used to set the TOF measurement range in mm.
    */
    inline void SetTOFMeasurementRangeInMm(FLTNB a_TOFMeasurementRangeInMm)
           {m_TOFMeasurementRangeInMm = a_TOFMeasurementRangeInMm;}
    /*!
      \fn      public inline int vProjector::GetTOFBinSizeInMm()
      \brief   This function is used to get the size in mm of a TOF bin.
      \return  TOF bin size in mm.
    */
    inline FLTNB GetTOFBinSizeInMm()
           {return m_TOFBinSizeInMm;}
    /*!
      \fn      public inline void vProjector::SetTOFBinSizeInMm()
      \param   FLTNB a_TOFBinSizeInMm
      \brief   This function is used to set the size of a TOF bin in mm.
    */
    inline void SetTOFBinSizeInMm(FLTNB a_TOFBinSizeInMm)
           {m_TOFBinSizeInMm = a_TOFBinSizeInMm;}
  // -------------------------------------------------------------------
  // Data members
  protected:
    // Few data from the oImageDimensionsAndQuantification to avoid getting them too often
    FLTNB mp_sizeVox[3];                   /*!< Local copy of the voxels' size */
    INTNB mp_nbVox[3];                     /*!< Local copy of the number of voxels */
    INTNB m_nbVoxXY;                       /*!< Local copy of the number of voxels per slice */
    FLTNB mp_halfFOV[3];                   /*!< Local copy of the FOV half dimensions */
    // Image dimensions
    oImageDimensionsAndQuantification* 
      mp_ImageDimensionsAndQuantification; /*!< Pointer to the oImageDimensionsAndQuantification object in use */
    // Scanner
    vScanner* mp_Scanner;                  /*!< Pointer to the vScanner object in use */
    // TOF related
    int m_TOFMethod;                       /*!< Integer tagging the type of TOF use */
    FLTNB m_TOFNbSigmas;                   /*!< The TOF distribution truncation factor (number of standard deviations) */
    INTNB m_TOFWeightingFcnNbSamples;         /*!< Length of the precomputed TOF weighting function */
    bool m_TOFWeightingFcnPrecomputedFlag; /*!< Choice of precomputed vs on the fly TOF weighting function */
    bool m_TOFBinProperProcessingFlag;     /*!< Take properly into account the TOF bin (the TOF Gaussian distribution is either convolved with or integrated over the TOF bin */
    HPFLTNB* mp_TOFWeightingFcn;           /*!< Precomputed TOF weighting function */
    FLTNB m_TOFBinSizeInMm;                /*!< The size in mm of a single TOF bin (either cumulative bin for histogrammed data or quantization bin for list-mode data) */
    FLTNB m_TOFMeasurementRangeInMm;       /*!< Total span of TOF measurements (in mm) */
    FLTNB m_TOFResolutionInMm;             /*!< Resolution (FWHM) of the time-of-flight Gaussian distribution (in mm) */
    FLTNB m_TOFGaussianNormCoef;           /*!< The normalization coefficient for the TOF Gaussian distribution (to ensure integral=1) */
    FLTNB m_TOFPrecomputedSamplingFactor;  /*!< Sampling factor for oversampling the precomputed TOF weighting function */
    // Voxel mask to remove voxels from projection
    bool* mp_mask;                         /*!< Mask for voxels: only true voxels will be taken into account for projector coefficients computation (currently 3D - XYZ)*/
    bool m_hasMask;                        /*!< Is there a mask for voxels */
    // Flag for POI
    bool m_applyPOI;                       /*!< Boolean that says if we apply POI info or not */
    // Flag for sensitivity computation
    bool m_sensitivityMode;                /*!< Boolean that says if we are computing the global sensitivity or not */
    // Flag that says if the projector is compatible with SPECT attenuation correction
    bool m_compatibleWithSPECTAttenuationCorrection; /*!< Boolean that says if the projector is compatible with SPECT attenuation correction */
    // Flag that says if the projection is compatible with compression
    // (i.e. when the crystal indices in the projection line will be -1 as they are not unique for a given event)
    bool m_compatibleWithCompression; /*!< Boolean that says if the projector is compatible with compression */
    // Verbosity
    int m_verbose;                         /*!< The verbose level */
    // Has been checked ?
    bool m_checked;                        /*!< Boolean that says if the parameters were checked or not */
    // Has been initialized ?
    bool m_initialized;                    /*!< Boolean that says if the projector was initialized or not */
};

// ----------------------------------------------------------------------
// Part of code that manages the auto declaration of children classes
// ----------------------------------------------------------------------

// Macro for the function that creates the object
#define FUNCTION_PROJECTOR(CLASS) \
  static vProjector *make_projector() { return new CLASS(); };

// Macro for the class that links the appropriate function to the map of objects
#define CLASS_PROJECTOR(NAME,CLASS)                                                           \
  class NAME##ProjectorCreator                                                                \
  {                                                                                           \
    public:                                                                                   \
      NAME##ProjectorCreator()                                                                \
        { sAddonManager::GetInstance()->mp_listOfProjectors[#NAME] = CLASS::make_projector; } \
  };                                                                                          \
  static NAME##ProjectorCreator ProjectorCreator##NAME;

#endif

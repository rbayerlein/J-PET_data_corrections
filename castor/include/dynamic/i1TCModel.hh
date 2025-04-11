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
  \brief    Declaration of class i1TCModel
*/


#ifndef I1TCMODEL_HH
#define I1TCMODEL_HH 1


/**
 * @defgroup 1CPT_METHOD 1cpt parameters optimization method
 *
 *    \brief Optimisation method (LS, NNLS) \n
 *           Defined in i1TCModel.hh
 * @{
 */
/** NNLS method (=0, default) */
#define METHOD_1CPT_NNLS 0
/** LS method (=1) */
#define METHOD_1CPT_LS 1
/** BF method (=2) */
#define METHOD_1CPT_BF 2
/** @} */


/**
 * @defgroup INTEGRATION_METHOD 
 *
 *    \brief Integration method for TACs (WPO, Trap) \n
 *           Defined in i1TCModel.hh
 * @{
 */
/** (W)eighed (P)arabola (O)verlapping method (=0, default) */
#define METHOD_INT_WPO 0
/** Trapezoid method (=1) */
#define METHOD_INT_TRAP 1
/** @} */


#include "vDynamicModel.hh"
#include "sAddonManager.hh"


/*!
  \class   i1TCModel
  \brief   This class implements a 1 compartiment model, to model kinetics of radiotracers such as radiowater
*/
class i1TCModel : public vDynamicModel
{
  // -----------------------------------------------------------------------------------------
  // Constructor & Destructor
  public:
    /*!
      \fn      i1TCModel::i1TCModel
      \brief   Constructor of i1TCModel. Simply set all data members to default values.
    */
    i1TCModel();
    /*!
      \fn      i1TCModel::~i1TCModel
      \brief   Destructor of i1TCModel
    */
    ~i1TCModel();


  // -----------------------------------------------------------------------------------------
  // Public member functions related to the initialization of the model
  public:
    // Function for automatic insertion (put the class name as the parameters and do not add semi-colon at the end of the line)
    FUNCTION_DYNAMICMODEL(i1TCModel)
    /*!
      \fn      i1TCModel::CheckSpecificParameters
      \brief   This function is used to check whether all member variables
               have been correctly initialized or not.
      \return  0 if success, positive value otherwise.
    */
    int CheckSpecificParameters();
    /*!
      \fn      i1TCModel::ReadAndCheckConfigurationFileSpecific
      \brief   This function is used to read options from a configuration file.
      \return  0 if success, other value otherwise.
    */
    int ReadAndCheckConfigurationFileSpecific();
    /*!
      \fn      i1TCModel::ReadAndCheckOptionsList
      \param   const string& a_optionsList : a list of parameters separated by commas
      \brief   This function is used to read parameters from a string.
      \return  0 if success, other value otherwise.
    */
    int ReadAndCheckOptionsList(string a_listOptions);
    /*!
      \fn      i1TCModel::InitializeSpecific
      \brief   This function is used to initialize parametric images and basis functions
      \todo    Read Interfile for parametric images initialization
      \return  0 if success, other value otherwise.
    */
    int InitializeSpecific();
    /*!
      \fn      i1TCModel::ShowHelp
      \brief   Print out specific help about the implementation of the model
               model and its initialization
    */
    void ShowHelpModelSpecific();


  // -------------------------------------------------------------------
  // Public member functions called by the main iterative algorithm class
    /*!
      \fn      i1TCModel::EstimateModelParameters
      \param   ap_ImageS : pointer to the ImageSpace
      \param   a_ite : index of the actual iteration (not used)
      \param   a_sset : index of the actual subset (not used)
      \brief   Estimate K1, k2, Va parametric images
      \return  0 if success, other value otherwise.
    */
    int EstimateModelParameters(oImageSpace* ap_Image, int a_ite, int a_sset);

    /*!
      \fn      i1TCModel::EstimateImageWithModel
      \param   ap_ImageS : pointer to the ImageSpace
      \param   a_ite : index of the actual iteration (not used)
      \param   a_sset : index of the actual subset (not used)
      \brief   Estimate image using model parametric images and basis functions
      \return  0 if success, other value otherwise.
    */
    int EstimateImageWithModel(oImageSpace* ap_Image, int a_ite, int a_sset);


  // -------------------------------------------------------------------
  // Private member functions 
  private:
    /*!
      \fn      i1TCModel::EstimateModelParametersWithLS
      \param   ap_ImageS : pointer to the ImageSpace
      \brief   Estimate K1, k2, Va parametric images using LS
      \return  0 if success, other value otherwise.
    */
    int EstimateModelParametersWithLS(oImageSpace* ap_ImageS);
    
    /*!
      \fn      i1TCModel::EstimateModelParametersWithNNLS
      \param   ap_ImageS : pointer to the ImageSpace
      \brief   Estimate K1, k2, Va parametric images using NNLS
      \return  0 if success, other value otherwise.
    */
    int EstimateModelParametersWithNNLS(oImageSpace* ap_ImageS);

    /*!
      \fn      i1TCModel::EstimateModelParametersWithBF
      \param   ap_ImageS : pointer to the ImageSpace
      \brief   Estimate K1, k2, Va parametric images using basis functions approach
      \return  0 if success, other value otherwise.
    */
    int EstimateModelParametersWithBF(oImageSpace* ap_ImageS);
    
    /*!
      \fn RRLS
      \param a_nP      : number of parameters (P)
      \param a_nT      : number of data samples (T)
      \param a2p_model : matrix containing the model functions (TxP elts)
      \param ap_data   : data vector (T elts)
      \param ap_w      : vector containing the weights (T elts)
      \param ap_result : vector recovering the estimated parameters (P elts)
      \brief   Non-linear least square estimator with Ridge-Regression
      \details This function will estimate a set of parameters (O) by 
               minimizing the weighted residual sum of squares (WRSS) 
               difference between the data and the model function
               Ô = [ X'*W*X + t.Rw ]-1 [ X'*W*y ]  
                 + [ X'*W*X + t.Rw ]-1 [ t.Rw*Rm ]
               y = data vector
               X = model matrix
               W = weights vector
               t = Ridge constant
               Rw= Ridge weights
               Rm= Ridge means
               
      \return 0 if success, other value otherwise.
    */
    int RRLS(uint16_t  a_nP,
             uint16_t  a_nT,
               FLTNB **a2p_model,
               FLTNB  *ap_data,
               FLTNB  *ap_w,
               FLTNB  *ap_result
            );
                    
    /*!
      \fn LS
      \param a_nP      : number of parameters (P)
      \param a_nT      : number of data samples (T)
      \param a2p_model : matrix containing the model functions (TxP elts)
      \param ap_data   : data vector (T elts)
      \param ap_w      : vector containing the weights (T elts)
      \param ap_result : vector recovering the estimated parameters (P elts)
      \brief   Non-linear least square estimator 
      \details This function will estimate a set of parameters (O) by 
               minimizing the weighted residual sum of squares (WRSS) 
               difference between the data and the model function
               Ô = [ X'WX ]-1 X'Wy 
               y = data vector
               X = model matrix
               W = weights vector
               
      \return 0 if success, other value otherwise.
    */
    int LS(uint16_t  a_nP,
           uint16_t  a_nT,
             FLTNB **a2p_model,
             FLTNB  *ap_data,
             FLTNB  *ap_w,
             FLTNB  *ap_result
          );
                        

    /*!
      \fn      WPOinc
      \param   a_time : time point for which the integral will be computed
      \param   tac : tac value at time t
      \param   b_tac : tac value at time t-1
      \param   bb_tac : tac value at time t-2
      \param   n_tac : tac value at time t+1
      \brief   Estimate the next integration value for a specific time point of a tac using WPO
      \return  The estimated value of the integral for this time point
    */
    FLTNB WPOinc(uint32_t a_time, FLTNB tac, FLTNB b_tac, FLTNB bb_tac, FLTNB n_tac);
    
    /*!
      \fn      IntegrateTAC
      \param   ap_tac : vector with nb time frames elements containing the time activity curve from which integral will be computed
      \param   ap_citac : (returned) vector with nb time frames elements recovering the cumulative integral for each time point t 
      \param   a_th : thread index
      \brief   Call one of the TAC integration method
      \return  0 if success, positive value otherwise
    */
    int IntegrateTAC(FLTNB* ap_tac, FLTNB* ap_citac, int a_th);
    
    /*!
      \fn      WPO
      \param   ap_tac : vector with nb time frames elements containing the time activity curve from which integral will be computed
      \param   ap_citac : (returned) vector with nb time frames elements recovering the cumulative integral for each time point t 
      \param   a_th : thread index
      \brief   return integral for each time point t (WPO_S), and cumulative integral (cumWPO_S )
      \return  0 if success, positive value otherwise
    */
    int WPO(FLTNB* ap_tac, FLTNB* ap_citac, int a_th);
    /*!
      \fn      Trapz
      \param   ap_tac : vector with nb time frames elements containing the time activity curve from which integral will be computed
      \param   ap_citac : (returned) vector with nb time frames elements recovering the cumulative integral for each time point t 
      \brief   return integral for each time point t (WPO_S), and cumulative integral (cumWPO_S )
      \return  0 if success, positive value otherwise
    */
    int Trapz(FLTNB* ap_tac, FLTNB* ap_citac);
    
    
    
  // -----------------------------------------------------------------------------------------
  // Data members
  protected:
    bool        m_savePImgFlag;           /*!<Flag indicating if parametric images should be written on disk (default=true) */
    int         m_OptimisationMethodFlag; /*!<Flag containing the optimization method */
    
    FLTNB**     m2p_ct;                   /*!<Vector containing voxel TACs (multithreaded) */
    FLTNB**     m2p_cti;                  /*!<Vector containing the integral of voxel TACs (multithreaded) */
    FLTNB*      mp_w;                     /*!<Vector containing the weights for NNLS estimation */
    FLTNB***    m3p_nnlsA;                /*!<2D coefficient matrix for NNLS estimation (multithreaded)*/
    FLTNB**     m2p_nnlsB;                /*!<1D vector for NNLS estimation, containing the solution (multithreaded)*/
    FLTNB**     m2p_nnlsMat;              /*!<1D working vector for NNLS estimation (multithreaded)*/
    uint16_t    m_nnlsN=3;                /*!<Number of parameters in NNLS/LS estimation */

    FLTNB*      mp_parUpperBounds;        /*!<Vector containing parameters upper bound */
    FLTNB*      mp_parLowerBounds;        /*!<Vector containing parameters lower bound */
    FLTNB       m_RRcst;                  /*!<Regularization parameter for ridge regression */
    bool        m_ridgeRegressionFlag;    /*!<Flag indicating if ridge regression is enabled (default=false) */
    
    FLTNB*      mp_DT2;                   /*!<Vector containing the half frame time duration (for integral computation) */
    // WPO
    FLTNB*      mp_wpoQ;
    FLTNB*      mp_wpoA;
    FLTNB**     m2p_wpoP;
    FLTNB**     m2p_wpoFD;
    FLTNB**     m2p_wpoBD;
    bool        m_intMethodFlag;          /*!<Flag indicating the method for TACs integral computation (default=true) */
        
    // Least-Square matrices for LS estimation
    oMatrix**   mp_Y;       /*!<Data  Matrix (threaded)*/
    oMatrix**   mp_X;       /*!<Model Matrix (threaded)*/
    oMatrix**   mp_Xt;      /*!<Matrix for computation (th)*/
    oMatrix**   mp_XtX;     /*!<Matrix for computation (th)*/
    oMatrix**   mp_LSnum;   /*!<Matrix for computation (th)*/
    oMatrix**   mp_LSden;   /*!<Matrix for computation (th)*/
    oMatrix**   mp_Theta;   /*!<Estimated parameters Matrix (th)*/

    // Least-Square with Ridge Regression variables
    oMatrix*    mp_RRm;     /*!<Ridge means   */
    oMatrix*    mp_RRw;     /*!<Ridge weights */
    oMatrix**   mp_RRnum;   /*!<Ridge numerator for computation (th)*/
    
    FLTNB*      mp_VaImage; /*!< Optional input image containing the blood volume */
};




// Class for automatic insertion (set here the visible dynamic model's name, put the class name as the parameters and do not add semi-colon at the end of the line)
CLASS_DYNAMICMODEL(_1TCM,i1TCModel)

#endif

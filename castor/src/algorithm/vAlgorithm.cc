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
  \ingroup algorithm
  \brief Implementation of class vAlgorithm
*/

#include "gVariables.hh"
#include "vAlgorithm.hh"
#include "iEventHistoCT.hh"

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

vAlgorithm::vAlgorithm()
{
  // Default number of iterations and subsets (1 for each)
  m_nbIterations = 1;
  mp_nbSubsets = (int*)malloc(1*sizeof(int));
  mp_nbSubsets[0] = 1;
  mp_outputIterations = NULL; 
  // set all members to default values
  m_verbose = -1;
  m_flagGPU = false;
  mp_ID = NULL;
  m2p_DataFile= NULL;
  mp_ProjectorManager = NULL;
  mp_OptimizerManager = NULL;
  mp_DeformationManager = NULL;
  mp_DynamicModelManager = NULL;
  mp_ImageSpace = NULL;
  mp_ImageConvolverManager = NULL;
  mp_ImageProcessingManager = NULL;
  m_nbBeds = -1;
  m_pathToInitialImg = "";
  m_pathToSensitivityImg = "";
  m_pathToMaskImg = "";
  m_saveSensitivityHistoFlag = false;
  m_saveImageAfterSubsets = false;
  m_saveDynamicBasisCoefficients = false;

}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

vAlgorithm::~vAlgorithm()
{
  // Delete some tables
  if (mp_nbSubsets!=NULL) free(mp_nbSubsets);
  if (mp_outputIterations!=NULL) free(mp_outputIterations);
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vAlgorithm::SetNbIterationsAndSubsets(const string& a_nbIterationsSubsets)
{
  // If the string is empty, then just keep the defaults
  if (a_nbIterationsSubsets=="") return 0;
  // Otherwise, reset the number of iterations to 0
  m_nbIterations = 0;

  // Copy the string
  string buf = a_nbIterationsSubsets;

  // Loop to search for commas
  size_t pos_comma;
  while ((pos_comma=buf.find_first_of(","))!=string::npos)
  {
    // Get the substring before the comma
    string sub_buf = buf.substr(0,pos_comma);
    // Search for columns
    size_t pos_column = sub_buf.find_first_of(":");
    if (pos_column==string::npos || pos_column==0)
    {
      Cerr("***** vAlgorithm::SetNbIterationsAndSubsets() -> Syntax problem in number of iterations and subsets !" << endl);
      return 1;
    }
    int iter = atoi(  sub_buf.substr(0,pos_column).c_str()  );
    int subs = atoi(  sub_buf.substr(pos_column+1).c_str()  );
    mp_nbSubsets = (int*)realloc(mp_nbSubsets,(m_nbIterations+iter)*sizeof(int));
    for (int it=0; it<iter; it++) mp_nbSubsets[m_nbIterations+it] = subs;
    m_nbIterations += iter;
    buf = buf.substr(pos_comma+1);
  }

  // Last couple of iterations:subsets
  size_t pos_column = buf.find_first_of(":");
  if (pos_column==string::npos || pos_column==0)
  {
    Cerr("***** vAlgorithm::SetNbIterationsAndSubsets() -> Syntax problem in number of iterations and subsets !" << endl);
    return 1;
  }
  int iter = atoi(  buf.substr(0,pos_column).c_str()  );
  int subs = atoi(  buf.substr(pos_column+1).c_str()  );
  mp_nbSubsets = (int*)realloc(mp_nbSubsets,(m_nbIterations+iter)*sizeof(int));
  for (int it=0; it<iter; it++) mp_nbSubsets[m_nbIterations+it] = subs;
  m_nbIterations += iter;

  if (m_verbose>=3) 
  {
    Cout("vAlgorithm::SetNbIterationsAndSubsets() ->  Selected numbers of subsets for each iteration:" << endl); 
    Cout("  Iteration / number of subsets "<< endl); 
    for (int it=0 ; it<m_nbIterations ; it++) Cout("    " << it+1 << "  /  " <<  mp_nbSubsets[it] << endl); 
  }
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vAlgorithm::SetOutputIterations(const string& a_outputIterations)
{
  if (m_verbose>=2) Cout("vAlgorithm::SetOutputIterations ..." << endl); 
    
  // Allocate the output iterations boolean table
  mp_outputIterations = (bool*)malloc(m_nbIterations*sizeof(bool));
  for (int it=0; it<m_nbIterations; it++) mp_outputIterations[it] = false;

  // If the list is empty, we save all iterations by default
  if (a_outputIterations=="")
  {
    for (int it=0; it<m_nbIterations; it++) mp_outputIterations[it] = true;
    return 0;
  }

  // Copy the string
  string buf = a_outputIterations;

  // Loop to search for commas
  size_t pos_comma;
  int current_iteration = -1;
  while ((pos_comma=buf.find_first_of(","))!=string::npos)
  {
    // Get the substring before the comma
    string sub_buf = buf.substr(0,pos_comma);
    // Search for columns
    size_t pos_column = sub_buf.find_first_of(":");
    if (pos_column==string::npos || pos_column==0)
    {
      Cerr("***** vAlgorithm::SetOutputIterations() -> Syntax problem in output iteration numbers !" << endl);
      return 1;
    }
    int step_iteration = atoi(  sub_buf.substr(0,pos_column).c_str()  );
    if (step_iteration<=0)
    {
      Cerr("***** vAlgorithm::SetOutputIterations() -> Iteration step must be strictly positive (here it is " << step_iteration << ") !" << endl);
      return 1;
    }
    int stop_iteration = atoi(  sub_buf.substr(pos_column+1).c_str()  );
    if (stop_iteration<=current_iteration)
    {
      Cerr("***** oIterationAlgorithm::SetOutputIterations() -> Output iteration number must be increasing through provided couples !" << endl);
      return 1;
    }
    for (int it=current_iteration+step_iteration; it<stop_iteration; it+=step_iteration) if (it<m_nbIterations) mp_outputIterations[it] = true;
    current_iteration = stop_iteration-1;
    buf = buf.substr(pos_comma+1);
  }

  // Last couple of iterations:subsets
  size_t pos_column = buf.find_first_of(":");
  if (pos_column==string::npos || pos_column==0)
  {
    // Special case if -1 is provided, it means we save the last iteration
    if (buf=="-1")
    {
      mp_outputIterations[m_nbIterations-1] = true;
      // We directly exist here
      return 0;
    }
    else
    {
      Cerr("***** vAlgorithm::SetOutputIterations() -> Syntax problem in output iteration numbers !" << endl);
      return 1;
    }
  }
  int step_iteration = atoi(  buf.substr(0,pos_column).c_str()  );
  if (step_iteration<=0)
  {
    Cerr("***** vAlgorithm::SetOutputIterations() -> Iteration step must be strictly positive (here it is " << step_iteration << ") !" << endl);
    return 1;
  }
  int stop_iteration = atoi(  buf.substr(pos_column+1).c_str()  );
  if (stop_iteration<=current_iteration)
  {
    Cerr("***** oIterationAlgorithm::SetOutputIterations() -> Output iteration number must be increasing through provided couples !" << endl);
    return 1;
  }
  for (int it=current_iteration+step_iteration; it<stop_iteration; it+=step_iteration) if (it<m_nbIterations) mp_outputIterations[it] = true;

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vAlgorithm::Run()
{
  if(m_verbose>=2) Cout("vAlgorithm::Iterate ..."<< endl); 
    
  #ifdef CASTOR_GPU
  if (m_flagGPU) 
    return RunGPU();
  else 
    return RunCPU();
  #else
  return RunCPU();
  #endif
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vAlgorithm::RunCPU()
{
// Here is the scheme of all functions call related to each step of the iterative reconstruction:
//   StepBeforeIterationLoop()
//   / Loop on iterations
//   | StepBeforeSubsetLoop(iteration)
//   |  / Loop on subsets
//   |  | StepPreProcessInsideSubsetLoop(iteration,subset)
//   |  | / Loop on bed positions
//   |  | | StepInnerLoopInsideSubsetLoop(iteration,subset,bed)
//   |  | StepPostProcessInsideSubsetLoop(iteration,subset)
//   |  StepAfterSubsetLoop(iteration)
//   StepAfterIterationLoop()

  // Synchronize MPI processes
  #ifdef CASTOR_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  #endif

  // Verbose
  if (m_verbose>=1) Cout("vAlgorithm::IterateCPU() -> Start algorithm for " << m_nbIterations << " iterations" << endl);

  // Initial clock for execution time
  clock_t clock_start_whole = clock();
  time_t time_start_whole = time(NULL);

  // Call before iteration function
  // Call before iteration function
  if (StepBeforeIterationLoop())
  {
    Cerr("***** vAlgorithm::IterateCPU() -> A problem occurred while calling StepBeforeIterationLoop() function !" << endl);
    return 1;
  }

  // Loop on iterations
  for (int iteration=0 ; iteration<m_nbIterations ; iteration++)
  {

    // Call before subset function
    if (StepBeforeSubsetLoop(iteration))
    {
      Cerr("***** vAlgorithm::IterateCPU() -> A problem occurred while calling StepBeforeSubsetLoop() function !" << endl);
      return 1;
    }

    // Loop on subsets
    for (int subset=0 ; subset<mp_nbSubsets[iteration] ; subset++)
    {
      // Verbose
      if (m_verbose>=1) 
      {
        Cout("vAlgorithm::IterateCPU() -> Start iteration " << iteration+1 << "/" << m_nbIterations 
          << " subset " << subset+1 << "/" << mp_nbSubsets[iteration] << endl);
      }
      
      // Clock start for subset execution time
      clock_t clock_start_subset = clock();
      time_t time_start_subset = time(NULL);

      // Call pre-process inside subset loop function
      if (StepPreProcessInsideSubsetLoop(iteration,subset))
      {
        Cerr("***** vAlgorithm::IterateCPU() -> A problem occurred while calling StepPreProcessInsideSubsetLoop() function !" << endl);
        return 1;
      }

      // Synchronize MPI processes
      #ifdef CASTOR_MPI
      MPI_Barrier(MPI_COMM_WORLD);
      #endif

      // Loop on bed positions
      for (int bed=0 ; bed<m_nbBeds ; bed++)
      {
        // Call the inner loop on events inside subset loop function
        if (StepInnerLoopInsideSubsetLoop(iteration,subset,bed))
        {
          Cerr("***** vAlgorithm::IterateCPU() -> A problem occurred while calling StepInnerLoopInsideSubsetLoop() function !" << endl);
          return 1;
        }
      } // End of beds loop

      // Synchronize MPI processes
      #ifdef CASTOR_MPI
      MPI_Barrier(MPI_COMM_WORLD);
      #endif

      // Call post-process inside subset loop function
      if (StepPostProcessInsideSubsetLoop(iteration,subset))
      {
        Cerr("***** vAlgorithm::IterateCPU() -> A problem occurred while calling StepPostProcessInsideSubsetLoop() function !" << endl);
        return 1;
      }

      // Clock stop for subset execution time
      clock_t clock_stop_subset = clock();
      time_t time_stop_subset = time(NULL);
      if (m_verbose>=2) Cout ("vAlgorithm::IterateCPU() -> Time spent for subset " << subset+1 << " | User: " << time_stop_subset-time_start_subset
                           << " sec | CPU: " << (clock_stop_subset-clock_start_subset)/((FLTNB)CLOCKS_PER_SEC) << " sec" << endl);


    } // End of subsets loop

              
    // Call after subset function
    if (StepAfterSubsetLoop(iteration))
    {
      Cerr("***** vAlgorithm::IterateCPU() -> A problem occurred while calling StepAfterSubsetLoop() function !" << endl);
      return 1;
    }


  } // End of iterations loop

  //  Call after iteration function
  if (StepAfterIterationLoop())
  {
    Cerr("***** vAlgorithm::IterateCPU() -> A problem occurred while calling StepAfterIterationLoop() function !" << endl);
    return 1;
  }

  // Final clock for execution time
  clock_t clock_stop_whole = clock();
  time_t time_stop_whole = time(NULL);
  if (m_verbose>=1) Cout("vAlgorithm::IterateCPU() -> Total time spent | User: " << time_stop_whole-time_start_whole 
                      << " sec | CPU: " << (clock_stop_whole-clock_start_whole)/((FLTNB)CLOCKS_PER_SEC) << " sec" << endl);
  

  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vAlgorithm::StepBeforeIterationLoop()
{
  if (m_verbose>=2) Cout("vAlgorithm::StepBeforeIterationLoop ... " << endl);

  // Check that the image space is already allocated
  if (mp_ImageSpace==NULL || !mp_ImageSpace->Checked())
  {
    Cerr("***** vAlgorithm::StepBeforeIterationLoop() -> Image space has not been created or was not checked !" << endl);
    return 1;
  }

  // Instantiate all the required images from the oImageSpace
  mp_ImageSpace->InstantiateImage();
  mp_ImageSpace->InstantiateForwardImage();
  mp_ImageSpace->InstantiateSensitivityImage(m_pathToSensitivityImg);
  mp_ImageSpace->InstantiateOutputImage();
  mp_ImageSpace->InstantiateVisitedVoxelsImage();

  if (mp_ImageSpace->InitAttenuationImage(m_pathToAtnImg) )
  {
    Cerr("***** vAlgorithm::StepBeforeIterationLoop()-> Error during attenuation image initialization !" << endl);  
    return 1;
  }

  if (mp_ImageSpace->InitSensitivityImage(m_pathToSensitivityImg))
  {
    Cerr("***** vAlgorithm::StepBeforeIterationLoop() -> An error occurred while initializing the sensitivity image !" << endl);
    return 1;
  }

  if (mp_ImageSpace->InitMultiModalImage(m_pathToMultiModalImg))
  {
    Cerr("***** vAlgorithm::StepBeforeIterationLoop()-> Error during multimodal image initialization !" << endl);
    return 1;
  }

  if (mp_ImageSpace->InitMaskImage(m_pathToMaskImg))
  {
    Cerr("***** vAlgorithm::StepBeforeIterationLoop()-> Error during mask image initialization !" << endl);
    return 1;
  }

  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vAlgorithm::StepBeforeSubsetLoop(int a_iteration)
{
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vAlgorithm::StepPreProcessInsideSubsetLoop(int a_iteration, int a_subset)
{
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vAlgorithm::StepPostProcessInsideSubsetLoop(int a_iteration, int a_subset)
{
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vAlgorithm::StepAfterSubsetLoop(int a_iteration)
{
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vAlgorithm::StepAfterIterationLoop()
{
  if (m_verbose>=2) Cout("vAlgorithm::StepAfterIterationLoop ... " << endl);
  
  // Deallocate everything
  mp_ImageSpace->DeallocateSensitivityImage();
  mp_ImageSpace->DeallocateMultiModalImage();
  mp_ImageSpace->DeallocateMaskImage();
  mp_ImageSpace->DeallocateForwardImage();
  mp_ImageSpace->DeallocateImage();
  mp_ImageSpace->DeallocateOutputImage();
  mp_ImageSpace->DeallocateVisitedVoxelsImage();

  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vAlgorithm::StepImageOutput(int a_iteration, int a_subset)
{
  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int vAlgorithm::InitSpecificOptions(string a_specificOptions)
{
  (void)a_specificOptions; // avoid 'unused parameter' compil. warnings
  // Normal end
  return 0;
}

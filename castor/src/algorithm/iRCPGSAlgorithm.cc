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

#include "iRCPGSAlgorithm.hh"
// TODO replace asserts with a better error management, propagate errors for a clean termination
#include <cassert>
#ifdef _WIN32
// Avoid compilation errors due to mix up between std::min()/max() and
// min max macros
#undef min
#undef max
#endif

iRCPGSAlgorithm::iRCPGSAlgorithm() : vAlgorithm()
{
  m_ddCRP = -1;
  m_backprojection = -1;
  m_neighbourhood = -1;
  m_gammaShape = -1.;
  m_gammaRate = -1.;
  m_ddcrpAlpha = -1.;
  m_ddcrpLogAlpha = 0.;
  m_multiModalLag = 0;
  m_ddcrpAlphaIncrement = 0.;
  m_meanClusterVolumeMin = -1.;
  m_meanClusterVolumeMax = -1.;
  m_currentMeanClusterVolume = 0.;
  mp_multiModalNoiseSigma = NULL;
  mp_multiModalParam = NULL;
  mp_voxelClusterMapping = NULL;
  mp_clusterValues = NULL;
  mp_clusterN = NULL;
  mp_clusterSensitivity = NULL;
  mpp_clusterMultiModal = NULL;
  mp_clusterCount = NULL;
  mp_nextLink = NULL;
  mp_listVoxelIndices = NULL;
  mp_listEventsIndices = NULL;
  m4p_EventsBackprojectionTrace = NULL;
  mp_listRelevantVoxelIndices = NULL;
  temp_count_multimodal = NULL;
  mp_permanentBackwardImage = NULL;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

iRCPGSAlgorithm::~iRCPGSAlgorithm()
{
  delete[] mp_clusterValues;
  mp_clusterValues = NULL;
  delete[] mp_voxelClusterMapping;
  mp_voxelClusterMapping = NULL;
  delete[] mp_nextLink;
  mp_nextLink = NULL;
  delete[] mp_listVoxelIndices;
  mp_listVoxelIndices = NULL;
  delete[] mp_clusterN;
  mp_clusterN = NULL;
  delete[] mp_clusterSensitivity;
  mp_clusterSensitivity = NULL;
  delete[] mpv_parentLinks;
  mpv_parentLinks = NULL;
  
  // optional variables, not always initialized and used
  if (mp_listEventsIndices!=NULL)
  {
    delete[] mp_listEventsIndices;
    mp_listEventsIndices = NULL;
  }
  if (mp_listRelevantVoxelIndices!=NULL)
  {
    delete[] mp_listRelevantVoxelIndices;
    mp_listRelevantVoxelIndices = NULL;
  }

  if (mpp_clusterMultiModal!=NULL)
  {
    for (int mmnb=0; mmnb<mp_ID->GetNbMultiModalImages(); mmnb++)
    {
      if (mpp_clusterMultiModal[mmnb]!=NULL) delete[] mpp_clusterMultiModal[mmnb];
    }
    delete[] mpp_clusterMultiModal;
    mpp_clusterMultiModal = NULL;
  }
  if (mp_clusterCount!=NULL)
  {
    delete[] mp_clusterCount;
    mp_clusterCount = NULL;
  }
  if (temp_count_multimodal!=NULL)
  {
    delete[] temp_count_multimodal;
    temp_count_multimodal = NULL;
  }
  if (mp_multiModalNoiseSigma!=NULL)
  {
    delete[] mp_multiModalNoiseSigma;
    mp_multiModalNoiseSigma = NULL;
  }
  if (mp_multiModalParam!=NULL)
  {
    delete[] mp_multiModalParam;
    mp_multiModalParam = NULL;
  }
  // delete helper variables used in backprojection modes 1 and 2
  if (m_backprojection==1 || m_backprojection==2)
  {
    for (int bed=0 ; bed<m_nbBeds ; bed++)
    {
      int64_t index;
      #pragma omp parallel for private(index) schedule(static,1)
      for (index=0 ; index<m2p_DataFile[bed]->GetSize() ; index++)
      {
        // Get the thread index
        int th = 0;
        #ifdef CASTOR_OMP
        th = omp_get_thread_num();
        #endif

        vEvent* event = m2p_DataFile[bed]->GetEvent(index, th);
        oProjectionLine *line = mp_ProjectorManager->ComputeProjectionLine(event, th);

        for (int b=0; b<line->GetNbTOFBins(); b++)
        {
          delete[] m4p_EventsBackprojectionTrace[bed][index][b];
        }
        delete[] m4p_EventsBackprojectionTrace[bed][index];
      }
      delete[] m4p_EventsBackprojectionTrace[bed];
    }
    delete[] m4p_EventsBackprojectionTrace;
    m4p_EventsBackprojectionTrace = NULL;
  }
}

// -----------------------------------------------------------------------------------------------------------------------------------------

int iRCPGSAlgorithm::InitSpecificOptions(string a_specificOptions)
{
  if (m_verbose>=2) Cout("iRCPGSAlgorithm::InitSpecificOptions() "<<a_specificOptions<<endl);

  // not great to do this allocation here TODO
  mp_multiModalNoiseSigma = new HPFLTNB[mp_ID->GetNbMultiModalImages()];
  for (int mmnb=0; mmnb<mp_ID->GetNbMultiModalImages(); mmnb++) mp_multiModalNoiseSigma[mmnb] = -1.;

  // TODO code copied from elsewhere, kind of a duplicate, so should be moved to more generic functions, to vAlgorithm
  // Search for a colon ":", this indicates that a configuration file is provided after the algorithm name
  size_t colon = a_specificOptions.find_first_of(":");
  size_t comma = a_specificOptions.find_first_of(",");
  
  string name = "";
  string list_options = "";
  string file_options = "";

  // Case 1: we have a colon : config file 
  if (colon!=string::npos)
  {
    // Get the name before the colon
    name = a_specificOptions.substr(0,colon);
    // Get the configuration file after the colon
    file_options = a_specificOptions.substr(colon+1);
    // List of options is empty
    list_options = "";
  }
  // Case 2: we have a comma : options list
  else if (comma!=string::npos)
  {
    // Get the name before the first comma
    name = a_specificOptions.substr(0,comma);
    // Get the list of options after the first comma
    list_options = a_specificOptions.substr(comma+1);
    // Configuration file is empty
    file_options = "";
  }
  // Case 3: no colon and no comma : default config file
  else
  {
    // Get the algorithm name
    name = a_specificOptions;
    // List of options is empty
    list_options = "";
    // Build the default configuration file
    sOutputManager* p_output_manager = sOutputManager::GetInstance();
    file_options = p_output_manager->GetPathToConfigDir() + "/algorithm/" + name + ".conf";
  }

  // TODO move this check elsewhere
  if (name!="RCPGS")
  {
    Cerr("***** iRCPGSAlgorithm::InitSpecificOptions() -> The given algorithm name is not RCPGS !" << endl);
    return 1;
  }
  
  // process options list, formatted as 'paramName=paramValue,paramName=paramValue,etc.'
  if(!list_options.empty())
  {
    vector<string> option_name;
    vector<string> option_value;

    // extract all parameter name-value pairs
    size_t pos_comma;
    while ((pos_comma=list_options.find_first_of(","))!=string::npos)
    {
      // Get the substring before the comma
      string sub_buf = list_options.substr(0,pos_comma);

      size_t pos_equal = sub_buf.find_first_of("=");
      if (pos_equal==string::npos || pos_equal==0)
      {
        Cerr("***** iRCPGSAlgorithm::InitSpecificOptions() -> Syntax problem in algorithm parameters !" << endl);
        return 1;
      }
      option_name.push_back(sub_buf.substr(0,pos_equal));
      option_value.push_back(sub_buf.substr(pos_equal+1));
      list_options = list_options.substr(pos_comma+1);
    }

    // last or single name-value pair
    if (list_options.length()>0)
    {
      size_t pos_equal = list_options.find_first_of("=");
      if (pos_equal==string::npos || pos_equal==0)
      {
        Cerr("***** iRCPGSAlgorithm::InitSpecificOptions() -> Syntax problem in algorithm parameters !" << endl);
        return 1;
      }
      option_name.push_back(list_options.substr(0,pos_equal));
      option_value.push_back(list_options.substr(pos_equal+1));
    }

    // temporary variable for noise std deviations for multimodal images
    vector<HPFLTNB> temp_multiModalNoiseSigma;

    // read extracted options
    for(size_t v=0;v<option_name.size();v++)
    {
      if (option_name.at(v)=="ddCRP")
      {
        int option[1];
        // Read them
        if (ReadStringOption(option_value.at(v), option, 1, ",", "ddCRP"))
        {
         Cerr("***** iRCPGSAlgorithm::InitSpecificOptions() -> Failed to read the list of options correctly!" << endl);
         return 1;
        }
        m_ddCRP = option[0];
      }
      else if (option_name.at(v)=="alpha")
      {
        HPFLTNB option[1];
        // Read them
        if (ReadStringOption(option_value.at(v), option, 1, ",", "alpha"))
        {
         Cerr("***** iRCPGSAlgorithm::InitSpecificOptions() -> Failed to read the list of options correctly!" << endl);
         return 1;
        }
        m_ddcrpAlpha = option[0];
      }
      else if (option_name.at(v)=="gammaShape")
      {
        HPFLTNB option[1];
        // Read them
        if (ReadStringOption(option_value.at(v), option, 1, ",", "gammaShape"))
        {
         Cerr("***** iRCPGSAlgorithm::InitSpecificOptions() -> Failed to read the list of options correctly!" << endl);
         return 1;
        }
        m_gammaShape = option[0];
      }
      else if (option_name.at(v)=="gammaRate")
      {
        HPFLTNB option[1];
        // Read them
        if (ReadStringOption(option_value.at(v), option, 1, ",", "gammaRate"))
        {
         Cerr("***** iRCPGSAlgorithm::InitSpecificOptions() -> Failed to read the list of options correctly!" << endl);
         return 1;
        }
        m_gammaRate = option[0];
      }
      else if (option_name.at(v)=="backprojection")
      {
        int option[1];
        // Read them
        if (ReadStringOption(option_value.at(v), option, 1, ",", "backprojection"))
        {
         Cerr("***** iRCPGSAlgorithm::InitSpecificOptions() -> Failed to read the list of options correctly!" << endl);
         return 1;
        }
        m_backprojection = option[0];
      }
      else if (option_name.at(v)=="multiModalNoiseSigma")
      {
        HPFLTNB option[1];
        // Read them
        if (ReadStringOption(option_value.at(v), option, 1, ",", "multiModalNoiseSigma"))
        {
         Cerr("***** iRCPGSAlgorithm::InitSpecificOptions() -> Failed to read the list of options correctly!" << endl);
         return 1;
        }
        temp_multiModalNoiseSigma.push_back(option[0]);
      }
      else if (option_name.at(v)=="multiModalLag")
      {
        INTNB option[1];
        // Read them
        if (ReadStringOption(option_value.at(v), option, 1, ",", "multiModalLag"))
        {
         Cerr("***** iRCPGSAlgorithm::InitSpecificOptions() -> Failed to read the list of options correctly!" << endl);
         return 1;
        }
        m_multiModalLag = option[0];
      }
      else if (option_name.at(v)=="meanClusterVolumeMin")
      {
        HPFLTNB option[1];
        // Read them
        if (ReadStringOption(option_value.at(v), option, 1, ",", "meanClusterVolumeMin"))
        {
         Cerr("***** iRCPGSAlgorithm::InitSpecificOptions() -> Failed to read the list of options correctly!" << endl);
         return 1;
        }
        m_meanClusterVolumeMin = option[0];
      }
      else if (option_name.at(v)=="meanClusterVolumeMax")
      {
        HPFLTNB option[1];
        // Read them
        if (ReadStringOption(option_value.at(v), option, 1, ",", "meanClusterVolumeMax"))
        {
         Cerr("***** iRCPGSAlgorithm::InitSpecificOptions() -> Failed to read the list of options correctly!" << endl);
         return 1;
        }
        m_meanClusterVolumeMax = option[0];
      }
      else if (option_name.at(v)=="alphaIncrement")
      {
        HPFLTNB option[1];
        // Read them
        if (ReadStringOption(option_value.at(v), option, 1, ",", "alphaIncrement"))
        {
         Cerr("***** iRCPGSAlgorithm::InitSpecificOptions() -> Failed to read the list of options correctly!" << endl);
         return 1;
        }
        m_ddcrpAlphaIncrement = option[0];
      }
    }
    if (temp_multiModalNoiseSigma.size()!=(size_t)mp_ID->GetNbMultiModalImages())
    {
      Cerr("***** iRCPGSAlgorithm::InitSpecificOptions() -> The number of provided noise standard deviations does not equal the number of provided multimodal images! " << endl);
      return 1;
    }
    for (int mmnb=0; mmnb<mp_ID->GetNbMultiModalImages(); mmnb++)
    {
      mp_multiModalNoiseSigma[mmnb] = temp_multiModalNoiseSigma.at(mmnb);
    }

    if (m_verbose>=2)
    {
      Cout("iRCPGSAlgorithm::InitSpecificOptions() -> Algorithm options read from command line "<<endl);
    }
  }
  // otherwise process config file
  else if(!file_options.empty())
  {
    ReadConfigurationFile(file_options);    
    if (m_verbose>=2)
    {
      Cout("iRCPGSAlgorithm::InitSpecificOptions() -> Algorithm options read from the configuration file "<<endl);
    }
  }

  // check parameters
  if (m_ddCRP>0 && m_ddcrpAlpha<0.)
  {
    Cerr("***** iRCPGSAlgorithm::InitSpecificOptions() -> The ddCRP alpha parameter must be non negative when ddCRP is used! " << endl);
    return 1;
  }
  if (m_ddCRP<0 || m_ddCRP>2)
  {
    Cerr("***** iRCPGSAlgorithm::InitSpecificOptions() -> The ddCRP type is not properly specified! " << endl);
    return 1;
  }
  if (m_backprojection<0 || m_backprojection>2)
  {
    Cerr("***** iRCPGSAlgorithm::InitSpecificOptions() -> The backprojection type is not properly specified! " << endl);
    return 1;
  }
  if (m_gammaShape<=0. || m_gammaRate<=0.)
  {
    Cerr("***** iRCPGSAlgorithm::InitSpecificOptions() -> Gamma prior parameters are not positive! " << endl);
    return 1;
  }
  for (int mmnb=0; mmnb<mp_ID->GetNbMultiModalImages(); mmnb++)
  {
    if (mp_multiModalNoiseSigma[mmnb]<=0.)
    {
      Cerr("***** iRCPGSAlgorithm::InitSpecificOptions() -> Wrong noise standard deviations for multimodal images! " << endl);
      return 1;
    }
  }

  if (m_meanClusterVolumeMin>0. || m_meanClusterVolumeMax>0.)
  {
    if (!(m_meanClusterVolumeMin>0. && m_meanClusterVolumeMax>0. && m_meanClusterVolumeMax>m_meanClusterVolumeMin))
    {
      Cerr("***** iRCPGSAlgorithm::InitSpecificOptions() -> Min/Max criteria for average cluster volume inconsistent! " << endl);
      return 1;
    }
    if (m_ddcrpAlphaIncrement<=1.)
    {
      Cerr("***** iRCPGSAlgorithm::InitSpecificOptions() -> Multiplicative update for ddCRP alpha must be >1. if min/max criteria for average cluster volume specified ! " << endl);
      return 1;
    }
    if (!(m_ddcrpAlpha>0.))
    {
      Cerr("***** iRCPGSAlgorithm::InitSpecificOptions() -> Options for adaptative alpha not compatible with alpha=0. ! " << endl);
      return 1;
    }
  }

  if (m_verbose>=2)
  {
    if (m_ddCRP>0)
    {
      Cout("iRCPGSAlgorithm::InitSpecificOptions() -> ddCRP parameters : type "<< m_ddCRP <<", alpha = "<<m_ddcrpAlpha<< ", multimodal lag = "<<m_multiModalLag<<" iterations "<< endl);
      if (!(m_ddcrpAlpha>0.))
      {
        Cout("                                                             as alpha=0, entering special mode when sampling the posterior clustering"<<endl);
      }
      if (m_meanClusterVolumeMax>0.)
      {
        Cout("                                                             alpha auto-tune with condition "<< m_meanClusterVolumeMin<<" < mean cluster volume < "<<m_meanClusterVolumeMax<< " mm3 and multiplicative increment = "<< m_ddcrpAlphaIncrement <<endl);
      }
    }
    else Cout("iRCPGSAlgorithm::InitSpecificOptions() -> ddCRP not used ! "<<endl);
    Cout("iRCPGSAlgorithm::InitSpecificOptions() -> Backprojection type = "<<m_backprojection<<endl);
  }

  return 0;
}

int iRCPGSAlgorithm::ReadConfigurationFile(const string& a_configurationFile)
{
  string key_word = "";

  key_word = "ddCRP type";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_ddCRP, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iRCPGSAlgorithm::ReadAndCheckConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }

  key_word = "ddCRP alpha";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_ddcrpAlpha, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iRCPGSAlgorithm::ReadAndCheckConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  
  key_word = "gamma prior shape";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_gammaShape, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iRCPGSAlgorithm::ReadAndCheckConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  
  key_word = "gamma prior rate";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_gammaRate, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iRCPGSAlgorithm::ReadAndCheckConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  
  key_word = "backprojection type";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_backprojection, 1, KEYWORD_MANDATORY))
  {
    Cerr("***** iRCPGSAlgorithm::ReadAndCheckConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }
  
  key_word = "multimodal noise sigma";
  if (ReadDataASCIIFile(a_configurationFile, key_word, mp_multiModalNoiseSigma, mp_ID->GetNbMultiModalImages(), KEYWORD_MANDATORY))
  {
    Cerr("***** iRCPGSAlgorithm::ReadAndCheckConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }

  key_word = "multimodal lag";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_multiModalLag, 1, KEYWORD_OPTIONAL))
  {
    Cerr("***** iRCPGSAlgorithm::ReadAndCheckConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }

  key_word = "mean cluster volume min";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_meanClusterVolumeMin, 1, KEYWORD_OPTIONAL))
  {
    Cerr("***** iRCPGSAlgorithm::ReadAndCheckConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }

  key_word = "mean cluster volume max";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_meanClusterVolumeMax, 1, KEYWORD_OPTIONAL))
  {
    Cerr("***** iRCPGSAlgorithm::ReadAndCheckConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }

  key_word = "ddCRP alpha increment";
  if (ReadDataASCIIFile(a_configurationFile, key_word, &m_ddcrpAlphaIncrement, 1, KEYWORD_OPTIONAL))
  {
    Cerr("***** iRCPGSAlgorithm::ReadAndCheckConfigurationFile() -> Failed to get the '" << key_word << "' keyword !" << endl);
    return 1;
  }

  // End
  return 0;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

int iRCPGSAlgorithm::InitializeHelperVar()
{
  if (m_verbose>=2) Cout("iRCPGSAlgorithm::InitializeHelperVar() "<<endl);

  int nbVoxels = mp_ID->GetNbVoxXYZ();
  bool multimodal_info = mp_ImageSpace->IsLoadedMultiModal();
  bool background_mask = mp_ImageSpace->IsLoadedMask();

  // number of neighbourhood voxels, currently fixed parameter, TODO implement different neighbourhoods
  m_neighbourhood = (mp_ID->GetNbVoxZ()>1)?6:4;
  if (m_ddCRP>0)
  {
    m_ddcrpLogAlpha = log(m_ddcrpAlpha);
    assert (!std::isnan(m_ddcrpLogAlpha));
  }

  mp_listRelevantVoxelIndices = new bool[nbVoxels];
  mp_voxelClusterMapping = new INTNB[nbVoxels];
  // Due to the current implementation, cluster indices <= nbVoxels, cluster array size thus being = nbVoxels+1
  mp_clusterValues = new HPFLTNB[nbVoxels+1];

  // the marginalized backprojection requires updates of mp_clusterN during parallel backprojection, so need for multithreaded mp_clusterN
  // approximative implementation, because modifications in cluster sums can't be shared between threads
  INTNB th_clusterN = (m_backprojection==2)?mp_ID->GetNbThreadsForProjection():1;
  mp_clusterN = new HPFLTNB*[th_clusterN];
  for (INTNB th=0; th<th_clusterN; th++) mp_clusterN[th] = new HPFLTNB[nbVoxels+1];
  
  mp_clusterSensitivity = new HPFLTNB[nbVoxels+1];

  // initialized only if multimodal info
  if (multimodal_info)
  {
    mpp_clusterMultiModal = new HPFLTNB*[mp_ID->GetNbMultiModalImages()];
    for (int mmnb=0; mmnb<mp_ID->GetNbMultiModalImages(); mmnb++)
    {
      mpp_clusterMultiModal[mmnb] = new HPFLTNB[nbVoxels+1];
    }
    mp_clusterCount = new INTNB[nbVoxels+1];
  }

  // description of ddCRP links
  mp_nextLink = new INTNB[nbVoxels];
  mpv_parentLinks = new vector<INTNB>[nbVoxels];
  mv_newClusters.assign(1,nbVoxels);
  mp_listVoxelIndices = new INTNB[nbVoxels];

  // permanent backward image for backprojection modes 1 and 2
  if (m_backprojection==1 || m_backprojection==2 ) mp_permanentBackwardImage = mp_ImageSpace->AllocateMiscellaneousImage();

  // helper variables for the backprojection modes 1 and 2
  // trace of voxels into which events were backprojected in the previous iteration 
  // initialize voxel indices at -1
  if (m_backprojection==1 || m_backprojection==2)
  {
    //mp_listEventsIndices = new INTNB*[m_nbBeds];
    m4p_EventsBackprojectionTrace = new INTNB***[m_nbBeds];
    for (int bed=0 ; bed<m_nbBeds ; bed++)
    {
      m4p_EventsBackprojectionTrace[bed] = new INTNB**[m2p_DataFile[bed]->GetSize()];
      //mp_listEventsIndices[bed] = new INTNB[m2p_DataFile[bed]->GetSize()];
    
      // Compute start and stop indices taking MPI into account (the vDataFile does that)
      int64_t index_start = 0;
      int64_t index_stop  = 0;
      m2p_DataFile[bed]->GetEventIndexStartAndStop(&index_start, &index_stop, 0, 1);

      int64_t index;
      // Keep the static scheduling with a chunk size at 1, it is important
      #pragma omp parallel for private(index) schedule(static, 1)
      for (index=0 ; index<m2p_DataFile[bed]->GetSize() ; index++)
      {
        // Get the thread index
        INTNB th = 0;
        #ifdef CASTOR_OMP
        th = omp_get_thread_num();
        #endif

        vEvent* event = m2p_DataFile[bed]->GetEvent(index, th);
        
        oProjectionLine *line = mp_ProjectorManager->ComputeProjectionLine(event, th);

        m4p_EventsBackprojectionTrace[bed][index] = new INTNB*[line->GetNbTOFBins()];

        for (int b=0; b<line->GetNbTOFBins(); b++)
        {
          m4p_EventsBackprojectionTrace[bed][index][b] = new INTNB[(INTNB)round(event->GetEventValue(b))];
          for (int e=0; e<(INTNB)round(event->GetEventValue(b)); e++)
          {
            m4p_EventsBackprojectionTrace[bed][index][b][e] = -1;
          }
        }
      }
    }
  }

  // voxel and clusters variables
  INTNB v;
  #pragma omp parallel for private(v) schedule(guided)
  for (v=0; v<nbVoxels; v++)
  {
    mp_listVoxelIndices[v] = v;

    mp_listRelevantVoxelIndices[v] = true;
    if (background_mask && (((INTNB)round(mp_ImageSpace->mp_maskImage[v])) == 0))
      mp_listRelevantVoxelIndices[v] = false;

    mp_voxelClusterMapping[v] = v;
    mp_clusterValues[v] = mp_ImageSpace->m4p_image[0][0][0][v];
    for (INTNB th=0; th<th_clusterN; th++) mp_clusterN[th][v] = 0.;
    mp_clusterSensitivity[v] = mp_ImageSpace->m5p_sensitivity[0][0][0][0][v];
    if (multimodal_info)
    {
      for (int mmnb=0; mmnb<mp_ID->GetNbMultiModalImages(); mmnb++)
      {
        mpp_clusterMultiModal[mmnb][v] = mp_ImageSpace->m2p_multiModalImage[mmnb][v];
      }
      mp_clusterCount[v] = 1;
    }
    mp_nextLink[v] = v;
    mpv_parentLinks[v].push_back(v);

    if (mp_permanentBackwardImage!=NULL) mp_permanentBackwardImage[v]=0.;
  }

  // parameter initialization for the free available cluster, will be recomputed anyway
  mp_clusterValues[nbVoxels] = 0.;
  for (INTNB th=0; th<th_clusterN; th++) mp_clusterN[th][nbVoxels] = 0.;
  mp_clusterSensitivity[nbVoxels] = 0.;
  if (multimodal_info)
  {
    for (int mmnb=0; mmnb<mp_ID->GetNbMultiModalImages(); mmnb++) mpp_clusterMultiModal[mmnb][nbVoxels] = 0.;
    mp_clusterCount[nbVoxels] = 0;
  }

  mp_multiModalParam = new HPFLTNB[mp_ID->GetNbMultiModalImages()];
  temp_count_multimodal = new HPFLTNB[mp_ID->GetNbMultiModalImages()];
  for (int mmnb=0; mmnb<mp_ID->GetNbMultiModalImages(); mmnb++)
  {
    temp_count_multimodal[mmnb] = 0.;
    mp_multiModalParam[mmnb] = 0.;
  }

  // End
  return 0;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

int iRCPGSAlgorithm::ProcessMultiModalInfo()
{
  if (m_verbose>=2) Cout("iRCPGSAlgorithm::ProcessMultiModalInfo() "<<endl);

  for (int mmnb=0; mmnb<mp_ID->GetNbMultiModalImages(); mmnb++)
  {
    FLTNB maxt = mp_ImageSpace->m2p_multiModalImage[mmnb][0];
    FLTNB mint = mp_ImageSpace->m2p_multiModalImage[mmnb][0];
    
    for (int v=1; v<mp_ID->GetNbVoxXYZ(); v++)
    {
      if (mp_ImageSpace->m2p_multiModalImage[mmnb][v]>maxt) maxt = mp_ImageSpace->m2p_multiModalImage[mmnb][v];
      if (mp_ImageSpace->m2p_multiModalImage[mmnb][v]<mint) mint = mp_ImageSpace->m2p_multiModalImage[mmnb][v];
    }

    // compute the prior parameter rho
    mp_multiModalParam[mmnb] = mp_multiModalNoiseSigma[mmnb]/maxt;

    if (m_verbose>=2) Cout("iRCPGSAlgorithm::StepBeforeIterationLoop() --> Multimodal image "<<mmnb<<" : noise std dev = "<<mp_multiModalNoiseSigma[mmnb]<<" , prior param rho = "<<mp_multiModalParam[mmnb]<<endl);
  }
  
  // End
  return 0;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

int iRCPGSAlgorithm::StepBeforeIterationLoop()
{
  if (vAlgorithm::StepBeforeIterationLoop())
  {
    Cerr("***** iRCPGSAlgorithm::StepBeforeIterationLoop() -> A problem occurred while calling StepBeforeIterationLoop() function !" << endl);
    return 1;
  }

  // Main image initialization
  if (mp_ImageSpace->InitImage(m_pathToInitialImg, 1.) )
  {
    Cerr("***** vAlgorithm::StepBeforeIterationLoop() -> An error occurred while reading the initialization image !" << endl);
    return 1;
  }

  // allocate backward images
  mp_ImageSpace->InstantiateBackwardImageFromDynamicBasis((m_backprojection==0)?1:2);

  // specific to this algorithm
  InitializeHelperVar();
  
  // currently the sensitivity image must be provided
  assert(mp_ImageSpace->IsLoadedSensitivity());

  // If subsets and OSEM-like backprojection (mode 0), then divide by the number of subsets because the acqusition duration is included in the sensitivity
  // No need for other backprojection modes because permanently dealing with the entire backprojected data for all the subsets
  if (mp_nbSubsets[0]>1 && m_backprojection==0)
  {
    if (m_verbose>=2)  Cout("iRCPGSAlgorithm::StepBeforeIterationLoop() --> Dividing the provided sensitivity image by the number of subsets "<<endl);
    for (INTNB v=0; v<mp_ID->GetNbVoxXYZ(); v++) mp_ImageSpace->m5p_sensitivity[0][0][0][0][v] /= ((FLTNB)mp_nbSubsets[0]);
  }

  // process additional (multimodal) images
  if (mp_ImageSpace->IsLoadedMultiModal()) ProcessMultiModalInfo();

  // End
  return 0;
}


// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iRCPGSAlgorithm::StepPreProcessInsideSubsetLoop(int a_iteration, int a_subset)
{
  if (m_verbose>=3) Cout("iRCPGSAlgorithm::StepPreProcessInsideSubsetLoop ... " << endl);

  // Initialize the correction backward image(s)
  mp_ImageSpace->InitBackwardImage();

  // Copy current image in forward-image buffer
  mp_ImageSpace->PrepareForwardImage();

  // End
  return 0;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

int iRCPGSAlgorithm::StepInnerLoopInsideSubsetLoop(int a_iteration, int a_subset, int a_bed)
{
  // Get the chrono manager singleton pointer
  sChronoManager* p_ChronoManager = sChronoManager::GetInstance();

  // Sample posterior complete data, backproject the acquired data, Gibbs Sampler Step 1
  p_ChronoManager->StartCustomStep(0,0);
  SampleConditionalCompleteData(a_iteration, a_subset, a_bed);  
  p_ChronoManager->StopCustomStep(0,0);

  // End
  return 0;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

int iRCPGSAlgorithm::StepPostProcessInsideSubsetLoop(int a_iteration, int a_subset)
{
  if (m_verbose>=3) Cout("iRCPGSAlgorithm::StepPostProcessInsideSubsetLoop ... " << endl);

  // Get the chrono manager singleton pointer
  sChronoManager* p_ChronoManager = sChronoManager::GetInstance();

  // Merge parallel results
  mp_ImageSpace->Reduce();

  // If mask provided, apply mask to sensitivity, the CleanNeverVisitedVoxels will take care of the images
  mp_ImageSpace->ApplyMaskToSensitivity();

  // Finalize processing of backward images if permanent backward image 
  if (mp_permanentBackwardImage!=NULL)
  {
    INTNB v;
    #pragma omp parallel for private(v) schedule(guided)
    for (v=0;v<mp_ID->GetNbVoxXYZ();v++)
    {
      mp_permanentBackwardImage[v] = mp_permanentBackwardImage[v] + mp_ImageSpace->m6p_backwardImage[1][0][0][0][0][v] + mp_ImageSpace->m6p_backwardImage[0][0][0][0][0][v];
    }
  }

  // update algorithm variables based on the new backprojection
  ComputeSumsPerClusters(a_iteration);

  // Sample posterior ddCRP (clustering), Gibbs Sampler Step 2
  if (m_ddCRP>0)
  {
    p_ChronoManager->StartCustomStep(0,1);

    SampleConditionalClustering(a_iteration);

    p_ChronoManager->StopCustomStep(0,1);
  }

  // Marginalized backprojection: as mp_clusterN threaded, and sampling step 2 (conditional clustering) modifies mp_clusterN only for the first thread,
  // copy the modifications to all the other threads for the next backprojection
  if (m_backprojection==2)
  {
    INTNB s;
    #pragma omp parallel for private(s) schedule(guided)
    for (s=0;s<mp_ID->GetNbVoxXYZ()+1;s++)
      for (INTNB th=1;th<mp_ID->GetNbThreadsForProjection();th++)
        mp_clusterN[th][s] = mp_clusterN[0][s];
  }

  // Sample posterior cluster intensities, Gibbs Sampler Step 3, once per iteration
  p_ChronoManager->StartCustomStep(0,2);
  SampleConditionalClusterIntensity();
  p_ChronoManager->StopCustomStep(0,2);

  // display performance info at each iteration
  p_ChronoManager->Display();

  // computed at each iteration because the variable is cleaned at each iteration after image clean-up
  UpdateVisitedVoxels();
  
  // Generate current posterior image sample from clustering/partition and cluster intensities
  GenerateCurrentImage();

  // if required, try to roughly update ddCRP alpha to obtain approximately the requested target average cluster volume
  // every 5 iterations and once per iteration
  if (m_meanClusterVolumeMax>0. && a_iteration>=5 && (a_iteration+1)%5==0 && a_subset==mp_nbSubsets[a_iteration]/2)
  {
    bool alpha_changed = false;
    // modify alpha if current average cluster volume less than m_meanClusterVolumeMax
    if (m_currentMeanClusterVolume>m_meanClusterVolumeMax)
    {
      m_ddcrpAlpha *= m_ddcrpAlphaIncrement;
      alpha_changed = true;
    }
    // modify alpha if current average cluster volume more than m_meanClusterVolumeMin
    else if (m_currentMeanClusterVolume<m_meanClusterVolumeMin)
    {
      m_ddcrpAlpha /= m_ddcrpAlphaIncrement;
      alpha_changed = true;
    }
    // Log the modification
    if (alpha_changed)
    {
      m_ddcrpLogAlpha = log(m_ddcrpAlpha);
      Cout("iRCPGSAlgorithm::StepPostProcessInsideSubsetLoop() -> Updated ddCRP alpha to "<<m_ddcrpAlpha<<endl);
    }
  }

  // Save the sensitivity image in histogram mode, if asked for
  if (mp_ID->GetMPIRank()==0 && mp_outputIterations[a_iteration] && m_saveSensitivityHistoFlag && m2p_DataFile[0]->GetDataMode()==MODE_HISTOGRAM)
  {
    // Get output manager to build the file name
    sOutputManager* p_output_manager = sOutputManager::GetInstance();
    // Build the file name
    string temp_sens = p_output_manager->GetPathName() + p_output_manager->GetBaseName();
    stringstream temp_it; temp_it << a_iteration + 1;
    stringstream temp_ss; temp_ss << a_subset + 1;
    temp_sens.append("_it").append(temp_it.str()).append("_ss").append(temp_ss.str()).append("_sensitivity");
    // Save sensitivity
    mp_ImageSpace->SaveSensitivityImage(temp_sens);
  }

  // Save the current image estimate if asked for
  if (mp_ID->GetMPIRank()==0 && mp_outputIterations[a_iteration] && m_saveImageAfterSubsets)
  {
    // Verbose
    if (m_verbose>=1) Cout("iRCPGSAlgorithm::StepPostProcessInsideSubsetLoop() -> Save image at iteration " << a_iteration+1 << " for subset " << a_subset+1 << endl);
    // Save image
    if (StepImageOutput(a_iteration,a_subset))
    {
      Cerr("***** iRCPGSAlgorithm::StepPostProcessInsideSubsetLoop() -> A problem occurred while saving images at iteration " << a_iteration+1 << " !" << endl);
      return 1;
    }
  }

  // End
  return 0;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

int iRCPGSAlgorithm::SampleConditionalCompleteData(int a_iteration, int a_subset, int a_bed)
{

  // Verbose
  if (m_verbose>=2)
  {
    if (m_nbBeds>1) Cout("iRCPGSAlgorithm::SampleConditionalCompleteData() -> Start loop over events for bed " << a_bed+1 << endl << flush);
    else Cout("iRCPGSAlgorithm::SampleConditionalCompleteData() -> Start loop over events" << endl << flush);
  }

  // Reinitialize 4D gate indexes
  mp_ID->ResetCurrentDynamicIndices();

  // Apply the bed offset for this bed position
  mp_ProjectorManager->ApplyBedOffset(a_bed);

  // Progression (increments of 2%)
  if (m_verbose>=VERBOSE_NORMAL && mp_ID->GetMPIRank()==0)
  {
    cout << "0   10%  20%  30%  40%  50%  60%  70%  80%  90%  100%" << endl;
    cout << "|----|----|----|----|----|----|----|----|----|----|" << endl;
    cout << "|" << flush;
  }

  int progression_percentage_old = 0;
  int progression_nb_bars = 0;
  uint64_t progression_printing_index = 0;

  // Compute start and stop indices taking MPI into account (the vDataFile does that)
  int64_t index_start = 0;
  int64_t index_stop  = 0;

  m2p_DataFile[a_bed]->GetEventIndexStartAndStop(&index_start, &index_stop, a_subset, mp_nbSubsets[a_iteration]);

  // Synchronize MPI processes
  #ifdef CASTOR_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  #endif

  // Set the number of threads for projections (right after this loop, we set back the number of threads for image computations)
  #ifdef CASTOR_OMP
  omp_set_num_threads(mp_ID->GetNbThreadsForProjection());
  #endif

  // This boolean is used to report any problem inside the parallel loop
  bool problem = false;

  // Get the chrono manager singleton pointer
  sChronoManager* p_ChronoManager = sChronoManager::GetInstance();

  // shuffle events
  if (m_backprojection==1)
  {
    //random_shuffle(mp_listEventsIndices[a_bed], mp_listEventsIndices[a_bed]+(index_stop-index_start));
  }

  // Launch the loop on events
  int64_t index_temp;
  // Keep the static scheduling with a chunk size at 1, it is important
  #pragma omp parallel for private(index_temp) schedule(static, 1)
  for ( index_temp = index_start  ;  index_temp < index_stop  ; index_temp += mp_nbSubsets[a_iteration])
  {
    // Get the thread index
    int th = 0;
    #ifdef CASTOR_OMP
    th = omp_get_thread_num();
    #endif

    // get random event
    INTNB index = index_temp;
    //if (m_backprojection==1) index = mp_listEventsIndices[a_bed][index_temp];

    // Print progression (do not log out with Cout here)
    if (m_verbose>=2 && th==0 && mp_ID->GetMPIRank()==0)
    {
      if (progression_printing_index%1000==0)
      {
        int progression_percentage_new = ((int)( (((float)(index-index_start+1))/((float)(index_stop-index_start)) ) * 100.));
        if (progression_percentage_new>=progression_percentage_old+2) // Increments of 2%
        {
          int nb_steps = (progression_percentage_new-progression_percentage_old)/2;
          for (int i=0; i<nb_steps; i++)
          {
            cout << "-" << flush;
            progression_nb_bars++;
          }
          progression_percentage_old += nb_steps*2;
        }
      }
      progression_printing_index++;
    }

    // Step 1: Get the current event for that thread index
    p_ChronoManager->StartIterativeDataUpdateStep1(th);

    vEvent* event = m2p_DataFile[a_bed]->GetEvent(index, th);

    if (event==NULL)
    {
      Cerr("***** iRCPGSAlgorithm::SampleConditionalCompleteData() -> An error occured while getting the event from index "
        << index << " (thread " << th << ") !" << endl);
      // Specify that there was a problem
      problem = true;
      // We must continue here because we are inside an OpenMP loop
      continue;
    }
    p_ChronoManager->StopIterativeDataUpdateStep1(th);

    // Step 2: Call the dynamic switch function that updates the current frame and gate numbers, and also detects involuntary patient motion
    p_ChronoManager->StartIterativeDataUpdateStep2(th);
    #ifdef CASTOR_DEBUG
    if (m_verbose>=4)
    {
      Cout("iRCPGSAlgorithm::SampleConditionalCompleteData() -> Step2: Check for Dynamic event (frame/gate switch, image-based deformation " << endl);
    }
    #endif
    int dynamic_switch_value = mp_ID->DynamicSwitch(index, event->GetTimeInMs(), a_bed, th);
    // If the DYNAMIC_SWITCH_CONTINUE is returned, then it means that we are not yet at the first frame
    if ( dynamic_switch_value == DYNAMIC_SWITCH_CONTINUE )
    {
      // Then we just skip this event
      continue;
    }

    p_ChronoManager->StopIterativeDataUpdateStep2(th);
      
    // Skip if empty event
    bool emptyEvent = true;
    for (int b=0; b<event->GetNbValueBins(); b++) 
      if (event->GetEventValue(b)>0.)
      {
        emptyEvent = false;
        break;
      }

    if (emptyEvent) continue;

    // Step 3: Compute the projection line
    p_ChronoManager->StartIterativeDataUpdateStep3(th);
    #ifdef CASTOR_DEBUG
    if (m_verbose>=4) Cout("iRCPGSAlgorithm::SampleConditionalCompleteData() -> Step3: Compute the projection line " << endl);
    #endif
    oProjectionLine *line = mp_ProjectorManager->ComputeProjectionLine(event, th);
    if (line==NULL)
    {
      Cerr("***** iRCPGSAlgorithm::SampleConditionalCompleteData() -> A problem occured while computing the projection line !" << endl);
      // Specify that there was a problem
      problem = true;
      // We must continue here because we are inside an OpenMP loop
      continue;
    }
    p_ChronoManager->StopIterativeDataUpdateStep3(th);

    if(line->NotEmptyLine()) 
    {
      // get time respiratory and cardiac frame indices
      int timeFrame = mp_ID->GetCurrentTimeFrame(th);
      int respGate = mp_ID->GetCurrentRespImage(th);
      int cardGate = mp_ID->GetCurrentCardImage(th);
      
      assert(timeFrame==0);assert(respGate==0);assert(cardGate==0);
      
      // Loop over time basis functions
      for (int tbf=0; tbf<mp_ID->GetNbTimeBasisFunctions(); tbf++)
      {
        // Get frame/basis coefficient and continue if null
        FLTNB time_basis_coef = mp_ID->GetTimeBasisCoefficient(tbf,timeFrame);
        if (time_basis_coef==0.) continue;
        // Loop over respiratory basis functions
        for (int rbf=0; rbf<mp_ID->GetNbRespBasisFunctions(); rbf++)
        {
          // Get resp_gate/basis coefficient and continue if null
          FLTNB resp_basis_coef = mp_ID->GetRespBasisCoefficient(rbf,respGate);
          if (resp_basis_coef==0.) continue;
          // Loop over cardiac basis functions
          for (int cbf=0; cbf<mp_ID->GetNbCardBasisFunctions(); cbf++)
          {
            // Get card_gate_basis coefficient and continue if null
            FLTNB card_basis_coef = mp_ID->GetCardBasisCoefficient(cbf,cardGate);
            if (card_basis_coef==0.) continue;
            // Compute global coefficient
            FLTNB global_basis_coef = time_basis_coef * resp_basis_coef * card_basis_coef;
            assert(global_basis_coef==1.);
            
            // Compute multiplicative correction for this LOR
            FLTNB multiplicative_corr = 1./(event->GetMultiplicativeCorrections() * mp_ID->GetQuantificationFactor(a_bed, timeFrame, respGate, cardGate));

            // Loop over TOF bins
            for (int b=0; b<line->GetNbTOFBins(); b++)
            {
            FLTNB temp_eventValue = event->GetEventValue(b);
            if (temp_eventValue>0.)
            {
              // multinomial backprojection step with permanent backward image, backprojection variables only updated
              if (m_backprojection==1 || m_backprojection==2)
              {
                INTNB numCounts = (INTNB)round(temp_eventValue);
                // the number of counts must be integer for this backprojection,
                // otherwise events tracking gets complicated
                assert(round(temp_eventValue)==temp_eventValue);

                // process acquired counts one by one
                for (int e=0; e<numCounts; e++)
                {
                  // remove the count backprojected at the previous iteration from all the variables
                  // if backprojected into random/scatter there is nothing to remove
                  if (a_iteration>0 && !(m4p_EventsBackprojectionTrace[a_bed][index][b][e]<0))
                  {
                    // voxel into which the count was backprojected at the previous iteration
                    INTNB previous_voxel = m4p_EventsBackprojectionTrace[a_bed][index][b][e];
                    INTNB previous_voxel_current_cluster = mp_voxelClusterMapping[previous_voxel];

                    // remove this count from the backprojected complete data
                    if (mp_permanentBackwardImage!=NULL)
                    {
                      mp_ImageSpace->m6p_backwardImage[1][th][tbf][rbf][cbf][previous_voxel] -= 1.;
                    }

                    if (m_backprojection==2)
                    {
                      // update cluster sums of backprojected counts, because this variable is used in marginalized backprojection,
                      // only the current thread is updated, modifications in other threads can't be taken into account,
                      // so approximative implementation when using more than one thread
                      mp_clusterN[th][previous_voxel_current_cluster] -= 1.;
                      assert(mp_clusterN[th][previous_voxel_current_cluster]>-0.001);
                    }
                  }

                  p_ChronoManager->StartCustomStep(th,3);
                  HPFLTNB temp_value = 0.;

                  // compute probabilities for backprojection
                  INTNB current_nb_voxels = line->GetCurrentNbVoxels(FORWARD,b);
                  vector<HPFLTNB> voxelProbabilities(current_nb_voxels+1,0.);
                  HPFLTNB probabilities_sum = 0.;

                  INTNB temp_voxelIndex = -1;
                  // voxels outcomes
                  if (a_iteration==0)
                  {
                    // first iteration : if the initial image = 1, the probability becomes aij
                    for (int vl=0; vl<current_nb_voxels; vl++)
                    {
                      temp_voxelIndex = line->GetVoxelIndex(FORWARD,b,vl);
                      temp_value = (!mp_listRelevantVoxelIndices[temp_voxelIndex]) ? 0. : line->GetVoxelWeights(FORWARD,b,vl) * multiplicative_corr * global_basis_coef;
                      
                      probabilities_sum += temp_value;
                      voxelProbabilities.at(vl+1) = temp_value;
                    }
                  }
                  else
                  {
                    // voxels outcomes
                    for (int vl=0; vl<current_nb_voxels; vl++)
                    {
                      temp_voxelIndex = line->GetVoxelIndex(FORWARD,b,vl);
                      // Aij (system matrix element) * expectation of conditional cluster intensity
                      temp_value = (!mp_listRelevantVoxelIndices[temp_voxelIndex]) ? 0. : line->GetVoxelWeights(FORWARD,b,vl) * multiplicative_corr * global_basis_coef;
                      if (m_backprojection==2)
                      { // multinomial backprojection marginalized over cluster intensity
                        temp_value *= ((mp_clusterN[th][mp_voxelClusterMapping[temp_voxelIndex]] + m_gammaShape)
                                      /(mp_clusterSensitivity[mp_voxelClusterMapping[temp_voxelIndex]] + m_gammaRate));
                      }
                      else
                      { // multinomial backprojection
                        temp_value *= mp_ImageSpace->m4p_forwardImage[tbf][rbf][cbf][temp_voxelIndex];
                      }

                      probabilities_sum += temp_value;
                      voxelProbabilities.at(vl) = temp_value;
                    }
                  }


                  // random and scatter counts as the last outcome
                  temp_value = event->GetAdditiveCorrections(b) * mp_ID->GetFrameDurationInSec(a_bed, timeFrame);
                  voxelProbabilities[current_nb_voxels] = temp_value;
                  probabilities_sum += temp_value;

                  p_ChronoManager->StopCustomStep(th,3);
                  assert(!(probabilities_sum<0.));
                  // if the probability distribution is not empty for this line
                  if (probabilities_sum>0.)
                  {
                    // get a random number from uniform distribution and multiply it by the sum of probabilities,
                    // faster than normalizing probabilities before sampling
                    //HPFLTNB current_random = gsl_rng_uniform(mpp_threadRandomGenerators[th])*probabilities_sum;
                    HPFLTNB current_random = sRandomNumberGenerator::GetInstance()->GenerateRdmNber()*probabilities_sum;
                    HPFLTNB cumulative_sum = 0.;
                    INTNB sampledIndex = 0;
                    for (sampledIndex=current_nb_voxels; sampledIndex>=0; sampledIndex--) 
                    {
                      cumulative_sum += voxelProbabilities[sampledIndex];
                      if (cumulative_sum>current_random) break;
                    }

                    // write the realization into the corresponding voxel and update other variables 
                    if (sampledIndex>=0 && sampledIndex<current_nb_voxels)
                    {
                      INTNB drawn_voxel_index = line->GetVoxelIndex(FORWARD,b,sampledIndex);
                      INTNB new_cluster = mp_voxelClusterMapping[drawn_voxel_index];

                      // add the count into the backprojected complete data
                      mp_ImageSpace->m6p_backwardImage[0][th][tbf][rbf][cbf][drawn_voxel_index] += 1.;

                      if (m_backprojection==2)
                      {
                        // update cluster sums of backprojected counts, because this variable is used in the backprojection
                        mp_clusterN[th][new_cluster] += 1.;
                      }

                      // update backprojection trace for the next iteration/subset
                      m4p_EventsBackprojectionTrace[a_bed][index][b][e] = drawn_voxel_index;
                    }
                    else // random/scatter outcome
                    {
                      // flag random/scatter outcome, not really used at present
                      m4p_EventsBackprojectionTrace[a_bed][index][b][e] = -2;
                    }
                  }
                }
              }  
              else // ordinary multinomial backprojection step, overwrites all the backprojection variables
              {
                p_ChronoManager->StartCustomStep(th,3);
                HPFLTNB temp_value = 0.;

                // Compute categorical probability distribution over voxels contributing to the line 
                // and the scatter/random source contributing to the line;
                // The probabilities are not normalized, the normalization factor is taken into 
                // account during sampling
                INTNB current_nb_voxels = line->GetCurrentNbVoxels(FORWARD,b);
                vector<HPFLTNB> voxelProbabilities(current_nb_voxels+1,0.);
                HPFLTNB probabilities_sum = 0.;

                // voxels
                for (int vl=0; vl<current_nb_voxels; vl++)
                {
                  INTNB temp_voxelIndex = line->GetVoxelIndex(FORWARD,b,vl);
                  temp_value = (!mp_listRelevantVoxelIndices[temp_voxelIndex]) ? 0. :line->GetVoxelWeights(FORWARD,b,vl) * multiplicative_corr * global_basis_coef
                  * mp_ImageSpace->m4p_forwardImage[tbf][rbf][cbf][temp_voxelIndex];
                  //* mp_clusterValues[mp_voxelClusterMapping[temp_voxelIndex]];

                  probabilities_sum += temp_value;
                  voxelProbabilities.at(vl) = temp_value;
                }

                // random and scatter counts as the last outcome
                temp_value = event->GetAdditiveCorrections(b) * mp_ID->GetFrameDurationInSec(a_bed, timeFrame);
                voxelProbabilities[current_nb_voxels] = temp_value;
                probabilities_sum += temp_value;

                assert(!(probabilities_sum<0.));
                p_ChronoManager->StopCustomStep(th,3);

                // if the probability distribution is not empty for this line
                if (probabilities_sum>0.)
                {
                  p_ChronoManager->StartCustomStep(th,4);
                  INTNB numCounts = (INTNB)round(temp_eventValue);
                  // if the number of counts is not integer (might happen for GE PET/MR sinogram data),
                  // pick randomly floor or ceil rounding
                  if (fabs(round(temp_eventValue)-temp_eventValue)>1.e-7)
                  {
                    if (sRandomNumberGenerator::GetInstance()->GenerateRdmNber()>0.5)
                    {
                      numCounts = (INTNB)ceil(temp_eventValue);
                    }
                    else
                    { 
                      numCounts = (INTNB)floor(temp_eventValue);
                    }
                  }
                  p_ChronoManager->StopCustomStep(th,4);
                  p_ChronoManager->StartCustomStep(th,5);
                  // the number of repetitions for the multinomial distribution = event value
                  for (int e=0; e<numCounts; e++)
                  {
                    // get a random number from uniform distribution and multiply it by the sum of probabilities,
                    // faster than normalizing probabilities before sampling
                    //HPFLTNB current_random = gsl_rng_uniform(mpp_threadRandomGenerators[th])*probabilities_sum;
                    HPFLTNB current_random = sRandomNumberGenerator::GetInstance()->GenerateRdmNber()*probabilities_sum;
                    HPFLTNB cumulative_sum = 0.;
                    INTNB sampledIndex = 0;
                    // faster if start from the random/scatter outcome, which probably has higher probability
                    for (sampledIndex=current_nb_voxels; sampledIndex>=0; sampledIndex--)
                    {
                      cumulative_sum += voxelProbabilities[sampledIndex];
                      if (cumulative_sum>current_random) break;
                    }
                    // write the realization into the corresponding voxel, and don't do anything if the sampled index refers to
                    // random and scatter source (the last probability in the list)
                    if (sampledIndex>=0 && sampledIndex<current_nb_voxels)
                    {
                      mp_ImageSpace->m6p_backwardImage[0][th][tbf][rbf][cbf][line->GetVoxelIndex(FORWARD,b,sampledIndex)] += 1.;
                    }
                    else if (sampledIndex<0)
                    {
                      // check : if nothing sampled throw a warning, though it should not happen
                      Cout("Warning: Nothing sampled, sampledIndex "<<sampledIndex<<", current_nb_voxels "<<current_nb_voxels<<endl);
                    }
                  }
                  p_ChronoManager->StopCustomStep(th,5);
                }
              }
            }
            }
          }
        }
      } 
    }
  } // End of events loop (OpenMP stops here)

  // Synchronize MPI processes
  #ifdef CASTOR_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  #endif

  // Set back the number of threads for image computation
  #ifdef CASTOR_OMP
  omp_set_num_threads(mp_ID->GetNbThreadsForImageComputation());
  #endif

  // End of progression printing (do not log out with Cout here)
  if (m_verbose>=2 && mp_ID->GetMPIRank()==0)
  {
    int progression_total_bars = 49;
    for (int i=0; i<progression_total_bars-progression_nb_bars; i++) cout << "-";
    cout << "|" << endl;
  }

  // If a problem was encountered, then report it here
  if (problem)
  {
    Cerr("***** iRCPGSAlgorithm::SampleConditionalCompleteData() -> A problem occured inside the parallel loop over events !" << endl);
    return 1;
  }

  // End
  return 0;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

int iRCPGSAlgorithm::ComputeSumsPerClusters(int a_iteration)
{
  if (m_verbose>=2) Cout("iRCPGSAlgorithm::ComputeSumsPerClusters()... " << endl);

  bool multimodal_info = mp_ImageSpace->IsLoadedMultiModal();

  // reinitialize and compute sums over clusters
  // sum of nij over clusters must be recomputed at each iteration, because nij have been resampled
  // recomputing sensitivity sum over clusters is not compulsory because the sensitivity is constant over iterations 
  // and the sums are updated at each clustering sampling, but it is done nevertheless, in order to avoid accumulations of numerical errors
  
  // marginalized backprojection needs threaded cluster sums
  INTNB th_clusterN = (m_backprojection==2)?mp_ID->GetNbThreadsForProjection():1;

  // reinitialize
  INTNB vp;
  #pragma omp parallel for private(vp) schedule(guided)
  for (vp=0; vp<(mp_ID->GetNbVoxXYZ()+1); vp++) 
  {
    mp_clusterN[0][vp] = 0.;
    for (INTNB th=1;th<th_clusterN;th++) mp_clusterN[th][vp] = 0.;
    mp_clusterSensitivity[vp] = 0.;
    if (multimodal_info)
    {
      for (int mmnb=0; mmnb<mp_ID->GetNbMultiModalImages(); mmnb++) mpp_clusterMultiModal[mmnb][vp] = 0.;
      mp_clusterCount[vp] = 0;
    }
  }

  // sum using voxel to cluster mapping
  for (int v=0; v<mp_ID->GetNbVoxXYZ(); v++)
  {
    // skip if background voxel
    if (!mp_listRelevantVoxelIndices[v]) continue;

    if (mp_permanentBackwardImage!=NULL)
    {
      assert(mp_permanentBackwardImage[v]>=-0.0001);
      mp_clusterN[0][mp_voxelClusterMapping[v]] += mp_permanentBackwardImage[v];
    }
    else
    {
      assert(mp_ImageSpace->m6p_backwardImage[0][0][0][0][0][v]>=0.);
      mp_clusterN[0][mp_voxelClusterMapping[v]] += mp_ImageSpace->m6p_backwardImage[0][0][0][0][0][v];
    }
    
    assert(mp_ImageSpace->m5p_sensitivity[0][0][0][0][v]>=0.);
    mp_clusterSensitivity[mp_voxelClusterMapping[v]] += mp_ImageSpace->m5p_sensitivity[0][0][0][0][v] ;
    if (multimodal_info)
    {
      for (int mmnb=0; mmnb<mp_ID->GetNbMultiModalImages(); mmnb++)
      {
        mpp_clusterMultiModal[mmnb][mp_voxelClusterMapping[v]] += mp_ImageSpace->m2p_multiModalImage[mmnb][v] ;
      }
      mp_clusterCount[mp_voxelClusterMapping[v]] += 1;
    }
  }

  // End
  return 0;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

int iRCPGSAlgorithm::UpdateVisitedVoxels()
{
  // Verbose
  if (m_verbose>=3) Cout("iRCPGSAlgorithm::UpdateVisitedVoxels() -> Tick visited voxels based on sensitivity" << endl);

  // zeros in the sensitivity map can be used to detect visited voxels, though not needed in this algorithm?
  int count=0;
  for (int v=0; v<mp_ID->GetNbVoxXYZ(); v++)
  {
    // Use the sensitivity image to find visited voxels (we just add 1 to discriminate between visited and not visited voxels)
    if (mp_ImageSpace->m5p_sensitivity[0][0][0][0][v]>0.) mp_ImageSpace->mp_visitedVoxelsImage[v] += 1.;
    else count++;
  }
  if (m_verbose>=2 && count>0) Cout("  --> Invisible voxels : "<<count<<" over "<<mp_ID->GetNbVoxXYZ()<<" !"<<endl);
  
  // End
  return 0;  
}

// -----------------------------------------------------------------------------------------------------------------------------------------

int iRCPGSAlgorithm::SampleConditionalClustering(int a_iteration)
{
  if (m_verbose>=2) Cout("iRCPGSAlgorithm::SampleConditionalClustering()... " << endl);
  // voxels eligible for assignment to the same cluster (here simple voxel neighbourhood)
  // variable size with maximum at m_neighbourhood+1

  //HPFLTNB mean_proba = 0., mean_proba_prev = 0., nb_proba = 0., std_proba = 0.;

  bool multimodal_info = mp_ImageSpace->IsLoadedMultiModal();

  vector<INTNB> fellow_voxels;
  fellow_voxels.reserve(m_neighbourhood+1);
  
  // non-normalized probabilities for the sampling of the next link
  // variable size with maximum at m_neighbourhood+1
  vector<HPFLTNB> probabilities;
  probabilities.reserve(m_neighbourhood+1);
  
  // voxels belonging to the same cluster as the current voxel,
  // after the next link of the current voxel has been cut (all parent voxels)
  vector<INTNB> voxels_after_cut;
  
  // dimensions
  INTNB dim_3D = mp_ID->GetNbVoxXYZ();

  // randomize voxel indices before sampling the ddCRP
  //random_shuffle(mp_listVoxelIndices, mp_listVoxelIndices+dim_3D);
  shuffle(mp_listVoxelIndices, mp_listVoxelIndices+dim_3D, sRandomNumberGenerator::GetInstance()->GetExtraGenerator(1));
  
  // resample voxels links
  for(int v=0; v<dim_3D; v++) 
  { 
    const int current_voxel = mp_listVoxelIndices[v];
    
    // skip if background voxel
    if (!mp_listRelevantVoxelIndices[current_voxel]) continue;

    const int previous_cluster = mp_voxelClusterMapping[current_voxel];

    // temporary new cluster
    int temp_new_cluster = mv_newClusters.at(mv_newClusters.size()-1);
    mv_newClusters.pop_back();
    assert( temp_new_cluster<=dim_3D );
    assert(temp_new_cluster!=previous_cluster);
    
    // get voxel indices in the neighbourhood
    ComputeFellowVoxelsList(fellow_voxels, current_voxel);

    // explore the parent links of the current voxel, check whether these links make a loop,
    // and build the list of voxels connected to the current voxel by parent links    
    bool loop = false;  
    voxels_after_cut.resize(0);

    // lambda function for optimization purposes
    function<void (int)> explorePreviousLinks= [this,&explorePreviousLinks,current_voxel,&voxels_after_cut,&loop] (int index) 
    {
      voxels_after_cut.push_back(index);
      for(int ip: mpv_parentLinks[index])
      {
        if (ip!=current_voxel) explorePreviousLinks(ip);
        else loop=true;
      }
    };  
    explorePreviousLinks(current_voxel);
    assert( voxels_after_cut.size()>0 );

    // cut the next (child) link of the current voxel, and update parent links
    size_t previous_next_link = mp_nextLink[current_voxel];
    for(size_t p=0; p<mpv_parentLinks[previous_next_link].size(); p++)
    {
      if (mpv_parentLinks[previous_next_link][p]==current_voxel)
      {
        mpv_parentLinks[previous_next_link].erase(mpv_parentLinks[previous_next_link].begin()+p);
        break;
      }
    }

    // assign the temporary new cluster to the current voxel and to all its parent voxels
    // compute also modifications of sums over the cluster
    HPFLTNB temp_count_N = 0.;
    HPFLTNB temp_count_sens = 0.;
    INTNB temp_count = 0;
    for (int mmnb=0; mmnb<mp_ID->GetNbMultiModalImages(); mmnb++) temp_count_multimodal[mmnb] = 0.;

    for(int nv : voxels_after_cut) 
    {
      if (mp_permanentBackwardImage!=NULL) temp_count_N += mp_permanentBackwardImage[nv];
      else temp_count_N += mp_ImageSpace->m6p_backwardImage[0][0][0][0][0][nv];
      
      temp_count_sens += mp_ImageSpace->m5p_sensitivity[0][0][0][0][nv];
      if (multimodal_info) 
      {
        temp_count++;
        for (int mmnb=0; mmnb<mp_ID->GetNbMultiModalImages(); mmnb++)
        {
          temp_count_multimodal[mmnb] += mp_ImageSpace->m2p_multiModalImage[mmnb][nv];
          assert (temp_count_multimodal[mmnb]>=0.);
        }
      }
		  mp_voxelClusterMapping[nv] = temp_new_cluster;
		}
    assert (temp_count_N>=0. && temp_count_sens>=0.);
    assert (temp_count>=0);
    
    // if this cluster contains the same voxels as the previous cluster (looped links), 
    // free the previous cluster, make it available for assignment
    // (all these voxels have already been assigned to the temporary new cluster)
    if (loop) mv_newClusters.push_back(previous_cluster);

    // update sums for the temporary new cluster and the remaining part of the previous_cluster (might be empty)
    mp_clusterN[0][previous_cluster] -= temp_count_N;
    mp_clusterSensitivity[previous_cluster] -= temp_count_sens;
    mp_clusterN[0][temp_new_cluster] = temp_count_N;
    mp_clusterSensitivity[temp_new_cluster] = temp_count_sens;

    assert(mp_clusterN[0][previous_cluster]>=0.);
    //if (mp_clusterSensitivity[previous_cluster]<-1.e-11) Cout("WARNING mp_clusterSensitivity[previous_cluster]<0! : "<<mp_clusterSensitivity[previous_cluster]<<" cluster idx "<<previous_cluster<<" temp_count_sens "<<temp_count_sens<<endl);
    assert(mp_clusterN[0][temp_new_cluster]>=0.);
    //if (mp_clusterSensitivity[temp_new_cluster]<-1.e-11) Cout("WARNING mp_clusterSensitivity[temp_new_cluster]<0! : "<<mp_clusterSensitivity[temp_new_cluster]<<" cluster idx "<<temp_new_cluster<<" temp_count_sens "<<temp_count_sens<<endl);
    
    if (multimodal_info) 
    {
      mp_clusterCount[previous_cluster] -= temp_count;
      mp_clusterCount[temp_new_cluster] = temp_count;
      assert(mp_clusterCount[previous_cluster]>=0);
      assert(mp_clusterCount[temp_new_cluster]>=0);
      for (int mmnb=0; mmnb<mp_ID->GetNbMultiModalImages(); mmnb++)
      {
        mpp_clusterMultiModal[mmnb][previous_cluster] -= temp_count_multimodal[mmnb];
        mpp_clusterMultiModal[mmnb][temp_new_cluster] = temp_count_multimodal[mmnb];
        assert(mpp_clusterMultiModal[mmnb][previous_cluster]>=0.);
        assert(mpp_clusterMultiModal[mmnb][temp_new_cluster]>=0.);
      }
    }

    // compute categorical probabilities for sampling the next link for the current voxel
    probabilities.resize(0);

    // add self-link if alpha not zero
    if (m_ddcrpAlpha>0.) probabilities.push_back(m_ddcrpLogAlpha);

    for(size_t fv=1; fv<fellow_voxels.size(); fv++) 
    {
      int fellow_cluster = mp_voxelClusterMapping[fellow_voxels.at(fv)];
      if(fellow_cluster!=temp_new_cluster)
      {
        // probability of a link that will cause the fusion of two clusters
        HPFLTNB n0 = mp_clusterN[0][temp_new_cluster];
        HPFLTNB n1 = mp_clusterN[0][fellow_cluster];
        HPFLTNB s0 = mp_clusterSensitivity[temp_new_cluster];
        HPFLTNB s1 = mp_clusterSensitivity[fellow_cluster];
        
        assert (!(n0<0. || n1<0. || s0<0. || s1<0.));

        //  ln (probability) to avoid numerical issues
        HPFLTNB temp = (n0+m_gammaShape)*log(s0+m_gammaRate)+(n1+m_gammaShape)*log(s1+m_gammaRate)-(n0+n1+m_gammaShape)*log(s0+s1+m_gammaRate)
          -lgamma(n0+m_gammaShape)-lgamma(n1+m_gammaShape)+lgamma(n0+n1+m_gammaShape)
          +lgamma(m_gammaShape)-m_gammaShape*log(m_gammaRate);

        // compute and add the multimodal influence only if beyond the specified number of iterations
        if (multimodal_info && a_iteration>=m_multiModalLag) 
        {
          HPFLTNB c0 = mp_clusterCount[temp_new_cluster];
          HPFLTNB c1 = mp_clusterCount[fellow_cluster];
          assert (c0>0. && c1>0.);          
          HPFLTNB help1 = c0*c1/(c0+c1);

          for (int mmnb=0; mmnb<mp_ID->GetNbMultiModalImages(); mmnb++)
          {
            HPFLTNB a0 = mpp_clusterMultiModal[mmnb][temp_new_cluster];
            HPFLTNB a1 = mpp_clusterMultiModal[mmnb][fellow_cluster];
            assert (a0>=0. && a1>=0.);
            HPFLTNB help2 = a0/c0-a1/c1;
            
            HPFLTNB temp_mm = (-log(mp_multiModalParam[mmnb]) - 0.5*( (help2*help2*help1)/(mp_multiModalNoiseSigma[mmnb]*mp_multiModalNoiseSigma[mmnb]) - log(help1) ));

            // if subsets, divide the MRI log influence by the number of subsets only if ordinary OSEM-like backprojection
            if (mp_nbSubsets[a_iteration]>1 && m_backprojection==0) temp_mm /= (HPFLTNB)mp_nbSubsets[a_iteration];

            // add the multimodal influence
            temp += temp_mm;

          }
        }
        assert (!std::isnan(temp));
        probabilities.push_back(temp);
        /*
        mean_proba += temp;
        nb_proba += 1.;
        std_proba += (temp-mean_proba_prev)*(temp-mean_proba_prev);
        */
      }
      else
      {
        // probability of a link which does not alter clusters, a link which loops into the same cluster (log(1))
        if (m_ddCRP==1)
          probabilities.push_back (0.);
        else if (m_ddCRP==2)
          probabilities.push_back(m_ddcrpLogAlpha);
      }
    }

    // apply exponential to log probabilities (tip to avoid numerical issues)
    HPFLTNB probabilities_sum=0.;
    for(size_t p=0; p<probabilities.size(); p++) 
    {
      probabilities[p] = exp(probabilities[p]);
      assert (!std::isnan(probabilities[p]));
      assert (!std::isinf(probabilities[p]));
      probabilities_sum += probabilities[p];
      assert (!std::isinf(probabilities_sum));
    }

    
    INTNB p=0;
    // sample the categorical distribution (multinomial distribution with a single repetition)
    // if at least one probability greater than the HPFLTNB minimum
    if (probabilities_sum>std::numeric_limits<HPFLTNB>::min())
    {
      HPFLTNB current_random = sRandomNumberGenerator::GetInstance()->GenerateExtraRdmNber(0)*probabilities_sum;
      HPFLTNB cumulative_sum = 0.;
      for(p=0; p<(INTNB)probabilities.size(); p++)
      {
        cumulative_sum += probabilities[p];
        if(cumulative_sum>current_random) break;
      }
    }
    // choose a self link in the (rare) cases when the probabilities are too low to be represented by HPFLTNB
    // or there are no or few fellow voxels (might occur with a mask)
    // this can happen only if alpha is zero and hasn't been added to the list of probabilities
    else
    {
      assert(!(m_ddcrpAlpha>0.));
      p = -1;
    }

    assert(p<(INTNB)probabilities.size());

    // if alpha is zero, the self link probability wasn't added in the list of probabilities,
    // so there is an index shift of 1 with respect to the list of fellow voxels
    if (!(m_ddcrpAlpha>0.)) p+=1;

    // take into account the new next link : update parent and next links
    int new_next_link = fellow_voxels.at(p);

    mpv_parentLinks[new_next_link].push_back(current_voxel);
    mp_nextLink[current_voxel] = new_next_link;

    // take into account the new next link : update clusters
    int new_next_link_cluster = mp_voxelClusterMapping[new_next_link];    
    // if the new link does not cause clusters fusion, there is nothing to do, the temporary new cluster
    // becomes the actual new cluster
    // if the new link causes clusters fusion, assign all the voxels to the cluster of the next link voxel
    // and release the temp_new_cluster as empty and free for assignment, and update sums per clusters
    if (new_next_link_cluster!=temp_new_cluster)
    {
      // update sums 
      mp_clusterN[0][new_next_link_cluster] += mp_clusterN[0][temp_new_cluster];
      mp_clusterSensitivity[new_next_link_cluster] += mp_clusterSensitivity[temp_new_cluster];
      
      assert(mp_clusterN[0][new_next_link_cluster]>=0.);
      assert(mp_clusterSensitivity[new_next_link_cluster]>=0.);
      
      if (multimodal_info)
      {
        mp_clusterCount[new_next_link_cluster] += mp_clusterCount[temp_new_cluster];
        assert(mp_clusterCount[new_next_link_cluster]>=0);
        for (int mmnb=0; mmnb<mp_ID->GetNbMultiModalImages(); mmnb++)
        {
          mpp_clusterMultiModal[mmnb][new_next_link_cluster] += mpp_clusterMultiModal[mmnb][temp_new_cluster];
          assert(mpp_clusterMultiModal[mmnb][new_next_link_cluster]>=0.);
        }
      }
      
      // release the temporary new cluster
      mv_newClusters.push_back(temp_new_cluster);
      // assign all the fusioned voxels to the cluster of the child voxel
      for(size_t fv=0; fv<voxels_after_cut.size(); fv++)
      {
        mp_voxelClusterMapping[voxels_after_cut.at(fv)] = new_next_link_cluster;
      }
    }
  }

  /*
  mean_proba/=nb_proba;
  mean_proba_prev = mean_proba;
  std_proba = sqrt(std_proba/nb_proba);
  Cout("iRCPGSAlgorithm::SampleConditionalClustering() -> fusion probabilities mean "<< mean_proba<<" std "<<std_proba<<endl);
  */

  // End
  return 0;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

int iRCPGSAlgorithm::SampleConditionalClusterIntensity()
{
  if (m_verbose>=2) Cout("iRCPGSAlgorithm::SampleConditionalClusterIntensity()... " << endl);
  
  // reinitialize clusters values to -1, 
  // trick for ensuring that the value of each cluster is sampled only once
  int vp;
  #pragma omp parallel for private(vp) schedule(guided)
  for (vp=0; vp<(mp_ID->GetNbVoxXYZ()+1); vp++)  mp_clusterValues[vp] = -1.;

  INTNB nb_clusters = 0; INTNB count_empty = 0; HPFLTNB maxN = 0.; INTNB nb_relevant_voxels = 0;
  // sample voxel intensity for each cluster, using the updated Gamma distribution
  for (int v=0; v<mp_ID->GetNbVoxXYZ(); v++)
  {
    // skip if background voxel
    if (!mp_listRelevantVoxelIndices[v]) continue;
    
    nb_relevant_voxels ++;
    
    int index = mp_voxelClusterMapping[v];
    
    if (m_ddCRP==0) assert(index==v);
    
    if (mp_clusterValues[index]<0.)
    {
      gamma_distribution<HPFLTNB> updated_gamma_distrib(mp_clusterN[0][index] + m_gammaShape, 1./(mp_clusterSensitivity[index] + m_gammaRate));
      mp_clusterValues[index] = updated_gamma_distrib(sRandomNumberGenerator::GetInstance()->GetExtraGenerator(1));
      
      assert(mp_clusterN[0][index]>=0. && mp_clusterSensitivity[index]>=0.);
      assert(!std::isnan(mp_clusterValues[index]) && !std::isinf(mp_clusterValues[index]) && mp_clusterValues[index]>=0.);
      
      nb_clusters++;
      if (mp_clusterValues[index]>maxN) maxN = mp_clusterValues[index];
      if (!(mp_clusterN[0][index]>0.)) count_empty++;
    }
  }
  
  // compute the current average cluster volume
  m_currentMeanClusterVolume = mp_ID->GetVoxSizeX() * mp_ID->GetVoxSizeY() * mp_ID->GetVoxSizeZ()* (FLTNB)nb_relevant_voxels / (FLTNB)nb_clusters;

  if (m_verbose>=2)
  {
    Cout("  --> Mean number of voxels per cluster = "<<((FLTNB)nb_relevant_voxels/(FLTNB)nb_clusters)<<endl);
    Cout("  --> Average cluster volume =  "<<m_currentMeanClusterVolume<<" mm3"<<endl); 
  }
  
  // End
  return 0;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

int iRCPGSAlgorithm::ComputeFellowVoxelsList(vector<INTNB>& a_fellow_voxels, int a_current_voxel)
{
  // X,Y,Z dimensions, ordered according to the array storage :
  // X is the innermost and Z is the outermost dimension
  INTNB dim_3D = mp_ID->GetNbVoxXYZ();
  INTNB dim_X = mp_ID->GetNbVoxX();
  INTNB dim_Z = mp_ID->GetNbVoxZ();
  INTNB dim_XY = mp_ID->GetNbVoxXY();
  
  // build the list of fellow voxels for the current voxel
  a_fellow_voxels.resize(0);
  // add first the current voxel itself
  a_fellow_voxels.push_back(a_current_voxel);
  
  if (dim_Z==1)
  {
    // build the list of 2D neighbours excluding diagonal voxels
    const int ix = a_current_voxel%dim_X;
    bool y_low = (a_current_voxel>=dim_X) && mp_listRelevantVoxelIndices[a_current_voxel-dim_X];
    bool y_high = (a_current_voxel<dim_XY-dim_X) && mp_listRelevantVoxelIndices[a_current_voxel+dim_X];
    bool x_low = (ix>0) && mp_listRelevantVoxelIndices[a_current_voxel-1];
    bool x_high = (ix<dim_X-1) && mp_listRelevantVoxelIndices[a_current_voxel+1];

    if (y_low) a_fellow_voxels.push_back(a_current_voxel-dim_X);
    if (y_high) a_fellow_voxels.push_back(a_current_voxel+dim_X);
    if (x_low) a_fellow_voxels.push_back(a_current_voxel-1);
    if (x_high) a_fellow_voxels.push_back(a_current_voxel+1);
    
    for(size_t i=0;i<a_fellow_voxels.size();i++) assert(a_fellow_voxels.at(i)>=0 && a_fellow_voxels.at(i)<dim_3D );
    assert(a_fellow_voxels.size()>0);
  }
  else
  {
    // build the list of 3D neighbours excluding diagonal voxels
    const int current_voxel_xy = a_current_voxel%(dim_XY);
    // perform the same processing as for 2D on current_voxel_xy
    const int ix = current_voxel_xy%dim_X;
    bool y_low = (current_voxel_xy>=dim_X) && mp_listRelevantVoxelIndices[a_current_voxel-dim_X];
    bool y_high = (current_voxel_xy<dim_XY-dim_X) && mp_listRelevantVoxelIndices[a_current_voxel+dim_X];
    bool x_low = (ix>0) && mp_listRelevantVoxelIndices[a_current_voxel-1];
    bool x_high = (ix<dim_X-1) && mp_listRelevantVoxelIndices[a_current_voxel+1];
    // check z dimension
    bool z_low = (a_current_voxel >= dim_XY) && mp_listRelevantVoxelIndices[a_current_voxel-dim_XY];
    bool z_high = (a_current_voxel < dim_3D-dim_XY) && mp_listRelevantVoxelIndices[a_current_voxel+dim_XY];

    if (y_low) a_fellow_voxels.push_back(a_current_voxel-dim_X);
    if (y_high) a_fellow_voxels.push_back(a_current_voxel+dim_X);
    if (x_low) a_fellow_voxels.push_back(a_current_voxel-1);
    if (x_high) a_fellow_voxels.push_back(a_current_voxel+1);
    if (z_low) a_fellow_voxels.push_back(a_current_voxel-dim_XY);
    if (z_high) a_fellow_voxels.push_back(a_current_voxel+dim_XY);
    
    for(size_t i=0;i<a_fellow_voxels.size();i++) assert(a_fellow_voxels.at(i)>=0 && a_fellow_voxels.at(i)<dim_3D );
    assert(a_fellow_voxels.size()>0);
  }
  
  // End
  return 0;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

int iRCPGSAlgorithm::GenerateCurrentImage()
{
  // Loop on voxels with OpenMP
  int vp;
  #pragma omp parallel for private(vp) schedule(guided)
  for (vp=0; vp<mp_ID->GetNbVoxXYZ(); vp++)
  {
    if (!mp_listRelevantVoxelIndices[vp])
    {
       mp_ImageSpace->m4p_image[0][0][0][vp] = 0.;
    }
    else
    {
      // Compute output image value
      mp_ImageSpace->m4p_image[0][0][0][vp] = mp_clusterValues[mp_voxelClusterMapping[vp]];
    }
  }
  
  // End
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iRCPGSAlgorithm::StepAfterIterationLoop()
{
  if (vAlgorithm::StepAfterIterationLoop())
  {
    Cerr("***** iRCPGSAlgorithm::StepAfterIterationLoop() -> A problem occurred while calling StepAfterIterationLoop() function !" << endl);
    return 1;
  }
  if (m_verbose>=2) Cout("iRCPGSAlgorithm::StepAfterIterationLoop ... " << endl);

  mp_ImageSpace->DeallocateBackwardImageFromDynamicBasis();

  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iRCPGSAlgorithm::StepBeforeSubsetLoop(int a_iteration)
{
  (void)a_iteration; // avoid 'unused parameter' compil. warnings
  if (m_verbose>=3) Cout("iRCPGSAlgorithm::StepBeforeSubsetLoop ... " << endl);

  // End
  return 0;
}

// -----------------------------------------------------------------------------------------------------------------------------------------

int iRCPGSAlgorithm::StepAfterSubsetLoop(int a_iteration)
{
  if (m_verbose>=3) Cout("iRCPGSAlgorithm::StepAfterSubsetLoop() -> Clean never visited voxels and save images if needed" << endl);
  // Clean never visited voxels
  mp_ImageSpace->CleanNeverVisitedVoxels();
  // Save the main image
  if (mp_ID->GetMPIRank()==0 && mp_outputIterations[a_iteration])
  {
    // Verbose
    if (m_verbose>=1) Cout("iRCPGSAlgorithm::StepAfterSubsetLoop() -> Save image at iteration " << a_iteration+1 << endl);
    // Save image
    if (StepImageOutput(a_iteration))
    {
      Cerr("***** iRCPGSAlgorithm::StepAfterSubsetLoop() -> A problem occurred while saving images at iteration " << a_iteration+1 << " !" << endl);
      return 1;
    }
  }

  // End
  return 0;
}


// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

int iRCPGSAlgorithm::StepImageOutput(int a_iteration, int a_subset)
{
  // =================================================================================================
  // Apply pre-processing steps
  // =================================================================================================
  // Copy the current image into the forward image
  mp_ImageSpace->PrepareForwardImage();

  // Apply output transaxial FOV masking
  if (mp_ImageSpace->ApplyOutputFOVMasking())
  {
    Cerr("***** iRCPGSAlgorithm::StepImageOutput() -> A problem occurred while applying output FOV masking !" << endl);
    return 1;
  }
  // Apply output flip
  if (mp_ImageSpace->ApplyOutputFlip())
  {
    Cerr("***** iRCPGSAlgorithm::StepImageOutput() -> A problem occurred while applying output flip !" << endl);
    return 1;
  }

  // =================================================================================================
  // Save frames/gates
  // =================================================================================================

  // Save output image (note that if no basis functions at all are in use, the the output image already points to the forward image)
  if (mp_ImageSpace->SaveOutputImage(a_iteration, a_subset))
  {
    Cerr("***** iRCPGSAlgorithm::StepImageOutput() -> A problem occurred while saving output image !" << endl);
    return 1;
  }

  if (((a_iteration+1)%50)==0)
  {
    // Get the output manager
    sOutputManager* p_output_manager = sOutputManager::GetInstance();

    // ----------------------
    // Build the file name
    // ----------------------

    string data_file = p_output_manager->GetPathName() + p_output_manager->GetBaseName()+"_intermediary";
    // Add a suffix for iteration
    if (a_iteration >= 0)
    {
      stringstream ss; ss << a_iteration + 1;
      data_file.append("_it").append(ss.str());
    }
    // Add extension
    data_file.append(".img");

    // ----------------------
    // Write image in file
    // ----------------------

    // Open file
    FILE* fout = fopen(data_file.c_str(),"wb");
    if (fout==NULL)
    {
      Cerr("***** iRCPGSAlgorithm::StepImageOutput() -> Failed to create output file '" << data_file << "' !" << endl);
      return 1;
    }

    // Write the file by converting the data to the output type, currently no specific output type, 
    // so the current data type is used 
    INTNB nb_data = 0;
    for (INTNB v=0; v<mp_ID->GetNbVoxXYZ(); v++)
    {
      INTNB buffer = -1;
      if (mp_listRelevantVoxelIndices[v]) buffer = ((INTNB)(mp_voxelClusterMapping[v]));
      nb_data += fwrite(&buffer,sizeof(INTNB),1,fout);
    }

    // Close file
    fclose(fout);

    // Check writing
    if (nb_data!=mp_ID->GetNbVoxXYZ())
    {
      Cerr("***** iRCPGSAlgorithm::StepImageOutput() -> Failed to write all data into the output file '" << data_file << "' !" << endl);
      return 1;
    }
  }

  // Normal end
  return 0;
}

// =====================================================================
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// =====================================================================

void iRCPGSAlgorithm::ShowHelpSpecific()
{
      cout << "" << endl;
      cout << "Random Clustering Prior - Gibbs Sampler : probabilistic Bayesian inference algorithm (no optimization)" << endl;
      cout << "M. Filipovic et al : PET reconstruction of the posterior image probability, including multimodal images." << endl;
      cout << "" << endl;
      cout << "The algorithm is launched using -prob, either with command line options (-prob RCPGS:option1name=option1value,option2name=option2value,...)" << endl;
      cout << "or with a config file (-prob RCPGS:config_file.conf), see also the default config file for examples of option values." << endl;
      cout << "" << endl;
      cout << "Required options :" << endl;
      cout << "  ddCRP = type of ddCRP prior (0 = no ddCRP, 1 = original ddCRP (recommended), 2 = modified ddCRP)" << endl;
      cout << "  alpha = ddCRP hyperparameter, the unnormalized probability of drawing a self link" << endl;
      cout << "  gammaShape = Gamma prior distribution shape parameter (0.5 recommended)" << endl;
      cout << "  gammaRate = Gamma prior distribution rate parameter (1.e-18 recommended)" << endl;
      cout << "  backprojection = type of multinomial backprojection of the current iteration/subset data (0 = the previous backprojection state is cleared, as for ML-EM (recommended), 1 = update of the previous backprojection state, 2 = backprojection marginalized over cluster intensity, implies update of the previous backprojection state )" << endl;
      cout << "  multiModalNoiseSigma = standard deviation of Gaussian noise in the provided additional image from another modality, must be repeated as many times as the -multimodal option" << endl;
      cout << "  multiModalLag = number of iterations after which the multimodal images start affecting voxels clustering" << endl;
      cout << "" << endl;
      cout << "Optional options (if one of these is specified, all the others must be specified as well):" << endl;
      cout << "  meanClusterVolumeMin = minimum threshold for average cluster volume (in mm3), used for tuning ddCRP alpha automatically through iterations" << endl;
      cout << "  meanClusterVolumeMax = maximum threshold for average cluster volume (in mm3), used for tuning ddCRP alpha automatically through iterations" << endl;
      cout << "  alphaIncrement = multiplicative increment for tuning ddCRP alpha through iterations" << endl;
      cout << "" << endl;
}

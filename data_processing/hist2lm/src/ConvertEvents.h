
#include <fstream>
#include <iostream>

/**
 * This class is responsible for converting events
 * from J-PET framework format into CASToR format.
 * It requires already generated CASToR LUT file
 * describing used geometry. Additionally it reqires
 * information about number of crystals in LUT file and
 * crystals size.
 *
 * As an input it takes uncategorized events (and it only
 * takes events that have exacly 2 hits).
 *
 * As an output it generates *.Cdf and .Cdh files
 * that can be used in CASToR for reconstruction.
 *
 * CASToR general documentation about LUT file and CASToR
 * format: https://castor-project.org/sites/default/files/2020-09/CASToR_general_documentation.pd**/

  std::ofstream fOutputStream;
  std::vector<std::tuple<float, float, float>> fCastorIDs;
  uint32_t fNumberOfEvents   = 0;
  uint32_t fNumberOfCrystals = 62400;
  float fCrystalSizeX = 0.f;
  float fCrystalSizeY = 0.f;
  float fCrystalSizeZ = 0.f;
  int fTOFFWHM             = 600;
  int fTOFMeasurementRange = 1500;
  std::string fInputFilePath  = "scanner";
  std::string fOutputFileBaseName;
  std::string fScannerName    = "scanner";
  const std::string kNumberOfCrystalsKey    = "62400";
  

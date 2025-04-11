#include "ConvertEventsTools.h"
#include "JPetLoggerInclude.h"
#include <tuple>
#include <cmath>
#include <limits>
#include <string>

bool ConvertEventsTools::readLUTFile(const std::string& path, int numberOfCrystals, std::vector<std::tuple<float, float, float>> &readedLutGeometry) {
  FILE* LUT_file = fopen(path.c_str(), "rb");
  if (LUT_file==NULL)
  {
    ERROR("Failed to open " << path << " file, aborting!");
    return false;
  }

  // Read data for each index
  int nb_data_read = 0;
  for (int i=0; i<numberOfCrystals; i++)
  {
    float x, y, z;
    float skip;
    // Read central crystal position X, then Y and Z
    nb_data_read += fread(&x,sizeof(x),1,LUT_file);
    nb_data_read += fread(&y,sizeof(y),1,LUT_file);
    nb_data_read += fread(&z,sizeof(z),1,LUT_file);
    // Read crystal orientation X, then Y and Z
    nb_data_read += fread(&skip,sizeof(skip),1,LUT_file);
    nb_data_read += fread(&skip,sizeof(skip),1,LUT_file);
    nb_data_read += fread(&skip,sizeof(skip),1,LUT_file);
    readedLutGeometry.push_back(std::make_tuple(x, y, z));
  }

  // Close file
  fclose(LUT_file);
  // Check reading
  if (nb_data_read!=numberOfCrystals*6) {
    ERROR("Failed to read correct number of crystals from lut file!");
    ERROR("Readed: " << nb_data_read << " Expected: " << numberOfCrystals * 6);
    return false;
  }
  return true;
}

/**
 * This function checks every crystal from LUT file and calculates
 * distance to hit position. This function assumes, that lowest distance
 * is crystal that we are looking for.
 *
 * TODO: add maximal allowed distance
 * TODO: add warning/error when distance is relative big
 * TODO: take orientation of crystal into account
 *
 * @param x X position (in mm) of hit
 * @param y Y position (in mm) of hit
 * @param z Z position (in mm) of hit
 * @param castorIDs vector containing LUT positions
 * @param scintillatorLength length of the scintillator on Z axis
 * @param crystalSizeZ size of crystal on Z axis
 *
 * @return ID of hit in CASToR LUT file
 */
uint32_t ConvertEventsTools::getCastorID(float x, float y, float z, const std::vector<std::tuple<float, float, float>> &castorIDs, const float scintillatorLength, const float crystalSizeZ) {
  // Calculate number of crystals on Z axis
  const int numberOfBinsInZ = std::round(scintillatorLength / crystalSizeZ);
  float lowestDistance2d = std::numeric_limits<float>::max();
  float lowestDistance = std::numeric_limits<float>::max();
  uint32_t nearestCrystal = 0;
  // First we want to find nearest position using only x and y position
  // In LUT we are saving data by first traverse along Z axis
  // so here we are checking only bins that have different x/y position
  // and skipping whole Z part
  for (int i = 0; i < castorIDs.size(); i = i + numberOfBinsInZ) {
    float diffX = std::get<0>(castorIDs[i]) - x;
    float diffY = std::get<1>(castorIDs[i]) - y;
    float distance = std::sqrt(diffX * diffX + diffY * diffY);
    if (distance < lowestDistance2d) {
      lowestDistance2d = distance;
      nearestCrystal = i;
    }
  }
  // After first initial loop
  // we are checking for each crystal along Z axis
  int currentCrystal = nearestCrystal;
  int endCrystal = nearestCrystal + numberOfBinsInZ;
  for (;currentCrystal < endCrystal; currentCrystal++) {
    float diffX = std::get<0>(castorIDs[currentCrystal]) - x;
    float diffY = std::get<1>(castorIDs[currentCrystal]) - y;
    float diffZ = std::get<2>(castorIDs[currentCrystal]) - z;
    float distance = std::sqrt(diffX * diffX + diffY * diffY + diffZ * diffZ);
    if (distance < lowestDistance) {
      lowestDistance = distance;
      nearestCrystal = currentCrystal;
    }
  }
  return nearestCrystal;
}

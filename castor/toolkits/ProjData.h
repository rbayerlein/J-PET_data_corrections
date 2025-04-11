#ifndef PROJDATA_H
#define PROJDATA_H

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>

#include "gDataConversionUtilities.hh"
#include "gOptions.hh"

#include "Segment.h"

class ProjData
{
  public:
    ProjData();
    ~ProjData();

    bool readSTIRProjData(const std::string& path);
    float getBinValue(uint32_t segment_id, uint32_t sinogram_id, uint32_t view_id, uint32_t tangential_id);
    float getScatterNumber(float posX1, float posY1, float posZ1, float posX2, float posY2, float posZ2);
    void printInfoAboutProjData();

  private:
    void readSTIRProjBinaryData();
    void getSublineAfterTheSeparator(std::string &line, std::string &sep);
    float get2DAngle(float posX1, float posY1, float posX2, float posY2);
    void calculatePointLineDistance(float x1, float y1, float x2, float y2, float &line_to_center_point_distance, float &directional_parameter);

    int segments_no;
    int views_no;
    int tangential_coord_no;
    int rings_no;
    float effective_bin_size; //in cm
    float distance_between_rings; //in cm
    float axial_fov; //in cm
    float tangential_fov; // in cm
    string header_path = "";
    string binary_path = "";

    std::vector<int> segments_sizes;
    std::vector<Segment> segments;



};
#endif

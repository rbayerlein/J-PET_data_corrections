#ifndef SEGMENT_H
#define SEGMENT_H

#include <vector>
#include <string>

#include "gDataConversionUtilities.hh"

class Segment

{
  public:
    Segment();
    Segment(std::vector<std::vector<std::vector<float>>> &seg_data);
    ~Segment();
    float getBinValue(uint32_t tangential_coord_id, uint32_t sinoram_id, uint32_t view_id);

  private:
    std::vector<std::vector<std::vector<float>>> segment_data;
    uint32_t sinograms_no, views_no, tangential_coord_no;

};
#endif

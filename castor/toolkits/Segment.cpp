#include "Segment.h"

Segment::Segment() {}

Segment::~Segment() {}

Segment::Segment(std::vector<std::vector<std::vector<float>>> &seg_data) 
{
  tangential_coord_no = seg_data.size();
  sinograms_no = seg_data[0].size();
  views_no = seg_data[0][0].size();
  segment_data = seg_data;
}

float Segment::getBinValue(uint32_t tangential_coord_id, uint32_t sinoram_id, uint32_t view_id)
{
  return segment_data[tangential_coord_id][sinoram_id][view_id];
}

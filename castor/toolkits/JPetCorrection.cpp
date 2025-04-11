#include "JPetCorrection.h"
#include <cmath>
#include "gVariables.hh"
#include "gOptions.hh"


JPetCorrection::JPetCorrection() {}

JPetCorrection::~JPetCorrection() {}

bool JPetCorrection::readCorrectionMatrixFile(const std::string& path, int number_of_correction_bins) { // Maybe better to move to constructor?

  if (path.empty()) {
    return false;
  }

  FILE* correction_matrix_file = fopen(path.c_str(), "rb");
  if (correction_matrix_file==NULL) {
    return false;
  }

  std::vector<float> correction_row;
  correction_row.reserve(number_of_correction_bins);
  int nb_data_read = 0;
  float correction_factor;
  for (int i = 0; i < number_of_correction_bins; i++) {
    for (int j = 0; j < number_of_correction_bins; j++) {
      nb_data_read += fread(&correction_factor, sizeof(correction_factor), 1, correction_matrix_file);
      correction_row.push_back(correction_factor);
    }
    correction_matrix->push_back(correction_row);
    correction_row.clear();
  }

  // Close file
  fclose(correction_matrix_file);
  correction_matrix_file = nullptr;
  // Check reading
  if (nb_data_read!=number_of_correction_bins*number_of_correction_bins) {
    return false;
  }
  return true;
}

float JPetCorrection::getCorrection(uint32_t castor_crystal_id1, uint32_t castor_crystal_id2) const {
  return correction_matrix->at(castor_crystal_id1).at(castor_crystal_id2);
}


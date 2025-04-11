#ifndef JPETCORRECTION_H
#define JPETCORRECTION_H

#include <string>
#include <vector>
#include <cassert>

class JPetCorrection
{
  public:
    JPetCorrection();
    ~JPetCorrection();
    bool readCorrectionMatrixFile(const std::string& path, int number_of_correction_bins);
    float getCorrection(uint32_t castorCrystalId1, uint32_t castorCrystalId2) const;
  private:
    uint32_t matrix_rows;
    int tof_bins;
    std::vector<std::vector<float>>* correction_matrix = new std::vector<std::vector<float>>();
};

struct CorrectionMatrix
{
  virtual uint64_t getCorrection(uint64_t castorID1,
                                 uint64_t castorID2) const = 0;

  virtual void incrementCorrection(uint64_t castorID1,
                                   uint64_t castorID2) = 0;

  virtual ~CorrectionMatrix() = default;
};

class GateCorrectionMatrix : public CorrectionMatrix
{
 public:
  GateCorrectionMatrix(uint32_t bins_elts, uint32_t nCrystalsTot): fBinsElts(bins_elts)
  {
      p_corr = new uint8_t*[bins_elts];
      for (size_t c = 0; c < bins_elts; c++)
      {
        p_corr[c] = new uint8_t[nCrystalsTot-c];

        for (size_t c2 = 0; c2 < nCrystalsTot - c; c2++)
          p_corr[c][c2] = 0;
      }
  }

  uint64_t getCorrection(uint64_t castorID1, uint64_t castorID2) const override
  {
    return p_corr[castorID1][castorID2];
  }

  void incrementCorrection(uint64_t castorID1, uint64_t castorID2) override
  {
    p_corr[castorID1][castorID2]++;
  }

  virtual ~GateCorrectionMatrix() {
    for(size_t b = 0; b < fBinsElts; b++)
      delete[] p_corr[b];

    delete[] p_corr;
  }

 private:
  uint8_t **p_corr = nullptr;
  uint32_t fBinsElts = 0;
};

class JPetCorrectionMatrix : public CorrectionMatrix
{
 public:
  JPetCorrectionMatrix(const std::string& path, const int numberOfCrystals): fCorrection_matrix()
  {
      fCorrection_matrix.readCorrectionMatrixFile(path, numberOfCrystals); //TODO: we should check if it returns false or not
  }

  uint64_t getCorrection(uint64_t moduleID1, uint64_t moduleID2) const override
  {
    return fCorrection_matrix.getCorrection(moduleID1, moduleID2);
  }

  uint64_t calculateModuleID(int crystalID, int rsectorID, int layerID) {
    uint64_t moduleID = crystalID * fRSectors * fPseudoUniqueLayerID +
                        rsectorID *             fPseudoUniqueLayerID +
                        int(layerID/(fStripsInModule * fPseudoCrystalsInZ)) * fStripsInModule +
                        layerID % fStripsInModule;
    return moduleID;
  }

  uint64_t getCorrection(int crystalID1, int rsectorID1, int layerID1,
                         int crystalID2, int rsectorID2, int layerID2) {
    return getCorrection(calculateModuleID(crystalID1, rsectorID1, layerID1),
                         calculateModuleID(crystalID2, rsectorID2, layerID2));
  }

  void incrementCorrection(uint64_t, uint64_t) override
  {
    return; //Do nothing
  }

 private:
  JPetCorrection fCorrection_matrix;
  int fRSectors = 24;
  int fStripsInModule = 16;
  int fPseudoCrystalsInZ = 286;
  int fPseudoUniqueLayerID = fStripsInModule * 2;
};

#endif

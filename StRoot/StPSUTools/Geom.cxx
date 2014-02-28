#include "Geom.h"

#include "St_base/StMessMgr.h"
#include "StFmsDbMaker/StFmsDbMaker.h"
#include "tables/St_fmsDetectorPosition_Table.h"
#include "tables/St_fmsChannelGeometry_Table.h"

using namespace std;
using namespace PSUGlobals;

ClassImp(Geom)

Geom::Geom() : FMSGeom(true) {
  InitDBGeom();
}

Geom::~Geom() { }

fmsDetectorPosition_st* Geom::find(int detectorId) {
  fmsDetectorPosition_st* positions(NULL);
  Table::iterator entry = mPositions.find(detectorId);
  if (entry != mPositions.end()) {
    positions = entry->second;
  }  // if
  return positions;
}

int Geom::ewNstbToDetectorId(int ew, int nstb) {
  int id(-1);
  // Only support FMS (ew == 2)
  // detector IDs are defined in the range [1, 14], with the 4 FMS subdetectors
  // being [9, 12], in the same order as nstb [1, 4]
  if (ew == 2 && nstb > 0 && nstb < 5) {
    id = nstb + 7;
  }  // if
  return id;
}

Float_t* Geom::ZFPD(Int_t ew, Int_t nstb) {
  fmsDetectorPosition_st* geometry = find(ewNstbToDetectorId(ew, nstb));
  if (geometry) {
    return &geometry->zoffset;
  }  // if
  return NULL;
}

Float_t* Geom::xOffset(Int_t ew, Int_t nstb) {
  fmsDetectorPosition_st* geometry = find(ewNstbToDetectorId(ew, nstb));
  if (geometry) {
    return &geometry->xoffset;
  }  // if
  return NULL;
}

Float_t* Geom::yOffset(Int_t ew, Int_t nstb) {
  fmsDetectorPosition_st* geometry = find(ewNstbToDetectorId(ew, nstb));
  if (geometry) {
    return &geometry->yoffset;
  }  // if
  return NULL;
}

Float_t* Geom::FpdTowWid(Int_t ew, Int_t nstb) {
  // I don't like this implementation, returning a pointer to access two floats
  // It relies on the data being aligned OK and seems dangerous. We should add
  // a more robust solution e.g. return a pair or 2-element vector.
  fmsDetectorPosition_st* geometry = find(ewNstbToDetectorId(ew, nstb));
  if (geometry) {
    return &geometry->xwidth;
  }  // if
  return NULL;
}

bool Geom::InitDBGeom() {
	StFmsDbMaker* fmsDbMaker = static_cast<StFmsDbMaker*>(
	  StMaker::GetChain()->GetMaker("fmsDb"));
  if (!fmsDbMaker) {
    LOG_ERROR << "Geom unable to locate an StFmsDbMaker - geometry will not "
      << "be initialised!" << endm;
    return false;
  }  // if
  fmsDetectorPosition_st* dbgeom = fmsDbMaker->DetectorPosition();
  if (dbgeom) {
    // Detector IDs count [0, N), so number of detectors is one greater than max
    const int nDetectors = fmsDbMaker->maxDetectorId() + 1;
    for (int i(0); i < nDetectors; ++i) {
      // The first detector has ID = 0. Subsequent detectors with no data will
      // also have detector ID = 0. This won't overwrite the first entry, as
      // map::insert won't insert if there is already an entry with that key.
      mPositions.insert(std::make_pair(dbgeom[i].detectorId, &dbgeom[i]));
    }  // for
    return true;
  }  // if
  return false;
}

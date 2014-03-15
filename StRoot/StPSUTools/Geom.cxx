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

const fmsDetectorPosition_st* Geom::find(Int_t detectorId) const {
  const fmsDetectorPosition_st* positions(NULL);
  Table::const_iterator entry = mPositions.find(detectorId);
  if (entry != mPositions.end()) {
    positions = entry->second;
  }  // if
  return positions;
}

const Float_t* Geom::ZFPD(Int_t detectorId) const {
  const fmsDetectorPosition_st* geometry = find(detectorId);
  if (geometry) {
    return &geometry->zoffset;
  }  // if
  return NULL;
}

const Float_t* Geom::xOffset(Int_t detectorId) const {
  const fmsDetectorPosition_st* geometry = find(detectorId);
  if (geometry) {
    return &geometry->xoffset;
  }  // if
  return NULL;
}

const Float_t* Geom::yOffset(Int_t detectorId) const {
  const fmsDetectorPosition_st* geometry = find(detectorId);
  if (geometry) {
    return &geometry->yoffset;
  }  // if
  return NULL;
}

const Float_t* Geom::FpdTowWid(Int_t detectorId) const {
  // I don't like this implementation, returning a pointer to access two floats
  // It relies on the data being aligned OK and seems dangerous. We should add
  // a more robust solution e.g. return a pair or 2-element vector.
  const fmsDetectorPosition_st* geometry = find(detectorId);
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

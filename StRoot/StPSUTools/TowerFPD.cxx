#include "StPSUTools/TowerFPD.h"

#include "StEvent/StFmsHit.h"
#include "StFmsDbMaker/StFmsDbMaker.h"

namespace PSUGlobals {
TowerFPD::TowerFPD()
    : mHit(NULL), mColumn(-1), mRow(-1), mCluster(-1) { }

TowerFPD::TowerFPD(const StFmsHit* fmsHit)
    : mHit(fmsHit), mColumn(-1), mRow(-1), mCluster(-1) { }

TowerFPD::~TowerFPD() { }

Bool_t TowerFPD::initialize(StFmsDbMaker* database) {
  if (!mHit || !database) {  // Check for invalid input
    return false;
  }  // if
  // Get row and column from the database
  mRow = database->getRowNumber(mHit->detectorId(), mHit->channel());
  mColumn = database->getColumnNumber(mHit->detectorId(), mHit->channel());
  return mRow > -1 && mColumn > -1;
}

Int_t TowerFPD::Compare(const TObject* tower) const {
  const TowerFPD* other = static_cast<const TowerFPD*>(tower);
  if (mHit->energy() < other->hit()->energy()) {
    return -1;
  } else if (mHit->energy() > other->hit()->energy()) {
    return 1;
  } else {
    return 0;
  }  // if
}

Bool_t TowerFPD::IsNeighbor(TowerFPD* other) {
  if (!other) {
    return false;
  }  // if
  return abs(mColumn - other->column()) + abs(mRow - other->row()) == 1;
}
}  // namespace PSUGlobals

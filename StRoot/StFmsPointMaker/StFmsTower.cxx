#include "StPSUTools/StFmsTower.h"

#include "StEvent/StFmsHit.h"
#include "StFmsDbMaker/StFmsDbMaker.h"

namespace FMSCluster {
StFmsTower::StFmsTower()
    : mHit(NULL), mColumn(-1), mRow(-1), mCluster(-1) { }

StFmsTower::StFmsTower(StFmsHit* fmsHit)
    : mHit(fmsHit), mColumn(-1), mRow(-1), mCluster(-1) { }

StFmsTower::~StFmsTower() { }

Bool_t StFmsTower::initialize(StFmsDbMaker* database) {
  if (!mHit || !database) {  // Check for invalid input
    return false;
  }  // if
  // Get row and column from the database
  mRow = database->getRowNumber(mHit->detectorId(), mHit->channel());
  mColumn = database->getColumnNumber(mHit->detectorId(), mHit->channel());
  return mRow > -1 && mColumn > -1;
}

Bool_t StFmsTower::isNeighbor(const StFmsTower& other) const {
  return abs(mColumn - other.column()) + abs(mRow - other.row()) == 1;
}
}  // namespace FMSCluster

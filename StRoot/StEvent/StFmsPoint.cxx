#include "StFmsPoint.h"

ClassImp(StFmsPoint)

StFmsPoint::StFmsPoint() : StObject() {
  mDetectorId = 0;
  mEnergy = -1.0;
  mX = -99.0;
  mY = -99.0;
  mId = -1;
}

StFmsPoint::~StFmsPoint() {}

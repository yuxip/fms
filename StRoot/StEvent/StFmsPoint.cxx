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

ostream& operator<<(ostream& os, const StFmsPoint& v) {
  return os << "StFmsPoint:\n\tenergy: " << v.energy()
                  << "\n\txpos:   " << v.x()
                  << "\n\typos:   " << v.y()
                  << "\n\tphId:   " << v.id()
                  << "\n";
}

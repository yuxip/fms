#include "StMessMgr.h"
#include "StFmsPoint.h"

ClassImp(StFmsPoint)

StFmsPoint::StFmsPoint() : StObject() {
  mEnergy = -1.0;
  mX = -99.0;
  mY = -99.0;
  mId = -1;
}

StFmsPoint::StFmsPoint(StFmsPoint& other) {
  LOG_DEBUG << "StFmsPoint copy constructor called " << endm;
  mEnergy = other.energy();
  mX = other.x();
  mY = other.y();
  mId = other.id();
  mFourMomentum = other.fourMomentum();
  mXYZLab = other.xyzLab();
}

StFmsPoint::~StFmsPoint() {}

void StFmsPoint::Print(const Option_t* opt) const {
  LOG_INFO << *this << endm;
}

ostream& operator<<(ostream& os, const StFmsPoint& v) {
  return os << "StFmsPoint:\n\tenergy: " << v.energy()
                  << "\n\txpos:   " << v.x()
                  << "\n\typos:   " << v.y()
                  << "\n\tphId:   " << v.id()
                  << "\n";
}

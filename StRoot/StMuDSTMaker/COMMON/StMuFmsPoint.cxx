// $Id$
//
// $Log$
/**
 \file      StMuFmsPoint.cxx
 \brief     Implementation of StMuFmsPoint, the MuDST FMS "point" class
 \author    Thomas Burton <tpb@bnl.gov>
 \date      2014
 \copyright Brookhaven National Lab
 */

#include "StMuFmsPoint.h"

#include <algorithm>  // For std::min
#include <cmath>

#include "StRoot/StEvent/StFmsPoint.h"
#include "StRoot/StMuDSTMaker/COMMON/StMuFmsCluster.h"

StMuFmsPoint::StMuFmsPoint(int detectorId, float energy,
                           float x, float y, float z)
    : mDetectorId(detectorId), mEnergy(energy), mX(x), mY(y), mZ(z) { }

StMuFmsPoint::StMuFmsPoint(const StFmsPoint& point) {
  set(point);
}

StMuFmsPoint::~StMuFmsPoint() { }

TVector3 StMuFmsPoint::momentum(float m) const {
  m = std::min(m, mEnergy);  // Prevent m > E
  TVector3 v(mX, mY, mZ);
  if (std::fabs(m) > 0.f) {
    v.SetMag(std::sqrt(std::pow(mEnergy, 2.f) - std::pow(m, 2.f)));
  } else {
    v.SetMag(mEnergy);
  }  // if
  return v;
}

TLorentzVector StMuFmsPoint::fourMomentum(float m) const {
  return TLorentzVector(momentum(m), mEnergy);
}

StMuFmsCluster* StMuFmsPoint::cluster() {
  return static_cast<StMuFmsCluster*>(mCluster.GetObject());
}

const StMuFmsCluster* StMuFmsPoint::cluster() const {
  return static_cast<const StMuFmsCluster*>(mCluster.GetObject());
}

void StMuFmsPoint::set(const StFmsPoint& point) {
  mDetectorId = point.detectorId();
  mEnergy = point.energy();
  mX = point.x();
  mY = point.y();
  // Calculate z coordinate from StFmsPoint 4-momentum as it doesn't store
  // z directly. z / x = pz / px, so...
  const TLorentzVector vec4 = point.fourMomentum();
  mZ = point.x() * vec4.Pz() / vec4.Px();
}

void StMuFmsPoint::setCluster(StMuFmsCluster* cluster) {
  mCluster = cluster;
}

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

#include <TLorentzVector.h>

#include "StRoot/StEvent/StFmsPoint.h"

StMuFmsPoint::StMuFmsPoint(int detectorId, float energy,
                           float x, float y, float z)
    : mDetectorId(detectorId), mEnergy(energy), mX(x), mY(y), mZ(z) { }

StMuFmsPoint::StMuFmsPoint(const StFmsPoint& point)
    : mDetectorId(point.detectorId()), mEnergy(point.energy()),
      mX(point.x()), mY(point.y()), mZ(0.) {
  // Calculate z coordinate from StFmsPoint 4-momentum as it doesn't store
  // z directly. z / x = pz / px, so...
  const TLorentzVector vec4 = point.fourMomentum();
  mZ = point.x() * vec4.Pz() / vec4.Px();
}

StMuFmsPoint::~StMuFmsPoint() { }

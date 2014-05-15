// $Id$
//
// $Log$
/**
 \file      StFmsPoint.cxx
 \brief     Implementation of StFmsPoint, the StEvent FMS photon structure
 \author    Yuxi Pan <yuxipan@physics.ucla.edu>
 \author    Thomas Burton <tpb@bnl.gov>
 \date      2014
 \copyright Brookhaven National Lab
 */
#include "StRoot/StEvent/StFmsPoint.h"

StFmsPoint::StFmsPoint() : StObject() {
  mDetectorId = 0;
  mEnergy = -1.0;
  mX = -99.0;
  mY = -99.0;
  mId = -1;
}

StFmsPoint::~StFmsPoint() {}

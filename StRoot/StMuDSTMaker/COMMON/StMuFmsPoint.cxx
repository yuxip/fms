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

StMuFmsPoint::StMuFmsPoint(int detectorId, float energy, float x, float y)
    : mDetectorId(detectorId), mEnergy(energy), mX(x), mY(y) { }

StMuFmsPoint::~StMuFmsPoint() { }

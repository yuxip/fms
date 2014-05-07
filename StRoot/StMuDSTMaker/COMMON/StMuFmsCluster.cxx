// $Id$
//
// $Log$
/**
 \file      StMuFmsCluster.cxx
 \brief     Implementation of StMuFmsCluster, the MuDST FMS cluster class
 \author    Thomas Burton <tpb@bnl.gov>
 \date      2014
 \copyright Brookhaven National Lab
 */

#include "StMuFmsCluster.h"

StMuFmsCluster::StMuFmsCluster(int detectorId, int category, float energy,
                               float x, float y)
    : mDetectorId(detectorId), mCategory(category), mEnergy(energy),
      mX(x), mY(y) { }

StMuFmsCluster::~StMuFmsCluster() { }

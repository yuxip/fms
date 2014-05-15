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

#include "StRoot/StEvent/StFmsCluster.h"

StMuFmsCluster::StMuFmsCluster(int detectorId, int category, float energy,
                               float x, float y)
    : mDetectorId(detectorId), mCategory(category), mEnergy(energy),
      mX(x), mY(y) { }

StMuFmsCluster::StMuFmsCluster(const StFmsCluster& cluster)
    : mDetectorId(cluster.detectorId()), mCategory(cluster.category()),
      mX(cluster.x()), mY(cluster.y()) { }

StMuFmsCluster::~StMuFmsCluster() { }

void StMuFmsCluster::Clear(Option_t* /* option */) {
  mHits.Clear();
  mPhotons.Clear();
}

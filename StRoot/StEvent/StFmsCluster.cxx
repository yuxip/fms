/*****************************************************************************
 *
 * $Id$
 *
 * Author: Thomas Burton, Yuxi Pan, 2014
 * ***************************************************************************
 *
 * Description: Implementation of StFmsCluster, a group of adjacent FMS towers
 *
 * ***************************************************************************
 *
 * $Log$
 *
 *****************************************************************************/
#include "StFmsCluster.h"

#include "StMessMgr.h"

StFmsCluster::StFmsCluster()
    : StObject(), mCategory(0), mNTowers(0), mEnergy(0.), mX(0.),
      mY(0.), mSigmaMin(0.), mSigmaMax(0.), mChi2Ndf1Photon(-1.),
      mChi2Ndf2Photon(-1.), mId(0) { }

StFmsCluster::~StFmsCluster() {
}

void StFmsCluster::Print(Option_t* /* not used */) const {
  LOG_INFO << "========StFmsCluster:\n\tcatag:\t" << category()
           << "\n\tnumberTower:\t" << nTowers()
           << "\n\tnPhoton:\t" << nPhotons()
           << "\n\tcluster energy:\t" << energy()
           << "\n\tx0:\t" << x()
           << "\n\ty0:\t" << y()
           << "\n\tsiggmaMax:\t" << sigmaMax()
           << "\n\tsigmaMax:\t" << sigmaMin()
           << "\n\tchi2NdfPh1:\t" << chi2Ndf1Photon()
           << "\n\tchi2NdfPh2:\t" << chi2Ndf2Photon()
           << "\n\tid:\t" << id()
           << "\n\tPhoton List:  \n" << endm;
}

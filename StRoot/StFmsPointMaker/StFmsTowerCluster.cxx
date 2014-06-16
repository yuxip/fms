// $Id$
//
// $Log$
/**
 \file      StFmsTowerCluster.cxx
 \brief     Implementation of StFmsTowerCluster, a cluster of FMS towers
 \author    Steven Heppelmann <steveheppelmann@gmail.com>
 \author    Yuxi Pan <yuxipan@physics.ucla.edu>
 \author    Thomas Burton <tpb@bnl.gov>
 \date      2014
 \copyright Brookhaven National Lab
 */
#include "StFmsPointMaker/StFmsTowerCluster.h"

#include <cmath>

#include <boost/foreach.hpp>

#include <TMath.h>
#include <TObjArray.h>
#include <TVector2.h>

#include "StEvent/StFmsCluster.h"
#include "StEvent/StFmsHit.h"

#include "StFmsPointMaker/StFmsTower.h"
#include "StFmsPointMaker/StFmsFittedPhoton.h"

namespace FMSCluster {
StFmsTowerCluster::StFmsTowerCluster(StFmsCluster* cluster)
    : mCluster(cluster), mEnergyCutoff(0.5) {
  Clear();
}

StFmsTowerCluster::~StFmsTowerCluster() {
}

void StFmsTowerCluster::Clear(const char* /* option */) { 
  mSigmaX = mSigmaY = mSigmaXY = mChiSquare = -1.;
  mThetaAxis = -10;
  for (Int_t i(0); i < kMaxPhotonsPerCluster; ++i) {
    mPhotons[i].Clear();
  }  // for
  mTowers.clear();
}

void StFmsTowerCluster::calculateClusterMoments(Float_t Ecoff) {
  mEnergyCutoff=Ecoff;
  Float_t w0, w1, mtmp, mx, my, sigx, sigy, sigXY;
  w0 = w1 = mtmp = mx = my = sigx = sigy = sigXY = 0;
  BOOST_FOREACH(const StFmsTower* tower, mTowers) {
    Float_t xxx, yyy;
    xxx = tower->column() - 0.5;
    yyy = tower->row() - 0.5;
    mtmp = log(tower->hit()->energy() + 1. - Ecoff) > 0 ?
           log(tower->hit()->energy() + 1. - Ecoff) : 0;
    w1 += mtmp;
    w0 += tower->hit()->energy();
    mx += mtmp * xxx;
    my += mtmp * yyy;
    sigx += mtmp * xxx * xxx;
    sigy += mtmp * yyy * yyy;
    sigXY += mtmp * xxx * yyy;
  }  // BOOST_FOREACH
  mCluster->setEnergy(w0);
  if (w1 > 0) {
    mCluster->setX(mx / w1);
    mCluster->setY(my / w1);
    mSigmaX = sqrt(fabs(sigx / w1 - std::pow(mCluster->x(), 2.)));
    mSigmaY = sqrt(fabs(sigy / w1 - std::pow(mCluster->y(), 2.)));
    mSigmaXY = sigXY / w1 - mCluster->x() * mCluster->y();
  } else {
    mCluster->setX(0.);
    mCluster->setY(0.);
    mSigmaX = 0;
    mSigmaY = 0;
    mSigmaXY = 0;
  }  // if
}

void StFmsTowerCluster::findClusterAxis() {
  Double_t dSigma2, aA, bB;
  dSigma2 = mSigmaX * mSigmaX - mSigmaY * mSigmaY;
  aA = sqrt(dSigma2 * dSigma2 + 4.0 * mSigmaXY * mSigmaXY) + dSigma2;
  bB = 2 * mSigmaXY;
	if (mSigmaXY < 1e-10) {
		if (aA < 1e-10) {
			bB = sqrt(dSigma2 * dSigma2 + 4.0 * mSigmaXY * mSigmaXY) - dSigma2;
			aA = 2 * mSigmaXY;
		}  // if
	}  // if
	mThetaAxis = atan2(bB, aA);
	Double_t myPi = TMath::Pi(); 
	while (mThetaAxis > (myPi / 2.0)) {
		mThetaAxis -= myPi;
	}  // while
	while (mThetaAxis < -(myPi / 2.0)) {
		mThetaAxis += myPi;
	}  // while
	mCluster->setSigmaMin(getSigma(mThetaAxis));
	mCluster->setSigmaMax(getSigma(mThetaAxis - TMath::Pi() / 2.0));
}

Double_t StFmsTowerCluster::getSigma(Double_t theta) const {
	Double_t sigma = 0;
	// 2-d vector vaxis define the axis
	TVector2 vaxis(cos(theta), sin(theta));
	// loop over all towers pointer in cluster
	StFmsTower* oneTower;
	float wnew =0;
	BOOST_FOREACH(const StFmsTower* tower, mTowers) {
		// the 2-d vector from the "center" of cluster to tower
		// "center" are at 0.5, 1.5, etc! Need shift of 0.5
		TVector2 v1(tower->column() - 0.5 - mCluster->x(),
		            tower->row() - 0.5 - mCluster->y());
		// perpendicular distance to the axis = length of the component of vector
		// "v1" that is norm to "vaxis"
		Double_t dis = (v1.Norm(vaxis)).Mod();
		// contribution to sigma
		//sigma += oneTower->energy * dis * dis;
		float wtmp = log(tower->hit()->energy() + 1. - mEnergyCutoff) > 0 ?
		             log(tower->hit()->energy() + 1. - mEnergyCutoff) : 0;
		wnew += wtmp;
		sigma += wtmp * dis * dis;
	}  // BOOST_FOREACH
	return wnew > 0 ? sqrt(sigma / wnew) : 0;
}
}  // namespace FMSCluster

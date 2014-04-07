#include "StPSUTools/HitCluster.h"

#include <cmath>

#include <TMath.h>
#include <TObjArray.h>
#include <TVector2.h>

#include "StEvent/StFmsCluster.h"
#include "StEvent/StFmsHit.h"

#include "StPSUTools/StFmsTower.h"
#include "StPSUTools/PhotonHitFPD.h"

namespace PSUGlobals {
HitCluster::HitCluster(StFmsCluster* cluster)
    : mTowers(NULL),  mCluster(cluster), Ecutoff(0.5) {
  Clear();
}

HitCluster::~HitCluster() {
  if (mTowers) {
    delete mTowers;
  }  // if
}

void HitCluster::Clear(const char* /* unused option inherited from TObject */) { 
  mSigmaX = mSigmaY = mSigmaXY = mChiSquare = -1.;
  mThetaAxis = -10;
  for (Int_t i(0); i < mMaxPhotonsPerCluster; ++i) {
    mPhotons[i].Clear();
  }  // for
  if (mTowers) {
    delete mTowers;
  }  // if
  mTowers = new TObjArray(50);
}

void HitCluster::FindClusterAxis() {
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
	mCluster->SetSigmaMin(GetSigma(mThetaAxis));
	mCluster->SetSigmaMax(GetSigma(mThetaAxis - TMath::Pi() / 2.0));
}

// Calculate sigma w.r.t the axis going through the "center" and of an angle
// "theta" in x-y plane
Double_t HitCluster::GetSigma(Double_t theta) {
	Double_t sigma = 0;
	// 2-d vector vaxis define the axis
	TVector2 vaxis(cos(theta), sin(theta));
	// loop over all towers pointer in cluster
	TowerFPD * oneTower;
	float wnew =0;
	for (Int_t it=0; it < mTowers->GetEntriesFast(); it++) {
		oneTower = (TowerFPD *)mTowers->At(it);
		// the 2-d vector from the "center" of cluster to tower
		// "center" are at 0.5, 1.5, etc! Need shift of 0.5
		TVector2 v1(oneTower->column() - 0.5 - mCluster->GetX0(),
		            oneTower->row() - 0.5 - mCluster->GetY0());
		// perpendicular distance to the axis = length of the component of vector
		// "v1" that is norm to "vaxis"
		Double_t dis = (v1.Norm(vaxis)).Mod();
		// contribution to sigma
		//sigma += oneTower->energy * dis * dis;
		float wtmp = log(oneTower->hit()->energy() + 1. - Ecutoff) > 0 ?
		             log(oneTower->hit()->energy() + 1. - Ecutoff) : 0;
		wnew += wtmp;
		sigma += wtmp * dis * dis;
	}  // for
	return wnew > 0 ? sqrt(sigma / wnew) : 0;
}

void HitCluster::CalClusterMoment(Float_t Ecoff) {
  Ecutoff=Ecoff;
  Float_t w0, w1, mtmp, mx, my, sigx, sigy, sigXY;
  w0 = w1 = mtmp = mx = my = sigx = sigy = sigXY = 0;
  TowerFPD * oneTow;
  TIter next(mTowers);
  while (oneTow=(TowerFPD*) next()) {
    Float_t xxx, yyy;
    xxx = oneTow->column() - 0.5;
    yyy = oneTow->row() - 0.5;
    mtmp = log(oneTow->hit()->energy() + 1. - Ecoff) > 0 ?
           log(oneTow->hit()->energy() + 1. - Ecoff) : 0;
    w1 += mtmp;
    w0 += oneTow->hit()->energy();
    mx += mtmp * xxx;
    my += mtmp * yyy;
    sigx += mtmp * xxx * xxx;
    sigy += mtmp * yyy * yyy;
    sigXY += mtmp * xxx * yyy;
  }  // while
  mCluster->SetClusterEnergy(w0);
  if (w1 > 0) {
    mCluster->SetX0(mx / w1);
    mCluster->SetY0(my / w1);
    mSigmaX = sqrt(fabs(sigx / w1 - std::pow(mCluster->GetX0(), 2.)));
    mSigmaY = sqrt(fabs(sigy / w1 - std::pow(mCluster->GetY0(), 2.)));
    mSigmaXY = sigXY / w1 - mCluster->GetX0() * mCluster->GetY0();
  } else {
    mCluster->SetX0(0.);
    mCluster->SetY0(0.);
    mSigmaX = 0;
    mSigmaY = 0;
    mSigmaXY = 0;
  }  // if
}
}  // namespace PSUGlobals

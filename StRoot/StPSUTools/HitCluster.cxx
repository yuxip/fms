#include "StPSUTools/HitCluster.h"

#include <TMath.h>
#include <TObjArray.h>
#include <TVector2.h>

#include "StEvent/StFmsHit.h"

#include "StPSUTools/TowerFPD.h"
#include "StPSUTools/PhotonHitFPD.h"
#include "StPSUTools/Yiqun.h"

namespace PSUGlobals {
HitCluster::HitCluster() : IsEUpdated(false), tow(0), Ecutoff(0.5) {
  Clear();
}

HitCluster::~HitCluster() {
  if (tow) {
    delete tow;
  }  // if
}

void HitCluster::Clear(const char* /* unused option inherited from TObject */) { 
  sigmaX = sigmaY = sigmaXY = chiSquare = -1;
  thetaAxis = -10;
  for (Int_t i(0); i < mMaxPhotonsPerCluster; ++i) {
    photon[i].Clear();
  }  // for
  if (tow) {
    delete tow;
  }  // if
  tow = new TObjArray(50);
}

void HitCluster::FindClusterAxis() {
  Double_t dSigma2, aA, bB;
  dSigma2 = sigmaX * sigmaX - sigmaY * sigmaY;
  aA = sqrt(dSigma2 * dSigma2 + 4.0 * sigmaXY * sigmaXY) + dSigma2;
  bB = 2 * sigmaXY;
	if (sigmaXY < 1e-10) {
		if (aA < 1e-10) {
			bB = sqrt(dSigma2 * dSigma2 + 4.0 * sigmaXY * sigmaXY) - dSigma2;
			aA = 2 * sigmaXY;
		}  // if
	}  // if
	thetaAxis = atan2(bB, aA);
	Double_t myPi = TMath::Pi(); 
	while (thetaAxis > (myPi / 2.0)) {
		thetaAxis -= myPi;
	}  // while
	while (thetaAxis < -(myPi / 2.0)) {
		thetaAxis += myPi;
	}  // while
	mSigmaMin = GetSigma(thetaAxis);
	mSigmaMax = GetSigma(thetaAxis - TMath::Pi() / 2.0);
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
	for (Int_t it=0; it<tow->GetEntriesFast(); it++) {
		oneTower = (TowerFPD *) tow->At(it);
		// the 2-d vector from the "center" of cluster to tower
		// "center" are at 0.5, 1.5, etc! Need shift of 0.5
		TVector2 v1(oneTower->column() - 0.5 - mX0, oneTower->row() - 0.5 - mY0);
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
  TIter next(tow);
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
  mEnergy = w0;
  if (w1 > 0) {
    mX0 = mx / w1;
    mY0 = my / w1;
    sigmaX = sqrt(fabs(sigx / w1 - mX0 * mX0));
    sigmaY = sqrt(fabs(sigy / w1 - mY0 * mY0));
    sigmaXY = sigXY / w1 - mX0 * mY0;
  } else {
    mX0 = 0;
    mY0 = 0;
    sigmaX = 0;
    sigmaY = 0;
    sigmaXY = 0;
  }  // if
}

void HitCluster::copyTo(StFmsCluster* cluster) const {
  cluster->SetCatag(mCatag);
  cluster->SetNumbTower(mNumbTower);
  cluster->SetNphoton(mNphoton);
  cluster->SetClusterEnergy(mEnergy);
  cluster->SetX0(mX0);
  cluster->SetY0(mY0);
  cluster->SetSigmaMax(mSigmaMax);
  cluster->SetSigmaMin(mSigmaMin);
  //	cluster->SetChi2NdfPh1(Chi2NdfPh1);
  //	--no such funtion in SH's package --Yuxi
  //	cluster->SetChi2NdfPh2(Chi2NdfPh2);
}
}  // namespace PSUGlobals

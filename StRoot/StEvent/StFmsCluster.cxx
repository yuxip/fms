#include "StFmsCluster.h"

#include "StMessMgr.h"

#define MAX_PHOTON_PER_CLUSTER 2

StFmsCluster::StFmsCluster()
    : StObject(), mCatag(0), mNumbTower(0), mNphoton(0), mEnergy(0.), mX0(0.),
      mY0(0.), mSigmaMax(0.), mSigmaMin(0.), mChi2NdfPh1(-1.), mChi2NdfPh2(-1.),
      mCluId(0) {}

StFmsCluster::~StFmsCluster() {
}

void StFmsCluster::Print(Option_t* /* not used */) const {
	LOG_INFO << "========StFmsCluster:\n\tcatag:\t" << GetCatag()
           << "\n\tnumberTower:\t" << GetNTower()
           << "\n\tnPhoton:\t" << GetNphoton()
           << "\n\tcluster energy:\t" << GetEnergy()
           << "\n\tx0:\t" << GetX0()
           << "\n\ty0:\t" << GetY0()
           << "\n\tsiggmaMax:\t" << GetSigmaMax()
           << "\n\tsigmaMin:\t" << GetSigmaMin()
           << "\n\tchi2NdfPh1:\t" << GetChi2NdfPh1()
           << "\n\tchi2NdfPh2:\t" << GetChi2NdfPh2()
           << "\n\tid:\t" << GetClusterId()
           << "\n\tPhoton List:  \n" << endm;
}

Bool_t StFmsCluster::SetNphoton(Int_t nPhoton) {
	if (nPhoton <= 0 || nPhoton > MAX_PHOTON_PER_CLUSTER) {
		LOG_ERROR << "StFmsCluster::SetNphoton() illegal nPhoton" << endm;
		return false;
	}  // if
	mNphoton = nPhoton;
	return true;
}

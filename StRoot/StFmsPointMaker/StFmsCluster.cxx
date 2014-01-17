#include "StFmsClHitCollection.h"
#include "StFmsPointCollection.h"
#include "StFmsCluster.h"
#include "StMessMgr.h"

#define MAX_PHOTON_PER_CLUSTER 2

ClassImp(StFmsCluster)

StFmsCluster::StFmsCluster():StObject() {

        mCatag           = 0;
        mNumbTower       = 0;
        mNphoton         = 0;
        mEnergy          = 0;
        mX0              = 0;
        mY0              = 0;
        mSigmaMax        = 0;
        mSigmaMin        = 0;
        mCluId           = 0;
        mChi2NdfPh1      = -1;
        mChi2NdfPh2      = -1;
	
	mPhotons = new StFmsPointCollection();
	mClhits = new StFmsClHitCollection();
	
}

StFmsCluster::~StFmsCluster() {
	
	LOG_DEBUG << " StFmsCluster destructor " << endm;
	if(mPhotons){
		mPhotons->Clear();
		delete mPhotons;
		mPhotons = 0;
	}
	if(mClhits){
		mClhits->Clear();
		delete mClhits;
		mClhits = 0;
	}
}

void StFmsCluster::Clear( const char* opt ) {

	if(mPhotons)mPhotons->Clear();
        if(mClhits)mClhits->Clear();
}

ostream& operator<<(ostream& os, const StFmsCluster& v){
	
	os << "========StFmsCluster:\n\tcatag:\t"<<v.GetCatag()
                   << "\n\tnumberTower:\t"<<v.GetNTower()
                   << "\n\tnPhoton:\t"<<v.GetNphoton()
                   << "\n\tcluster energy:\t"<<v.GetEnergy()
                   << "\n\tx0:\t"<<v.GetX0()
                   << "\n\ty0:\t"<<v.GetY0()
                   << "\n\tsiggmaMax:\t"<<v.GetSigmaMax()
                   << "\n\tsigmaMin:\t"<<v.GetSigmaMin()
                   << "\n\tchi2NdfPh1:\t"<<v.GetChi2NdfPh1()
                   << "\n\tchi2NdfPh2:\t"<<v.GetChi2NdfPh2()
                   << "\n\tid:\t"<<v.GetClusterId()
                   << "\n\tPhoton List:  \n";
	return os;
}

void StFmsCluster::Print(Option_t *option) const {

	LOG_INFO << *this << endm;
	mPhotons->Print();
	mClhits->Print();
}

Bool_t StFmsCluster::SetNphoton( Int_t nPhoton ){

	if(nPhoton<=0||nPhoton>MAX_PHOTON_PER_CLUSTER){
		LOG_ERROR << "StFmsCluster::SetNphoton() illegal nPhoton " << endm;
		return false;
	}
	
	mNphoton = nPhoton;
	return true;
}

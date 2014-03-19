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
}

StFmsCluster::~StFmsCluster() {
	
	LOG_DEBUG << " StFmsCluster destructor " << endm;
	for (unsigned i(0); i < mPhotons.size(); ++i) {
	  if (mPhotons.at(i)) {
	    delete mPhotons.at(i);
	    mPhotons.at(i) = NULL;
	  }  // if
	}  // for
}

void StFmsCluster::Clear( const char* opt ) {

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
}

Bool_t StFmsCluster::SetNphoton( Int_t nPhoton ){

	if(nPhoton<=0||nPhoton>MAX_PHOTON_PER_CLUSTER){
		LOG_ERROR << "StFmsCluster::SetNphoton() illegal nPhoton " << endm;
		return false;
	}
	
	mNphoton = nPhoton;
	return true;
}

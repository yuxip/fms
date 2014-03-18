#include "StMessMgr.h"
#include "StFmsClusterCollection.h"

ClassImp(StFmsClusterCollection)

StFmsClusterCollection::StFmsClusterCollection():TDataSet("StFmsClusterCollection")  { /* no op */ }

StFmsClusterCollection::~StFmsClusterCollection() {
	
	LOG_DEBUG << "StFmsClusterCollection destructor " << endm;
	Int_t n = mClusters.size();
	for(Int_t i=0; i<n; i++) { 
		delete mClusters[i]; 
		mClusters[i] = 0;
	}

	mClusters.clear();
}

void StFmsClusterCollection::Clear( const char* opt ) {
	
	Int_t n = mClusters.size();
        for (Int_t i=0; i<n; i++) {
                delete mClusters[i];
                mClusters[i] = 0;
        }
	LOG_INFO<<"  StFmsClusterCollection::Clear() "<<endm;
        mClusters.clear();
}

void StFmsClusterCollection::Print( Option_t *option ) const {

        LOG_INFO << "+++++ StFmsClusterCollection--NumberOfClusters: " << this->NumberOfClusters() << endm;

        for(StPtrVecFmsClusterConstIterator iclu = mClusters.begin(); iclu != mClusters.end(); iclu++){
                (*iclu)->Print();
        }

        LOG_INFO << "+++++ END StFmsClusterCollection " <<endm;
}

unsigned int StFmsClusterCollection::NumberOfClusters() const { return mClusters.size(); }

void StFmsClusterCollection::AddCluster(StFmsCluster* cluster) { mClusters.push_back(cluster); }

StPtrVecFmsCluster& StFmsClusterCollection::clusters() { return mClusters; }

const StPtrVecFmsCluster& StFmsClusterCollection::clusters() const { return mClusters; }


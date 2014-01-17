#include "StMessMgr.h"
#include "StFmsClHitCollection.h"

ClassImp(StFmsClHitCollection)

StFmsClHitCollection::StFmsClHitCollection():TDataSet("StFmsClHitCollection") { /* no op */ }

StFmsClHitCollection::~StFmsClHitCollection() {
	
	LOG_DEBUG << "StFmsClHitCollection destructor " << endm;
	int n = mHits.size();
	for (int i=0; i<n; i++) { delete mHits[i]; }
    	mHits.clear();
}

void StFmsClHitCollection::Clear( const char* opt ) {
	
	int n = mHits.size();
	for (int i=0; i<n; i++) {
		delete mHits[i];
		mHits[i] = 0;
	}
	
	mHits.clear();
}

void StFmsClHitCollection::Print( Option_t *option ) const {
	
	LOG_INFO << "===== StFmsClHitCollection--NumberOfHits: " << this->NumberOfHits() << endm; 

	for(StFmsClHitConstIterator ihit = mHits.begin(); ihit != mHits.end(); ihit++){
		(*ihit)->Print();
	}

	LOG_INFO << "===== END StFmsClHitCollection " <<endm;
}	

unsigned int StFmsClHitCollection::NumberOfHits() const { return mHits.size(); }

void StFmsClHitCollection::AddHit(StFmsClHit* hit) { mHits.push_back(hit); }

StPtrVecFmsClHit& StFmsClHitCollection::hits() { return mHits; }

const StPtrVecFmsClHit& StFmsClHitCollection::hits() const { return mHits; }


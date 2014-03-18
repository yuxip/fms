#include "StMessMgr.h"
#include "StFmsPointCollection.h"

ClassImp(StFmsPointCollection)

StFmsPointCollection::StFmsPointCollection():TDataSet("StFmsPointCollection") { /* no op */ }

StFmsPointCollection::~StFmsPointCollection() {
	
	LOG_DEBUG << "StFmsPointCollection destructor " << endm;
	Int_t n = mPoints.size();
	for(Int_t i=0; i<n; i++) { delete mPoints[i]; }
	mPoints.clear();
}

void StFmsPointCollection::Clear( const char* opt ) {
	
	Int_t n = mPoints.size();
        for (Int_t i=0; i<n; i++) {
                delete mPoints[i];
                mPoints[i] = 0;
        }

        mPoints.clear();
}

void StFmsPointCollection::Print( Option_t *option ) const {

        LOG_INFO << "----- StFmsPointCollection--NumberOfPoints: " << this->NumberOfPoints() << endm;

        for(StPtrVecFmsPointConstIterator ipts = mPoints.begin(); ipts != mPoints.end(); ipts++){
                (*ipts)->Print();
        }

        LOG_INFO << "----- END StFmsPointCollection " <<endm;
}

unsigned int StFmsPointCollection::NumberOfPoints() const { return mPoints.size(); }

void StFmsPointCollection::AddPoint(StFmsPoint* point) { mPoints.push_back(point); }

StPtrVecFmsPoint& StFmsPointCollection::points() { return mPoints; }

const StPtrVecFmsPoint& StFmsPointCollection::points() const { return mPoints; }


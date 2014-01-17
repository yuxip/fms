#include "StMessMgr.h"
#include "StFmsPoint.h"

ClassImp(StFmsPoint)

StFmsPoint::StFmsPoint():StObject(){
	
	mEnergy = -1.0;
	mXpos	= -99.0;
	mYpos	= -99.0;
	mPhid	= -1;
	
}

StFmsPoint::StFmsPoint ( StFmsPoint& other ) {
	
	LOG_DEBUG << "StFmsPoint copy constructor called " <<endm;
	mEnergy = other.GetEnergy();
	mXpos	= other.GetXpos();
	mYpos	= other.GetYpos();
	mPhid	= other.GetPhid();
	mFourMomentum = other.fourMomentum();
	mPointXYZLab = other.pointXYZLab();
}

StFmsPoint::~StFmsPoint() {}

void StFmsPoint::Print( const Option_t* opt ) const { LOG_INFO << *this << endm; }

ostream& operator<<( ostream& os, const StFmsPoint& v ) {
	
	return os <<"StFmsPoint:\n\tenergy: "<<v.GetEnergy()
                  <<"\n\txpos:   "<<v.GetXpos()
                  <<"\n\typos:   "<<v.GetYpos()
                  <<"\n\tphId:   "<<v.GetPhid()
                  <<"\n";
}

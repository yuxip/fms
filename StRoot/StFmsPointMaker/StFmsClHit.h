#ifndef StFmsClHit_HH
#define StFmsClHit_HH

#include"StObject.h"
#include "Stiostream.h"

#include <vector>

#ifndef ST_NO_NAMESPACES
using std::vector;
#endif

//class that represents a hit on fms cell, to be associated with recontructed clusters and
//Yuxi Pan 03/31/2013


class StFmsClHit : public StObject {

public:
		StFmsClHit();
		~StFmsClHit();
		
		StFmsClHit(Int_t nstb, Float_t row0, Float_t col0, UInt_t adc, Float_t energy, UChar_t status);

	Int_t	GetNstb()	const	{ return mNstb;   };
	Float_t	GetRow()	const	{ return mRow;	  };
	Float_t	GetCol()	const	{ return mCol;	  };
	Float_t	GetEnergy()	const	{ return mEnergy; };
	UInt_t	GetAdc()	const	{ return mAdc;    };
	UChar_t	GetStatus()	const	{ return mStatus; };
	

	Int_t	SetEnergy( Float_t energy )	{ mEnergy	= energy; return 1; };
	Int_t	SetAdc( UInt_t adc)		{ mAdc		= adc; return 1; };	
	Int_t	SetStatus( UChar_t status)	{ mStatus	= status; return 1; };
	Bool_t	Legal( Int_t iew, Int_t nstb, Float_t row0, Float_t col0 );

	void Print( Option_t *option = "" )	const;

private:
	Int_t	mNstb;
	Float_t	mRow;
	Float_t	mCol;
	UInt_t	mAdc;
	Float_t	mEnergy;
	UChar_t	mStatus;
	

	ClassDef(StFmsClHit,1)
};

ostream& operator<<(ostream&, const StFmsClHit&);

typedef vector<StFmsClHit*> StPtrVecFmsClHit;
typedef StPtrVecFmsClHit::iterator StFmsClHitIterator;
typedef StPtrVecFmsClHit::const_iterator StFmsClHitConstIterator;

#endif					

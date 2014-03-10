//container of StFmsClHit
//Yuxi Pan --03/31/2013
#ifndef StFmsClHitCollection_HH
#define StFmsClHitCollection_HH

#include "StObject.h"
#include "TDataSet.h"

#include "StFmsClHit.h"

class StFmsClHitCollection : public TDataSet {

public:
	StFmsClHitCollection();
	~StFmsClHitCollection();
	void Clear(const char* opt="");
	void Print(Option_t *option="") const;

	void AddHit(StFmsClHit*);
	unsigned int NumberOfHits() const;
	
	StPtrVecFmsClHit&		hits();
	const StPtrVecFmsClHit&		hits() const;
	
private:
	StPtrVecFmsClHit mHits;	//StFmsClHitCollection owns the hits
	
	ClassDef(StFmsClHitCollection,1)
};
#endif

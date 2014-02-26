//Call Yiqun class to find FMS clusters and
//fit clusters with photon hypothesis (shower fit)
//adapted from PSU code by Yuxi Pan --03/31/2013

#ifndef StFmsPointMaker_HH
#define StFmsPointMaker_HH

#include <vector>

#include <TMatrix.h>

#ifndef StMaker_H
#define StMaker_H
#include "StMaker.h"
#endif
#include "StPSUTools/Yiqun.h"
#include "StPSUTools/TowerFPD.h"
#include "StPSUTools/Geom.h"
using namespace std;
using namespace PSUGlobals;

class StFmsDbMaker;
class StFmsHitMaker;
class StFmsClusterCollection;
class StFmsPointCollection;
class StFmsPointMaker : public StMaker {

public:
	StFmsPointMaker(const char* name);
	~StFmsPointMaker();
	
	void  Clear(const char* opt="");
	Int_t Init();
	Int_t InitRun(Int_t runNumber);    //called by StMaker when switch to a new run#
	Int_t Make();
	Int_t Finish();


private:

	Int_t FindPoint();			//  --interface to the actual photon reconstruction
	Bool_t initialiseEnergyMatrices();
  Bool_t Legal(Int_t iew, Int_t nstb, Int_t row0, Int_t col0);

  StFmsDbMaker* mFmsDbMaker;
	std::vector<TMatrix> mEnergyMatrices;
	StFmsClusterCollection* mFmsClColl;	//! --clusters (and points within cluster) to be added to TDataSet
	//StFmsPointCollection*   mFmsPtsColl;	//! --all the points (photons) extracted from clusters
	Geom* fmsgeom;
	
	ClassDef(StFmsPointMaker,0)
};
#endif

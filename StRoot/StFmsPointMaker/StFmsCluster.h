#ifndef StStFmsCluster_HH
#define StStFmsCluster_HH

#include "StObject.h"
#include "Stiostream.h"
#include "TLorentzVector.h"

#include <vector>

#ifndef ST_NO_NAMESPACES
using std::vector;
#endif

#define MAX_PHOTON_PER_CLUSTER 2

//class that represents a tower cluster
//Yuxi Pan 03/31/2013

class StFmsClHitCollection;
class StFmsPointCollection;

class StFmsCluster : public StObject {
	
public: 
	StFmsCluster();
	~StFmsCluster();
	
	void Print( Option_t *option = "" ) const;
	void Clear( const char* opt = "" );	

	Int_t GetNstb()		const { return mNstb; }
	Int_t GetCatag()	const { return mCatag; }
	Int_t GetNTower()	const { return mNumbTower; }
	Int_t GetNphoton()	const { return mNphoton; }
	Float_t GetEnergy()	const { return mEnergy; }
	Float_t GetX0()		const { return mX0; }
	Float_t GetY0()		const { return mY0; }
	Float_t GetSigmaMax()	const { return mSigmaMax; }
	Float_t GetSigmaMin()	const { return mSigmaMin; }
	Float_t GetChi2NdfPh1()	const { return mChi2NdfPh1; }
	Float_t GetChi2NdfPh2()	const { return mChi2NdfPh2; }
	Int_t GetClusterId()	const { return mCluId; }
	TLorentzVector GetFourMomentum() const { return mFourMomentum; }

	StFmsPointCollection* GetPointCollection() { return mPhotons; }
	StFmsClHitCollection* GetClHitCollection() { return mClhits; }	
	
	
	void SetNstb( Int_t nstb )		{ mNstb = nstb; }	
	void SetCatag( Int_t catag )		{ mCatag = catag; }
	void SetNumbTower( Int_t numbTower )	{ mNumbTower = numbTower; }
	Bool_t SetNphoton( Int_t nPhoton );
	void SetClusterEnergy( Float_t energy ) { mEnergy = energy; }
	void SetX0 ( Float_t x0 )		{ mX0 = x0; }
	void SetY0 ( Float_t y0 )		{ mY0 = y0; }
	void SetSigmaMax ( Float_t sigmaMax )	{ mSigmaMax = sigmaMax; }
	void SetSigmaMin ( Float_t sigmaMin )	{ mSigmaMin = sigmaMin; }
	void SetChi2NdfPh1(Float_t chi2ndfph1)	{ mChi2NdfPh1 = chi2ndfph1; }
	void SetChi2NdfPh2(Float_t chi2ndfph2)	{ mChi2NdfPh2 = chi2ndfph2; }
	void SetClusterId ( Float_t cluid )	{ mCluId = cluid; }
	void SetFourMomentum ( TLorentzVector p4 ) { mFourMomentum = p4; }

private:
	
	Int_t		mNstb;			// Nstb starts from 1
	Int_t		mCatag  ;               // catagory of cluster (1: 1-photon,  2: 2-photon,  0: could be either 1- or 2-photon)
	Int_t		mNumbTower;             // number of non_zero towers in the cluster
	Int_t		mNphoton;               // number of photons contained in this cluster: ( nPhoton <= MAX_PHOTON_PER_CLUSTER ! )
	Float_t		mEnergy ;               // total energy  contained in this cluster (0th moment)
	Float_t		mX0     ;               // mean x ("center of gravity") in local grid coordinate (1st moment)
	Float_t		mY0     ;               // mean y ("center of gravity") in local grid coordinate (1st moment)     
	Float_t		mSigmaMax;              // maximum 2nd moment (along major axis)
	Float_t		mSigmaMin;              // minimum 2nd moment
	Float_t		mChi2NdfPh1;            // chi2 / ndf for 1-photon fit
	Float_t		mChi2NdfPh2;            // chi2 / ndf for 2-photon fit
	Int_t		mCluId     ;            // eventwise cluster id, also include nstb info.

	TLorentzVector  mFourMomentum;		// cluster four momentum;

	StFmsPointCollection*	mPhotons; 	//->
						//fitted points (photons) in the cluster
	StFmsClHitCollection*	mClhits;	//->

						//an array of tower hits of the current cluster
	ClassDef(StFmsCluster,1)
};

ostream& operator<<(ostream&, const StFmsCluster&);

typedef vector<StFmsCluster*> StPtrVecFmsCluster;
typedef StPtrVecFmsCluster::iterator StFmsClusterIterator;
typedef StPtrVecFmsCluster::const_iterator StFmsClusterConstIterator;

#endif

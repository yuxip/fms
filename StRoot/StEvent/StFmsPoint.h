#ifndef StStFmsPoint_HH
#define StStFmsPoint_HH

#include "StObject.h"
#include "Stiostream.h"

#include "TLorentzVector.h"

#include <vector>

#ifndef ST_NO_NAMESPACES
using std::vector;
#endif


//class that represents a fitted photon
//Yuxi Pan 03/31/2013

class StFmsPoint : public StObject {
	
public:
	StFmsPoint();
	StFmsPoint( StFmsPoint& other );
	~StFmsPoint();
	
	void Print(const Option_t* opt = "") const;

	Float_t energy()	const { return mEnergy; }
	Float_t GetXpos()	const { return mXpos;	}	//in cm
	Float_t GetYpos()	const { return mYpos;	}	//in cm
	Int_t 	GetPhid()	const { return mPhid;	}
	Int_t   GetParentCluId()  const { return mCluid;  }
	Int_t	GetParentNclPh()  const { return mNclph;	}
	TLorentzVector	fourMomentum()	const { return mFourMomentum; }
	TVector3	pointXYZLab() const { return mPointXYZLab; }
	
	void SetEnergy	( Float_t energy )	{ mEnergy = energy; }
	void SetXpos  	( Float_t xpos   )	{ mXpos = xpos; }
	void SetYpos  	( Float_t ypos   )	{ mYpos = ypos; }
	void SetPhotonId( Int_t	  phid	 )	{ mPhid = phid; }
	void SetParentCluId( Int_t cluid  )	{ mCluid = cluid; }
	void SetParentNclPh( Int_t nclph  )	{ mNclph = nclph; }

	void setFourMomentum ( TLorentzVector p4 ) { mFourMomentum = p4; }
	void SetPointXYZLab ( TVector3 phpos3 ) { mPointXYZLab = phpos3; }


private:

	Float_t mEnergy;          // fitted energy
        Float_t mXpos  ;          // fitted (relative) x-position
        Float_t mYpos  ;          // fitted (relative) y-position
        Int_t   mPhid;            // photon id within event, also include det12 info
	Int_t   mCluid;		  // id of the parent cluster
	Int_t 	mNclph;		  //# of photons produced by the parent cluster (1 or 2)	

	TLorentzVector mFourMomentum;
	TVector3 mPointXYZLab;	  //photon coordinate in lab fraom
	
	ClassDef(StFmsPoint,1)
};

ostream& operator<<(ostream&, const StFmsPoint&);

#endif

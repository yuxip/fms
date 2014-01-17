////////////////////////////////////////////////////////////////////
//
//  definition of a FPD tower (position and energy deposit)
//
////////////////////////////////////////////////////////////////////

#ifndef FPD_TOWER_H
#define FPD_TOWER_H

#include "CalibStr.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>

#include "TMath.h"
#include "TVector2.h"

#include "TObjArray.h"
#include "TMinuit.h"
#include "TVector3.h"
namespace PSUGlobals {//$NMSPC
class TowerFPD : public TObject {

 public:

	Float_t energy;
	Int_t   col;         // start from 0: count the columns, moving horizontally (STAR x-coord)
	Int_t   row;         // start from 0: count the rows,    moving vertically   (STAR y-coord)
	Int_t   cluster;
	Int_t   cluster2;
	Int_t   adc_over_ped;

	TowerFPD();

	TowerFPD(Float_t ene, Int_t towX, Int_t towY, Int_t clu);

	~TowerFPD() 
	  {
	    if(Lnk_LRUD)delete Lnk_LRUD;
	  };
 
	Bool_t IsEqual(const TObject *obj) const {
		return ( ( col == ((TowerFPD *) obj)->col ) && ( row == ((TowerFPD *) obj)->row ) ) ;
	}

	Bool_t IsSortable() const {return kTRUE;}

	Int_t Compare(const TObject *obj) const ;

	Bool_t IsNeighbor(TowerFPD *a);

  // New Structure Aug 2009
  
	Bool_t LinkDefined; //true after SetContext called
	Bool_t CalibSet;//true after SetContext called
	Bool_t Live; //true if cell has calibration constants
	Bool_t LocalHighTower();
	Bool_t IsolatedDeadCell();
	void ClearClusterInfo(TObjArray* to);
	Int_t CreateContiguous();
	Int_t AddToContiguous(Int_t set,Int_t hitcounter=0);

	//Get the NCl'th list (Must be Deleted);
	TObjArray* CreateClusterList(Int_t NCl);
	Int_t ClusterCategory(TObjArray*);
	Bool_t TowerFilter();
	Bool_t NeighborsAlive(Float_t dcell );
  
	Bool_t SetContext(TObjArray* towers,CalibStr* gain,CalibStr* gaincorr,Int_t IEW,Int_t NSTB);

	void Print(void);

	TowerFPD& operator=(const TowerFPD& rhs) {

		// TowerFPD assignment operator
		//
		if (this != &rhs) {
			energy  = rhs.energy  ;
			col     = rhs.col     ;
			row     = rhs.row     ;
			cluster = rhs.cluster ;
			adc_over_ped = rhs.adc_over_ped;
			LinkDefined = rhs.LinkDefined;
			CalibSet = rhs.CalibSet;
			Live = rhs.Live;
			if(rhs.TowerSet)TowerSet=new TObjArray(*(rhs.TowerSet));
			Gain=rhs.Gain;
			GainCorr=rhs.GainCorr;
		}
		return *this;
	};
	Int_t ClaimedID;
	TObjArray* Lnk_LRUD;
	CalibStr* Gain;
	CalibStr* GainCorr;
	Int_t IEW;
	Int_t NSTB;
	TObjArray* TowerSet;
	
 private:
	ClassDef (TowerFPD, 7);
	
};
}

#endif

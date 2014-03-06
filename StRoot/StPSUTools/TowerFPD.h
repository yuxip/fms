////////////////////////////////////////////////////////////////////
//
//  definition of a FPD tower (position and energy deposit)
//
////////////////////////////////////////////////////////////////////

#ifndef FPD_TOWER_H
#define FPD_TOWER_H

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

class StFmsHit;

namespace PSUGlobals {//$NMSPC
class TowerFPD : public TObject {

 public:

  const StFmsHit* hit;
	Int_t   col;         // start from 0: count the columns, moving horizontally (STAR x-coord)
	Int_t   row;         // start from 0: count the rows,    moving vertically   (STAR y-coord)
	Int_t   cluster;
	Int_t   adc_over_ped;

	TowerFPD();

	TowerFPD(const StFmsHit* fmsHit, Int_t towX, Int_t towY, Int_t clu);

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

	Bool_t SetContext(TObjArray* towers,Int_t IEW,Int_t NSTB);

	void Print(void);

	TowerFPD& operator=(const TowerFPD& rhs) {
		if (this != &rhs) {
			col     = rhs.col     ;
			row     = rhs.row     ;
			cluster = rhs.cluster ;
			adc_over_ped = rhs.adc_over_ped;
		}
		return *this;
	};
	TObjArray* Lnk_LRUD;
	Int_t IEW;
	Int_t NSTB;
 private:
	ClassDef (TowerFPD, 7);
};
}

#endif

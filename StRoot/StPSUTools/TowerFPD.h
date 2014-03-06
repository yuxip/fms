////////////////////////////////////////////////////////////////////
//
//  definition of a FPD tower (position and energy deposit)
//
////////////////////////////////////////////////////////////////////

#ifndef FPD_TOWER_H
#define FPD_TOWER_H

#include <TObjArray.h>

class StFmsHit;
class StFmsDbMaker;

namespace PSUGlobals {//$NMSPC
class TowerFPD : public TObject {
 public:
  const StFmsHit* hit;
	Int_t   col;         // start from 0: count the columns, moving horizontally (STAR x-coord)
	Int_t   row;         // start from 0: count the rows,    moving vertically   (STAR y-coord)
	Int_t   cluster;
	TowerFPD();
	TowerFPD(const StFmsHit* fmsHit);
	~TowerFPD();
	/**
	 Initialize additional tower geometry information from the database.
	 
	 Return true upon successful initialization, false if something goes wrong.
	 Important: an uninitialized tower should NOT be used!
	 */
	Bool_t initialize(StFmsDbMaker*);
	Bool_t IsEqual(const TObject *obj) const;
	Bool_t IsSortable() const;
	Int_t Compare(const TObject *obj) const ;
	Bool_t IsNeighbor(TowerFPD *a);
	Bool_t SetContext(TObjArray* towers);
	void Print(void);
	TowerFPD& operator=(const TowerFPD& rhs);
	TObjArray* Lnk_LRUD;
	ClassDef (TowerFPD, 7)
};

inline Bool_t TowerFPD::IsEqual(const TObject* obj) const {
  const TowerFPD* other = static_cast<const TowerFPD*>(obj);
  return col == other->col && row == other->row;
}

inline Bool_t TowerFPD::IsSortable() const { return kTRUE; }
}  // namespace PSUGlobals
#endif

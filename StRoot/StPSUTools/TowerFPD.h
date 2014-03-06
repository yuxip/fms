////////////////////////////////////////////////////////////////////
//
//  definition of a FPD tower (position and energy deposit)
//
////////////////////////////////////////////////////////////////////

#ifndef FPD_TOWER_H
#define FPD_TOWER_H

#include <TObjArray.h>

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
	~TowerFPD();
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

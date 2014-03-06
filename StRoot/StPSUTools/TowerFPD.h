#ifndef FPD_TOWER_H
#define FPD_TOWER_H

#include <TObjArray.h>

class StFmsHit;
class StFmsDbMaker;

namespace PSUGlobals {//$NMSPC
/*
 Lightweight wrapper around an StFmsHit for use in tower clustering.
 
 Clustering requires identifying tower neighbors, which is most easily done
 via row and column number. As these aren't stored in StFmsHit, we store an
 StFmsHit pointer here with its row and column number. This removes the need
 to recalculate row and column (which requires database access) each time they
 are needed. We also store the index of the cluster (if any) that the tower
 becomes associated with during clustering.
 
 The TowerFPD does not own the StFmsHit; it merely references it. Therefore it
 is vital that the StFmsHit have a longer lifetime than the TowerFPD i.e.
 do not clear your StFmsCollection until after you have finished clustering!
 
 Inherits from TObject so it can be stored in a ROOT container.
 */
class TowerFPD : public TObject {
 public:
  /** Default constructor */
	TowerFPD();
	/**
	 Constructor.
	 
	 Initialize with an StFmsHit, which the TowerFPD does *now* own. It should
	 therefore have longer lifetime that the TowerFPD.
	 */
	TowerFPD(const StFmsHit* fmsHit);
	/** Assignment operator */
	TowerFPD& operator=(const TowerFPD& rhs);
	/** Destructor */
	~TowerFPD();
	/**
	 Initialize tower row and column information from the database.
	 
	 Return true upon successful initialization, false if something goes wrong.
	 Important: an uninitialized tower should NOT be used!
	 */
	Bool_t initialize(StFmsDbMaker*);
	Bool_t IsEqual(const TObject* tower) const;
	/** Returns true, as TowerFPD can be sorted in a ROOT container */
	Bool_t IsSortable() const;
	/**
	 Function for sorting TowerFPD in order of ascending energy
	 
	 See TObject::Compare() for the convention of return values
	 */
	Int_t Compare(const TObject* tower) const;
	/**
	 Test if another TowerFPD is a neighbor of this tower
	 
	 A neighbor is the tower immediately above, below, left or right of this one
	 i.e. NOT diagonally adjacent towers
	 */
	Bool_t IsNeighbor(TowerFPD* tower);
	/**
	 Set this tower's neighbors from a list of towers
	 
	 Loop over all towers in the list and store any that are immediately
	 adjacent to this one
	 */
	Bool_t SetContext(TObjArray* towers);
  const StFmsHit* hit;  // Not owned by TowerFPD
	Int_t   col;  // Column number, starts at 0, moves horizontally (STAR x-coord)
	Int_t   row;  // Row number, starts at 0, moves vertically (STAR y-coord)
	Int_t   cluster;  // Index of cluster the tower is associated with
	TObjArray* Lnk_LRUD;  // List of neighbor towers
	ClassDef (TowerFPD, 7)
};

inline Bool_t TowerFPD::IsEqual(const TObject* obj) const {
  const TowerFPD* other = static_cast<const TowerFPD*>(obj);
  return col == other->col && row == other->row;
}

inline Bool_t TowerFPD::IsSortable() const { return kTRUE; }
}  // namespace PSUGlobals
#endif

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
	TowerFPD(StFmsHit* fmsHit);
	/** Destructor */
	~TowerFPD();
	/**
	 Initialize tower row and column information from the database.
	 
	 Return true upon successful initialization, false if something goes wrong.
	 Important: an uninitialized tower should NOT be used!
	 */
	Bool_t initialize(StFmsDbMaker*);
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
	/** Returns the hit information for this tower (NULL if unknown) */
	StFmsHit* hit();
	/** Returns the hit information for this tower (NULL if unknown) */
	const StFmsHit* hit() const;
	/** Returns the row of this tower (-1 if unknown) */
	Int_t column() const;
	/** Returns the column of this tower (-1 if unknown) */
	Int_t row() const;
	/** Returns the cluster index of this tower (-1 if unassociated) */
	Int_t cluster() const;
	/** Sets the cluster index and returns the new index */
	Int_t setCluster(Int_t cluster);

 protected:
  StFmsHit* mHit;  //!< Hit information, not owned by TowerFPD
	Int_t mColumn;  ///< Column number, starts at 0, horizontal (STAR x-coord)
	Int_t mRow;  ///< Row number, starts at 0, vertical (STAR y-coord)
	Int_t mCluster;  ///< Index of cluster the tower is associated with
	ClassDef(TowerFPD, 7)
};

inline Bool_t TowerFPD::IsSortable() const { return kTRUE; }

inline StFmsHit* TowerFPD::hit() { return mHit; }

inline const StFmsHit* TowerFPD::hit() const { return mHit; }

inline Int_t TowerFPD::column() const { return mColumn; }

inline Int_t TowerFPD::row() const { return mRow; }

inline Int_t TowerFPD::cluster() const { return mCluster; }

inline Int_t TowerFPD::setCluster(Int_t cluster) {
  mCluster = cluster;
  return cluster;
}
}  // namespace PSUGlobals
#endif

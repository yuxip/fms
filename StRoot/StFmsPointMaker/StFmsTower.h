#ifndef FPD_TOWER_H
#define FPD_TOWER_H

#include <TObjArray.h>

class StFmsHit;
class StFmsDbMaker;

namespace FMSCluster {  // $NMSPC
/*
 Lightweight wrapper around an StFmsHit for use in tower clustering.
 
 Clustering requires identifying tower neighbors, which is most easily done
 via row and column number. As these aren't stored in StFmsHit, we store an
 StFmsHit pointer here with its row and column number. This removes the need
 to recalculate row and column (which requires database access) each time they
 are needed. We also store the index of the cluster (if any) that the tower
 becomes associated with during clustering.
 
 The StFmsTower does not own the StFmsHit; it merely references it. Therefore it
 is vital that the StFmsHit have a longer lifetime than the StFmsTower i.e.
 do not clear your StFmsCollection until after you have finished clustering!
 
 Inherits from TObject so it can be stored in a ROOT container.
 */
class StFmsTower : public TObject {
 public:
  /** Default constructor */
	StFmsTower();
	/**
	 Constructor.
	 
	 Initialize with an StFmsHit, which the StFmsTower does *now* own. It should
	 therefore have longer lifetime that the StFmsTower.
	 */
	StFmsTower(StFmsHit* fmsHit);
  // Use default copy constructor and assignment operator
	/** Destructor */
	~StFmsTower();
	/**
	 Initialize tower row and column information from the database.
	 
	 Return true upon successful initialization, false if something goes wrong.
	 Important: an uninitialized tower should NOT be used!
	 */
	Bool_t initialize(StFmsDbMaker*);
	/** Returns true, as StFmsTower can be sorted in a ROOT container */
	Bool_t IsSortable() const { return kTRUE; }
	/**
	 Test if another StFmsTower is a neighbor of this tower
	 
	 A neighbor is the tower immediately above, below, left or right of this one
	 i.e. NOT diagonally adjacent towers
	 */
	Bool_t isNeighbor(const StFmsTower& tower) const;
	/** Returns the hit information for this tower (NULL if unknown) */
	const StFmsHit* hit() const { return mHit; }
	/** Returns the row of this tower (-1 if unknown) */
	Int_t column() const { return mColumn; }
	/** Returns the column of this tower (-1 if unknown) */
	Int_t row() const { return mRow; }
	/** Returns the cluster index of this tower (-1 if unassociated) */
	Int_t cluster() const { return mCluster; }
	/** Sets the cluster index and returns the new index */
	void setCluster(Int_t cluster) { mCluster = cluster; }

 protected:
  const StFmsHit* mHit;  //!< Hit information, not owned by StFmsTower
	Int_t mColumn;  ///< Column number, starts at 0, horizontal (STAR x-coord)
	Int_t mRow;  ///< Row number, starts at 0, vertical (STAR y-coord)
	Int_t mCluster;  ///< Index of cluster the tower is associated with
	ClassDef(StFmsTower, 7)
};
}  // namespace FMSCluster
#endif
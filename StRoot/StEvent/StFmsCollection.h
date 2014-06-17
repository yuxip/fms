/***************************************************************************
 *
 * $Id: StFmsCollection.h,v 2.1 2010/01/08 22:42:31 ullrich Exp $
 *
 * Author: Jingguo Ma, Dec 2009
 ***************************************************************************
 *
 * Description:
 *
 ***************************************************************************
 *
 * $Log: StFmsCollection.h,v $
 * Revision 2.1  2010/01/08 22:42:31  ullrich
 * Initial Revision.
 *
 **************************************************************************/
#ifndef StFmsCollection_hh
#define StFmsCollection_hh

#include "St_base/StObject.h"
#include "StEvent/StContainers.h"  // StSPtrVecFmsHit/Cluster/Point definitions

class StFmsHit;
class StFmsCluster;
class StFmsPoint;

/**
 Collection of all hits (towers), clusters and points (photons) in the FMS.
 
 This collection owns all these objects, and is itself owned by StEvent.
 It is therefore vital to *not* delete any of the objects stored in this
 container yourself - the collection will handle freeing memory.
 Similarly, any object added to the collection via an add() method must be
 allocated with new, and not be owned anywhere else.
 */
class StFmsCollection : public StObject {
 public:
  /** Constructor */
  StFmsCollection();
  /** Destructor. Deletes all hits, clusters and points */
  ~StFmsCollection();
  /** Add a hit to the collection */
  void addHit(StFmsHit*);
  /** Add a cluster to the collection */
  void addCluster(StFmsCluster*);
  /** Add a point to the collection */
  void addPoint(StFmsPoint*);
  /** Return the number of hits in the collection */
  unsigned int numberOfHits() const;
  /** Return the number of clusters in the collection */
  unsigned int numberOfClusters() const;
  /** Return the number of points in the collection */
  unsigned int numberOfPoints() const;
  /** Return the hit list */
  StSPtrVecFmsHit& hits();
  /** \overload */
  const StSPtrVecFmsHit& hits() const;
  /** Return the cluster list */
  StSPtrVecFmsCluster& clusters();
  /** \overload */
  const StSPtrVecFmsCluster& clusters() const;
  /** Return the cluster list */
  StSPtrVecFmsPoint& points();
  /** \overload */
  const StSPtrVecFmsPoint& points() const;
    
 private:
  StSPtrVecFmsHit mHits;  ///< Owns all FMS hits
  StSPtrVecFmsCluster mClusters;  ///< Owns all FMS clusters
  StSPtrVecFmsPoint mPoints;  ///< Owns all FMS points (photons)
  ClassDef(StFmsCollection, 2)
};

#endif

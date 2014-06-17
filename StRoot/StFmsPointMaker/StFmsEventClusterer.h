// $Id$
//
// $Log$
/**
 \file      StFmsEventClusterer.h
 \brief     Declaration of StFmsEventClusterer, manager class for clustering
 \author    Steven Heppelmann <steveheppelmann@gmail.com>
 \author    Yuxi Pan <yuxipan@physics.ucla.edu>
 \author    Thomas Burton <tpb@bnl.gov>
 \date      2014
 \copyright Brookhaven National Lab
 */
#ifndef STROOT_STFMSPOINTMAKER_STFMSEVENTCLUSTERER_H_
#define STROOT_STFMSPOINTMAKER_STFMSEVENTCLUSTERER_H_

#include <vector>

#include "TObject.h"

#include "StFmsPointMaker/StFmsClusterFinder.h"

namespace FMSCluster {  // $NMSPC
class StFmsClusterFitter;
class StFmsFittedPhoton;
class StFmsGeometry;
class StFmsEventClusterer: public TObject {
 public:
  /**
   Constructor.

   Clustering is done per sub-detector, so initialize with the detector ID and
   FMS geometry information.

   See FMSCluster::StFmsDetectorId for valid detector IDs.
   */
  StFmsEventClusterer(const StFmsGeometry* geometry, Int_t detectorId);
  /** Destructor. */
  ~StFmsEventClusterer();
  /**
   Returns the ID of the detector for which clustering is being performed.

   See FMSCluster::StFmsDetectorId for valid detector IDs.
   */
  int detector() const { return mDetectorId; }
  /**
   Perform cluster-finding and photon-fitting on a list of towers.
   
   Returns true if photon fits to all clusters succeed, or false if one or more
   clusters have a photon fit with a chi-square exceeding the maximum allowed
   value.
   */
  Bool_t cluster(std::vector<FMSCluster::StFmsTower>* towers);
#ifndef __CINT__  // Hide ClusterList from CINT as it uses C++11
  /** Returns the list of clusters in this detector for the event. */
  ClusterList& clusters() { return mClusters; }
  /** \overload */
  const ClusterList& clusters() const { return mClusters; }
#endif  // __CINT__

 private:
#ifndef __CINT__  // Hide ClusterList from CINT as it uses C++11
  /** ClusterList is defined in StFmsClusterFinder.h */
  typedef ClusterList::iterator ClusterIter;
  /** ClusterList is defined in StFmsClusterFinder.h */
  typedef ClusterList::const_iterator ClusterConstIter;
#endif  // __CINT__
  /**
   Disallow copy construction.

   With the various collections of towers, clusters and photons, many of which
   are dynamically allocated, it becomes too complicated to easily implement
   copying, and there is little reason to need it. Therefore prevent it from
   happening.
   */
  StFmsEventClusterer(const StFmsEventClusterer&);
  /** Disallow assignment. */
  StFmsEventClusterer& operator=(const StFmsEventClusterer&);
  /**
   Perform all clustering and photon-fitting for the current event.

   Returns the final status for the event: true for "all good", false for
   "something was bad".
   \todo Change return type to bool.
   */
  Int_t fitEvent();
  /**
   The energy deposit in a cluster by a photon.

   This is the sum over all towers making the cluster.
   */
  Double_t photonEnergyInCluster(Double_t towerWidth,
                                 const StFmsTowerCluster* cluster,
                                 const StFmsFittedPhoton* photon) const;
  /**
   The energy deposit in a tower by a photon.

   Calculated using the shower-shape function.
   */
  Double_t photonEnergyInTower(Double_t towerWidth, const StFmsTower* tower,
                               const StFmsFittedPhoton* photon) const;
  /**
   Perform a 1-photon fit on a single cluster.

   Returns the &chi;<sup>2</sup> of the fit.
   */
  Float_t fitOnePhoton(StFmsTowerCluster* cluster);
#ifndef __CINT__  // Hide Cluster(Const)Iter from CINT as it uses C++11
  /*
   Perform a global fit of all photons in an event.

   Update the (x, y) positions end energies of the photons in each cluster based
   on a global fit including all photons.
   Only makes sense when there is more than one photon in the event.
   Arguments:
    - nPh, number of photons in the event
    - nCl, number of clusters containing those photons
    - first, iterator to the first cluster

   Returns the &chi;<sup>2</sup> of the fit.
   */
  Float_t globalFit(const Int_t, const Int_t, ClusterIter first);
  /*
   Special 2-photon fit for a single cluster.

   Cluster moments must have been calculated first

   Returns the &chi;<sup>2</sup> of the fit.
   */
  Float_t fit2PhotonClust(ClusterIter cluster);
  /*
   Run tests on the lower-energy photon in a 2-photon cluster.

   Return true if the photon passes tests, in which case it is a real photon.
   Return false if it fails, in which case it is a bogus photon due to some
   problem in reconstruction - the cluster is actually a 1-photon cluster.

   Arguments:
    - clusterIndex: index of cluster to test
    - nRealClusters: total number of clusters in the event
   */
  bool validate2ndPhoton(ClusterConstIter cluster) const;
  ClusterList mClusters;  ///< List of clusters in this sub-detector/event
#endif  // __CINT__
  StFmsClusterFinder mClusterFinder;   ///< Cluster-finding routine
  const StFmsGeometry* mGeometry;   ///< FMS geometry for current run
  Int_t mDetectorId;   ///< ID of this FMS sub-detector
  std::vector<FMSCluster::StFmsTower>* mTowers;   ///< Towers to cluster
  StFmsClusterFitter* mFitter;   ///< Routine for fitting photons to clusters
  std::vector<Float_t> mTowerWidthXY;   ///< Geometry for this sub-detector (cm)
  ClassDef(StFmsEventClusterer, 0)
};
}  // namespace FMSCluster
#endif  // STROOT_STFMSPOINTMAKER_STFMSEVENTCLUSTERER_H_

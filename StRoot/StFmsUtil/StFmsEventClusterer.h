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

#include <memory>  // For std::unique_ptr
#include <vector>

#include "TObject.h"

#include "StFmsUtil/StFmsClusterFinder.h"

namespace FMSCluster {  // $NMSPC
class StFmsClusterFitter;
class StFmsFittedPhoton;
class StFmsGeometry;
/**
 Performs clustering and photon-fitting for a single detector and event.

 Clusters towers in a single sub-detector and performs photon shower-shape
 fitting on all the formed clusters.
 */
class StFmsEventClusterer: public TObject {
 public:
  /**
   Constructor.

   Clustering is done per sub-detector, so initialize with the detector ID and
   FMS geometry information. The geometry object is not owned by
   StFmsEventClusterer.

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
   Perform clustering for the current event.

   Returns the number of found clusters.
   */
  Int_t findClusters();
  /**
   Fit all found clusters.

   Returns true if all clusters are fit successfully, false if any fail.
   */
  Bool_t fitClusters();
  /**
   Globally refit all photons from all clusters.

   Returns true if all clusters are fit successfully, false if any fail.
   */
  Bool_t refitClusters();
  /**
   The energy deposit in a cluster by a photon.

   This is the sum over all towers making the cluster.
   */
  Double_t photonEnergyInCluster(const StFmsTowerCluster* cluster,
                                 const StFmsFittedPhoton* photon) const;
  /**
   The energy deposit in a tower by a photon.

   Calculated using the shower-shape function.
   */
  Double_t photonEnergyInTower(const StFmsTower* tower,
                               const StFmsFittedPhoton* photon) const;
  /**
   Perform a 1-photon fit on a single cluster.

   Returns the &chi;<sup>2</sup> of the fit.
   */
  Double_t fit1PhotonCluster(StFmsTowerCluster* cluster);
#ifndef __CINT__  // Hide Cluster(Const)Iter from CINT as it uses C++11
  /*
   Special 2-photon fit for a single cluster.

   Cluster moments must have been calculated first

   Returns the &chi;<sup>2</sup> of the fit.
   */
  Double_t fit2PhotonCluster(ClusterIter cluster);
  /*
   Fit an ambiguous cluster (one that isn't obviously 1- or 2-photon).

   First tries a 1-photon fit. If that fit is good enough, it is set as a
   1-photon cluster.
   Otherwise tries a 2-photon fit and chooses the better result.

   Returns the category of the cluster (EFmsClusterCategory in StFmsCluster).
   */
  Int_t fitAmbiguousCluster(ClusterIter cluster);
  /*
   Perform a global fit of all photons in an event.

   Update the (x, y) positions end energies of the photons in each cluster based
   on a global fit including all photons.
   Only makes sense when there is more than one photon in the event.
   Also, should only be done after prior single- or N-photon fits of the
   clusters concerned, in order to provide starting values for the photon
   positions and energies.
   Arguments:
    - nPhotons, number of photons in the event
    - nClusters, number of clusters containing those photons
    - first, iterator to the first cluster

   Returns the &chi;<sup>2</sup> of the fit.
   */
  Double_t fitGlobalClusters(unsigned int nPhotons, unsigned int nClusters,
                            ClusterIter first);
  /*
   Run tests on the lower-energy photon in a 2-photon cluster.

   Return true if the photon passes tests, in which case it is a real photon.
   Return false if it fails, in which case it is a bogus photon due to some
   problem in reconstruction - the cluster is actually a 1-photon cluster.

   Arguments:
  	const iterator of ClusterList
   */
  bool validate2ndPhoton(ClusterConstIter cluster) const;
  ClusterList mClusters;  ///< List of clusters in this sub-detector/event
  std::unique_ptr<StFmsClusterFitter> mFitter;   ///< Performs photon fits
#endif  // __CINT__
  StFmsClusterFinder mClusterFinder;   ///< Cluster-finding routine
  const StFmsGeometry* mGeometry;   ///< FMS geometry for current run
  Int_t mDetectorId;   ///< ID of this FMS sub-detector
  std::vector<FMSCluster::StFmsTower>* mTowers;   ///< Towers to cluster
  std::vector<Double_t> mTowerWidthXY;   ///< Geometry for this sub-detector (cm)
  ClassDef(StFmsEventClusterer, 0)
};
}  // namespace FMSCluster
#endif  // STROOT_STFMSPOINTMAKER_STFMSEVENTCLUSTERER_H_

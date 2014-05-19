#ifndef STROOT_STFMSPOINTMAKER_STFMSEVENTCLUSTERER_H_
#define STROOT_STFMSPOINTMAKER_STFMSEVENTCLUSTERER_H_

#include <vector>

#ifndef __CINT__
// http://www.boost.org/doc/libs/1_55_0/libs/ptr_container/doc/ptr_container.html
#include <boost/ptr_container/ptr_list.hpp>
#endif  // __CINT__

#include "TF2.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TRandom.h"
#include "TVector3.h"

#include "StFmsGeometry.h"
#include "StFmsClusterFitter.h"
#include "StFmsTowerCluster.h"
#include "StFmsFittedPhoton.h"
#include "StFmsClusterFinder.h"

namespace FMSCluster {  // $NMSPC
class StFmsEventClusterer: public TObject {
 public:
  /** Constructor */
  StFmsEventClusterer(const StFmsGeometry* geometry, Int_t detectorId);
  /** Destructor */
  ~StFmsEventClusterer();
  /** Return the ID of the detector for which clustering is being performed */
  int detector() const { return mDetectorId; }
  /**
   Perform cluster finding and photon fitting on a list of towers
   
   Return true if photon fits to all clusters succeeds, or false if one or more
   clusters have a photon fit with a chi-square exceeding the maximum allowed
   value.
   */
  Bool_t cluster(std::vector<FMSCluster::StFmsTower>* towers);
#ifndef __CINT__  // Hide ClusterList from CINT as it uses Boost
  /** Return the list of clusters in this detector for the event */
  ClusterList& clusters() { return mClusters; }
  /** \overload */
  const ClusterList& clusters() const { return mClusters; }
#endif  // __CINT__

 private:
#ifndef __CINT__  // Hide ClusterList from CINT as it uses Boost
  // ClusterList is defined in StFmsClusterFinder.h
  typedef ClusterList::iterator ClusterIter;
  typedef ClusterList::const_iterator ClusterConstIter;
#endif  // __CINT__
  /** Disallow copy construction */
  StFmsEventClusterer(const StFmsEventClusterer&);
  /** Disallow assignment */
  StFmsEventClusterer& operator=(const StFmsEventClusterer&);
  Int_t fitEvent();
  Double_t photonEnergyInCluster(Double_t towerWidth,
                                 const StFmsTowerCluster* cluster,
                                 const StFmsFittedPhoton* photon) const;
  Double_t photonEnergyInTower(Double_t towerWidth, const StFmsTower* tower,
                               const StFmsFittedPhoton* photon) const;
  Float_t fitOnePhoton(StFmsTowerCluster* cluster);
#ifndef __CINT__  // Hide ClusterList from CINT as it uses Boost
  Float_t globalFit(const Int_t, const Int_t, ClusterIter);
  Float_t fit2PhotonClust(ClusterIter cluster);
  bool validate2ndPhoton(ClusterConstIter cluster) const;
  ClusterList mClusters;
#endif  // __CINT__
  /** Fit clusters for all towers in this detector for the event */
  StFmsClusterFinder mClusterFinder;
  const StFmsGeometry* mGeometry;
  Int_t mDetectorId;
  std::vector<FMSCluster::StFmsTower>* mTowers;
  StFmsClusterFitter* mFitter;
  std::vector<Float_t> mTowerWidthXY;
  ClassDef(StFmsEventClusterer, 7)
};
}  // namespace FMSCluster
#endif  // STROOT_STFMSPOINTMAKER_STFMSEVENTCLUSTERER_H_

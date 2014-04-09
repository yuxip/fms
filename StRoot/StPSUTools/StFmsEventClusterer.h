#ifndef STFMSEVENTCLUSTERER_H
#define STFMSEVENTCLUSTERER_H

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

namespace PSUGlobals {//$NMSPC
class StFmsEventClusterer: public TObject {
 public:
  /** Constructor */
  StFmsEventClusterer(StFmsGeometry* pgeom, Int_t detectorId);
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
  Bool_t cluster(std::vector<PSUGlobals::StFmsTower>* towers);
#ifndef __CINT__
  /** Return the list of clusters in this detector for the event */
  ClusterList& clusters() { return mClusters; }
  /** \overload */
  const ClusterList& clusters() const { return mClusters; }
#endif  // __CINT__

 private:
  Int_t fitEvent();
  Double_t photonEnergyInCluster(Double_t towerWidth,
                                 StFmsTowerCluster* cluster,
                                 StFmsFittedPhoton* photon);
  Double_t photonEnergyInTower(Double_t towerWidth, StFmsTower* tower,
                               StFmsFittedPhoton* photon);
#ifndef __CINT__
  Float_t fitOnePhoton(StFmsTowerCluster*);
  // ClusterList is defined in StFmsClusterFinder.h
  typedef ClusterList::iterator ClusterIter;
  Float_t globalFit(const Int_t, const Int_t, ClusterIter);
  Float_t fit2PhotonClust(ClusterIter);
  bool validate2ndPhoton(ClusterIter cluster);
  ClusterList mClusters;
#endif  // __CINT__
  /** Fit clusters for all towers in this detector for the event */
  StFmsClusterFinder mClusterFinder;
  StFmsGeometry* mGeometry;
  Int_t mDetectorId;
  std::vector<PSUGlobals::StFmsTower>* mTowers;
  StFmsClusterFitter* mFitter;
  std::vector<Float_t> mTowerWidthXY;
  
  static TF1* mEDepCorrection;
  /**
   Response function for nonlinear energy correction, based on cerenkov studies.
   */
  static TF1* GetEDepCorrection();
  ClassDef(StFmsEventClusterer,7);
};
}  // namespace PSUGlobals
#endif

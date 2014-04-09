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

#define MAX_NUMER_CLUSTERS 6

namespace PSUGlobals {//$NMSPC
class ToweFPD;
class StFmsEventClusterer: public TObject {
 public:
  /** Constructor */
  StFmsEventClusterer(StFmsGeometry* pgeom, Int_t detectorId);
  /** Destructor */
  ~StFmsEventClusterer();
  /** Return the ID of the detector for which clustering is being performed */
  int detector() const { return mDetectorId; }
  typedef std::vector<StFmsTower> TowerList;
  /**
   Perform cluster finding and photon fitting on a list of towers
   
   Return true if photon fits to all clusters succeeds, or false if one or more
   clusters have a photon fit with a chi-square exceeding the maximum allowed
   value.
   */
  Bool_t cluster(TowerList* towers);
#ifndef __CINT__
  /** Return the list of clusters in this detector for the event */
  ClusterList& clusters() { return mClusters; }
  /** \overload */
  const ClusterList& clusters() const { return mClusters; }
#endif  // __CINT__

 private:
#ifndef __CINT__
  ClusterList mClusters;
  Float_t FitOnePhoton(StFmsTowerCluster*);
  // ClusterList is defined in StFmsClusterFinder.h
  typedef ClusterList::iterator ClusterIter;
  Float_t GlobalFit(const Int_t, const Int_t, ClusterIter);
  Float_t Fit2PhotonClust(ClusterIter);
  bool validate2ndPhoton(ClusterIter cluster);
#endif  // __CINT__
  /** Fit clusters for all towers in this detector for the event */
  Int_t FitEvent();
  Double_t EnergyInClusterByPhoton(Double_t widthLG, StFmsTowerCluster*, StFmsFittedPhoton*);
  Double_t EnergyInTowerByPhoton(Double_t, StFmsTower* , StFmsFittedPhoton* );
  StFmsClusterFinder mClusterFinder;
  StFmsGeometry* p_geom;
  Int_t mDetectorId;
  TowerList* towers;
  StFmsClusterFitter* fitter;
  Double_t step[3*StFmsClusterFitter::MAX_NUMB_PHOTONS+1];
  std::vector<Float_t> widLG;
  static TF1* EDepCorrection;
  /**
   Response function for nonlinear energy correction, based on cerenkov studies.
   */
  static TF1* GetEDepCorrection();
  ClassDef(StFmsEventClusterer,7);
};
}  // namespace PSUGlobals
#endif

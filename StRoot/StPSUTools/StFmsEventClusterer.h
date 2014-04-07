#ifndef YIQUN_H
#define YIQUN_H

#include <iostream>
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
#include "PhotonHitFPD.h"
#include "StFmsClusterFinder.h"

#define MAX_NUMER_CLUSTERS 6

namespace PSUGlobals {//$NMSPC
class ToweFPD;
class Yiqun: public TObject {
 public:
  TowerUtil* pTowerUtil;
  Int_t nrows;
  Int_t ncols;
#ifndef __CINT__
  ClusterList mClusters;
  Float_t FitOnePhoton(HitCluster*);
  // ClusterList is defined in StFmsClusterFinder.h
  typedef ClusterList::iterator ClusterIter;
  Float_t GlobalFit(const Int_t, const Int_t, ClusterIter);
  Float_t Fit2PhotonClust(ClusterIter);
  bool validate2ndPhoton(ClusterIter cluster);
  ClusterList& clusters() { return mClusters; }
  const ClusterList& clusters() const { return mClusters; }
#endif  // __CINT__
  Int_t FitEvent(Int_t nTows, Int_t &nClusts, Int_t &nRealClusts, Bool_t &junkyEvent);
  Double_t EnergyInClusterByPhoton(Double_t widthLG, HitCluster*, PhotonHitFPD*);
  Double_t EnergyInTowerByPhoton(Double_t, TowerFPD* , PhotonHitFPD* );
  Yiqun(Geom* pgeom, Int_t detectorId);
  typedef std::vector<TowerFPD> TowerList;
  /**
   Perform cluster finding and photon fitting on a list of towers
   
   Return true if photon fits to all clusters succeeds, or false if one or more
   clusters have a photon fit with a chi-square exceeding the maximum allowed
   value.
   */
  Bool_t cluster(TowerList* towers);
  ~Yiqun();
  Geom* p_geom;
  Int_t mDetectorId;
  Int_t NTower;
  TowerList* towers;
  TObjArray* tow_Arr;
  FitTower* fitter;
  Int_t NPh; 
  Int_t NClusts;
  Int_t NRealClusts;
  Float_t posDif_2PC;    // in unit of Lead-glass
  Float_t eneRat_2PC;
  Float_t dggPara[6];
  Float_t thetaPara;
  Float_t posDif_Gl;    // in unit of "cm"
  Float_t eneRat_Gl;
  Float_t maxGood1PhChi2NDF;
  Float_t minHTEneOverPhoton;
  Float_t maxHTEneOverPhoton;
  Float_t maxRatioSpill;
  Float_t MaxChi2Catag2;
  Float_t minRealClusterEne;
  Int_t maxHitsInRealCluster;
  Double_t step[3*FitTower::MAX_NUMB_PHOTONS+1];
  std::vector<Float_t> widLG;
  static TF1* EDepCorrection;
  /**
   Response function for nonlinear energy correction, based on cerenkov studies.
   */
  static TF1* GetEDepCorrection();
  ClassDef(Yiqun,7);
};
}  // namespace PSUGlobals
#endif

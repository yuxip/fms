#ifndef YIQUN_H
#define YIQUN_H
#define MAX_NUMER_CLUSTERS 6

#include "Geom.h"
#include "TLorentzVector.h"
#include "TObjArray.h"
#include "TF2.h"
#include "TVector3.h"
#include "FitTower.h"
#include "TowerUtil.h"
#include "HitCluster.h"
#include "PhotonHitFPD.h"

#include "TMath.h"
#include "TRandom.h"
#include <iostream>

namespace PSUGlobals {//$NMSPC
class Yiqun: public TObject
{
 private:
  
  ClassDef(Yiqun,7);
  
 public:
  TowerUtil* pTowerUtil;
  Int_t nrows;
  Int_t ncols;
  
  Float_t FitOnePhoton(HitCluster*);
  Float_t GlobalFit(const Int_t, const Int_t, HitCluster*);
  Float_t Fit2PhotonClust(HitCluster*);
  Int_t FitEvent(Int_t nTows, Int_t &nClusts, Int_t &nRealClusts, Bool_t &junkyEvent);
  Double_t EnergyInClusterByPhoton(Double_t widthLG, HitCluster*, PhotonHitFPD*);
  Double_t EnergyInTowerByPhoton(Double_t, TowerFPD* , PhotonHitFPD* );
  typedef std::vector<TowerFPD> TowerList;
  Yiqun(TowerList* pEm,Geom* pgeom,Int_t iew ,Int_t nstb);
  void Y(TowerList*);
  ~Yiqun();
  Geom* p_geom;
  Int_t EW;
  Int_t NSTB;
  Int_t NTower;
  TowerList* towers;
  TObjArray* tow_Arr;
  FitTower* fitter;
  HitCluster clust[MAX_NUMER_CLUSTERS];
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
  Double_t widLG[2];
  bool validate2ndPhoton(int clusterIndex, int nRealClusters);
  static TF1* EDepCorrection;
  /**
   Response function for nonlinear energy correction, based on cerenkov studies.
   */
  static TF1* GetEDepCorrection();
};
}
#endif

#ifndef TOWER_UTIL_H
#define TOWER_UTIL_H

#include "TH1.h"
#include "TObjArray.h"
#include "TowerFPD.h"
#include "HitCluster.h"
#include "TPaveText.h"
#include "TH2.h"
#include "TF2.h"

namespace PSUGlobals {//$NMSPC
class TowerUtil {
 public:
  TowerUtil();
  ~TowerUtil();
  Int_t FindTowerCluster(TObjArray *arrTow, HitCluster *clust);
  void CalClusterMoment(HitCluster *clust);
  Int_t CatagBySigmXY(HitCluster *clust);
  void PrintTowers(const TowerFPD *tows);
  void SetMomentEcutoff(Float_t ecoff=0.5){Ecutoff=ecoff;};
  Float_t GetMomentEcutoff(){return Ecutoff;};
 private:
  Float_t Ecutoff;  
  TObjArray* neighbor;
  TObjArray* arrValley;
  ClassDef(TowerUtil,3);
};
}  // namespace PSUGlobals
#endif

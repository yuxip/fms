#ifndef TOWER_UTIL_H
#define TOWER_UTIL_H

#include <list>

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
  typedef std::list<TowerFPD*> TowerList;
 private:
  Float_t Ecutoff;  
//  static const Int_t NS_ETA_NUM = 7;
//  static const Int_t NS_PHI_NUM = 7;
//  static const Int_t TB_ETA_NUM = 5;
//  static const Int_t TB_PHI_NUM = 5;
  static const Int_t maxNClusters = 6;
  static const Int_t nNSTow = 49;//NS_ETA_NUM * NS_PHI_NUM;
  static const Int_t nTBTow = 25;//TB_ETA_NUM * TB_PHI_NUM ;
  // number of "peaks" that has the same shortest distance to a "valley" tower
  Int_t nClusts;
  Int_t nPeakSameDist[nNSTow];
  // store which "peaks" have the same shortest distance to a "valley" tower
  Int_t peaksToValley[nNSTow][maxNClusters];
  TObjArray* neighbor;
  TObjArray* arrValley;
  unsigned associateTowersWithClusters(TowerList& neighbor,
                                                HitCluster *clust,
                                                TObjArray* arrValley);
  unsigned associateResidualTowersWithClusters(TowerList& neighbor,
                                                HitCluster *clust);
  ClassDef(TowerUtil,3);
};
}  // namespace PSUGlobals
#endif

#ifndef TOWER_UTIL_H
#define TOWER_UTIL_H

#include <list>

#include <Rtypes.h>

class TObjArray;

namespace PSUGlobals {//$NMSPC
class HitCluster;
class TowerFPD;
class TowerUtil {
 public:
  TowerUtil();
  ~TowerUtil();
  Int_t FindTowerCluster(TObjArray* towers, HitCluster* clusters);
  void CalClusterMoment(HitCluster* cluster);
  Int_t CatagBySigmXY(HitCluster* cluster);
  void SetMomentEcutoff(Float_t ecoff=0.5) { Ecutoff=ecoff; }
  Float_t GetMomentEcutoff() { return Ecutoff; }
  typedef std::list<TowerFPD*> TowerList;
 private:
  static const Int_t maxNClusters = 6;
  static const Int_t nNSTow = 49;  //NS_ETA_NUM * NS_PHI_NUM;
  static const Int_t nTBTow = 25;  //TB_ETA_NUM * TB_PHI_NUM ;
  Float_t Ecutoff;  
  // number of "peaks" that has the same shortest distance to a "valley" tower
  Int_t nClusts;
  unsigned associateTowersWithClusters(TowerList& neighbors,
                                       HitCluster* clusters,
                                       TObjArray* valleys);
  unsigned associateValleyTowersWithClusters(TowerList& neighbors,
                                             HitCluster* clusters,
                                             TObjArray* valleys);
  unsigned associateResidualTowersWithClusters(TowerList& neighbors,
                                               HitCluster* clusters);
  unsigned associateSubThresholdTowersWithClusters(TowerList& towers,
                                                   HitCluster* clusters);
  ClassDef(TowerUtil,3);
};
}  // namespace PSUGlobals
#endif

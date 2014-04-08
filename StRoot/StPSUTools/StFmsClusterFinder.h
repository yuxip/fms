#ifndef TOWER_UTIL_H
#define TOWER_UTIL_H

#include <list>

#ifndef __CINT__
// http://www.boost.org/doc/libs/1_55_0/libs/ptr_container/doc/ptr_container.html
#include <boost/ptr_container/ptr_list.hpp>
#endif  // __CINT__

#include <Rtypes.h>

class TObjArray;

namespace PSUGlobals {//$NMSPC
class StFmsTowerCluster;
class TowerFPD;
#ifndef __CINT__
typedef boost::ptr_list<StFmsTowerCluster> ClusterList;
#endif  // __CINT__
class TowerUtil {
 public:
  typedef std::list<TowerFPD*> TowerList;
  TowerUtil();
  ~TowerUtil();
#ifndef __CINT__
  Int_t FindTowerCluster(TowerList* towers, ClusterList* clusters);
#endif  // __CINT__
  void CalClusterMoment(StFmsTowerCluster* cluster);
  Int_t CatagBySigmXY(StFmsTowerCluster* cluster);
  void SetMomentEcutoff(Float_t ecoff=0.5) { Ecutoff=ecoff; }
  Float_t GetMomentEcutoff() { return Ecutoff; }
 private:
  static const Int_t maxNClusters = 6;
  static const Int_t nNSTow = 49;  //NS_ETA_NUM * NS_PHI_NUM;
  static const Int_t nTBTow = 25;  //TB_ETA_NUM * TB_PHI_NUM ;
  Float_t Ecutoff;  
  // number of "peaks" that has the same shortest distance to a "valley" tower
  Int_t nClusts;
#ifndef __CINT__
  unsigned locateClusterSeeds(TowerList* towers, TowerList* neighbors,
                              ClusterList* clusters);
  unsigned associateTowersWithClusters(TowerList* neighbors,
                                       ClusterList* clusters,
                                       TObjArray* valleys);
  unsigned associateValleyTowersWithClusters(TowerList* neighbors,
                                             ClusterList* clusters,
                                             TObjArray* valleys);
  unsigned associateResidualTowersWithClusters(TowerList* neighbors,
                                               ClusterList* clusters);
  unsigned associateSubThresholdTowersWithClusters(TowerList* towers,
                                                   ClusterList* clusters);
#endif  // __CINT__
  ClassDef(TowerUtil,3);
};
}  // namespace PSUGlobals
#endif

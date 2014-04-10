#ifndef STFMSCLUSTERFINDER_H
#define STFMSCLUSTERFINDER_H

#include <list>

#ifndef __CINT__  // Hide Boost headers from CINT as they will confuse it
// http://www.boost.org/doc/libs/1_55_0/libs/ptr_container/doc/ptr_container.html
#include <boost/ptr_container/ptr_list.hpp>
#endif  // __CINT__

#include <Rtypes.h>  // Provides ROOT ClassDef macro

class TObjArray;

namespace FMSCluster {  // $NMSPC
class StFmsTowerCluster;
class StFmsTower;
// Typedef a pointer list of tower clusters for convenience
// Hide it from CINT as it won't know how to handle the Boost header
#ifndef __CINT__
typedef boost::ptr_list<StFmsTowerCluster> ClusterList;
#endif  // __CINT__
/** Form clusters adjacent FMS towers */
class StFmsClusterFinder {
 public:
  // Typedef StFmsTower pointer list for convenience
  typedef std::list<FMSCluster::StFmsTower*> TowerList;
  /** Constructor */
  StFmsClusterFinder();
  /** Destructor */
  ~StFmsClusterFinder();
  /**
   Calculate moments (means, sigmas of tower positions) for a cluster
   
   Also updates cluster with the current number of towers in it.
   */
  void calculateClusterMoments(StFmsTowerCluster* cluster) const;
  /**
   Categorise a cluster based on its energy and tower distribution
   
   Set the cluster's category field and return the category
   See EFmsClusterCategory in StEvent/StFmsCluster.h for valid categories
   */
  int categorise(StFmsTowerCluster* cluster);
  /** Set energy cutoff on towers for when calculating cluster moments */
  void setMomentEnergyCutoff(float cutoff = 0.5) { mEnergyCutoff = cutoff; }
  /** Return energy cutoff on towers for when calculating cluster moments */
  float getMomentEnergyCutoff() const { return mEnergyCutoff; }
#ifndef __CINT__  // Hide ClusterList from CINT
  /**
   Find clusters from a collection of input towers
   
   Populate the cluster list with the found clusters
   Return the number of found clusters
   Arguments
    - towers: input tower list. Note it may be modified (elements erased or
              reordered) so do not rely on its contents after calling
              findClusters()
    - clusters: output cluster list. Newly found clusters are appended to the
                list.
   */
  int findClusters(TowerList* towers, ClusterList* clusters);
#endif  // __CINT__

 private:
#ifndef __CINT__  // Hide ClusterList from CINT
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
  static const Int_t kMaxNClusters = 6;
  Float_t mEnergyCutoff;  
  Int_t mNClusts;
  ClassDef(StFmsClusterFinder, 3)
};  // class StFmsClusterFinder
}  // namespace FMSCluster
#endif

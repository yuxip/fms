#include "TowerUtil.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <list>
#include <memory>

#include <TCollection.h>
#include <TObjArray.h>

#include "StPSUTools/TowerFPD.h"
#include "StPSUTools/HitCluster.h"

namespace {
/*
 Cluster-finding constants used within this file.
 Some are only used once, but it is easier to keep track of their values, and
 more self-documenting, if they are all named, and set in the same location.
 */
const Float_t maxDistanceFromPeak = 0.3;
const Int_t minTowerCatag02 = 5;
const Float_t cutEcSigma[2][2] = {{2.1, 7.0}, {2.1, 2.0}};
const Float_t minEcSigma2Ph = 35.;
const Float_t maxEcSigma1Ph = 10.;
const Float_t minTowerEnergy = 0.01;
const Float_t minRatioPeakTower = 1.6;
// Extreme distance between towers (no distance can be this large!)
const Float_t ExtremelyFaraway = 99999 ;

typedef PSUGlobals::TowerUtil::TowerList TowerList;
typedef TowerList::iterator TowerIter;
typedef TowerList::reverse_iterator TowerRIter;

using PSUGlobals::TowerFPD;

/*
 Test for a tower that can be a cluster peak.
 
 Returns true if a tower can be a peak tower, given a global population of
 known non-peak towers and a minimum ratio between the energy of peak towers and
 non-peak towers. Returns false if the tower can not possibly be a peak. Note
 that returning true does not mean the tower *is* a peak, merely that it *can*
 be (i.e. it is consistent with that hypothesis given this input).
 */
Bool_t couldBePeakTower(TowerFPD* tower, TowerList* nonPeakTowers) {
  Bool_t couldBePeak(true);
  for (TowerIter i = nonPeakTowers->begin(); i != nonPeakTowers->end(); ++i) {
    // Compare this tower's energy with that of its immediate neighbours
    if (tower->IsNeighbor(*i)) {
      if (tower->energy < minRatioPeakTower * (*i)->energy) {
        couldBePeak = false;
        break;
      }  // if
    }  // if
  } // end of for loop over non-peak towers
  return couldBePeak;
}

/**
 Populate an STL container of pointers from a ROOT collection.
 
 Returns the size of the filled container.
 No objects are duplicated; the STL container points to the same objects (and
 in the same order) as the ROOT collection.
 */
template<typename StlContainer>
typename StlContainer::size_type fillStlContainerFromRootCollection(
    const TCollection& collection, StlContainer* container) {
  TIter next(&collection);
  typedef typename StlContainer::value_type Pointer;
  Pointer element(NULL);
  while ((element = static_cast<Pointer>(next()))) {
    container->push_back(element);
  }  // while
  return container->size();
};

/**
 Comparison function to sort towers in order of ascending energy.
 */
struct AscendingTowerEnergySorter {
  bool operator()(const TowerFPD* a, const TowerFPD* b) const {
    return a->Compare(b) < 0;
  }
};

/**
 Comparison function to sort towers in order of descending energy.
 */
struct DescendingTowerEnergySorter {
  bool operator()(const TowerFPD* a, const TowerFPD* b) const {
    return a->energy > b->energy;
  }
};

/**
 Predicate determining if a test tower is a neighbour of a reference tower.
 
 A neighbor is defined as the tower immediately above, below or to the side
 i.e not diagonally adjacent cells.
 If the test tower has energy below the global cutoff, always return false,
 even if it is physically a neighbour of the reference tower.
 */
struct TowerIsNeighbor
    : public std::binary_function<TowerFPD*, TowerFPD*, bool> {
  bool operator()(TowerFPD* test, TowerFPD* reference) const {
    if(test->energy < minTowerEnergy) {
      return false;
    }  // if
    return test->IsNeighbor(reference);
  }
};  // End of class TowerIsNeighbor

/*
 Predicate testing for tower energy above the global cutoff
 */
struct TowerEnergyIsAboveThreshold {
  bool operator()(const TowerFPD* tower) const {
    return !(tower->energy < minTowerEnergy);
  }
};

/*
 Filter out towers below the minimum energy threshold from a list.
 
 Returns a pointer list of below-threshold towers and erases those towers from
 the input list. The order of towers after filtering is not guaranteed.
 */
TowerList filterTowersBelowEnergyThreshold(TowerList* towers) {
  // Move towers above threshold to the front and those below to the end
  // newEnd marks the end of the above-threshold towers and the beginning of the
  // below-threshold towers
  TowerIter newEnd = std::partition(towers->begin(), towers->end(),
                                    TowerEnergyIsAboveThreshold());
  // Store the below-threshold towers in a new list
  TowerList belowThreshold(newEnd, towers->end());
  // Remove the below-threshold towers from the input list
  towers->erase(newEnd, towers->end());
  return belowThreshold;
}

/* There are different ways of calculating a tower-to-cluster distance */
enum ETowerClusterDistance {
  kPeakTower,  // Distance from tower to peak tower in cluster
  kClusterCenter  // Distance from tower to calculated center of cluster
};

/**
 Sort the towers in an array of clusters in order of ascending energy.
 
 \todo Replace cluster list with an STL list, which will make sorting and
       reversing much simpler than with a ROOT container
 */
void sortTowersEnergyAscending(PSUGlobals::HitCluster* cluster, int nClusters) {
  for (Int_t i(0); i < nClusters; ++i) {
    (cluster[i].tow)->UnSort();
    (cluster[i].tow)->Sort();  // Sorts by *descending* energy
    // Do an exchange of towers: 0<-->N-1, 1<-->N-2... to get ascending
    Int_t n = (cluster[i].tow)->GetEntriesFast();
    for(Int_t j(0); j < n / 2; ++j) {
      TowerFPD* t1 = static_cast<TowerFPD*>((cluster[i].tow)->RemoveAt(j));
      TowerFPD* t2 = static_cast<TowerFPD*>(
        (cluster[i].tow)->RemoveAt(n - 1 - j));
      (cluster[i].tow)->AddAt(t1, n - 1 - j);
      (cluster[i].tow)->AddAt(t2, j);
    }  // for
  }  // for
}
}  // unnamed namespace

namespace PSUGlobals {
/*
 Association information between a tower and clusters.
 
 Store a tower, and a list of clusters with which it could potentially be
 associated. This information is used in determining the true cluster with
 which the tower is actually associated.
 
 Inherits from TObject to allow it to be placed in a ROOT container.
 */ 
class TowerClusterAssociation : public TObject {
 public:
  /*
   Constructor.
   
   Initialise with the tower of interest.
   */
  TowerClusterAssociation(TowerFPD* tower) : mTower(tower) { }
  /* Returns this tower */
  TowerFPD* tower() { return mTower; }
  /* Returns the list of potential associate clusters */
  std::list<HitCluster*>* clusters() { return &mClusters; }
  /*
   Calculate the separation between this tower and another.
   
   The separation is in local row-column coordinates i.e. the distance is in
   number of towers, not in cm.
   */
  double separation(TowerFPD* tower) {
    return sqrt(pow(tower->col - mTower->col, 2.) +
                pow(tower->row - mTower->row, 2.));
  }
  /*
   Calculate the separation between this tower and a cluster.

   Distances are calculated as follows:
    - distance=kPeakTower: distance from tower to cluster peak tower
    - distance=kClusterCenter: distance from tower to cluster center, based on
                               cluster moment calculation
   The separation is in local row-column coordinates i.e. the distance is in
   number of towers, not in cm.
   */
  double separation(HitCluster* cluster, const ETowerClusterDistance distance) {
    if (kPeakTower == distance) {
      TowerFPD* peak = static_cast<TowerFPD*>(cluster->tow->First());
      return separation(peak);
    } else {
      // Use calculated cluster center (x0, y0)
      return sqrt(pow(cluster->x0 - (mTower->col - 0.5), 2.) +
                  pow(cluster->y0 - (mTower->row - 0.5), 2.));
    }  // if
  }
  /*
   Returns true if this tower can be associated with a cluster
   
   Association of this tower with a cluster is defined as:
    - this tower energy less than the cluster peak tower energy
    - physically adjacent to at least one tower in the cluster, so long as...
    - this energy less than minRatioPeakTower * adjacent tower energy
      i.e. this tower cannot fulfil the criterion for being a peak wrt to
      the adjacent tower
   This is really a "potential" association; we use the information in this
   class to determine the (single) cluster that this tower is actually part of,
   in the case that there is more than one potential associate.
   */
  bool canAssociate(HitCluster* cluster) {
    // The peak tower in a cluster is always the first
    TowerFPD* peak = static_cast<TowerFPD*>(cluster->tow->First());
    // Make sure that this tower has lower energy than the peak, but be careful;
    // because of digitization, it is possible that the "neighbor" tower
    // has the exact same energy as the peak tower, not just less
    if (peak->energy < mTower->energy) {
      return false;
    }  // if
    // Loop over all towers in this cluster to see if this tower is
    // physically adjacent to any of them.
    for(Int_t i(0); i < cluster->tow->GetEntriesFast(); ++i) {
      TowerFPD* clusterTower = static_cast<TowerFPD*>(cluster->tow->At(i));
      // Place an energy selection when determining adjacent towers, as a
      // neighbor cannot exceed an adjacent tower by a factor more than
      // minRatioPeakTower, otherwise it will be considered a peak itself.
      if (mTower->IsNeighbor(clusterTower) &&
          mTower->energy < minRatioPeakTower * clusterTower->energy) {
        return true;  // Stop looping once we find any match
      }  // if
    }  // for loop over all towers in a cluster
    return false;
  }
  /*
   Attempt to add a potential associate cluster with this tower.
   
   Add the cluster to the list of potential associates if this tower can
   associate with it (see canAssociate()). However, if there is already one or
   more clusters in the list:
   - if the new cluster is closer to the tower, replace the existing cluster(s)
   - if the new cluster if further away, do not add it
   - if the new cluster is *exactly* the same separation as the existing
     cluster(s), add it to the list but keep the existing ones
   i.e. at any time there can only be clusters of the same (minimal) separation
   from the tower, but there can be multiple clusters of identical separation.

   See separation(HitCluster*) for the meaning of the distance argument.
   
   Return true if the new cluster is added, false if not.
   */
  bool add(HitCluster* cluster, const ETowerClusterDistance distance) {
    bool inserted(false);
    if (canAssociate(cluster)) {
      if (mClusters.empty()) {
        // There is nothing in the list yet, add the cluster, simples!
        mClusters.push_back(cluster);
        mTower->cluster = cluster->index;
        inserted = true;
      } else {
        // Cluster(s) are already present, so only add the new one if it is
        // not further away. If it is closer, remove the existing cluster.
        double distNew = separation(cluster, distance);
        double distOld = separation(mClusters.front(), distance);
        /** \todo I don't like using simple float comparison here, look into a
                  more robust method */
        // If the new cluster is closer, remove the old ones
        if (distNew < distOld) {
          mClusters.clear();
        } // if
        // Add the new cluster if it is not further away than existing ones
        if (distNew <= distOld) {
          mClusters.push_back(cluster);
          mTower->cluster = cluster->index;
          inserted = true;
        }  // if
      }  // if
    }  // if
    return inserted;
  }
  /*
   Calculate the nearest cluster out of the list of potential associates.
   
   The distance is that between this tower and the cluster centre, (x0, y0),
   therefore HitCluster::CalClusterMoment() must have been called before doing
   this in order to calculate x0 and y0 of the cluster.
   
   Returns NULL if there are no clusters in the list.
   */
  HitCluster* nearestCluster() {
    HitCluster* nearest(NULL);
    double minDist = ExtremelyFaraway;
    std::list<HitCluster*>::iterator i;
    for (i = mClusters.begin(); i != mClusters.end(); ++i) {
      float distance = separation(*i, kClusterCenter);
      // Check if the distance to the "center" of this cluster is smaller
      if (distance < minDist) {
        minDist = distance;
        nearest = *i;
      }  // if
    }  // for
    return nearest;
  }
  TowerFPD* mTower;
  std::list<HitCluster*> mClusters;
};

/**
 Associate tower with clusters
 
 Go through a list of unassociated neighbor towers and try to associate each one
 with a cluster. If a neighbor is associated with a single cluster, add it to
 that cluster and remove it from the neighbor list.
 
 If a neighbor could associate with more than one cluster based on currently
 available information, remove it from the neighbor list and add it to the
 valley list. We will work out the association of the valley towers later.
 
 Return the number of neighbors either associated with clusters or placed in the
 valley i.e. the number removed from the neighbor list.
 */
unsigned TowerUtil::associateTowersWithClusters(TowerList& neighbor,
                                                HitCluster *clust,
                                                TObjArray* valleys) {
  TowerList associated;  // Store neighbors we associate
  // Towers are sorted in ascending energy, so use reverse iterator to go from
  // highest to lowest energy
  TowerRIter tower;
  for (tower = neighbor.rbegin(); tower != neighbor.rend(); ++tower) {
    // Populate association information of this tower with each cluster
    std::auto_ptr<TowerClusterAssociation> association(
      new TowerClusterAssociation(*tower));
    for (Int_t i(0); i < nClusts; ++i) {
      association->add(&clust[i], kPeakTower);
    }  // loop over all clusters
    // Attempt to move the tower to the appropriate cluster
    if (association->clusters()->size() == 1) {
      // Only one peak is closest to the tower; the tower belongs to this peak
      association->clusters()->front()->tow->Add(*tower);
      associated.push_back(*tower);
    } else if (association->clusters()->size() > 1) {
      // Multiple potential clusters, need to do something more sophisticated
      // Add this association to the "valley" array so we can use it later
      valleys->Add(association.release());
      associated.push_back(*tower);
    }  // if
  } // loop over TObjArray "neighbor"
  // Remove associated neighbors from the neighbor list
  for (TowerIter i = associated.begin(); i != associated.end(); ++i) {
    neighbor.remove(*i);
  }  // for
  return associated.size();
}

/**
 Associate tower with clusters
 
 These are towers that were left over after the first round of association.
 The clusters have now had towers associated with them, so the cluster moment
 can be calculated to give a better measure of the cluster centre. This is then
 used to calculate the tower-cluster separation.
 
 Return the number of neighbors associated with clusters.
 */
unsigned TowerUtil::associateResidualTowersWithClusters(TowerList& neighbor,
                                                        HitCluster* clust) {
  TowerList associated;
  TowerRIter tower;
  for (tower = neighbor.rbegin(); tower != neighbor.rend(); ++tower) {
    // Populate tower-cluster association information
    TowerClusterAssociation association(*tower);
    for (Int_t i(0); i < nClusts; ++i) {
      // There are already some towers in the cluster so we can use a computed
      // cluster center to give a better estimate of tower-cluster separation
      CalClusterMoment(&clust[i]);
      association.add(&clust[i], kClusterCenter);
    }  // loop over all clusters
    if (!association.clusters()->empty()) {
      HitCluster* cluster = association.clusters()->front();
      (*tower)->cluster = cluster->index;
      cluster->tow->Add(*tower);
      associated.push_back(*tower);
    }  // if
  } // loop over TObjArray "neighbor"
  for (TowerIter i = associated.begin(); i != associated.end(); ++i) {
    neighbor.remove(*i);
  }  // for
  return associated.size();
}

unsigned TowerUtil::associateValleyTowersWithClusters(TowerList& neighbors,
                                                      HitCluster* clusters,
                                                      TObjArray* valleys) {
  unsigned size = neighbors.size();
  for (Int_t i(0); i < valleys->GetEntriesFast(); ++i) {
    TowerClusterAssociation* association =
      static_cast<TowerClusterAssociation*>(valleys->At(i));
    HitCluster* cluster = association->nearestCluster();
    if (cluster) {
      // Move the tower to the appropriate cluster
      association->tower()->cluster = cluster->index;
      neighbors.remove(association->tower());
      cluster->tow->Add(association->tower());
    } else {
      std::cout << "Something is wrong! The following \"Valley\" tower does "
        << "not belong to any cluster! Error!" << std::endl;
      association->tower()->Print();
      std::cout << "!!!!!!!!\n" << std::endl;
    }  // if (cluster)
  }  // end of for loop over valley towers
  return size - neighbors.size();
}

TowerUtil::TowerUtil() : nClusts(0) {
  SetMomentEcutoff();
}

TowerUtil::~TowerUtil() {
}

Int_t TowerUtil::FindTowerCluster(TObjArray *inputTow, HitCluster *clust) {
  TowerList arrTow;
  fillStlContainerFromRootCollection(*inputTow, &arrTow);
  // Remove towers below energy threshold, but save them for later use
  TowerList belowThreshold = filterTowersBelowEnergyThreshold(&arrTow);
  // the neighbor TObjArray
  TowerList neighbor;
  // the "valley" TObjArray
  TObjArray valleys(16);
  valleys.SetOwner(true);
  arrTow.sort(DescendingTowerEnergySorter());
  // "TObjArray::Sort()" sorts the array from small to big ( [0]<=[1]<=...<=[48] )
  // need to take care of that
  // have to go last object first!
  // The algoriths is such, first get the highest tower (which is ALWAYS the last one),
  // and it and all its neighbors to the next cluster. Then repeat the process over the
  // remaining towers.
  while (!arrTow.empty() && nClusts < maxNClusters) {
    // By design, this tower is the highest tower in "arrTow", but it could be lower
    // than a tower in "neighbor"
    TowerFPD* high = arrTow.front();
    arrTow.pop_front();
		// 2003-08-15
		// Fix a logical loop hole in deciding if a tower is
		//    a peak; Need to first compare the highest tower
		//    remained in "arrTow" to all towers in "neighbor" which is its neighbor,
		//    and if it is lower than any of those, it is
		//    a neighbor. Move it to "neighbor" and continue to
		//    the next tower in "arrTow".
    if (couldBePeakTower(high, &neighbor)) {
      // Add "high" to cluster and move towers neighboring "high" to "neighbor"
      high->cluster = nClusts;
      clust[nClusts].index = nClusts;
      (clust[nClusts].tow)->Add(high);
      nClusts++ ;
      // Partition the remaining towers so that neighbours of the high tower are
      // placed at the beginning, and non-neighbours placed at the end. Use
      // stable_partition so we don't alter the energy ordering.
      TowerIter neighborEnd =
        std::stable_partition(arrTow.begin(), arrTow.end(),
                              std::bind2nd(TowerIsNeighbor(), high));
      // Copy neighbors to the neighbor list, erase them from the tower list
      neighbor.insert(neighbor.end(), arrTow.begin(), neighborEnd);
      arrTow.erase(arrTow.begin(), neighborEnd);
    } else {  // Not a peak, add it to the neighbor collection
      neighbor.push_back(high);
    }  // when "high" is a "peak"
    // A tower separated from neighbors only by towers of the same energy will
    // become a peak by the above logic. To close this loophole, loop again
    // over towers and move any with energy <= any of its neighbors to the
    // neighbor list.
    TowerIter towerIter = arrTow.begin();
    while (towerIter != arrTow.end()) {
      // Need to remove list items whilst iterating, so be careful to increment
      // the iterator before erasing items to avoid iterator invalidation
      if (!couldBePeakTower(*towerIter, &neighbor)) {
        neighbor.push_back(*towerIter);
        arrTow.erase(towerIter++);  // Increment will evaluate before erase()
      } else {
        ++towerIter;
      }  // if
    }  // while
  }  // End of for loop over "arrTow"
  // We have now found all peaks. Now decide the affiliation of neighbor towers
  // i.e. which peak each neighbor is associated with in a cluster.
  // First, we need to sort the "neighbor" TObjArray, because we want to
  // consider the "neighbor" towers from higher towers to lower towers
  neighbor.sort(AscendingTowerEnergySorter());
  // Keep trying to make tower-cluster associations until we make an entire loop
  // through all neighbors without successfully associating anything. Then stop,
  // otherwise we end up in an infinite loop when we can't associate all the
  // neighbors with a cluster (which we usually can't).
  unsigned nAssociations(0);
  do {
    nAssociations = associateTowersWithClusters(neighbor, clust, &valleys);
  } while (nAssociations > 0);
  // Calculate the moments of clusters. We need to do this before calling
  // TowerClusterAssociation::nearestCluster, which uses the cluster moment
  // to determine tower-cluster separations for the valley towers.
  for (Int_t i(0); i < nClusts; ++i) {
    CalClusterMoment(&clust[i]);
  }  // for
  // Ambiguous "valley" towers that were equally spaced between clusters can
  // now be associated
  associateValleyTowersWithClusters(neighbor, clust, &valleys);
  // If there are still towers left in "neighbor", distribute them to clusters
  do {
    nAssociations = associateResidualTowersWithClusters(neighbor, clust);
  } while (nAssociations > 0);
  sortTowersEnergyAscending(clust, nClusts);
  // Recalculate various moment of clusters
  for(Int_t i(0); i < nClusts; ++i) {
    CalClusterMoment(&clust[i]);
  }  // for
  // now add those "zero" towers to the clusters
  // those towers serve the purpose of preventing the creation of bogus peak
  // (peak where there is no energy deposited at the tower)
  TowerList toRemove;
  TowerIter tower;
  for (tower = belowThreshold.begin(); tower != belowThreshold.end(); ++tower) {
    TowerClusterAssociation association(*tower);
    // loop over all clusters
    for(Int_t i(0); i < nClusts; ++i) {
      association.add(&clust[i], kPeakTower);
    }  // for
    HitCluster* cluster = association.nearestCluster();
    if (cluster &&
        association.separation(cluster, kClusterCenter) < maxDistanceFromPeak) {
      (*tower)->cluster = cluster->index;
      toRemove.push_back(*tower);
      cluster->tow->Add(*tower);
    }  // if
  }  // for
  return nClusts;
}

// cluster moment calculation is now a seperate function
void TowerUtil::CalClusterMoment(HitCluster *clust) {
  if(clust) {
    clust->CalClusterMoment(Ecutoff);
  }
  clust->numbTower = clust->tow->GetEntriesFast() ;
}

// 2003-08-28
// a new way of catagorize clusters
// 2003-09-13
// new parameters
// lines of seperation:
// 1. y=1.825(x-6), below this line, (and x>19) is almost certainly 2-photon (1 exception)
// 2. y=1.825(x-0), above this, (and x<25) is almost certainly 1-photon (some exceptions, but probably very weak 2nd photon)
// catagozie cluster (1-photon, 2-photon, or need-to-fit-for-both-and-see)
Int_t TowerUtil::CatagBySigmXY(HitCluster *clust) {
  // only need to be called once "ReadParamters()"
  // if the number of towers in a cluster is less than "minTowerCatag02"
  //    consider the cluster a catag-1 cluster (with only 1 gamma)
  if( clust->numbTower < minTowerCatag02 ) {
    clust->catag = 1 ;
    return 1;
  }
  Float_t sMaxEc = clust->sigmaMax * clust->energy;
  if( clust->energy < cutEcSigma[0][0] * ( sMaxEc - cutEcSigma[0][1] ) ) {
    if( sMaxEc > minEcSigma2Ph ) {
      clust->catag = 2 ;
      return 2;
    } else {
      clust->catag = 0 ;
      return 0;
    }
  } else if( clust->energy > cutEcSigma[1][0] * ( sMaxEc - cutEcSigma[1][1] ) ) {
    if( sMaxEc < maxEcSigma1Ph ) {
      clust->catag = 1 ;
      return 1;
    } else {
      clust->catag = 0 ;
      return 0;
    }
  } else {
    clust->catag = 0 ;
  }
  return clust->catag ;
}
}  // namespace PSUGlobals

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>

#include <algorithm>
#include <list>

#include <TCollection.h>
#include "TText.h"
#include "TCutG.h"
#include "TCanvas.h"

#include "TowerUtil.h"
using namespace std;
using namespace PSUGlobals;

ClassImp(TowerUtil);

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

typedef TowerUtil::TowerList TowerList;
typedef TowerList::iterator TowerIter;
typedef TowerList::reverse_iterator TowerRIter;

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
  while((element = static_cast<Pointer>(next()))) {
    container->push_back(element);
  }  // while
  return container->size();
}

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
 Predicate for testing if a test tower is a neighbour of a reference tower.
 
 If the test tower
 has energy below the global cutoff, always return false, even if it is a neighbour
 of the reference tower.
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

struct TowerEnergyIsAboveThreshold {
  bool operator()(const TowerFPD* tower) const {
    return !(tower->energy < minTowerEnergy);
  }
};

/*
 Filter out towers below the minimum energy threshold.
 
 Returns a pointer list of below-threshold towers.
 Erases the below-threshold towers from the input list.
 The order of towers after filtering is not guaranteed.
 */
TowerList filterTowersBelowEnergyThreshold(TowerList* towers) {
  // Move towers above threshold to the front, those below to the end
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
}  // unnamed namespace

unsigned TowerUtil::associateTowersWithClusters(TowerList& neighbor,
                                                HitCluster *clust,
                                                TObjArray* arrValley) {
  // distance to peak of cluster
  Float_t distToClust[maxNClusters] ;
  TowerList associated;
  for (TowerRIter i = neighbor.rbegin(); i != neighbor.rend(); ++i) {
    TowerFPD* nbT = *i;
    // which cluster should this tower belong?
    Int_t whichCluster(-1);
    // The smallest distance to a peak tower
    Float_t minDist = ExtremelyFaraway ;
    for (Int_t ic=0; ic<nClusts; ic++) {
      // first set the distance to the peak of a cluster to be unreasonably high
      distToClust[ic] = ExtremelyFaraway ;
      // peak-tower is always the first one in cluster's array
      TowerFPD* pkT = static_cast<TowerFPD*>((clust[ic].tow)->First());
      // Make sure the neighbor has lower energy than the peak, but be careful;
      // because of digitization, it is possible that the "neighbor" tower
      // has the exact same energy as the peak tower, not just less.
      if (pkT->energy < nbT->energy) {
        continue;
      }  // if
      // Loop over all towers in this cluster to see if the neighbor tower is
      // adjacent to any of them.
      for(Int_t jt=0; jt<(clust[ic].tow)->GetEntriesFast(); jt++) {
        TowerFPD* towInClust = static_cast<TowerFPD*>((clust[ic].tow)->At(jt));
        // Place an energy selection when determining adjacent towers, as a
        // neighbor cannot exceed an adjacent tower by a factor more than
        // minRatioPeakTower, otherwise it will be considered a peak itself.
        if (nbT->IsNeighbor(towInClust) &&
            nbT->energy < minRatioPeakTower * towInClust->energy) {
          // Calculate distance to the "peak" of the cluster
          // Distance = sqrt(row-separation^2 + column-separation^2)
          distToClust[ic] = sqrt(pow(pkT->col - nbT->col, 2.) +
                                 pow(pkT->row - nbT->row, 2.));
          // Once a tower is found to be a neighbor of any tower in a cluster
          // it is not necessary to contine
          break;
        }  // if
      }  // for loop over all towers in a cluster
    	// Check if the distance to the peak of this cluster is smaller
    	if (distToClust[ic] < minDist) {
  	    minDist = distToClust[ic] ;
	      whichCluster = ic ;
	    }  // if
    } // loop over all clusters
    // Attempt to move the tower to the appropriate cluster
    // If there is only one cluster at the minimal distance, associate it with
    // that cluster.
    // If there are multiple clusters at the same distance we need to make a
    // more sophisticated association. This will be done later, so save the
    // necessary information now.
    if (minDist < ExtremelyFaraway) {
      // Loop over all clusters, and count the number of "peaks" that have the
      // same, minimal, distance to the neighbor. If this is one, we can
      // unambiguously assign the tower to that cluster now. If it is more than
      // one, we save some information and will assign it later. We will refer
      // to these ambiguous towers as "valley" towers here.
      Int_t valleyIndex = arrValley->GetEntriesFast();
      nPeakSameDist[valleyIndex] = 0;  // Number of peaks at minimal dist
      for (Int_t clusterIndex = 0; clusterIndex < nClusts; clusterIndex++) {
        if (distToClust[clusterIndex] == minDist) {
          int peakIndex = nPeakSameDist[valleyIndex];
          // The cluster index for the mth peak at this separation from the
          // nth valley tower
          peaksToValley[valleyIndex][peakIndex] = clusterIndex;
          nPeakSameDist[valleyIndex]++;
        }  // if
      }  // for
      if (nPeakSameDist[valleyIndex] == 1) {
        // Only one "peak" is closest to "nbT". "nbT" belongs to this "peak"!
        nbT->cluster = whichCluster;
        associated.push_back(nbT);
        (clust[whichCluster].tow)->Add(nbT);
      } else if(nPeakSameDist[valleyIndex] > 1) {
        // Add this tower to the "valley" array so we can associate it with a
        // cluster later
        associated.push_back(nbT);
        arrValley->Add(nbT);
      } else {
        cout << "Something wrong in your logic! nPeakSameDist = " << nPeakSameDist << "! Error!" << endl;
      }  // if
    } else {
      cout << "tpbdebug didn't find cluster with dist < ExtremelyFaraway" << endl;
    }  // if (minDist < ExtremelyFaraway)
  } // loop over TObjArray "neighbor"
  // Remove associated neighbors from the neighbor list
  for (TowerIter i = associated.begin(); i != associated.end(); ++i) {
    neighbor.remove(*i);
  }  // for
  return associated.size();
}

unsigned TowerUtil::associateResidualTowersWithClusters(TowerList& neighbor,
                                                        HitCluster *clust) {
  Int_t jjn = neighbor.size() - 1 ;
  // distance to peak of cluster
  Float_t distToClust[maxNClusters] ;
  TowerFPD *nbT;
  TowerFPD *pkT;
  TowerList associated;
  std::cout << "tpbdebug handling " << neighbor.size() << " residual towers" << std::endl;
  for (TowerRIter i = neighbor.rbegin(); i != neighbor.rend(); ++i) {
    nbT = *i;
    Int_t numbTowBefore = neighbor.size();
    // which cluster should this tower belong?
    Int_t whichCluster;
    // the smallest distance to the peaks
    // 2003-09-30
    // now the smallest distance to clusters
    Float_t minDist;
    minDist = ExtremelyFaraway ;
    for(Int_t ic=0; ic<nClusts; ic++) {
	  // first set the distance to the peak of a cluster to be unreasonably high
  	  distToClust[ic] = ExtremelyFaraway ;
      // peak-tower is always the first one in cluster's array
      pkT = (TowerFPD *) (clust[ic].tow)->First();
      // do not consider if the peak is lower than the "neighbor" tower
      if( pkT->energy < nbT->energy )
        continue;
      // loop over all towers in this cluster
      TowerFPD *towInClust;
      for(Int_t jt=0; jt<(clust[ic].tow)->GetEntriesFast(); jt++) {
	      // check if "nbT" is neigboring any tower in this cluster
	      towInClust = (TowerFPD *) (clust[ic].tow)->At(jt) ;
	      if( nbT->IsNeighbor( towInClust ) ) {
          // make sure that "nbT" is not more than "minRatioPeakTower" times of "towInClust"
          if( nbT->energy > (minRatioPeakTower * towInClust->energy) ) 
            continue;
          // calculate distance to the "peak" of the cluster
          Float_t delc, delr;
          CalClusterMoment(&clust[ic]);
          delc = clust[ic].x0 - (nbT->col - 0.5) ;
          delr = clust[ic].y0 - (nbT->row - 0.5) ;
          distToClust[ic] = sqrt( delc * delc + delr * delr ) ;
          // once if the tower is found to be a neighbor of (and does not have higher enough energy than) one tower
          //   in the cluster, it is not necessary to contine
          break;
		    }
	    } // loop over all towers in a cluster
      // check if the distance to the peak of this cluster is smaller
      if( distToClust[ic] < minDist ) {
        minDist = distToClust[ic] ;
        whichCluster = ic ;
      }
    } // loop over all clusters
    // move the tower to the appropriate cluster
    if( minDist < ExtremelyFaraway ) {
  	  nbT->cluster = whichCluster ;
      associated.push_back(nbT);
	    (clust[whichCluster].tow)->Add(nbT);
	    std::cout << "tpbdebug associating residual neighbor " << jjn << " with cluster " << whichCluster << std::endl;
	  }
    // move forward on to the next tower
    jjn-- ;
    if( jjn == -1 ) {
      // Counts number of towers in "neighbor". If the next iteration does not move any tower to a cluster
      //    (same number of towers), we have an infinite loop. Break out and print out error message!
      Int_t numbTowBefore;
      numbTowBefore = neighbor.size();
      jjn = neighbor.size() - 1;
  	  if( numbTowBefore > 0 && (jjn + 1) == numbTowBefore ) {
	      cout << "Infinite loop! The following towers are not claimed by any cluster!" << endl;
	      for (TowerIter i = neighbor.begin(); i != neighbor.end(); ++i) {
	        cout << "\t";
	        (*i)->Print();
	      }  // for
	      cout << "\n  Algorithm could not deal with the case when \"neighbor\" is higher than the ";
	      cout << "closest \"peak\", but towers surrounding it all belongs to that cluster!\n" << endl;
	      cout << "Check the event! No (or bad) pedestal subtraction!?? \n" << endl;
//	      break;
	    }
  	}
  } // loop over TObjArray "neighbor"
  for (TowerIter i = associated.begin(); i != associated.end(); ++i) {
    neighbor.remove(*i);
  }  // for
  return associated.size();
}

TowerUtil::TowerUtil() : nClusts(0) {
  arrValley=NULL;
  neighbor=NULL;
  SetMomentEcutoff();
};

TowerUtil::~TowerUtil() {
  if(arrValley)delete arrValley;
  if(neighbor)delete neighbor;
};

Int_t TowerUtil::FindTowerCluster(TObjArray *inputTow, HitCluster *clust) {
  TowerList arrTow;
  fillStlContainerFromRootCollection(*inputTow, &arrTow);
  // Remove towers below energy threshold, but save them for later use
  TowerList belowThreshold = filterTowersBelowEnergyThreshold(&arrTow);
  // the neighbor TObjArray
  TowerList neighbor;
  // the "valley" TObjArray
  arrValley = new TObjArray(16);
  arrValley->Clear();
  arrTow.sort(DescendingTowerEnergySorter());
  // "TObjArray::Sort()" sorts the array from small to big ( [0]<=[1]<=...<=[48] )
  // need to take care of that
  // have to go last object first!
  // The algoriths is such, first get the highest tower (which is ALWAYS the last one),
  // and it and all its neighbors to the next cluster. Then repeat the process over the
  // remaining towers.
  while(!arrTow.empty() && nClusts < maxNClusters) {
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
  // Each loop of association goes through the entire remaining neighbor list
  // and tries to associate each neighbor with a cluster. If a neighbor is
  // associated with a cluster add it to that cluster and remove it from the
  // neighbor list.
  // Keep doing this until we make an entire loop through all neighbors without
  // successfully associating anything (nSuccessfulAssociations == 0). Then we
  // stop, lest we end up in an infinite loop when we can't associate all the
  // neighbors with a cluster.
  unsigned nSuccessfulAssociations(0);
  do {
    nSuccessfulAssociations =
      associateTowersWithClusters(neighbor, clust, arrValley);
  } while(nSuccessfulAssociations > 0);
  // distance to peak of cluster
  Float_t distToClust[maxNClusters] ;
  TowerFPD *nbT;
  TowerFPD *pkT;
  // All towers in "neighbor" must belong to a cluster.
  // Loop over them, check if it is bordering a cluster.
  //   If yes, move it to the nearest cluster, and move on to the next tower.
  //   If no, move on to the next tower.
  //   When reach the end of the array, start from the beginning.
  Int_t jjn = neighbor.size() - 1;
  // 2003-10-12
  // calculate the moment of clusters, then decide where the "valley" towers should go
  for(Int_t ic=0; ic<nClusts; ic++) {
    CalClusterMoment(&clust[ic]);
  }  
  Int_t numbVal = arrValley->GetEntriesFast();
  std::cout << "tpbdebug handling " << numbVal << " valley towers" << std::endl;
  for(Int_t iVal = 0; iVal<numbVal; iVal++) {
    nbT = (TowerFPD *) arrValley->At(iVal);
    // which cluster should this tower belong?
    Int_t whichCluster = -1;
    Float_t minDist;
    minDist = ExtremelyFaraway ;
    for(Int_t lc=0; lc<nPeakSameDist[iVal]; lc++) {
      // which cluster is this?
      Int_t jkc;
      jkc = peaksToValley[iVal][lc] ;
      Float_t delc, delr;
      delc = clust[jkc].x0 - (nbT->col - 0.5) ;
      delr = clust[jkc].y0 - (nbT->row - 0.5) ;
      distToClust[lc] = sqrt( delc * delc + delr * delr ) ;
      // check if the distance to the "center" of this cluster is smaller
      if( distToClust[lc] < minDist ) {
        minDist = distToClust[lc] ;
        whichCluster = jkc ;
      }
    }
    // move the tower to the appropriate cluster
    if( minDist < ExtremelyFaraway ) {
      nbT->cluster = whichCluster ;
      neighbor.remove(nbT);
      (clust[whichCluster].tow)->Add(nbT);
      std::cout << "tpbdebug associate valley tower " << iVal << " with cluster " << whichCluster << std::endl;
    } else {
      cout << "Something is wrong! The following \"Valley\" tower does not belong to any cluster! Error!" << endl;
      nbT->Print();
      cout << "!!!!!!!!\n" << endl;
    }
  }  // end of for loop over valley towers
  // 2003-10-12
  // If there are still towers left in "neighbor", distribute them to clusters
  do {
    nSuccessfulAssociations =
      associateResidualTowersWithClusters(neighbor, clust);
  } while(nSuccessfulAssociations > 0);
  // 2003-08-26
  // Sort towers by energy (descending, higher energy towers first)
  for(Int_t jc=0; jc<nClusts; jc++) {
    (clust[jc].tow)->UnSort();
    (clust[jc].tow)->Sort();
    // TObjArray sort ascending! Not what I want!
    // Do an exchange of towers: 0<-->nTT-1, 1<-->nTT-2, etc.
    TowerFPD *tmp1;
    TowerFPD *tmp2;
    Int_t nTT;
    nTT = (clust[jc].tow)->GetEntriesFast();
    for(Int_t itt=0; itt<nTT/2; itt++) {
      tmp1 = (TowerFPD *) (clust[jc].tow)->RemoveAt(itt) ; 
      tmp2 = (TowerFPD *) (clust[jc].tow)->RemoveAt(nTT-1-itt) ;
      (clust[jc].tow)->AddAt(tmp1, nTT-1-itt);
      (clust[jc].tow)->AddAt(tmp2, itt);
    }
  }
  // 2003-08-30
  // put center of towers at 0.5 lgd, because this is the more natural way.
  // calculate various moment of clusters
  for(Int_t ic=0; ic<nClusts; ic++) {
    CalClusterMoment(&clust[ic]);
    // set initial nPhoton to 0 & catag to be -1
    clust[ic].nPhoton =  0 ;
    clust[ic].catag   = -1 ;
  }
  // 2003-09-14
  // now add those "zero" towers to the clusters
  // those towers serve the purpose of preventing the creation of bogus peak
  //   (peak where there is no energy deposited at the tower)
  TowerList toRemove;
  for (TowerIter i = belowThreshold.begin(); i != belowThreshold.end(); ++i) {
    nbT = *i;
    // which cluster should this tower belong?
    Int_t whichCluster = 0;
    // the smallest distance to the peaks
    Float_t minDist;
    minDist = ExtremelyFaraway ;
    // loop over all clusters
    for(Int_t ijc=0; ijc<nClusts; ijc++) {
      Float_t dist, delc, delr;
      // peak-tower is always the first one in cluster's array
      pkT = (TowerFPD *) (clust[ijc].tow)->First();
      // distance to this peak
      delc = nbT->col - pkT->col;
      delr = nbT->row - pkT->row;
      // 2003-10-11
      // distance to "center" of cluster is used
      // 			delc = clust[ijc].x0 - (nbT->col - 0.5) ;
      // 			delr = clust[ijc].y0 - (nbT->row - 0.5) ;
      dist = sqrt( delc*delc + delr*delr ) ;
      // since the higher-peak cluster is considered first, when "dist" is the same, favor the
      // higher-peak cluster (naturally)
      if( dist < minDist ) {
        minDist = dist ;
        whichCluster = ijc;
      }
    }
    // if the distance is smaller than "maxDistanceFromPeak"
    //    move the tower to the appropriate cluster
    // Do not want to add too many "zero" towers to a cluster!
    if( minDist < maxDistanceFromPeak ) {
      nbT->cluster = whichCluster ;
      toRemove.push_back(nbT);
      (clust[whichCluster].tow)->Add(nbT);
    }
  }
  for (TowerIter i = toRemove.begin(); i != toRemove.end(); ++i) {
    arrTow.remove(*i);
  }  // for
  neighbor.clear();
  arrValley->Clear();
  delete arrValley;
  arrValley=NULL;
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

void TowerUtil::PrintTowers(const TowerFPD *tows) {
  for(Int_t j=0; j<nNSTow; j++) {
    if( j%7 == 0 ) printf("\n");
    printf("%3d", tows[j].cluster+1);
  }
  for(Int_t j=0; j<nNSTow; j++) {
    if( j%7 == 0 ) printf("\n");
    printf("%10.4f", tows[j].energy);
  }
  printf("\n\n");
}

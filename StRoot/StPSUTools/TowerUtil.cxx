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

typedef TowerUtil::TowerList TowerList;
typedef TowerList::iterator TowerIter;

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
}  // unnamed namespace

unsigned TowerUtil::associateTowersWithClusters(TowerList& neighbor,
                                                HitCluster *clust,
                                                TObjArray* arrValley) {
  Int_t jjn = neighbor.size() - 1 ;
  const Float_t ExtremelyFaraway = 99999 ;
  // distance to peak of cluster
  Float_t distToClust[maxNClusters] ;
  TowerFPD *nbT;
  TowerFPD *pkT;
  TowerList associated;
  typedef TowerList::reverse_iterator TowerRIter;
  for (TowerRIter i = neighbor.rbegin(); i != neighbor.rend(); ++i) {
    nbT = *i;
    Int_t numbTowBefore = neighbor.size();
    // towers in "neighbor" should NEVER be lower than "minTowerEnergy"
    if( nbT->energy < minTowerEnergy ) {
      cout << "Something is wrong! A tower in \"neighbor\" has energy " << nbT->energy;
      cout << ". Lower than " << minTowerEnergy << ".\n" << endl;
    }
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
      // 2003-08-15
      // Because of digitization, it is possible that the "neighbor" tower
      // has the exact energy as the peak tower!!!
      // Change the condition from "if( pkT->energy <= nbT->energy )"
      // to "if( pkT->energy < nbT->energy )", so that this "neighbor" tower
      // would not be hung dry!
      // 2003-09-27
      // On rare occasions, this requirement that a "neighbor" must not be higher than
      //    the "peak" could leave the "neighbor" hang on dry (infinite loop). It happens
      //    (run 4126039, EN, for example. Pedestals subtraction not properly done?!) when
      //    "neighbor" is higher than the closest "peak", but towers surrounding it all belongs
      //    to that cluster.
      // The requirement is reasonable. So we have to find a way to break out of the infinite
      //    loop and print out an error message. (Counts the number of towers remaining in
      //    "neighbor". If it is the same as the last iteration, we have an infinite loop!)
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
          // 2003-10-03
          // Revert to using distance to "peak" tower.
          // 2003-09-30
          // Calculate cluster moment (only need center position) while adding
          //    "neighbor" towers to clusters. Use the center of cluster to
          //    ecide the distance from a "neighbor" tower to a cluster (previously
          //    use the "peak" tower position). This way, a tower just in between
          //    of two "peaks" will probably go to the cluster of lower "peak".
          //    Thus, the fitting program will less likely to try to find a 2nd
          //    peak in the direction of another (lower) cluster.
          // 					CalClusterMoment(&clust[ic]);
          //  					delc = clust[ic].x0 - (nbT->col - 0.5) ;
          // 					delr = clust[ic].y0 - (nbT->row - 0.5) ;
          delc = pkT->col - nbT->col ;
          delr = pkT->row - nbT->row ;
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
      // loop over all clusters, and count the number of "peaks" that have the same distance to "nbT"
      Int_t numbValleyTower;
      numbValleyTower = arrValley->GetEntriesFast();
      nPeakSameDist[numbValleyTower] = 0 ;
      for(Int_t llc=0; llc<nClusts; llc++) {
        if( distToClust[llc] == minDist ) {
          int npeaks = nPeakSameDist[numbValleyTower];
          peaksToValley[numbValleyTower][npeaks] = llc ;
          nPeakSameDist[numbValleyTower] ++ ;
        }
      }
      if( nPeakSameDist[numbValleyTower] == 1 ) {
        // Only one "peak" is closest to "nbT". "nbT" belongs to this "peak"!
        nbT->cluster = whichCluster ;
        associated.push_back(nbT);
        (clust[whichCluster].tow)->Add(nbT);
        std::cout << "tpbdebug associate neighbour " << jjn << " with cluster " << whichCluster << std::endl;
      }
      else if( nPeakSameDist[numbValleyTower] > 1 ) {
        associated.push_back(nbT);
        arrValley->Add(nbT);
        std::cout << "tpbdebug add neighbour " << jjn << " to valley array" << std::endl;
      }
      else {
        cout << "Something wrong in your logic! nPeakSameDist = " << nPeakSameDist << "! Error!" << endl;
      }
    } else {
      cout << "tpbdebug didn't find cluster with dist < ExtremelyFaraway" << endl;
    }
    // move forward on to the next tower
    jjn-- ;
    if( jjn == -1 ) {  // Means we've been through the entire list again
      // Counts number of towers in "neighbor". If the next iteration does not move any tower to a cluster
      //    (same number of towers), we have an infinite loop. Break out and print out error message!
      jjn = neighbor.size() - 1;
      if( numbTowBefore > 0 && (jjn + 1) == numbTowBefore ) {
        std::cout << "tpbdebug <deprecated> breaking out of neighbour association with " << neighbor.size() << " remaining" << std::endl;
//  	    break;
	    }  // if
    }
  } // loop over TObjArray "neighbor"
  for (TowerIter i = associated.begin(); i != associated.end(); ++i) {
    neighbor.remove(*i);
  }  // for
  return associated.size();
}

unsigned TowerUtil::associateResidualTowersWithClusters(TowerList& neighbor,
                                                HitCluster *clust,
                                                TObjArray* arrValley) {
  Int_t jjn = neighbor.size() - 1 ;
  const Float_t ExtremelyFaraway = 99999 ;
  // distance to peak of cluster
  Float_t distToClust[maxNClusters] ;
  TowerFPD *nbT;
  TowerFPD *pkT;
  TowerList associated;
  std::cout << "tpbdebug handling " << neighbor.size() << " residual towers" << std::endl;
  typedef TowerList::reverse_iterator TowerRIter;
  for (TowerRIter i = neighbor.rbegin(); i != neighbor.rend(); ++i) {
    nbT = *i;
    Int_t numbTowBefore = neighbor.size();
    // towers in "neighbor" should NEVER be lower than "minTowerEnergy"
    if (nbT->energy < minTowerEnergy) {
      cout << "Something is wrong! A tower in \"neighbor\" has energy " << nbT->energy;
      cout << ". Lower than " << minTowerEnergy << ".\n" << endl;
    }  // if
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
  // the neighbor TObjArray
  TowerList neighbor;
  // the "valley" TObjArray
  arrValley = new TObjArray(16);
  arrValley->Clear();
  // We assume that "TObjArray *inputTow" is already sorted, however, to be safe,
  // sort tower by energy, if not done already
  arrTow.sort(AscendingTowerEnergySorter());
  TowerFPD *high;
  // "TObjArray::Sort()" sorts the array from small to big ( [0]<=[1]<=...<=[48] )
  // need to take care of that
  // have to go last object first!
  // The algoriths is such, first get the highest tower (which is ALWAYS the last one),
  // and it and all its neighbors to the next cluster. Then repeat the process over the
  // remaining towers.
  while(!arrTow.empty()) {
    // By design, this tower is the highest tower in "arrTow", but it could be lower
    // than a tower in "neighbor"
    TowerFPD* high = arrTow.back();
    // when tower energy is less than minTowerEnergy, break out the loop!
    if(high->energy < minTowerEnergy) {
      break;
    }  // if
		// 2003-08-15
		// Fix a logical loop hole in deciding if a tower is
		//    a peak; Need to first compare the highest tower
		//    remained in "arrTow" to all towers in "neighbor" which is its neighbor,
		//    and if it is lower than any of those, it is
		//    a neighbor. Move it to "neighbor" and continue to
		//    the next tower in "arrTow".
    TowerFPD *resTow;
    Bool_t isPeak = couldBePeakTower(high, &neighbor);
    // if "high" is not a peak, move it to "neighbor"
    if (!isPeak) {
      std::cout << "tpbdebug At arrTow size " << arrTow.size() << " got is not peak" << std::endl;
      arrTow.remove(high);
      neighbor.push_back(high);
    } else {
      std::cout << "tpbdebug At arrTow size " << arrTow.size() << " got IS peak" << std::endl;
      // else if "high" is a peak
      // remove the high tower from the original TObjArray, and add it to the next cluster
      arrTow.remove(high);
      std::cout << "tpbdebug Adding tower with E = " << high->energy << " to cluster " << nClusts << std::endl;
      high->cluster = nClusts;
      (clust[nClusts].tow)->Add(high);
      // loop over rest of original TObjArray, and move any towers neighboring "high"
      // to "neighbor"
      // Find the first tower above the minimum energy threshold, as we will want
      // to exclude these in many cases (we may want to create two separate lists,
      // with above- and below-threshold towers, but let's try to mimic Steve's
      // implementation for now).
      TowerIter aboveThreshold = std::find_if(arrTow.begin(), arrTow.end(),
                                              TowerEnergyIsAboveThreshold());
      // Partition the remaining hits so that neighbours of the high tower are
      // placed at the end, and non-neighbours placed at the start. Use
      // stable_partition so we don't alter the energy ordering.
      TowerIter newEnd = std::stable_partition(aboveThreshold, arrTow.end(),
          std::not1(std::bind2nd(TowerIsNeighbor(), high)));
      // Now copy the neighbors to the neighbor list, and erase them from the
      // master tower list
      unsigned initialSize = arrTow.size();
      neighbor.insert(neighbor.end(), newEnd, arrTow.end());
      arrTow.erase(newEnd, arrTow.end());
      arrTow.sort(AscendingTowerEnergySorter());
      std::cout << "tpbdebug removed " << initialSize - arrTow.size() << " neighbours" << std::endl;
    } // when "high" is a "peak"
    // 2003-09-27
    // Do the follow, no matter "high" is a "peak" or "neighbor"!
    // To close the logic loop hole that a tower which is seperated (by towers of the same energy) from
    //   the "neighbor" towers becomes a peak, just because it happens to be in front of those towers
    //   of the same energy (since the sorting can not distinguish that).
    // We need to again loop over "arrTow", move any tower (that is neighboring any of the "neighbor"
    //   towers and also has lower (or equal) energy than that "neighbor" tower) to "neighbor"
    // Every time we remove a tower from "arrTow" (except we just simply go over all items in TObjArray
    //   sequentially without worrying the relative order), we need to compress "arrTow"!
    //   Since we assume that it has no emty slot!
    TowerIter aboveThreshold = std::find_if(arrTow.begin(), arrTow.end(),
                                            TowerEnergyIsAboveThreshold());
    unsigned initialSize = arrTow.size();
    TowerList above;
    std::copy(aboveThreshold, arrTow.end(), std::back_inserter(above));
    above.reverse();
    for (TowerIter i = above.begin(); i != above.end(); ++i) {
      if (!couldBePeakTower(*i, &neighbor)) {
        neighbor.push_back(*i);
        arrTow.remove(*i);
      }  // if
    }  // for
    std::cout << "tpbdebug removed " << initialSize - arrTow.size() <<
      " adjacent towers that could not be peak towers" << std::endl;
    // increment "nClusts" when we find a "peak"
    if (isPeak) {
      nClusts++ ;
      // when reaching maximum number of clusters, break out the loop!
      if (nClusts >= maxNClusters) {
        break;
      }  // if
    }  // if
  }  // End of for loop over "arrTow"
  // now that we know all the peaks, let's decide the affiliation of
  // those remote neighbor-towers in "neighbor" TObjArray
  // First, we need to sort the "neighbor" TObjArray, because we want to consider the "neighbor" towers
  //    from higher towers to lower towers
  neighbor.sort(AscendingTowerEnergySorter());
  std::cout << "tpbdebug " << neighbor.size() << " neighbours need to be affiliated with clusters" << std::endl;
  // extremely faraway distance (no distance between towers can be this large!)
  const Float_t ExtremelyFaraway = 99999 ;
  // distance to peak of cluster
  Float_t distToClust[maxNClusters] ;
  TowerFPD *nbT;
  TowerFPD *pkT;
  // Each loop of association goes through the entire remaining neighbor list
  // and tries to associate each neighbor with a cluster.
  // Keep doing this until we make an entire loop through all neighbors without
  // successfully associating anything (nSuccessfulAssociations == 0). Then we stop, lest
  // we end up in an infinite loop when we can't associate all neighbors.
  unsigned nSuccessfulAssociations(0);
  do {
    nSuccessfulAssociations =
      associateTowersWithClusters(neighbor, clust, arrValley);
  } while(nSuccessfulAssociations > 0);
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
      associateResidualTowersWithClusters(neighbor, clust, arrValley);
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
  for (TowerIter i = arrTow.begin(); i != arrTow.end(); ++i) {
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

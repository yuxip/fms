// $Id$
//
// $Log$
/**
 \file      StFmsEventClusterer.cxx
 \brief     Implementation of StFmsEventClusterer, manager class for clustering
 \author    Steven Heppelmann <steveheppelmann@gmail.com>
 \author    Yuxi Pan <yuxipan@physics.ucla.edu>
 \author    Thomas Burton <tpb@bnl.gov>
 \date      2014
 \copyright Brookhaven National Lab
 */
#include "StFmsPointMaker/StFmsEventClusterer.h"

#include <algorithm>
#include <functional>
#include <iterator>
#include <list>
#include <numeric>

#include "TF2.h"  // To use shower-shape function
#include "TMath.h"
#include "TRandom.h"  // For ROOT global random generator, gRandom

#include "StRoot/St_base/StMessMgr.h"
#include "StEvent/StFmsCluster.h"
#include "StEvent/StFmsHit.h"

#include "StFmsPointMaker/StFmsClusterFitter.h"
#include "StFmsPointMaker/StFmsFittedPhoton.h"
#include "StFmsPointMaker/StFmsGeometry.h"
#include "StFmsPointMaker/StFmsTower.h"
#include "StFmsPointMaker/StFmsTowerCluster.h"

using namespace std;
using namespace FMSCluster;

namespace {
// We use the tower list defined in StFmsTowerCluster throughout this file.
// Define some typedefs for convenience.
typedef StFmsTowerCluster::Towers Towers;

/* Helper function to add numbers of photons using std::accumulate */
int accumulatePhotons(int nPhotons, const ClusterList::value_type& cluster) {
  return nPhotons + cluster->cluster()->nPhotons();
}

/* Sum the total number of photons in a list of clusters */
template<class Iterator>
int sumPhotonsOverClusters(Iterator start, Iterator end) {
  return std::accumulate(start, end, 0, accumulatePhotons);
}

/* Sum the total number of photons in a list of clusters */
template<class Container>
int sumPhotonsOverClusters(const Container& clusters) {
  return sumPhotonsOverClusters(clusters.begin(), clusters.end());
}

/* Unary predicate for selecting bad clusters. */
struct IsBadCluster
    : public std::unary_function<const ClusterList::value_type&, bool> {
  // Set minimum allowed cluster energy and maximum number of towers
  IsBadCluster(double minEnergy, unsigned maxTowers)
      : energy(minEnergy), towers(maxTowers) { }
  bool operator()(const ClusterList::value_type& cluster) const {
    return cluster->cluster()->energy() <= energy ||
           cluster->towers().size() > towers;
  }
  double energy;
  unsigned towers;
};

/*
 Returns a pointer to the lowest energy photon in a cluster
 
 Assumes the cluster is either 1- or 2-photon
 Returns nullptr if there is no photon in the cluster
 */
const StFmsFittedPhoton* findLowestEnergyPhoton(
    const StFmsTowerCluster* cluster) {
  const StFmsFittedPhoton* photon = nullptr;
  switch (cluster->cluster()->nPhotons()) {
    case 1:
      photon = &(cluster->photons()[0]);
      break;
    case 2:
      if (cluster->photons()[0].energy < cluster->photons()[1].energy) {
        photon = &(cluster->photons()[0]);
      } else {
        photon = &(cluster->photons()[1]);
      }  // if
    default:
      break;  // photon remains nullptr
  }  // switch
  return photon;
}

/* Functor returning true if a tower matches (row, column) */
struct HasRowColumn {
  int row, column;
  HasRowColumn(int r, int c) : row(r), column(c) { }
  bool operator()(const StFmsTower* tower) const {
    return tower->row() == row && tower->column() == column;
  }
};

/*
 Search towers in a cluster for one matching a row and column number

 Return a pointer to the matching tower if one is found, nullptr otherwise.
 */
const StFmsTower* searchClusterTowers(int row, int column,
                                      const StFmsTowerCluster& cluster) {
  auto found = std::find_if(cluster.towers().begin(), cluster.towers().end(),
                            HasRowColumn(row, column));
  if (found != cluster.towers().end()) {
    return *found;
  }  // if
  return nullptr;
}

/*
 Yields fit parameter start points and limits for 2-photon fit.

 See SFmsClusterFitter::fit2PhotonCluster() for parameter meanings
 */
void setup2PhotonFitParameters(std::vector<double>& start,
                               std::vector<double>& steps,
                               std::vector<double>& lower,
                               std::vector<double>& upper,
                               double x, double y,
                               const StFmsTowerCluster* towerCluster) {
  const auto cluster = towerCluster->cluster();
  start = {
    2,
    x * cluster->x(),
    y * cluster->y(),
    2.2 * x * cluster->sigmaMax(),
    towerCluster->thetaAxis(),
    gRandom->Uniform(-0.1, 0.1),
    cluster->energy(),
  };
  steps = {0, 0.02, 0.02, 0.01, 0.01, 0.01, 0.1};
  const double sigmaMaxE = cluster->sigmaMax() * cluster->energy();
  double maxTheta = cluster->sigmaMin() / cluster->sigmaMax() / 2.8;
  maxTheta = std::min(maxTheta, TMath::PiOver2());
  lower = {
    1.5,
    start.at(1) - 0.2 * x,
    start.at(2) - 0.2 * y,
    std::max(18. / pow(sigmaMaxE, 0.8), 0.5) * x,
    start.at(4) - maxTheta,
    -1.,
    start.at(6) * 0.95
  };
  upper = {
    2.5,
    start.at(1) + 0.2 * x,
    start.at(2) + 0.2 * y,
    std::min(0.085 * (60. - sigmaMaxE), 3.5) * x,
    start.at(4) + maxTheta,
    1.,
    start.at(6) * 1.05
  };
  // With the above approach the limits on parameter 3 can sometimes go beyond
  // sensible values, so limit them.
  lower.at(3) = std::min(lower.at(3), start.at(3) * 0.9);
  upper.at(3) = std::max(upper.at(3), start.at(3) * 1.1);
}
}  // unnamed namespace

StFmsEventClusterer::StFmsEventClusterer(const StFmsGeometry* geometry,
                                         Int_t detectorId)
    : mClusterFinder(0.5), mGeometry(geometry), mDetectorId(detectorId) { }

StFmsEventClusterer::~StFmsEventClusterer() {
  if (mFitter) {
    delete mFitter;
  }  // if
}

Bool_t StFmsEventClusterer::cluster(std::vector<StFmsTower>* towerList) {
  mTowers = towerList;
  mTowerWidthXY = mGeometry->towerWidths(mDetectorId);
  /** \todo Test of number of towers should be detector-dependent */
  if (mTowers->size() > 578) {
    LOG_ERROR << "Too many towers for Fit" << endm;
    return false;
  }  // if
  mFitter = new StFmsClusterFitter(mGeometry, mDetectorId);
  return fitEvent();  // Return true for success
}

Int_t StFmsEventClusterer::fitEvent() {
  // Possible alternative clusters for 1-photon fit: for catagory 0
  StFmsClusterFinder::TowerList towerList;
  for (auto i = mTowers->begin(); i != mTowers->end(); ++i) {
    towerList.push_back(&(*i));
  }  // for
  mClusterFinder.findClusters(&towerList, &mClusters);
  // Cluster energy should be at least 2 GeV (parameter "minRealClusterEne")
  /** \todo Use detector ID enums here */
  if (mDetectorId == 8 || mDetectorId == 9) {
    mClusters.remove_if(IsBadCluster(0.75, 25));
  } else {  // Different cuts for small cell
    mClusters.remove_if(IsBadCluster(2.0, 49));
  }  // if
  // Must do moment analysis before catagorization
  for (auto i = mClusters.begin(); i != mClusters.end(); ++i) {
    (*i)->findClusterAxis(mClusterFinder.momentEnergyCutoff());
  }  // for
  // Loop over clusters, catagorize, guess the photon locations for cat 0 or 2
  // clusters then fit, compare, and choose the best fit
  bool badEvent = false;
  const double max2PhotonFitChi2 = 10.;
  for (auto cluster = mClusters.begin(); cluster != mClusters.end();
       ++cluster) {
    Int_t clustCatag = mClusterFinder.categorise(cluster->get());
    // point to the real TObjArray that contains the towers to be fitted
    // it is the same tower array for the cluster or all alternative clusters
    mFitter->setTowers(&(*cluster)->towers());
    // Number of Degree of Freedom for the fit
    if (clustCatag == k1PhotonCluster) {
      // Do 1-photon fit
      fitOnePhoton(cluster->get());
    } else if (clustCatag == k2PhotonCluster) {
      // Do 2-photon fit
      fit2PhotonClust(cluster);
      badEvent = (*cluster)->chiSquare() > max2PhotonFitChi2;
    } else if (clustCatag == kAmbiguousCluster) {
      // for catagory-0 cluster, first try 1-photon fit!
      // If the fit is good enough, it is 1-photon. Else also
      // try 2-photon fit, and find the best fit (including 1-photon fit).
      Bool_t is2Photon = true;
      double chiSq1 = fitOnePhoton(cluster->get());
      const StFmsFittedPhoton photon = (*cluster)->photons()[0];  // Cache
      double chiSq2(NAN);  // Only set if do 2-photon fit
      // Decide if this 1-photon fit is good enough
      if (chiSq1 < 5.) {
        is2Photon = false;
      } else {
        // The 1-photon fit isn't good enough, so try 2-photon fit
        chiSq2 = fit2PhotonClust(cluster);
        // Check that the 2-photon fit didn't result in a bogus 2nd photon
        if (validate2ndPhoton(cluster)) {
          // If the 2nd photon in the 2-photon cluster is real, select 1- or 2-
          // photon based on fit chi2/ndf.
          is2Photon = chiSq2 <= chiSq1;
        } else {
          // The 2nd photon was bogus, it must really be a 1-photon cluster,
          // even though the chi2 is poor
          is2Photon = false;
        }  // if (validate2ndPhoton...)
      }  // if (chiSq1...)
      // Now fill in the fit result, either 1- or 2-photon
      if (is2Photon) {
        // 2-photon fit is better
        (*cluster)->cluster()->setNPhotons(2);
        (*cluster)->setChiSquare(chiSq2);
        // Flag the event as bad if the fit chi2/ndf is too bad
        badEvent = (*cluster)->chiSquare() > max2PhotonFitChi2;
      } else {
        // 1-photon fit is better
        (*cluster)->cluster()->setNPhotons(1);
        (*cluster)->setChiSquare(chiSq1);
        (*cluster)->photons()[0] = photon;
      }  // if (is2Photon)
    } else {  // Invalid cluster category
      // should not happen!
      LOG_ERROR << "Your logic of catagory is wrong! Something impossible " <<
        "happens! This a catagory-" << clustCatag <<
        " clusters! Don't know how to fit it!" << endm;
    }  // if (clustCatag...)
  }  // Loop over all real clusters
  const int nPh = sumPhotonsOverClusters(mClusters);
  if (nPh > StFmsClusterFitter::maxNFittedPhotons()) {
    // myFitter can only do up to "maxNFittedPhotons()"-photon fit
    LOG_WARN << "Can not fit " << nPh << " (more than " <<
      StFmsClusterFitter::maxNFittedPhotons() << " photons!" << endm;
    return nPh;
  }  // if
  // For global fit, add all towers from all clusters
  Towers allTow;
  for (auto cluster = mClusters.begin(); cluster != mClusters.end();
       ++cluster) {
    allTow.insert(allTow.end(), (*cluster)->towers().begin(),
                  (*cluster)->towers().end());
  }  // for
  mFitter->setTowers(&allTow);
  // Only do global fit for 2 or more clusters (2-photon fit for one cluster
  // already has global fit)
  if (mClusters.size() > 1) {
    globalFit(nPh, mClusters.size(), mClusters.begin());
    // Check for errors in the global fit - the number of photons returned by
    // the global fit should equal the sum of photons in the fitted clusters
    const int iph = sumPhotonsOverClusters(mClusters);
    if (iph != nPh) {
      LOG_ERROR << "total nPh=" << nPh << " iPh=" << iph << endm;
    }  // if
  }  // if (mClusters.size() > 1)
  return !badEvent;
}

Double_t StFmsEventClusterer::photonEnergyInCluster(
    const StFmsTowerCluster* cluster,
    const StFmsFittedPhoton* photon) const {
  double energy = 0.;
  const Towers& towers = cluster->towers();
  for (auto tower = towers.begin(); tower != towers.end(); ++tower) {
    energy += photonEnergyInTower(*tower, photon);
  }  // for
  return energy;
}

Double_t StFmsEventClusterer::photonEnergyInTower(
    const StFmsTower* tower,
    const StFmsFittedPhoton* photon) const {
  double x = (tower->column() - 0.5) * mTowerWidthXY.at(0) - photon->xPos;
  double y = (tower->row() - 0.5) * mTowerWidthXY.at(1) - photon->yPos;
  return photon->energy * mFitter->showerShapeFunction()->Eval(x, y);
}

/* 1-photon fitting function */
Float_t StFmsEventClusterer::fitOnePhoton(StFmsTowerCluster* towerCluster) {
  auto cluster = towerCluster->cluster();
  // 4 parameters: nPhotons, cluster x, cluster y  and cluster energy.
  const std::vector<double> start = {
    1.0,
    mTowerWidthXY.at(0) * cluster->x(),
    mTowerWidthXY.at(1) * cluster->y(),
    cluster->energy()
  };
  const std::vector<double> delta = {
    0.5,
    mTowerWidthXY.at(0) * 0.5,
    mTowerWidthXY.at(1) * 0.5,
    cluster->energy() * 0.15
  };
  std::vector<double> lower, upper;
  for (unsigned i(0); i < start.size(); ++i) {
    lower.push_back(start.at(i) - delta.at(i));
    upper.push_back(start.at(i) + delta.at(i));
  }  // for
  PhotonList photons;
  Double_t chiSquare = mFitter->fit(start, std::vector<double>(),
                                    lower, upper, &photons);
  if (photons.empty()) {  // check return status in case of a bad fit
    LOG_ERROR << "1-photon Minuit fit found no photons" << endm;
  } else {
    towerCluster->photons()[0] = photons.back();
  }  // if
  cluster->setNPhotons(photons.size());
  const int nDegreesOfFreedom =
    std::max(int(towerCluster->towers().size()) - 3, 1);
  towerCluster->setChiSquare(chiSquare / nDegreesOfFreedom);
  return towerCluster->chiSquare();
}

/* Global fitting function, fitting photons across all clusters */
Float_t StFmsEventClusterer::globalFit(unsigned nPhotons,
                                       const unsigned nClusters,
                                       ClusterIter first) {
  ClusterIter end = first;
  std::advance(end, nClusters);  // Marks end point for cluster iteration
  const unsigned totalPhotons = sumPhotonsOverClusters(first, end);
  if (totalPhotons != nPhotons) {
    LOG_WARN << "Global fit called for " << nPhotons << " but found " <<
      totalPhotons << "... will proceed with " << totalPhotons << endm;
    nPhotons = totalPhotons;
  }  // if
  if (int(nPhotons) > StFmsClusterFitter::maxNFittedPhotons() || nPhotons < 2) {
    LOG_ERROR << "Global fit cannot fit " << nPhotons << " photons" << endm;
    return -9999;
  }  // if
  // Fit has 1 parameter for the number of photons plus 3 per photon (x, y, E)
  std::vector<double> start(1, nPhotons);
  std::vector<double> lower(1, 0.5);
  std::vector<double> upper(1, StFmsClusterFitter::maxNFittedPhotons() + 0.5);
  for (ClusterIter cluster = first; cluster != end; ++cluster) {
    for (int i = 0; i < (*cluster)->cluster()->nPhotons(); i++) {
      start.push_back((*cluster)->photons()[i].xPos);
      lower.push_back(start.back() - 1.25);
      upper.push_back(start.back() + 1.25);
      start.push_back((*cluster)->photons()[i].yPos);
      lower.push_back(start.back() - 1.25);
      upper.push_back(start.back() + 1.25);
      start.push_back((*cluster)->photons()[i].energy);
      lower.push_back(start.back() * (1 - 0.3));  // Limit to +/- 30% energy
      upper.push_back(start.back() * (1 + 0.3));
    }  // for
  }  // for
  PhotonList photons;
  Double_t chiSquare = mFitter->fit(start, std::vector<double>(),
                                    lower, upper, &photons);
  if (photons.size() == nPhotons) {
    // Put the fit result back in the clusters
    PhotonList::const_iterator photon = photons.begin();
    for (ClusterIter cluster = first; cluster != end; ++cluster) {
      for (int i = 0; i < (*cluster)->cluster()->nPhotons(); ++i, ++photon) {
        (*cluster)->photons()[i] = *photon;
      }  // for loop over photons
    }  // for loop over clusters
  } else {
    LOG_WARN << "Global Minuit fit found " << photons.size() <<
      " photons but expected " << nPhotons << endm;
  }  // if
  return chiSquare;
}

/* 2-photon fitting function */
Float_t StFmsEventClusterer::fit2PhotonClust(ClusterIter towerCluster) {
  std::vector<double> start, steps, lower, upper;
  setup2PhotonFitParameters(start, steps, lower, upper,
                            mTowerWidthXY.at(0), mTowerWidthXY.at(1),
                            towerCluster->get());
  PhotonList photons;
  Double_t chiSquare = mFitter->fit2PhotonCluster(start, steps, lower, upper,
                                                  &photons);
  if (photons.empty()) {
    LOG_WARN << "Minuit fit returns error!" << endm;
  }  // if
  // Do a global fit, using result of 1st fit as starting point
  (*towerCluster)->photons()[0] = photons.front();
  (*towerCluster)->photons()[1] = photons.back();
  (*towerCluster)->cluster()->setNPhotons(photons.size());
  chiSquare = globalFit(2, 1, towerCluster);
  int nDegreesOfFreedom = (*towerCluster)->towers().size() - 6;
  if (nDegreesOfFreedom < 1) {
    nDegreesOfFreedom = 1;
  }  // if
  (*towerCluster)->setChiSquare(chiSquare / nDegreesOfFreedom);
  return (*towerCluster)->chiSquare();
}

/*
 Further information:
 If one photon peak lies on top of a low (compared to photon energy) or even
 zero tower, this photon is definitely bogus. This could happen if a nearby
 cluster took away towers that might make the fit have a large Chi-square.
 So first check that the fitted photon of the lower-energy photon (we assume the
 higher energy photon should be fine) is over one of the non-zero towers of the
 cluster. First of all, this ensures that we don't have an "outside of cluster"
 bogus photon, i.e. a bogus photon that could be the result of minimizing the
 chi-square over towers that do not include the supposed peak tower. 
*/
bool StFmsEventClusterer::validate2ndPhoton(ClusterConstIter cluster) const {
  // Select the lower-energy of the two photons
  const StFmsFittedPhoton* photon = findLowestEnergyPhoton(cluster->get());
  // Tower row and column where the fitted photon of lower energy should hit
  int column = 1 + (Int_t)(photon->xPos / mTowerWidthXY.at(0));
  int row = 1 + (Int_t)(photon->yPos / mTowerWidthXY.at(1));
  // Now check whether this tower is one of the non-zero towers of the cluster
  // The temporary StFmsTower only needs row and column set for the test
  const StFmsTower* tower = searchClusterTowers(row, column, **cluster);
  // If tower is non-nullptr, the photon does hit in a tower in this cluster.
  if (!tower) {
    return false;
  }  // if
  // Now test the photon and tower properties.
  // Check if the fitted energy is too large compared to the energy of the tower
  if (tower->hit()->energy() < 0.25 * photon->energy) {
    return false;
  }  // if
  // Check if the 2nd photon's "High-Tower" enery is too large compared to its
  // fitted energy. If so, it is probably splitting one photon into two
  Double_t eSS = photonEnergyInTower(tower, photon);
  if (tower->hit()->energy() > 1.5 * eSS) {
    return false;
  }  // if
  // Check that the 2nd photon is not near the edge of another cluster
  // Namely, we check what would be the energy deposited in other clusters by
  // this photon vs. energy deposited in its own cluster
  // If the ratio is too high, this fitted photon is probably a bogus one
  Double_t energyInOwnCluster = photonEnergyInCluster(cluster->get(), photon);
  // Loop over all clusters except its own
  for (ClusterConstIter i = mClusters.begin(); i != mClusters.end(); ++i) {
    if (i != cluster) {  // Skip the photon's own cluster
      if (photonEnergyInCluster(i->get(), photon) >
          (0.2 * energyInOwnCluster)) {
        return false;  // Stop as soon as we fail for one cluster
      }  // if
    }  // if
  }  // for
  return true;  // The photon passed all tests; it's real
}

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
 Gives fit parameter start points and limits for 1-photon fit.

 See SFmsClusterFitter::fit() for parameter meanings
 */
struct OnePhotonFitParameters {
  std::vector<double> start, lower, upper;
  OnePhotonFitParameters(const std::vector<float>& xyWidth,
                         const StFmsCluster* cluster) {
    const double x = xyWidth.at(0);
    const double y = xyWidth.at(1);
    start = {
      1.0,
      x * cluster->x(),
      y * cluster->y(),
      cluster->energy()
    };
    const std::vector<double> delta = {
      0.5,
      x * 0.5,
      y * 0.5,
      cluster->energy() * 0.15
    };
    for (unsigned i(0); i < start.size(); ++i) {
      lower.push_back(start.at(i) - delta.at(i));
      upper.push_back(start.at(i) + delta.at(i));
    }  // for
  }
};

/*
 Gives fit parameter start points and limits for 2-photon fit.

 See SFmsClusterFitter::fit2PhotonCluster() for parameter meanings
 */
struct TwoPhotonFitParameters {
  std::vector<double> start, steps, lower, upper;
  TwoPhotonFitParameters(const std::vector<float>& xyWidth,
                         const StFmsTowerCluster* towerCluster) {
    const double x = xyWidth.at(0);
    const double y = xyWidth.at(1);
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
};

typedef ClusterList::iterator ClusterIter;
/* Gives fit parameters for global photon fit */
struct GlobalPhotonFitParameters {
  std::vector<double> start, lower, upper;
  GlobalPhotonFitParameters(unsigned nPhotons,
                            ClusterIter first, ClusterIter end)
    // Initialise N-photons parameters as the first element
    : start(1, nPhotons), lower(1, 0.5),
      upper(1, StFmsClusterFitter::maxNFittedPhotons() + 0.5) {
    // Append (x, y, E) fit parameters for each photon
    for (auto cluster = first; cluster != end; ++cluster) {
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
  }
};
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
  if (findClusters()) {
    return fitClusters() && refitClusters();
  }  // if
  return false;
}

Int_t StFmsEventClusterer::findClusters() {
  StFmsClusterFinder::TowerList towerList;
  for (auto i = mTowers->begin(); i != mTowers->end(); ++i) {
    towerList.push_back(&(*i));
  }  // for
  mClusterFinder.findClusters(&towerList, &mClusters);
  switch (mDetectorId) {
    case kFmsNorthLarge:  // Deliberate fall-through
    case kFmsSouthLarge:
      mClusters.remove_if(IsBadCluster(0.75, 25));
      break;
    case kFmsNorthSmall:  // Deliberate fall-through
    case kFmsSouthSmall:
      mClusters.remove_if(IsBadCluster(2.0, 49));
      break;
    default:
      break;
  }  // switch
  // Must do moment analysis before catagorization
  for (auto i = mClusters.begin(); i != mClusters.end(); ++i) {
    (*i)->findClusterAxis(mClusterFinder.momentEnergyCutoff());
  }  // for
  return mClusters.size();
}

Bool_t StFmsEventClusterer::fitClusters() {
  // Loop over clusters, catagorize, guess the photon locations for cat 0 or 2
  // clusters then fit, compare, and choose the best fit
  bool badFit = false;
  for (auto iter = mClusters.begin(); iter != mClusters.end(); ++iter) {
    int category = mClusterFinder.categorise(iter->get());
    mFitter->setTowers(&(*iter)->towers());
    switch (category) {
      case k1PhotonCluster:
        fitOnePhoton(iter->get());
        break;
      case k2PhotonCluster:
        fit2PhotonClust(iter);
        break;
      case kAmbiguousCluster:
        category = fitAmbiguousCluster(iter);
        break;
      default:
        LOG_ERROR << "The logic of cluster catagory is wrong and something "
          << "impossible has happened! This a catagory-" << category <<
          " cluster! Do not know how to fit it!" << endm;
        break;
    }  // switch
    if (category == k2PhotonCluster && (*iter)->chiSquare() > 10.) {
      badFit = true;
    }  // if
  }  // Loop over all real clusters
  return !badFit;
}

Bool_t StFmsEventClusterer::refitClusters() {
  // Only do a global fit for 2 or more clusters (2-photon fit for one cluster
  // already performs a global fit as part of its normal procedure)
  if (mClusters.size() < 2) {
    return true;
  }  // if
  Towers towers;
  for (auto i = mClusters.begin(); i != mClusters.end(); ++i) {
    towers.insert(towers.end(), (*i)->towers().begin(), (*i)->towers().end());
  }  // for
  mFitter->setTowers(&towers);
  const int nPhotons = sumPhotonsOverClusters(mClusters);
  globalFit(nPhotons, mClusters.size(), mClusters.begin());
  return nPhotons == sumPhotonsOverClusters(mClusters);  // Shouldn't change
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
  OnePhotonFitParameters parameters(mTowerWidthXY, towerCluster->cluster());
  PhotonList photons;
  double chiSquare = mFitter->fit(parameters.start, std::vector<double>(),
                                  parameters.lower, parameters.upper, &photons);
  if (photons.empty()) {  // check return status in case of a bad fit
    LOG_ERROR << "1-photon Minuit fit found no photons" << endm;
  } else {
    towerCluster->photons()[0] = photons.back();
  }  // if
  towerCluster->cluster()->setNPhotons(photons.size());
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
  GlobalPhotonFitParameters parameters(nPhotons, first, end);
  PhotonList photons;
  Double_t chiSquare = mFitter->fit(parameters.start, std::vector<double>(),
                                    parameters.lower, parameters.upper,
                                    &photons);
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
  TwoPhotonFitParameters parameters(mTowerWidthXY, towerCluster->get());
  PhotonList photons;
  double chiSquare =
    mFitter->fit2PhotonCluster(parameters.start, parameters.steps,
                               parameters.lower, parameters.upper, &photons);
  if (photons.size() == 2) {
    (*towerCluster)->photons()[0] = photons.front();
    (*towerCluster)->photons()[1] = photons.back();
  } else {
    LOG_WARN << "2-photon Minuit fit found " << photons.size() << " photons"
      << endm;
  }  // if
  (*towerCluster)->cluster()->setNPhotons(photons.size());
  chiSquare = globalFit(2, 1, towerCluster);
  const int nDegreesOfFreedom = std::max(1,
    int((*towerCluster)->towers().size() - 6));
  (*towerCluster)->setChiSquare(chiSquare / nDegreesOfFreedom);
  return (*towerCluster)->chiSquare();
}

/* Distinguish an ambiguous cluster as either 1- or 2-photon */
Int_t StFmsEventClusterer::fitAmbiguousCluster(ClusterIter towerCluster) {
  const double chiSquare1Photon = fitOnePhoton(towerCluster->get());
  const StFmsFittedPhoton photon = (*towerCluster)->photons()[0];  // Cache
  // Decide if this 1-photon fit is good enough, if not try 2-photon fit
  int category = k1PhotonCluster;
  if (chiSquare1Photon >= 5.) {
    if (fit2PhotonClust(towerCluster) <= chiSquare1Photon &&
        validate2ndPhoton(towerCluster)) {
      category = k2PhotonCluster;
    }  // if
  }  // if
  if (category == k2PhotonCluster) {  // 2-photon fit is better
    (*towerCluster)->cluster()->setNPhotons(2);
  } else {  // 1-photon fit was better, restore it's properties
    (*towerCluster)->setChiSquare(chiSquare1Photon);
    (*towerCluster)->photons()[0] = photon;
    (*towerCluster)->cluster()->setNPhotons(1);
  }  // if
  return category;
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
  // Find the tower hit by the lowest energy photon in a cluster
  const StFmsFittedPhoton* photon = findLowestEnergyPhoton(cluster->get());
  int column = 1 + int(photon->xPos / mTowerWidthXY.at(0));
  int row = 1 + int(photon->yPos / mTowerWidthXY.at(1));
  const StFmsTower* tower = searchClusterTowers(row, column, **cluster);
  // If tower is nullptr, the photon doesn't hit in a tower in this cluster.
  if (!tower) {
    return false;
  }  // if
  // Check if the fitted energy is too large compared to the energy of the tower
  if (tower->hit()->energy() < 0.25 * photon->energy) {
    return false;
  }  // if
  // Check if the 2nd photon's "high-hower" enery is too large compared to its
  // fitted energy. If so, it is probably splitting one photon into two
  if (tower->hit()->energy() > 1.5 * photonEnergyInTower(tower, photon)) {
    return false;
  }  // if
  // Check that the 2nd photon is not near the edge of another cluster
  const double energyInOwnCluster =
    photonEnergyInCluster(cluster->get(), photon);
  for (ClusterConstIter i = mClusters.begin(); i != mClusters.end(); ++i) {
    if (i != cluster) {  // Skip the photon's own cluster
      if (photonEnergyInCluster(i->get(), photon) > 0.2 * energyInOwnCluster) {
        return false;  // Stop as soon as we fail for one cluster
      }  // if
    }  // if
  }  // for
  return true;  // The photon passed all tests; it's real
}

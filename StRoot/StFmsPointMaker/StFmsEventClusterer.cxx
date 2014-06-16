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
#include "StFmsEventClusterer.h"

#include <algorithm>
#include <functional>
#include <iterator>
#include <list>
#include <numeric>

#include <boost/foreach.hpp>

#include <TRandom.h>  // For ROOT global random generator, gRandom

#include "StRoot/St_base/StMessMgr.h"
#include "StEvent/StFmsCluster.h"
#include "StEvent/StFmsHit.h"

#include "StFmsPointMaker/StFmsTower.h"

using namespace std;
using namespace FMSCluster;

namespace {
/* Helper function to add numbers of photons using std::accumulate */
int accumulatePhotons(int nPhotons, const StFmsTowerCluster& cluster) {
  return nPhotons + cluster.cluster()->nPhotons();
}

/* Unary predicate for selecting bad clusters. */
struct IsBadCluster
    : public std::unary_function<const StFmsTowerCluster&, bool> {
  // Set minimum allowed cluster energy and maximum number of towers
  IsBadCluster(double minEnergy, unsigned maxTowers)
      : energy(minEnergy), towers(maxTowers) { }
  bool operator()(const StFmsTowerCluster& cluster) const {
    return cluster.cluster()->energy() <= energy ||
           cluster.towers().size() > towers;
  }
  double energy;
  unsigned towers;
};

/*
 Returns a pointer to the lowest energy photon in a cluster
 
 Assumes the cluster is either 1- or 2-photon
 Returns NULL if there is no photon in the cluster
 */
const StFmsFittedPhoton* findLowestEnergyPhoton(
    const StFmsTowerCluster* cluster) {
  const StFmsFittedPhoton* photon(NULL);
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
      break;  // photon remains NULL
  }  // switch
  return photon;
}

/*
 Search towers in a cluster for one matching a row and column number
 
 Return a pointer to the matching tower if one is found, NULL otherwise.
 */
const StFmsTower* searchClusterTowers(int row, int column,
                                      const StFmsTowerCluster& cluster) {
  const StFmsTower* match(NULL);
  BOOST_FOREACH(const StFmsTower* tower, cluster.towers()) {
    if (tower->row() == row && tower->column() == column) {
      match = tower;
      break;
    }  // if
  }  // BOOST_FOREACH
  return match;
}
}  // unnamed namespace

StFmsEventClusterer::StFmsEventClusterer(const StFmsGeometry* geometry,
                                         Int_t detectorId) {
  mGeometry = geometry;
  mDetectorId = detectorId;
}

StFmsEventClusterer::~StFmsEventClusterer() {
  if (mFitter) {
    delete mFitter;
  }  // if
}

Bool_t StFmsEventClusterer::cluster(std::vector<StFmsTower>* towerList) {
  mTowers = towerList;
  mClusterFinder.setMomentEnergyCutoff(.5);  
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
  std::vector<StFmsTower>::iterator towerIter;
  for (towerIter = mTowers->begin(); towerIter != mTowers->end(); ++towerIter) {
    towerList.push_back(&(*towerIter));
  }  // for
  mClusterFinder.findClusters(&towerList, &mClusters);
  // Cluster energy should be at least 2 GeV (parameter "minRealClusterEne")
  if (mDetectorId == 8 || mDetectorId == 9) {
    mClusters.erase_if(IsBadCluster(0.75, 25));
  } else {
    mClusters.erase_if(IsBadCluster(2.0, 49));  // Different cuts for small cell
  }  // if
  // Must do moment analysis before catagorization
  for (ClusterIter i = mClusters.begin(); i != mClusters.end(); ++i) {
    i->findClusterAxis(mClusterFinder.momentEnergyCutoff());
  }  // for
  // Loop over clusters, catagorize, guess the photon locations for cat 0 or 2
  // clusters then fit, compare, and choose the best fit
  bool badEvent = false;
  const double max2PhotonFitChi2 = 10.;
  for (ClusterIter cluster = mClusters.begin(); cluster != mClusters.end();
       ++cluster) {
    Int_t clustCatag = mClusterFinder.categorise(&(*cluster));
    // point to the real TObjArray that contains the towers to be fitted
    // it is the same tower array for the cluster or all alternative clusters
    mFitter->setTowers(&cluster->towers());
    // Number of Degree of Freedom for the fit
    if (clustCatag == k1PhotonCluster) {
      // Do 1-photon fit
      fitOnePhoton(&(*cluster));
    } else if (clustCatag == k2PhotonCluster) {
      // Do 2-photon fit
      fit2PhotonClust(cluster);
      badEvent = cluster->chiSquare() > max2PhotonFitChi2;
    } else if (clustCatag == kAmbiguousCluster) {
      // for catagory-0 cluster, first try 1-photon fit!
      // If the fit is good enough, it is 1-photon. Else also
      // try 2-photon fit, and find the best fit (including 1-photon fit).
      Bool_t is2Photon = true;
      double chiSq1 = fitOnePhoton(&(*cluster));
      const StFmsFittedPhoton photon = cluster->photons()[0];  // Cache photon
      double chiSq2(NAN);  // Only set if do 2-photon fit
      // Decide if this 1-photon fit is good enough
      if (chiSq1 < 5.) {
        is2Photon = false ;
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
        cluster->cluster()->setNPhotons(2);
        cluster->setChiSquare(chiSq2);
        // Flag the event as bad if the fit chi2/ndf is too bad
        badEvent = cluster->chiSquare() > max2PhotonFitChi2;
      } else {
        // 1-photon fit is better
        cluster->cluster()->setNPhotons(1);
        cluster->setChiSquare(chiSq1);
        cluster->photons()[0] = photon;
      }  // if (is2Photon)
    } else {  // Invalid cluster category
      // should not happen!
      LOG_ERROR << "Your logic of catagory is wrong! Something impossible " <<
        "happens! This a catagory-" << clustCatag <<
        " clusters! Don't know how to fit it!" << endm;
    }  // if (clustCatag...)
  }  // Loop over all real clusters
  Int_t nPh = std::accumulate(mClusters.begin(), mClusters.end(), 0,
                              accumulatePhotons);
  if(nPh > StFmsClusterFitter::maxNFittedPhotons()) {
    // myFitter can only do up to "maxNFittedPhotons()"-photon fit
    LOG_WARN << "Can not fit " << nPh << " (more than " <<
      StFmsClusterFitter::maxNFittedPhotons() << " photons!" << endm;
    return nPh;
  }  // if
  // For global fit, add all towers from all clusters
  std::list<StFmsTower*> allTow;
  Int_t ndfg = 0 ;
  for (ClusterIter cluster = mClusters.begin(); cluster != mClusters.end();
       ++cluster) {
    allTow.insert(allTow.end(), cluster->towers().begin(),
                  cluster->towers().end());
    ndfg += (cluster->towers().size() -
             3 * cluster->cluster()->nPhotons());
  }  // for
  if(ndfg <= 0) {
    ndfg = 1;
  }  // if
  mFitter->setTowers(&allTow);
  // Only do global fit for 2 or more clusters (2-photon fit for one cluster
  // already has global fit)
  if (mClusters.size() > 1) {
    double chiSqG = globalFit(nPh, mClusters.size(), mClusters.begin()) / ndfg;
    // Check for errors in the global fit - the number of photons returned by
    // the global fit should equal the sum of photons in the fitted clusters
    Int_t iph = std::accumulate(mClusters.begin(), mClusters.end(), 0,
                                accumulatePhotons);
    if(iph != nPh) {
      LOG_ERROR << "total nPh=" << nPh << " iPh=" << iph << endm;
    }  // if
  } else if (mClusters.size() == 1) {
      double chiSqG = mClusters.front().chiSquare();
  } else {
    double chiSqG = -1 ;
  }  // if (mClusters.size() > 1)
  return !badEvent;
}

Double_t StFmsEventClusterer::photonEnergyInCluster(
    Double_t widthLG,
    const StFmsTowerCluster *p_clust,
    const StFmsFittedPhoton *p_photon) const {
  Double_t eSS = 0;
  // Sum depositions by the photon in all towers of this cluster
  BOOST_FOREACH(const StFmsTower* tower, p_clust->towers()) {
    eSS += photonEnergyInTower(widthLG, tower, p_photon);
  }  // BOOST_FOREACH
  return eSS;
}

Double_t StFmsEventClusterer::photonEnergyInTower(
    Double_t widthLG,
    const StFmsTower *p_tower,
    const StFmsFittedPhoton* p_photon) const {
  Double_t xx = ((Double_t)p_tower->column() - 0.5) *
                 mTowerWidthXY[0] - p_photon->xPos;
  Double_t yy = ((Double_t)p_tower->row() - 0.5) *
                 mTowerWidthXY[1] - p_photon->yPos;
  Double_t eSS = p_photon->energy *
                 mFitter->showerShapeFunction()->Eval(xx, yy, 0);
  return eSS;
}

Float_t StFmsEventClusterer::fitOnePhoton(StFmsTowerCluster* p_clust) {
  // 4 parameters are passed to the fitting routine: nPhotons, cluster x
  // position, cluster y position and cluster energy. Set the starting points
  // for the fitting routine, plus lower and upper bounds on allowed values.
  // - set starting points for the fit parameters:
  const Double_t start[4] = {
    1.0, mTowerWidthXY[0] * p_clust->cluster()->x(),
    mTowerWidthXY[1] * p_clust->cluster()->y(),
    p_clust->cluster()->energy()};
  // - maximum deviations from the start points during fit:
  const Double_t delta[4] = {
    0.5, 0.5 * mTowerWidthXY[0], 0.5 * mTowerWidthXY[1],
    0.15 * p_clust->cluster()->energy()};
  // - set lower and upper physical limits of fit parameters = start +/- delta
  //   The parameters will stay within these ranges during the fit
  Double_t lowLim[4], upLim[4];
  for (int i(0); i < 4; ++i) {
    lowLim[i] = start[i] - delta[i];
    upLim[i] = start[i] + delta[i];
  }  // for
  PhotonList photons;
  Double_t chiSq = mFitter->fit(start, NULL, lowLim, upLim, &photons);
  if (photons.empty()) {  // check return status in case of a bad fit
    LOG_ERROR << "1-photon Minuit fit returns error!" << endm;
  }  // if
  p_clust->photons()[0] = photons.back();
  p_clust->cluster()->setNPhotons(photons.size());
  int ndf = p_clust->towers().size() - 3;
  if (ndf <= 0) {
    ndf = 1;
  }  // if
  p_clust->setChiSquare(chiSq / ndf);
  return p_clust->chiSquare();
}

Float_t StFmsEventClusterer::globalFit(const Int_t nPh, const Int_t nCl,
                                       ClusterIter first) {
  // By design, we can only fit up to "maxNFittedPhotons()" photons
  if (nPh > StFmsClusterFitter::maxNFittedPhotons() || nPh < 2) {
    LOG_ERROR << "Global fit! Can not fit " << nPh << " photons!" << endm;
    return -9999;
  }  // if
  // Check that there is at least one cluster
  if (nCl < 1) {
    LOG_ERROR << nCl << " clusters! Global fit will NOT work!" << endm;
    return -9999;
  }  // if
  // Fit has 3 parameters per photon (x, y, E), plus 1 for the number of photons
  const Int_t nParam = 3 * StFmsClusterFitter::maxNFittedPhotons() + 1;
  // Starting position, lower and upper limit of parameters
  Double_t start[nParam], lowLim[nParam], upLim[nParam];
  // The positions (e.b. cluster->photons()[jp].xPos) are already in unit of cm
  // Clusters have already had all their fields properly filled
  // (for example cluster[].photons()[0] should NOT be NULL!)
  Int_t totPh = 0;
  // Loop over all clusters
  /** \todo Improve this implementation? This approach is necessary because the
            original code uses a pointer to a StFmsTowerCluster in the fitting
            routines as both a pointer to a single cluster, and an array of
            clusters .*/
  ClusterIter end = first;
  std::advance(end, nCl);
  for (ClusterIter cluster = first; cluster != end; ++cluster) {
    // Loop over all photons in cluster
    for (Int_t jp = 0; jp < cluster->cluster()->nPhotons(); jp++) {
      if (totPh > StFmsClusterFitter::maxNFittedPhotons()) {
        LOG_ERROR << "Total # of photons in " << nCl << " clusters is at least "
          << totPh << "! I can NOT do fit!" << endm;
        return -9999;
      }  // if
      // Note positions are in centimetres, not tower unites
      Int_t kpar = 3 * totPh + 1;
      start[kpar] = cluster->photons()[jp].xPos;
      lowLim[kpar] = start[kpar] - 1.25;
      upLim[kpar] = start[kpar] + 1.25;
      kpar++;
      start[kpar] = cluster->photons()[jp].yPos;
      lowLim[kpar] = start[kpar] - 1.25;
      upLim[kpar] = start[kpar] + 1.25;
      kpar++;
      start[kpar] = cluster->photons()[jp].energy;
      lowLim[kpar] = start[kpar] * (1 - 0.3);  // Limit to +/- 30% energy
      upLim[kpar] = start[kpar] * (1 + 0.3);
      totPh++ ;
    }  // for
  }  // for
  if (totPh != nPh) {
    LOG_WARN << "WARNING! Total # of photons in " << nCl <<
      " clusters is at least " << totPh << "! Not the same as the nPh = "
      << nPh << "! I will try " << totPh << " instead!" << endm;
  }  // if
  // Set the number-of-photons fit parameter
  start[0] = totPh;
  lowLim[0] = 0.5;
  upLim[0] = StFmsClusterFitter::maxNFittedPhotons() + 0.5 ;
  PhotonList photons;
  Double_t chiSq = mFitter->fit(start, NULL, lowLim, upLim, &photons);
  if (photons.empty()) {
    LOG_WARN << "Global Minuit fit returns error!" << endm;
  }  // if
  // Put the fit result back in the clusters
  // Loop over all clusters
  PhotonList::const_iterator photonIter = photons.begin();
  for (ClusterIter cluster = first; cluster != end; ++cluster) {
    // Loop over all photons in cluster
    for (Int_t jp = 0; jp < cluster->cluster()->nPhotons();
         jp++, ++photonIter) {
      cluster->photons()[jp] = *photonIter;
    }  // for loop over photons
  }  // for loop over clusters
  // Evaluate the Chi-square function and return it
  return chiSq;
}

Float_t StFmsEventClusterer::fit2PhotonClust(ClusterIter p_clust) {
  const Double_t step2[7] = {0, 0.02, 0.02, 0.01, 0.01, 0.01, 0.1};
  Double_t ratioSigma = p_clust->cluster()->sigmaMin() /
                        p_clust->cluster()->sigmaMax();
  Double_t maxTheta = ratioSigma / 2.8;
  if (maxTheta > (TMath::Pi() / 2.0)) {
    maxTheta = TMath::Pi() / 2.0;
  }  // if
  // Use for restricting d_gg
  Double_t EcSigmaMax = p_clust->cluster()->energy() *
                        p_clust->cluster()->sigmaMax();
  // Starting position, lower and upper limit of parameters
  Double_t start[7], lowLim[7], upLim[7];
  // First parameter is the number of photons, which is constant = 2 photons
  start[0] = 2;
  lowLim[0] = 1.5;
  upLim[0] = 2.5;
  // Parameter starting points and limits are determined by looking at cluster
  // information
  //  - xPi and yPi: rarely do they go beyond 0.3 unit of lgd
  //  - theta:       have a narrow theta range (for r=sigmaMax/sigmaMax,
  //                 |theta|<0.5*r/0.65 when r<0.65, and linear increase from
  //                 0.5 to Pi/2 for 0.65<r<1)
  //  - E_gg:        given by Ec (+/- 20% or less)
  //  - z_gg:        should just let it vary from -1 to 1.
  //  - d_gg:        a lower bound is given by r=sqrt(sigmaX^2+sigmaY^2). 
  //                 d_gg > Max( 2.5*(r-0.6), 0.5 )
  start[1]  = mTowerWidthXY[0] * p_clust->cluster()->x();
  start[2]  = mTowerWidthXY[1] * p_clust->cluster()->y();
  start[6]  = p_clust->cluster()->energy();
  start[4]  = p_clust->thetaAxis();
  const float dggPara[6] = {18.0, 2.2, 0.5, 60.0, 0.085, 3.5};
  start[3] = dggPara[1] * mTowerWidthXY[0] * p_clust->cluster()->sigmaMax();
  // Randomize the starting point of Z_gg (from -0.1 to 0.1)
  start[5]  = 0.1 * (2 * gRandom->Rndm() - 1);
  lowLim[1] = start[1] - 0.2 * mTowerWidthXY[0];
  lowLim[2] = start[2] - 0.2 * mTowerWidthXY[1];
  lowLim[6] = start[6] * (1. - 0.05);
  upLim[1]  = start[1] + 0.2 * mTowerWidthXY[0];
  upLim[2]  = start[2] + 0.2 * mTowerWidthXY[1];
  upLim[6]  = start[6] * (1. + 0.05);
  lowLim[4] = start[4] - maxTheta;
  lowLim[5] = - 1.0;
  lowLim[3] = dggPara[0] / pow(EcSigmaMax, 0.8);
  if (lowLim[3] < dggPara[2]) {
    lowLim[3] = dggPara[2];
  }  // if
  lowLim[3] *= mTowerWidthXY[0];
  if (lowLim[3] >= start[3]) {
    lowLim[3] = start[3] * 0.9;
  }  // if
  upLim[3] = dggPara[4] * (dggPara[3] - EcSigmaMax);
  if (upLim[3] > dggPara[5]) {
    upLim[3] = dggPara[5];
  }  // if
  upLim[3] *= mTowerWidthXY[0] ;
  if (upLim[3] <= start[3]) {
    upLim[3] = start[3] * 1.1;
  }  // if
  upLim[4] = start[4] + maxTheta;
  upLim[5] = 1.0;
  // Call special 2-photon-cluster mFitter
  PhotonList photons;
  Double_t chiSq = mFitter->fit2PhotonCluster(start, step2, lowLim, upLim,
                                              &photons);
  if (photons.empty()) {
    LOG_WARN << "Minuit fit returns error!" << endm;
  }  // if
  // Do a global fit, using result of 1st fit as starting point
  // Need to set "nPhoton" before calling "globalFit(..)"
  p_clust->photons()[0] = photons.front();
  p_clust->photons()[1] = photons.back();
  p_clust->cluster()->setNPhotons(photons.size());
  chiSq = globalFit(2, 1, p_clust);
  int ndf = p_clust->towers().size() - 6;
  if (ndf <= 0) {
    ndf = 1;
  }  // if
  p_clust->setChiSquare(chiSq / ndf);
  return p_clust->chiSquare();
};

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
  const StFmsFittedPhoton* photon = findLowestEnergyPhoton(&(*cluster));
  // Tower row and column where the fitted photon of lower energy should hit
  int column = 1 + (Int_t)(photon->xPos / mTowerWidthXY[0]);
  int row = 1 + (Int_t)(photon->yPos / mTowerWidthXY[1]);
  // Now check whether this tower is one of the non-zero towers of the cluster
  // The temporary StFmsTower only needs row and column set for the test
  const StFmsTower* tower = searchClusterTowers(row, column, *cluster);
  // If tower is non-NULL, the photon does hit in a tower in this cluster.
  if (!tower) {
    return false;
  }  // if
  // Now test the photon and tower properties.
  // Check if the fitted energy is too large compared to the energy of the tower
  if(tower->hit()->energy() < 0.25 * photon->energy) {
    return false;
  }  // if
  // Check if the 2nd photon's "High-Tower" enery is too large compared to its
  // fitted energy. If so, it is probably splitting one photon into two
  Double_t eSS = photonEnergyInTower(mTowerWidthXY[0], tower, photon);
  if(tower->hit()->energy() > 1.5 * eSS) {
    return false;
  }  // if
  // Check that the 2nd photon is not near the edge of another cluster
  // Namely, we check what would be the energy deposited in other clusters by
  // this photon vs. energy deposited in its own cluster
  // If the ratio is too high, this fitted photon is probably a bogus one
  Double_t energyInOwnCluster =
    photonEnergyInCluster(mTowerWidthXY[0], &(*cluster), photon);
  // Loop over all clusters except its own
  for (ClusterConstIter i = mClusters.begin(); i != mClusters.end(); ++i) {
    if (i != cluster) {  // Skip the photon's own cluster
      if (photonEnergyInCluster(mTowerWidthXY[0], &(*i), photon) > 
          (0.2 * energyInOwnCluster)) {
        return false;  // Stop as soon as we fail for one cluster
      }  // if
    }  // if
  }  // for
  return true;  // The photon passed all tests; it's real
}

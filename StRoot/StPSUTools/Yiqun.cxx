#include "Yiqun.h"

using namespace std;
using namespace PSUGlobals;

#include <TCanvas.h>
#include <TRandom.h>  // For ROOT global random generator, gRandom

#include <algorithm>
#include <functional>
#include <list>
#include <numeric>

#include "StEvent/StFmsHit.h"

namespace {
/* Helper function to add numbers of photons using std::accumulate */
int accumulatePhotons(int nPhotons, const HitCluster& cluster) {
  return nPhotons + cluster.nPhoton;
}

/*
 Lack of c++0x in STAR-standard gcc as of writing so we provide our own
 implementation of copy_if - http://en.cppreference.com/w/cpp/algorithm/copy
 */
template<class InputIterator, class OutputIterator, class UnaryPredicate>
OutputIterator copy_if(InputIterator start, InputIterator end, 
                       OutputIterator outStart, UnaryPredicate predicate) {
  return std::remove_copy_if(start, end, outStart, std::not1(predicate));
}

/* Unary predicate for selecting good clusters. */
struct IsGoodCluster : public std::unary_function<const HitCluster*, bool> {
  // Initialise with minimum energy and maximum number of towers
  IsGoodCluster(double minEnergy, int maxTowers)
    : energy(minEnergy), towers(maxTowers) { }
  // Returns true if the cluster has both energy > min energy and number of
  // towers <= maximum number of towers. Returns false otherwise.
  // Use a pointer argument to avoid reference-to-reference problems when
  // binding the predicate.
  bool operator()(const HitCluster* cluster) const {
    return cluster->energy > energy && cluster->tow->GetEntries() <= towers;
  }
  double energy;
  int towers;
};

/*
 Returns a pointer to the lowest energy photon in a cluster
 
 Assumes the cluster is either 1- or 2-photon
 Returns NULL if there is no photon in the cluster
 */
PhotonHitFPD* findLowestEnergyPhoton(HitCluster* cluster) {
  PhotonHitFPD* photon(NULL);
  switch (cluster->nPhoton) {
    case 1:
      photon = &(cluster->photon[0]);
      break;
    case 2:
      if (cluster->photon[0].energy < cluster->photon[1].energy) {
        photon = &(cluster->photon[0]);
      } else {
        photon = &(cluster->photon[1]);
      }  // if
    default:
      break;  // photon remains NULL
  }  // switch
  return photon;
}

/*
 Search towers in a cluster for one matching a test tower
 
 Matching is done via TowerFPD::IsEqual
 Return a pointer to the matching tower if one is found, NULL otherwise.
 */
TowerFPD* searchClusterTowers(const TowerFPD& test, const HitCluster& cluster) {
  TowerFPD* match(NULL);
  for (Int_t i(0); i < cluster.numbTower; ++i) {
    TowerFPD* tower = static_cast<TowerFPD*>(cluster.tow->At(i));
    if (test.IsEqual(tower)) {
      match = tower;
      break;
    }  // if
  }  // for
  return match;
}
}  // unnamed namespace

// Static members
TF1* Yiqun::EDepCorrection(NULL);

Float_t Yiqun::FitOnePhoton(HitCluster* p_clust) {
  // 4 parameters are passed to the fitting routine: nPhotons, cluster x
  // position, cluster y position and cluster energy. Set the starting points
  // for the fitting routine, plus lower and upper bounds on allowed values.
  // - set starting points for the fit parameters:
  const Double_t start[4] = {
    1.0, widLG[0] * p_clust->x0, widLG[1] * p_clust->y0, p_clust->energy};
  // - maximum deviations from the start points during fit:
  const Double_t delta[4] = {
    0.5, 0.5 * widLG[0], 0.5 * widLG[1], 0.15 * p_clust->energy};
  // - set lower and upper physical limits of fit parameters = start +/- delta
  //   The parameters will stay within these ranges during the fit
  Double_t lowLim[4], upLim[4];
  for (int i(0); i < 4; ++i) {
    lowLim[i] = start[i] - delta[i];
    upLim[i] = start[i] + delta[i];
  }  // for
  Int_t status = fitter->Fit(start, fitter->step, lowLim, upLim);
  if (status != 0) {  // check return status in case of a bad fit
    std::cerr << "Minuit fit returns " << status << "!" << std::endl;
  }  // if
  // Get the fit results for starting positions and errors
  Double_t param[4];
  Double_t error[4];
  fitter->fMn->GetParameter(0, param[0], error[0]);
  // There are 3 parameters per photon, plus the 1st parameter, which gives the
  // number of photons
  Int_t nPar = 3 * ((Int_t)param[0]) + 1;
  for (Int_t i(1); i < nPar; ++i) {
    fitter->fMn->GetParameter(i, param[i], error[i]);
  }  // for
  // put the fit result back in "clust"
  p_clust->photon[0].xPos    = param[1];
  p_clust->photon[0].errXPos = error[1];
  p_clust->photon[0].yPos    = param[2];
  p_clust->photon[0].errYPos = error[2];
  p_clust->photon[0].energy  = param[3];
  p_clust->photon[0].errEne  = error[3];
  // evaluate the Chi-square function
  Double_t gradient[4];
  Double_t chiSq;
  Int_t iflag = 1;
  fitter->fMn->Eval(4, gradient, chiSq, param, iflag);
  p_clust->nPhoton = 1;
  int ndf = p_clust->tow->GetEntriesFast() - 3;
  if (ndf <= 0) {
    ndf = 1;
  }  // if
  p_clust->chiSquare = chiSq / ndf;
  return p_clust->chiSquare;
}

/*
 Perform a global fit of all photons in an event.
 
 Update the (x, y) positions end energies of the photons in each cluster based
 on a global fit including all photons.
 Only makes sense when there is more than one photon in the event.
 Arguments:
  nPh - number of photons in the event
  nCl - number of clusters containing those photons
  p_clust - cluster array
 */
Float_t Yiqun::GlobalFit(const Int_t nPh, const Int_t nCl,
                         HitCluster *p_clust) {
  // By design, we can only fit up to "MAX_NUMB_PHOTONS" (currently 4) photons
  if (nPh > FitTower::MAX_NUMB_PHOTONS || nPh < 2) {
    std::cout << "Global fit! Can not fit " << nPh << " photons! ERROR!" << "\n";
    return -9999;
  }  // if
  // Check that there is at least one cluster
  if (nCl < 1) {
    std::cout << nCl << " clusters! Global fit will NOT work! ERROR!" << "\n";
    return -9999;
  }  // if
  // Fit has 3 parameters per photon (x, y, E), plus 1 for the number of photons
  const Int_t nParam = 3 * FitTower::MAX_NUMB_PHOTONS + 1;
  // Fit parameters, errors, and gradients of function
  Double_t param[nParam];
  Double_t error[nParam];
  Double_t gradient[nParam];
  // Starting position, lower and upper limit of parameters
  Double_t start[nParam], lowLim[nParam], upLim[nParam];
  // The positions (e.b. clust[ic].photon[jp].xPos) are already in unit of cm
  // Clusters have already had all their fields properly filled
  // (for example cluster[].photon[0] should NOT be NULL!)
  Int_t totPh = 0;
  // Loop over all clusters
  for (Int_t ic=0; ic< nCl; ic++) {
    // Loop over all photons in cluster
    for (Int_t jp = 0; jp < p_clust[ic].nPhoton; jp++) {
      if (totPh > FitTower::MAX_NUMB_PHOTONS) {
        std::cout << "Total # of photons in " << nCl << " clusters is at least "
          << totPh << "! I can NOT do fit! ERROR!" << "\n";
        return -9999;
      }  // if
      // Note positions are in centimetres, not tower unites
      Int_t kpar = 3 * totPh + 1;
      start[kpar] = p_clust[ic].photon[jp].xPos;
      lowLim[kpar] = start[kpar] - posDif_Gl;
      upLim[kpar] = start[kpar] + posDif_Gl;
      kpar++;
      start[kpar] = p_clust[ic].photon[jp].yPos;
      lowLim[kpar] = start[kpar] - posDif_Gl;
      upLim[kpar] = start[kpar] + posDif_Gl;
      kpar++;
      start[kpar] = p_clust[ic].photon[jp].energy;
      lowLim[kpar] = start[kpar] * (1 - eneRat_Gl);
      upLim[kpar] = start[kpar] * (1 + eneRat_Gl);
      totPh++ ;
    }  // for
  }  // for
  if (totPh != nPh) {
    std::cout << "WARNING! Total # of photons in " << nCl <<
      " clusters is at least " << totPh << "! Not the same as the nPh = ";
    std::cout  << nPh << "! I will try " << totPh << " instead!" << "\n";
  }  // if
  // Set the number-of-photons fit parameter
  start[0] = totPh;
  lowLim[0] = 0.5;
  lowLim[0] = FitTower::MAX_NUMB_PHOTONS + 0.5 ;
  // Fit status, and flag needed by fitter
  Int_t status, iflag=1;
  status = fitter->Fit(start, fitter->step, lowLim, upLim);
  if (status != 0) {
    std::cout << "Minuit fit returns " << status << "!" << "\n";
  }  // if
  // Get the fit results
  fitter->fMn->GetParameter(0, param[0], error[0]);
  Int_t nPar = 3 * ((Int_t)param[0]) + 1;
  for(Int_t i = 1; i < nPar; ++i) {
    fitter->fMn->GetParameter(i, param[i], error[i]);
  }  // for
  // Put the fit result back in the clusters
  // Loop over all clusters
  Int_t tPh = 0 ;
  for (Int_t ic = 0; ic < nCl; ic++) {
    // Loop over all photons in cluster
    for (Int_t jp = 0; jp < p_clust[ic].nPhoton; jp++) {
      Int_t kpar = 3 * tPh + 1 ;
      p_clust[ic].photon[jp].xPos    = param[kpar] ;
      p_clust[ic].photon[jp].errXPos = error[kpar] ;
      p_clust[ic].photon[jp].yPos    = param[kpar+1] ;
      p_clust[ic].photon[jp].errYPos = error[kpar+1] ;
      p_clust[ic].photon[jp].energy  = param[kpar+2] ;
      p_clust[ic].photon[jp].errEne  = error[kpar+2] ;
      tPh++;
    }  // for loop over photons
  }  // for loop over clusters
  // Evaluate the Chi-square function and return it
  Double_t chiSq;
  fitter->fMn->Eval(7, gradient, chiSq, param, iflag);
  return chiSq;
}

/*
 Special 2-photon cluster fitting
 
 Cluster moments must have been calculated first
 */
Float_t Yiqun::Fit2PhotonClust(HitCluster* p_clust) {
  const Double_t step2[7] = {0, 0.02, 0.02, 0.01, 0.01, 0.01, 0.1} ;
  Double_t ratioSigma = p_clust->sigmaMin / p_clust->sigmaMax ;
  Double_t maxTheta = ratioSigma / thetaPara ;
  if (maxTheta > (TMath::Pi() / 2.0)) {
    maxTheta = TMath::Pi() / 2.0;
  }  // if
  // Use for restricting d_gg
  Double_t EcSigmaMax = p_clust->energy * p_clust->sigmaMax ;
  Double_t chiSq;
  // Fit parameters (starting positions), errors, and gradients of function
  Double_t param[7];
  Double_t error[7];
  Double_t gradient[7];
  // Starting position, lower and upper limit of parameters
  Double_t start[7], lowLim[7], upLim[7];
  // Fit status, and flag
  Int_t status, iflag=1;
  // First parameter is the number of photons, which is constant = 2 photons
  start[0] = 2;
  lowLim[0] = 1.5;
  lowLim[0] = 2.5;
  // Parameter starting points and limits are determined by looking at cluster
  // information
  //  - xPi and yPi: rarely do they go beyond 0.3 unit of lgd
  //  - theta:       have a narrow theta range (for r=sigmaMin/sigmaMax,
  //                 |theta|<0.5*r/0.65 when r<0.65, and linear increase from
  //                 0.5 to Pi/2 for 0.65<r<1)
  //  - E_gg:        given by Ec (+/- 20% or less)
  //  - z_gg:        should just let it vary from -1 to 1.
  //  - d_gg:        a lower bound is given by r=sqrt(sigmaX^2+sigmaY^2). 
  //                 d_gg > Max( 2.5*(r-0.6), 0.5 )
  start[1]  = widLG[0] * p_clust->x0;
  start[2]  = widLG[1] * p_clust->y0;
  start[6]  = p_clust->energy;
  start[4]  = p_clust->thetaAxis;
  start[3] = dggPara[1] * widLG[0] * p_clust->sigmaMax ;
  // Randomize the starting point of Z_gg (from -0.1 to 0.1)
  start[5]  = 0.1 * (2 * gRandom->Rndm() - 1);
  lowLim[1] = start[1] - posDif_2PC * widLG[0];
  lowLim[2] = start[2] - posDif_2PC * widLG[1];
  lowLim[6] = start[6] * (1 - eneRat_2PC);
  upLim[1]  = start[1] + posDif_2PC * widLG[0];
  upLim[2]  = start[2] + posDif_2PC * widLG[1];
  upLim[6]  = start[6] * (1 + eneRat_2PC);
  lowLim[4] = start[4] - maxTheta;
  lowLim[5] = - 1.0;
  lowLim[3] = dggPara[0] / pow(EcSigmaMax, 0.8);
  if (lowLim[3] < dggPara[2]) {
    lowLim[3] = dggPara[2];
  }  // if
  lowLim[3] *= widLG[0];
  if (lowLim[3] >= start[3]) {
    lowLim[3] = start[3] * 0.9;
  }  // if
  upLim[3] = dggPara[4] * (dggPara[3] - EcSigmaMax);
  if (upLim[3] > dggPara[5]) {
    upLim[3] = dggPara[5];
  }  // if
  upLim[3] *= widLG[0] ;
  if (upLim[3] <= start[3]) {
    upLim[3] = start[3] * 1.1;
  }  // if
  upLim[4] = start[4] + maxTheta;
  upLim[5] = 1.0;
  // Call special 2-photon-cluster fitter
  status = fitter->Fit2Pin1Clust(start, step2, lowLim, upLim);
  if (status != 0) {
    std::cout << "Minuit fit returns " << status << "!" << "\n";
  }  // if
  // Get the fit result
  fitter->fMn->GetParameter(0, param[0], error[0]);
  Int_t nPar = 3*((Int_t) param[0])+1;
  for (Int_t ipar = 1; ipar < nPar; ipar++) {
    fitter->fMn->GetParameter(ipar, param[ipar], error[ipar]);
  }  // for
  // Evaluate the Chi-square function
  fitter->fMn->Eval(7, gradient, chiSq, param, iflag);
  if (status != 0) {
    std::cout << "Minuit fit returns " << status << "!" << "\n";
  }  // if
  // Put the fit result back in "clust"
  p_clust->photon[0].xPos    = param[1] + cos(param[4]) * param[3] * (1 - param[5]) / 2.0 ;
  p_clust->photon[0].errXPos = error[1] + (cos(param[4])*error[3]-error[4]*sin(param[4])*param[3])*(1-param[5])/2 - cos(param[4])*param[3]*error[5]/2.0;
  p_clust->photon[0].yPos    = param[2] + sin(param[4]) * param[3] * (1 - param[5]) / 2.0 ;
  p_clust->photon[0].errYPos = error[2] + (sin(param[4])*error[3]+error[4]*cos(param[4])*param[3])*(1-param[5])/2 - sin(param[4])*param[3]*error[5]/2.0;
  p_clust->photon[0].energy  = param[6] * (1 + param[5]) / 2.0 ;
  p_clust->photon[0].errEne  = error[6]*(1+param[5])/2.0 + param[6]*error[5]/2.0;
  // Second photon
  p_clust->photon[1].xPos    = param[1] - cos(param[4]) * param[3] * (1 + param[5]) / 2.0 ;
  p_clust->photon[1].errXPos = error[1] + (-cos(param[4])*error[3]+error[4]*sin(param[4])*param[3])*(1+param[5])/2 - cos(param[4])*param[3]*error[5]/2.0;
  p_clust->photon[1].yPos    = param[2] - sin(param[4]) * param[3] * (1 + param[5]) / 2.0 ;
  p_clust->photon[1].errYPos = error[2] + (sin(param[4])*error[3]-error[4]*cos(param[4])*param[3])*(1+param[5])/2 - sin(param[4])*param[3]*error[5]/2.0;
  p_clust->photon[1].energy  = param[6] * (1 - param[5]) / 2.0 ;
  p_clust->photon[1].errEne  = error[6]*(1-param[5])/2.0 - param[6]*error[5]/2.0;
  // Do a global fit, using result of 1st fit as starting point
  // Need to set "nPhoton" before calling "GlobalFit(..)"
  p_clust->nPhoton = 2;
  chiSq = GlobalFit(2, 1, p_clust);
  int ndf = p_clust->tow->GetEntriesFast() - 6;
  if (ndf <= 0) {
    ndf = 1;
  }  // if
  p_clust->chiSquare = chiSq / ndf ;
  return p_clust->chiSquare;
};

/*
 Run tests on the lower-energy photon in a 2-photon cluster
 
 Return true if the photon passes tests, in which case it is a real photon
 
 Return false if it fails, in which case it is a bogus photon due to some
 problem in reconstruction - the cluster is actually a 1-photon cluster

 Further information:
 If one photon peak lies on top of a low (compared to photon energy) or even
 zero tower, this photon is definitely bogus. This could happen if a nearby
 cluster took away towers that might make the fit have a large Chi-square.
 So first check that the fitted photon of the lower-energy photon (we assume the
 higher energy photon should be fine) is over one of the non-zero towers of the
 cluster. First of all, this ensures that we don't have an "outside of cluster"
 bogus photon, i.e. a bogus photon that could be the result of minimizing the
 chi-square over towers that do not include the supposed peak tower. 
 
 Arguments:
  - clusterIndex: index of cluster to test
  - nRealClusters: total number of clusters in the event
 */
bool Yiqun::validate2ndPhoton(int clusterIndex, int nRealClusters) {
  HitCluster& cluster = clust[clusterIndex];  // The cluster of interest
  // Select the lower-energy of the two photons
  PhotonHitFPD* photon = findLowestEnergyPhoton(&cluster);
  // Tower row and column where the fitted photon of lower energy should hit
  int column = 1 + (Int_t)(photon->xPos / widLG[0]);
  int row = 1 + (Int_t)(photon->yPos / widLG[1]);
  // Now check whether this tower is one of the non-zero towers of the cluster
  // The temporary TowerFPD only needs row and column set for the test
  TowerFPD* tower =
    searchClusterTowers(TowerFPD(NULL, column, row, 0), cluster);
  // If tower is non-NULL, the photon does hit in a tower in this cluster.
  if (!tower) {
    return false;
  }  // if
  // Now test the photon and tower properties.
  // Check if the fitted energy is too large compared to the energy of the tower
  if(tower->hit->energy() < minHTEneOverPhoton * photon->energy) {
    return false;
  }  // if
  // Check if the 2nd photon's "High-Tower" enery is too large compared to its
  // fitted energy. If so, it is probably splitting one photon into two
  Double_t eSS = EnergyInTowerByPhoton(widLG[0], tower, photon);
  if(tower->hit->energy() > maxHTEneOverPhoton * eSS) {
    return false;
  }  // if
  // Check that the 2nd photon is not near the edge of another cluster
  // Namely, we check what would be the energy deposited in other clusters by
  // this photon vs. energy deposited in its own cluster
  // If the ratio is too high, this fitted photon is probably a bogus one
  Double_t energyInOwnCluster =
    EnergyInClusterByPhoton(widLG[0], &cluster, photon);
  // Loop over all clusters except its own
  for (Int_t i = 0; i < nRealClusters; ++i) {
    if (i != clusterIndex) {  // Skip the photon's own cluster
      if (EnergyInClusterByPhoton(widLG[0], &clust[i], photon) > 
          (maxRatioSpill * energyInOwnCluster)) {
        return false;  // Stop as soon as we fail for one cluster
      }  // if
    }  // if
  }  // for
  return true;  // The photon passed all tests; it's real
}

Int_t Yiqun::FitEvent(Int_t nTows, Int_t &nClusts, Int_t &nRealClusts,
                      Bool_t &junkyEvent) {
  // Possible alternative clusters for 1-photon fit: for catagory 0
  HitCluster altClu ;
  nClusts = 0 ;
  for (Int_t ic = 0; ic < MAX_NUMER_CLUSTERS; ic++) {
    if (clust[ic].tow) {
      clust[ic].tow->Clear();
    }  // if
    clust[ic].Clear();
  }  // for
  TowerUtil::TowerList towerList;
  std::vector<TowerFPD>::iterator towerIter;
  for (towerIter = towers->begin(); towerIter != towers->end(); ++towerIter) {
    towerList.push_back(&(*towerIter));
  }  // for
  nClusts = pTowerUtil->FindTowerCluster(&towerList, clust);
  std::list<HitCluster*> pClusters;
  for (int i(0); i < nClusts; ++i) {
    pClusters.push_back(&clust[i]);
  }  // for  
  std::list<HitCluster*> myClusters;
  copy_if(pClusters.begin(), pClusters.end(), std::back_inserter(myClusters),
          IsGoodCluster(minRealClusterEne, maxHitsInRealCluster));
  // Cluster energy should be at least 2 GeV (parameter "minRealClusterEne")
  // Move all clusters above "minRealClusterEne" to front, and
  // clusters below "minRealClusterEne" to the end. This way, even if a low
  // energy cluster were before a high energy one, we would not miscount!
  // "ifront" movers from the 1st clusters, "jend" moves from the last cluster
  Int_t ifront = 0;
  Int_t jend = nClusts - 1;
  while(1) {
    // Break out the loop when "ifront" meets "jend"
    if (ifront > jend) {
      break;
    }  // if
    if (clust[ifront].energy >= minRealClusterEne 
        && clust[ifront].tow->GetEntries()<= maxHitsInRealCluster) {
      ifront++;
    } else {
      // Search from "jend" and forward, until find a clusters that is higher
      // than "minRealClusterEne"
      while (jend >= ifront && (clust[jend].energy < minRealClusterEne
             || clust[jend].tow->GetEntries() > maxHitsInRealCluster)) {
        jend--;
      }  // while
      // Check that "jend" is in front of "iend", we have exhausted the array
      if (jend < ifront) {
        break;
      } else {
        // Exchange the clusters: "ifront"<-->"jend", then increase "ifront" by
        // 1 and decrease "jend" by 1;
        altClu        = clust[ifront] ;
        clust[ifront] = clust[jend]   ;
        clust[jend]   = altClu        ;
        ifront ++;
        jend   --;
      }  // if (jend < ifront)
    }  // if (clust[ifront].energy...)
  }  // while
  // At the end of the above code, "ifront" counts the number of clusters above
  // "minRealClusterEne"
  nRealClusts = ifront ;
  if (nRealClusts != int(myClusters.size())) {
    std::cerr << "Ya blew it!" << std::endl;
    exit(-1);
  } else if (nRealClusts != nClusts) {
    std::cerr << "huzzah! Steve's algorithm gives " << nRealClusts <<
      ", mine gives " << myClusters.size() << " (out of " << nClusts <<
      " total clusters)" << std::endl;
  }  // if
  // Must do moment analysis before catagorization
  for (Int_t iic = 0; iic < nRealClusts; iic++) {
    clust[iic].FindClusterAxis(pTowerUtil->GetMomentEcutoff());
  }  // for
  // Loop over clusters, catagorize, guess the photon locations for cat 0 or 2
  // clusters then fit, compare, and choose the best fit
  junkyEvent = false;  // Bad event?
  for(Int_t icc=0; icc<nRealClusts; icc++) {
    Int_t clustCatag = pTowerUtil->CatagBySigmXY(&clust[icc]);
    // point to the real TObjArray that contains the towers to be fitted
    // it is the same tower array for the cluster or all alternative clusters
    fitter->tow2Fit = clust[icc].tow;
    // Number of Degree of Freedom for the fit
    if (clustCatag == k1PhotonCluster) {
      // Do 1-photon fit
      FitOnePhoton(&clust[icc]);
    } else if (clustCatag == k2PhotonCluster) {
      // Do 2-photon fit
      Fit2PhotonClust(&clust[icc]);
      junkyEvent = clust[icc].chiSquare > MaxChi2Catag2;
    } else if (clustCatag == kAmbiguousCluster) {
      // for catagory-0 cluster, first try 1-photon fit!
      // If the fit is good enough, it is 1-photon. Else also
      // try 2-photon fit, and find the best fit (including 1-photon fit).
      Bool_t is2Photon = true;
      altClu = clust[icc];
      double chiSq1 = FitOnePhoton(&altClu);
      double chiSq2(NAN);  // Only set if do 2-photon fit
      // Decide if this 1-photon fit is good enough
      if (chiSq1 < maxGood1PhChi2NDF) {
        is2Photon = false ;
      } else {
        // The 1-photon fit isn't good enough, so try 2-photon fit
        chiSq2 = Fit2PhotonClust(&clust[icc]);
        // Check that the 2-photon fit didn't result in a bogus 2nd photon
        if (validate2ndPhoton(icc, nRealClusts)) {
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
        clust[icc].nPhoton = 2;
        clust[icc].chiSquare = chiSq2;
        // Flag the event as bad if the fit chi2/ndf is too bad
        junkyEvent = clust[icc].chiSquare > MaxChi2Catag2;
      } else {
        // 1-photon fit is better
        clust[icc].nPhoton = 1;
        clust[icc].chiSquare = chiSq1;
        clust[icc].photon[0] = altClu.photon[0];
      }  // if (is2Photon)
    } else {  // Invalid cluster category
      // should not happen!
      std::cerr << "Your logic of catagory is wrong! Something impossible " <<
        "happens! This a catagory-" << clustCatag;
      std::cerr << " clusters! Don't know how to fit it!\n" << "\n";
    }  // if (clustCatag...)
  }  // Loop over all real clusters
  Int_t nPh = std::accumulate(clust, clust + nRealClusts, 0, accumulatePhotons);
  if(nPh > FitTower::MAX_NUMB_PHOTONS) {
    // myFitter can only do up to "MAX_NUMB_PHOTONS"-photon fit
    std::cout << "Can not fit " << nPh << " (more than " <<
      FitTower::MAX_NUMB_PHOTONS << " photons!" << "\n";
    return nPh;
  }  // if
  // For global fit, add all towers from all clusters
  TObjArray allTow(nTows);
  Int_t ndfg = 0 ;
  for(Int_t jjc = 0; jjc < nRealClusts; jjc++) {
    allTow.AddAll(clust[jjc].tow);
    ndfg += (clust[jjc].tow->GetEntriesFast() - 3 * clust[jjc].nPhoton);
  }  // for
  if(ndfg <= 0) {
    ndfg = 1;
  }  // if
  fitter->tow2Fit = &allTow;
  // Only do global fit for 2 or more clusters (2-photon fit for one cluster
  // already has global fit)
  if (nRealClusts > 1) {
    double chiSqG = GlobalFit(nPh, nRealClusts, clust) / ndfg;
    // Check for errors in the global fit - the number of photons returned by
    // the global fit should equal the sum of photons in the fitted clusters
    Int_t iph = std::accumulate(clust, clust + nRealClusts, 0,
                                accumulatePhotons);
    if(iph != nPh) {
      std::cerr << "ERROR total nPh=" << nPh << " iPh=" << iph << std::endl;
    }  // if
  } else if (nRealClusts == 1) {
      double chiSqG = clust[0].chiSquare ;
  } else {
    double chiSqG = -1 ;
  }  // if (nRealClusts > 1)
  return nPh;
}

/*
 Calculate the energy deposit in a cluster by a photon from shower-shape
 function. In this case, we only need to consider non-zero towers.
 */
Double_t Yiqun::EnergyInClusterByPhoton(Double_t widthLG, HitCluster *p_clust,
                                        PhotonHitFPD *p_photon) {
  Double_t eSS = 0;
  // Sum depositions by the photon in all towers of this cluster
  for(Int_t it=0; it<p_clust->numbTower; it++) {
    eSS += EnergyInTowerByPhoton(widthLG, (TowerFPD*)p_clust->tow->At(it),
                                 p_photon);
  };
  return eSS;
}

/*
 Calculate the energy deposit in a tower by a photon from shower-shape
 function.
 */
Double_t Yiqun::EnergyInTowerByPhoton(Double_t widthLG, TowerFPD *p_tower,
                                      PhotonHitFPD* p_photon) {
  Double_t xx = ((Double_t)p_tower->col - 0.5) * widLG[0] - p_photon->xPos;
  Double_t yy = ((Double_t)p_tower->row - 0.5) * widLG[1] - p_photon->yPos;
  Double_t eSS = p_photon->energy *
                 fitter->GetFunctShowShape()->Eval( xx, yy, 0 );
  return eSS;
}

Yiqun::Yiqun(TowerList* pEm,Geom* pgeom,Int_t iew,Int_t nstb) {
  p_geom=pgeom;
  EW=iew;
  NSTB=nstb;
  Y(pEm);
}

void Yiqun::Y(TowerList* pEm) {   
  towers = pEm;   
  NPh=0;
  NTower=0;
  NClusts=0;
  NRealClusts=0;
  pTowerUtil=new TowerUtil();
  pTowerUtil->SetMomentEcutoff(.5);  
  tow_Arr=0;
  widLG[0]=widLG[1]=(p_geom->FpdTowWid(EW,NSTB))[0];
  if (p_geom->FMSGeom) {
    widLG[1] = (p_geom->FpdTowWid(EW, NSTB))[1];
  }  // if
  NTower=pEm->size();
  if (NTower > 578) {
    printf("Too many towers for Fit\n");
    /** \todo Need to handle this more gracefully. We CANNOT exit during a
              production run */
    exit(-1);
  }  // if
  Int_t cnt = -1;
  tow_Arr = new TObjArray(NTower);
  for (TowerList::iterator i = pEm->begin(); i != pEm->end(); ++i) {
    if (i->hit->energy() > 0.001) {
      tow_Arr->Add(&(*i));
    }  // if
  }  // for
  TIter next(tow_Arr);
  while (TowerFPD* tow = (TowerFPD*)next()){
    tow->SetContext(tow_Arr, EW, NSTB);
  }  // while
  posDif_2PC=0.2;
  eneRat_2PC=0.05;
  dggPara[0]=18.0;
  dggPara[1]=2.2;
  dggPara[2]=0.5;
  dggPara[3]=60.0;
  dggPara[4]=0.085;
  dggPara[5]=3.5;
  thetaPara=2.8;
  // Global fit (more than 1 cluster)
  posDif_Gl=1.25;
  eneRat_Gl=0.3;
  maxGood1PhChi2NDF=5;
  minHTEneOverPhoton=0.25;
  maxHTEneOverPhoton=1.5;
  maxRatioSpill=0.2;
  MaxChi2Catag2=10.;
  minRealClusterEne=2.0;
  maxHitsInRealCluster=49;
  if(EW==2 && (NSTB==1 || NSTB==2)) {
    maxHitsInRealCluster=10;
    minRealClusterEne=.75;
    maxHitsInRealCluster=25;
  }  // if
  fitter = new FitTower(p_geom,EW,NSTB);
  fitter->SetTWidthCM(*(p_geom->FpdTowWid(EW, NSTB)));
  if (p_geom->FMSGeom) {
    fitter->SetXYTWidthCM(p_geom->FpdTowWid(EW, NSTB));
  }  // if
  Bool_t badEvent(true);
  NPh = FitEvent(NTower, NClusts, NRealClusts, badEvent);
}

Yiqun::~Yiqun() {
  fitter->GetFunctShowShape()->Draw("colz");
  gPad->SetLogz(true);
  gPad->Print("showerShape.eps");
  if (fitter) {
    delete fitter;
  }  // if
  if (pTowerUtil) {
    delete pTowerUtil;
  }  // if
  if (tow_Arr) {
    delete tow_Arr;
  }  // if
}

TF1* Yiqun::GetEDepCorrection() {
  if (!Yiqun::EDepCorrection) {
    // This instantiation would seemingly leak memory, as it is not freed
    // anywhere. However, ROOT stores all functions in a global function list,
    // accessible via gROOT->GetListOfFunctions(). Therefore the ROOT
    // environment should take care of deleting it.
    Yiqun::EDepCorrection = new TF1("EDepCorrection",
                                    "(1.3-.15*exp(-(x)/[0])-.6*exp(-(x)/[1]))",
                                    1, 250);
    Yiqun::EDepCorrection->SetParameter(0, 10.);
    Yiqun::EDepCorrection->SetParameter(1, 70.);
  }  // if
  return Yiqun::EDepCorrection;
}

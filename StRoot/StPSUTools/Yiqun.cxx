#include "Yiqun.h"

#include <TCanvas.h>
#include <TRandom.h>  // For ROOT global random generator, gRandom

#include <algorithm>
#include <functional>
#include <iterator>
#include <list>
#include <numeric>

#include "StEvent/StFmsHit.h"

#include "StPSUTools/TowerFPD.h"

using namespace std;
using namespace PSUGlobals;

namespace {
/* Helper function to add numbers of photons using std::accumulate */
int accumulatePhotons(int nPhotons, const HitCluster& cluster) {
  return nPhotons + cluster.GetNphoton();
}

/* Unary predicate for selecting bad clusters. */
struct IsBadCluster : public std::unary_function<const HitCluster&, bool> {
  // Set minimum allowed cluster energy and maximum number of towers
  IsBadCluster(double minEnergy, int maxTowers)
      : energy(minEnergy), towers(maxTowers) { }
  bool operator()(const HitCluster& cluster) const {
    return cluster.GetEnergy() <= energy ||
           cluster.towers()->GetEntries() > towers;
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
  switch (cluster->GetNphoton()) {
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
TowerFPD* searchClusterTowers(int row, int column, const HitCluster& cluster) {
  TowerFPD* match(NULL);
  for (Int_t i(0); i < cluster.GetNTower(); ++i) {
    TowerFPD* tower = static_cast<TowerFPD*>(cluster.towers()->At(i));
    if (tower->row() == row && tower->column() == column) {
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
    1.0, widLG[0] * p_clust->GetX0(), widLG[1] * p_clust->GetY0(),
    p_clust->GetEnergy()};
  // - maximum deviations from the start points during fit:
  const Double_t delta[4] = {
    0.5, 0.5 * widLG[0], 0.5 * widLG[1], 0.15 * p_clust->GetEnergy()};
  // - set lower and upper physical limits of fit parameters = start +/- delta
  //   The parameters will stay within these ranges during the fit
  Double_t lowLim[4], upLim[4];
  for (int i(0); i < 4; ++i) {
    lowLim[i] = start[i] - delta[i];
    upLim[i] = start[i] + delta[i];
  }  // for
  PhotonList photons;
  Double_t chiSq = fitter->Fit(start, fitter->step, lowLim, upLim, &photons);
  if (photons.empty()) {  // check return status in case of a bad fit
    std::cerr << "1-photon Minuit fit returns error!" << std::endl;
  }  // if
  p_clust->photons()[0] = photons.back();
  p_clust->SetNphoton(photons.size());
  int ndf = p_clust->towers()->GetEntriesFast() - 3;
  if (ndf <= 0) {
    ndf = 1;
  }  // if
  p_clust->setChiSquare(chiSq / ndf);
  return p_clust->chiSquare();
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
                         ClusterIter first) {
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
  // The positions (e.b. cluster->photons()[jp].xPos) are already in unit of cm
  // Clusters have already had all their fields properly filled
  // (for example cluster[].photons()[0] should NOT be NULL!)
  Int_t totPh = 0;
  // Loop over all clusters
  /** \todo Improve this implementation? This approach is necessary because the
            original code uses a pointer to a HitCluster in the fitting routines
            as both a pointer to a single cluster, and an array of clusters .*/
  ClusterIter end = first;
  std::advance(end, nCl);
  for (ClusterIter cluster = first; cluster != end; ++cluster) {
    // Loop over all photons in cluster
    for (Int_t jp = 0; jp < cluster->GetNphoton(); jp++) {
      if (totPh > FitTower::MAX_NUMB_PHOTONS) {
        std::cout << "Total # of photons in " << nCl << " clusters is at least "
          << totPh << "! I can NOT do fit! ERROR!" << "\n";
        return -9999;
      }  // if
      // Note positions are in centimetres, not tower unites
      Int_t kpar = 3 * totPh + 1;
      start[kpar] = cluster->photons()[jp].xPos;
      lowLim[kpar] = start[kpar] - posDif_Gl;
      upLim[kpar] = start[kpar] + posDif_Gl;
      kpar++;
      start[kpar] = cluster->photons()[jp].yPos;
      lowLim[kpar] = start[kpar] - posDif_Gl;
      upLim[kpar] = start[kpar] + posDif_Gl;
      kpar++;
      start[kpar] = cluster->photons()[jp].energy;
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
  upLim[0] = FitTower::MAX_NUMB_PHOTONS + 0.5 ;
  // Fit status, and flag needed by fitter
  Int_t status, iflag=1;
  PhotonList photons;
  Double_t chiSq = fitter->Fit(start, fitter->step, lowLim, upLim, &photons);
  if (photons.empty()) {
    std::cout << "Global Minuit fit returns error!" << "\n";
  }  // if
  // Put the fit result back in the clusters
  // Loop over all clusters
  PhotonList::const_iterator photonIter = photons.begin();
  for (ClusterIter cluster = first; cluster != end; ++cluster) {
    // Loop over all photons in cluster
    for (Int_t jp = 0; jp < cluster->GetNphoton(); jp++, ++photonIter) {
      cluster->photons()[jp] = *photonIter;
    }  // for loop over photons
  }  // for loop over clusters
  // Evaluate the Chi-square function and return it
  return chiSq;
}

/*
 Special 2-photon cluster fitting
 
 Cluster moments must have been calculated first
 */
Float_t Yiqun::Fit2PhotonClust(ClusterIter p_clust) {
  const Double_t step2[7] = {0, 0.02, 0.02, 0.01, 0.01, 0.01, 0.1} ;
  Double_t ratioSigma = p_clust->GetSigmaMin() / p_clust->GetSigmaMax();
  Double_t maxTheta = ratioSigma / thetaPara ;
  if (maxTheta > (TMath::Pi() / 2.0)) {
    maxTheta = TMath::Pi() / 2.0;
  }  // if
  // Use for restricting d_gg
  Double_t EcSigmaMax = p_clust->GetEnergy() * p_clust->GetSigmaMax();
  // Starting position, lower and upper limit of parameters
  Double_t start[7], lowLim[7], upLim[7];
  // Fit status, and flag
  Int_t status, iflag=1;
  // First parameter is the number of photons, which is constant = 2 photons
  start[0] = 2;
  lowLim[0] = 1.5;
  upLim[0] = 2.5;
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
  start[1]  = widLG[0] * p_clust->GetX0();
  start[2]  = widLG[1] * p_clust->GetY0();
  start[6]  = p_clust->GetEnergy();
  start[4]  = p_clust->thetaAxis();
  start[3] = dggPara[1] * widLG[0] * p_clust->GetSigmaMax();
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
  PhotonList photons;
  Double_t chiSq = fitter->Fit2Pin1Clust(start, step2, lowLim, upLim, &photons);
  if (photons.empty()) {
    std::cout << "Minuit fit returns error!" << "\n";
  }  // if
  // Do a global fit, using result of 1st fit as starting point
  // Need to set "nPhoton" before calling "GlobalFit(..)"
  p_clust->photons()[0] = photons.front();
  p_clust->photons()[1] = photons.back();
  p_clust->SetNphoton(photons.size());
  chiSq = GlobalFit(2, 1, p_clust);
  int ndf = p_clust->towers()->GetEntriesFast() - 6;
  if (ndf <= 0) {
    ndf = 1;
  }  // if
  p_clust->setChiSquare(chiSq / ndf);
  return p_clust->chiSquare();
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
bool Yiqun::validate2ndPhoton(ClusterIter cluster) {
  // Select the lower-energy of the two photons
  PhotonHitFPD* photon = findLowestEnergyPhoton(&(*cluster));
  // Tower row and column where the fitted photon of lower energy should hit
  int column = 1 + (Int_t)(photon->xPos / widLG[0]);
  int row = 1 + (Int_t)(photon->yPos / widLG[1]);
  // Now check whether this tower is one of the non-zero towers of the cluster
  // The temporary TowerFPD only needs row and column set for the test
  TowerFPD* tower = searchClusterTowers(row, column, *cluster);
  // If tower is non-NULL, the photon does hit in a tower in this cluster.
  if (!tower) {
    return false;
  }  // if
  // Now test the photon and tower properties.
  // Check if the fitted energy is too large compared to the energy of the tower
  if(tower->hit()->energy() < minHTEneOverPhoton * photon->energy) {
    return false;
  }  // if
  // Check if the 2nd photon's "High-Tower" enery is too large compared to its
  // fitted energy. If so, it is probably splitting one photon into two
  Double_t eSS = EnergyInTowerByPhoton(widLG[0], tower, photon);
  if(tower->hit()->energy() > maxHTEneOverPhoton * eSS) {
    return false;
  }  // if
  // Check that the 2nd photon is not near the edge of another cluster
  // Namely, we check what would be the energy deposited in other clusters by
  // this photon vs. energy deposited in its own cluster
  // If the ratio is too high, this fitted photon is probably a bogus one
  Double_t energyInOwnCluster =
    EnergyInClusterByPhoton(widLG[0], &(*cluster), photon);
  // Loop over all clusters except its own
  for (ClusterIter i = mClusters.begin(); i != mClusters.end(); ++i) {
    if (i != cluster) {  // Skip the photon's own cluster
      if (EnergyInClusterByPhoton(widLG[0], &(*i), photon) > 
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
  nClusts = 0 ;
  TowerUtil::TowerList towerList;
  std::vector<TowerFPD>::iterator towerIter;
  for (towerIter = towers->begin(); towerIter != towers->end(); ++towerIter) {
    towerList.push_back(&(*towerIter));
  }  // for
  nClusts = pTowerUtil->FindTowerCluster(&towerList, &mClusters);
  // Cluster energy should be at least 2 GeV (parameter "minRealClusterEne")
  mClusters.erase_if(IsBadCluster(minRealClusterEne, maxHitsInRealCluster));
  // Must do moment analysis before catagorization
  for (ClusterIter i = mClusters.begin(); i != mClusters.end(); ++i) {
    i->FindClusterAxis(pTowerUtil->GetMomentEcutoff());
  }  // for
  // Loop over clusters, catagorize, guess the photon locations for cat 0 or 2
  // clusters then fit, compare, and choose the best fit
  junkyEvent = false;  // Bad event?
  for (ClusterIter cluster = mClusters.begin(); cluster != mClusters.end();
       ++cluster) {
    Int_t clustCatag = pTowerUtil->CatagBySigmXY(&(*cluster));
    // point to the real TObjArray that contains the towers to be fitted
    // it is the same tower array for the cluster or all alternative clusters
    fitter->tow2Fit = cluster->towers();
    // Number of Degree of Freedom for the fit
    if (clustCatag == k1PhotonCluster) {
      // Do 1-photon fit
      FitOnePhoton(&(*cluster));
    } else if (clustCatag == k2PhotonCluster) {
      // Do 2-photon fit
      Fit2PhotonClust(cluster);
      junkyEvent = cluster->chiSquare() > MaxChi2Catag2;
    } else if (clustCatag == kAmbiguousCluster) {
      // for catagory-0 cluster, first try 1-photon fit!
      // If the fit is good enough, it is 1-photon. Else also
      // try 2-photon fit, and find the best fit (including 1-photon fit).
      Bool_t is2Photon = true;
      double chiSq1 = FitOnePhoton(&(*cluster));
      const PhotonHitFPD photon = cluster->photons()[0];  // Cache the photon
      double chiSq2(NAN);  // Only set if do 2-photon fit
      // Decide if this 1-photon fit is good enough
      if (chiSq1 < maxGood1PhChi2NDF) {
        is2Photon = false ;
      } else {
        // The 1-photon fit isn't good enough, so try 2-photon fit
        chiSq2 = Fit2PhotonClust(cluster);
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
        cluster->SetNphoton(2);
        cluster->setChiSquare(chiSq2);
        // Flag the event as bad if the fit chi2/ndf is too bad
        junkyEvent = cluster->chiSquare() > MaxChi2Catag2;
      } else {
        // 1-photon fit is better
        cluster->SetNphoton(1);
        cluster->setChiSquare(chiSq1);
        cluster->photons()[0] = photon;
      }  // if (is2Photon)
    } else {  // Invalid cluster category
      // should not happen!
      std::cerr << "Your logic of catagory is wrong! Something impossible " <<
        "happens! This a catagory-" << clustCatag;
      std::cerr << " clusters! Don't know how to fit it!\n" << "\n";
    }  // if (clustCatag...)
  }  // Loop over all real clusters
  Int_t nPh = std::accumulate(mClusters.begin(), mClusters.end(), 0,
                              accumulatePhotons);
  if(nPh > FitTower::MAX_NUMB_PHOTONS) {
    // myFitter can only do up to "MAX_NUMB_PHOTONS"-photon fit
    std::cout << "Can not fit " << nPh << " (more than " <<
      FitTower::MAX_NUMB_PHOTONS << " photons!" << "\n";
    return nPh;
  }  // if
  // For global fit, add all towers from all clusters
  TObjArray allTow(nTows);
  Int_t ndfg = 0 ;
  for (ClusterIter cluster = mClusters.begin(); cluster != mClusters.end();
       ++cluster) {
    allTow.AddAll(cluster->towers());
    ndfg += (cluster->towers()->GetEntriesFast() - 3 * cluster->GetNphoton());
  }  // for
  if(ndfg <= 0) {
    ndfg = 1;
  }  // if
  fitter->tow2Fit = &allTow;
  // Only do global fit for 2 or more clusters (2-photon fit for one cluster
  // already has global fit)
  if (mClusters.size() > 1) {
    double chiSqG = GlobalFit(nPh, mClusters.size(), mClusters.begin()) / ndfg;
    // Check for errors in the global fit - the number of photons returned by
    // the global fit should equal the sum of photons in the fitted clusters
    Int_t iph = std::accumulate(mClusters.begin(), mClusters.end(), 0,
                                accumulatePhotons);
    if(iph != nPh) {
      std::cerr << "ERROR total nPh=" << nPh << " iPh=" << iph << std::endl;
    }  // if
  } else if (mClusters.size() == 1) {
      double chiSqG = mClusters.front().chiSquare();
  } else {
    double chiSqG = -1 ;
  }  // if (nRealClusts > 1)
  nRealClusts = mClusters.size();
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
  for(Int_t it=0; it<p_clust->GetNTower(); it++) {
    eSS += EnergyInTowerByPhoton(widthLG, (TowerFPD*)p_clust->towers()->At(it),
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
  Double_t xx = ((Double_t)p_tower->column() - 0.5) * widLG[0] - p_photon->xPos;
  Double_t yy = ((Double_t)p_tower->row() - 0.5) * widLG[1] - p_photon->yPos;
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
    if (i->hit()->energy() > 0.001) {
      tow_Arr->Add(&(*i));
    }  // if
  }  // for
  TIter next(tow_Arr);
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

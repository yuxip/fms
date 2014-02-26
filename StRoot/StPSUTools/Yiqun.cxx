#include "Yiqun.h"
#include <TRandom.h>  // For ROOT global random generator, gRandom
using namespace std;
using namespace PSUGlobals;
#include <TCanvas.h>
#include <algorithm>
#include <functional>
namespace {
// Lack of c++0x in STAR-standard gcc as of writing so we provide our own
// implementation of copy_if - http://en.cppreference.com/w/cpp/algorithm/copy
template<class InputIterator, class OutputIterator, class UnaryPredicate>
OutputIterator copy_if(InputIterator start, InputIterator end, 
                       OutputIterator outStart, UnaryPredicate predicate) {
  return std::remove_copy_if(start, end, outStart, std::not1(predicate));
}
// Unary predicate for selecting good clusters.
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
}  // unnamed namespace
#include <list>
ClassImp(Yiqun)


// the real pointers are declared here
//
//-----------------------------static

TF1* Yiqun::EDepCorrection(NULL);

Float_t Yiqun::FitOnePhoton(HitCluster* p_clust)
{
  // 4 parameters are passed to the fitting routine: nPhotons, x, y, energy
  // Starting points for the fit parameters:
  const Double_t start[4] = {
    1.0, widLG[0] * p_clust->x0, widLG[1] * p_clust->y0, p_clust->energy};
  // Maximum allowed deviations from the start points:
  const Double_t delta[4] = {
    0.5, 0.5 * widLG[0], 0.5 * widLG[1], 0.15 * p_clust->energy};
  // Lower and upper physical limits of fit parameters = start +/- delta
  Double_t lowLim[4], upLim[4];
  for (int i(0); i < 4; ++i) {
    lowLim[i] = start[i] - delta[i];
    upLim[i] = start[i] + delta[i];
  }  // for
	Int_t status = fitter->Fit(start, fitter->step, lowLim, upLim);
	// check return status
	if (status != 0) {
		std::cerr << "Minuit fit returns " << status << "!" << std::endl;
	}  // if
	// fit parameters(starting positions), errors, and gradients of function
	Double_t param[4];
	Double_t error[4];
	// get the fit result
	fitter->fMn->GetParameter(0, param[0], error[0]);
	// 3 parameters per photon, and the 1st parameter giving the number of photons
	Int_t nPar = 3*((Int_t) param[0])+1;
	for (Int_t ipar = 1; ipar < nPar; ipar++) {
		fitter->fMn->GetParameter(ipar, param[ipar], error[ipar]);
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
	return chiSq;
}

// 2003-10-10
// Global fit, return x (or y) positions (or errors) in "cm"!
//
// fit "nPh" of photons, using "*fitter" on "nCl" of clusters "*clust"
//
// Check "nCl>=1 && nPh>=2", this is a global fit. Restrict the range accordingly.
//

Float_t Yiqun::GlobalFit(const Int_t nPh, const Int_t nCl, HitCluster *p_clust)
{


	// by design, we can only fit up to "MAX_NUMB_PHOTONS" (currently 4) photons
	//
	if( nPh > FitTower::MAX_NUMB_PHOTONS || nPh < 2 ) {
	  std::cout << "Global fit! Can not fit " << nPh << " photons! ERROR!" << "\n";
		return -9999;
	}

	// check that there is at least one cluster
	//
	if( nCl < 1 ) {
		std::cout << nCl << " clusters! Global fit will NOT work! ERROR!" << "\n";
		return -9999;
	}


	// fit parameters(starting positions), errors, and gradients of function
	//
	Double_t param[3*FitTower::MAX_NUMB_PHOTONS+1];
	Double_t error[3*FitTower::MAX_NUMB_PHOTONS+1];
	Double_t gradient[3*FitTower::MAX_NUMB_PHOTONS+1];

	// starting position, lower and upper limit of parameters
	//
	Double_t start[3*FitTower::MAX_NUMB_PHOTONS+1], lowLim[3*FitTower::MAX_NUMB_PHOTONS+1], upLim[3*FitTower::MAX_NUMB_PHOTONS+1];


	// at least 2 photons: global fit!
	// The positions (such as "clust[ic].photon[jp].xPos") are already in unit of "cm"!!!
	//
	// here all clusters have all their fields properly filled (for example cluster[].photon[0] should NOT be NULL!)
	//
	Int_t totPh;
	totPh = 0;

	// loop over all clusters
	//
	for(Int_t ic=0; ic< nCl; ic++) {

		// loop over all photons in cluster
		//
		for(Int_t jp=0; jp<p_clust[ic].nPhoton; jp++)
		  {
		    
		    if( totPh > FitTower::MAX_NUMB_PHOTONS ) {
		      std::cout << "Total # of photons in " << nCl << " clusters is at least " << totPh << "! I can NOT do fit! ERROR!" << "\n";
		      return -9999;
		    }
		    
		    Int_t kpar;
		    kpar = 3 * totPh + 1 ;
		    
		    start[ kpar] = p_clust[ic].photon[jp].xPos;
		    lowLim[kpar] = start[kpar] - posDif_Gl ;        // in "cm", not in unit of lgd!
		    upLim[ kpar] = start[kpar] + posDif_Gl ;        // in "cm", not in unit of lgd!
		    
		    kpar++ ;
		    start[ kpar] = p_clust[ic].photon[jp].yPos;
		    lowLim[kpar] = start[kpar] - posDif_Gl ;        // in "cm", not in unit of lgd!
		    upLim[ kpar] = start[kpar] + posDif_Gl ;        // in "cm", not in unit of lgd!
		    
		    kpar++ ;
		    start[ kpar] = p_clust[ic].photon[jp].energy;
		    lowLim[kpar] = start[kpar] * (1 - eneRat_Gl) ;
		    upLim[ kpar] = start[kpar] * (1 + eneRat_Gl) ;

		    totPh++ ;
		  }			
	}
	
	if( totPh != nPh ) {
		std::cout << "WARNING! Total # of photons in " << nCl << " clusters is at least " << totPh << "! Not the same as the nPh = ";
		std::cout  << nPh << "! I will try " << totPh << " instead!" << "\n";
	}
	start[0] = totPh ;

	lowLim[0] = 0.5;
	lowLim[0] = FitTower::MAX_NUMB_PHOTONS + 0.5 ;
	// fit status, and flag
	//
	Int_t status, iflag=1;

	// call fitter->Fit(...), return status
	//
	status = fitter->Fit(start, fitter->step, lowLim, upLim);


	// 2003-10-09
	// check return status
	//
	if( status != 0 ) {
		std::cout << "Minuit fit returns " << status << "!" << "\n";
	}


	// get the fit result
	//
	fitter->fMn->GetParameter(0, param[0], error[0]);
	Int_t nPar = 3*((Int_t) param[0])+1;
	for(Int_t ipar=1; ipar<nPar; ipar++) {
		fitter->fMn->GetParameter(ipar, param[ipar], error[ipar]);
	}


	// put the fit result back in "clust"
	//
	// loop over all clusters
	//
	Int_t tPh = 0 ;
	for(Int_t ic=0; ic< nCl; ic++) {

		// loop over all photons in cluster
		//
		for(Int_t jp=0; jp<p_clust[ic].nPhoton; jp++) {

			Int_t kpar;
			kpar = 3 * tPh + 1 ;

			p_clust[ic].photon[jp].xPos    = param[kpar] ;
			p_clust[ic].photon[jp].errXPos = error[kpar] ;

			p_clust[ic].photon[jp].yPos    = param[kpar+1] ;
			p_clust[ic].photon[jp].errYPos = error[kpar+1] ;

			p_clust[ic].photon[jp].energy  = param[kpar+2] ;
			p_clust[ic].photon[jp].errEne  = error[kpar+2] ;

			tPh ++;
		}

	}

	// evaluate the Chi-square function
	//
	Double_t chiSq;
	fitter->fMn->Eval(7, gradient, chiSq, param, iflag);

	return chiSq;
	
}



// special initial step for 2-photon-cluster
//

// 2003-08-30
// special 2-photon cluster fitting
//
Float_t Yiqun::Fit2PhotonClust(HitCluster* p_clust)
{
  const Double_t step2[7] = {0, 0.02, 0.02, 0.01, 0.01, 0.01, 0.1} ;
  // sigma of cluster
  //
  // 	Double_t sigma;
  // 	sigma  = clust->sigmaX * clust->sigmaX ;
  // 	sigma += clust->sigmaY * clust->sigmaY ;
  // 	sigma = sqrt(sigma);
  // ratio of sigmaMin/sigmaMax
  //
  Double_t ratioSigma = p_clust->sigmaMin / p_clust->sigmaMax ;
  
  // possible maximum angle
  //
  // 	Double_t maxTheta = asin( thetaPara * ratioSigma) ;
  //
  // 2003-10-11
  // tigher limits
  //
  Double_t maxTheta = ratioSigma / thetaPara ;
  if( maxTheta > (TMath::Pi()/2.0) ) {
    maxTheta = TMath::Pi() / 2.0 ;
  }
  // use for restrict d_gg
  //
  Double_t EcSigmaMax = p_clust->energy * p_clust->sigmaMax ;
  
  
  Double_t chiSq;
  
  // fit parameters(starting positions), errors, and gradients of function
  //
  Double_t param[7];
  Double_t error[7];
  Double_t gradient[7];
  
  // starting position, lower and upper limit of parameters
  //
  Double_t start[7], lowLim[7], upLim[7];
  
  // fit status, and flag
  //
  Int_t status, iflag=1;
  
  
  // constant parameter: 2 photon
  //
  start[0] = 2;
  lowLim[0] = 1.5;
  lowLim[0] = 2.5;
  
  
  // parameter starting points and limits are determined by looking at cluster information
  //
  //  xPi and yPi:   rarely do they go beyond 0.3 unit of lgd
  //  theta:         have a narrow theta range (for r=sigmaMin/sigmaMax, |theta|<0.5*r/0.65
  //                      when r<0.65, and linear increase from 0.5 to Pi/2 for 0.65<r<1)
  //  E_gg:          given by Ec (+/- 20% or less)
  //  z_gg:          should just let it vary from -1 to 1.
  //  d_gg:          a lower bound is given by r=sqrt(sigmaX^2+sigmaY^2). 
  //                      d_gg > Max( 2.5*(r-0.6), 0.5 )
  //
  start[1]  = widLG[0] * p_clust->x0;
  start[2]  = widLG[1] * p_clust->y0;
  start[6]  = p_clust->energy;
  
  start[4]  = p_clust->thetaAxis;
  
  // 2003-10-10
  // d_gg should be more closely related to sigma_max!
  //
  // 2003-10-09
  // Change start[3] (related to  2 * sigma)
  //
  // 	Double_t dggPi0Peak;
  // 	dggPi0Peak = towerWidthCM * sigma * 2.0;
  //
  start[3] = dggPara[1] * widLG[0] * p_clust->sigmaMax ;
  
  
  // randomize the starting point of Z_gg (from -0.1 to 0.1)
  //
  start[5]  = 0.1 * (2 * gRandom->Rndm() - 1) ;
  
    lowLim[1] = start[1] - posDif_2PC * widLG[0] ;
    lowLim[2] = start[2] - posDif_2PC * widLG[1] ;
    lowLim[6] = start[6] * (1 - eneRat_2PC) ;
    
    upLim[1]  = start[1] + posDif_2PC * widLG[0] ;
    upLim[2]  = start[2] + posDif_2PC * widLG[1] ;
    upLim[6]  = start[6] * (1 + eneRat_2PC) ;
    
    lowLim[4] = start[4] - maxTheta ;
    lowLim[5] = - 1.0 ;
    
    // 		// old limit
    // 		//
    // 		lowLim[3] = dggPara[1] * (dggPara[0] - EcSigmaMax ) ;
    //
    // 2003-10-12
    // new limit
    //
    lowLim[3] = dggPara[0] / pow(EcSigmaMax, 0.8) ;
    if( lowLim[3] < dggPara[2] ) {
      lowLim[3] = dggPara[2] ;
    }
    lowLim[3] *= widLG[0];
    if( lowLim[3] >= start[3] ) {
      lowLim[3] = start[3] * 0.9 ;
    }
    
    upLim[3] = dggPara[4] * (dggPara[3] - EcSigmaMax ) ;
    if( upLim[3] > dggPara[5] ) {
      upLim[3] = dggPara[5] ;
    }
    upLim[3] *= widLG[0] ;
    if( upLim[3]  <= start[3] ) {
      upLim[3] = start[3] * 1.1 ;
    }
    
    upLim[4]  = start[4] + maxTheta ;
    upLim[5]  =   1.0 ;
  
  
  // 	// 2003-10-13
  // 	// fix E_total and theta
  // 	//
  // 	fitter->fMn->FixParameter(6);
  // 	fitter->fMn->FixParameter(4);
  
  
  // call special 2-photon-cluster fitter
  //
  status = fitter->Fit2Pin1Clust(start, step2, lowLim, upLim);
  
  // 2003-10-09
  // check return status
  //
  if( status != 0 ) {
    std::cout << "Minuit fit returns " << status << "!" << "\n";
  }
  
  
  // get the fit result
  //
  fitter->fMn->GetParameter(0, param[0], error[0]);
  Int_t nPar = 3*((Int_t) param[0])+1;
  for(Int_t ipar=1; ipar<nPar; ipar++) {
    fitter->fMn->GetParameter(ipar, param[ipar], error[ipar]);
  }
  //  float mm=param[6]*sqrt(1-param[5]*param[5])*param[3]/730./2.;
  // printf("mm=%f fitx=%f fity=%f dgg=%f zgg=%f egg=%f \n",mm,
  //	 param[1],param[2],param[3],param[5],param[6]);

  // evaluate the Chi-square function
  //
  fitter->fMn->Eval(7, gradient, chiSq, param, iflag);

  // 	// 2003-10-13
  // 	// 2nd fit: use result of 1st fit as starting point, and free all fixed parameter
  // 	//
  // 	for(Int_t ipar=1; ipar<nPar; ipar++) {
  // 		start[ipar] = param[ipar];
  // 	}
  // 	fitter->fMn->mnfree(0);
  
  // 	// call special 2-photon-cluster fitter
  // 	//
  // 	status = fitter->Fit2Pin1Clust(start, step2, lowLim, upLim);
  
  
  // 2003-10-09
  // check return status
  //
  if( status != 0 ) {
    std::cout << "Minuit fit returns " << status << "!" << "\n";
  }
  
  
  // 	// get the fit result
  // 	//
  // 	fitter->fMn->GetParameter(0, param[0], error[0]);
  // 	for(Int_t ipar=1; ipar<nPar; ipar++) {
  // 		fitter->fMn->GetParameter(ipar, param[ipar], error[ipar]);
  // 	}
  
  // 	// evaluate the Chi-square function
  // 	//
  // 	fitter->fMn->Eval(7, gradient, chiSq, param, iflag);
  
  
  // put the fit result back in "clust"
  //
  p_clust->photon[0].xPos    = param[1] + cos(param[4]) * param[3] * (1 - param[5]) / 2.0 ;
  p_clust->photon[0].errXPos = error[1] + (cos(param[4])*error[3]-error[4]*sin(param[4])*param[3])*(1-param[5])/2 - cos(param[4])*param[3]*error[5]/2.0;
  p_clust->photon[0].yPos    = param[2] + sin(param[4]) * param[3] * (1 - param[5]) / 2.0 ;
  p_clust->photon[0].errYPos = error[2] + (sin(param[4])*error[3]+error[4]*cos(param[4])*param[3])*(1-param[5])/2 - sin(param[4])*param[3]*error[5]/2.0;
  p_clust->photon[0].energy  = param[6] * (1 + param[5]) / 2.0 ;
  p_clust->photon[0].errEne  = error[6]*(1+param[5])/2.0 + param[6]*error[5]/2.0;
  
  p_clust->photon[1].xPos    = param[1] - cos(param[4]) * param[3] * (1 + param[5]) / 2.0 ;
  p_clust->photon[1].errXPos = error[1] + (-cos(param[4])*error[3]+error[4]*sin(param[4])*param[3])*(1+param[5])/2 - cos(param[4])*param[3]*error[5]/2.0;
  p_clust->photon[1].yPos    = param[2] - sin(param[4]) * param[3] * (1 + param[5]) / 2.0 ;
  p_clust->photon[1].errYPos = error[2] + (sin(param[4])*error[3]-error[4]*cos(param[4])*param[3])*(1+param[5])/2 - sin(param[4])*param[3]*error[5]/2.0;
  p_clust->photon[1].energy  = param[6] * (1 - param[5]) / 2.0 ;
  p_clust->photon[1].errEne  = error[6]*(1-param[5])/2.0 - param[6]*error[5]/2.0;
  
  
  
  
  // 2003-10-13
  // do a global fit, using result of 1st fit as starting point
  //
  // need to set "nPhoton" before calling "GlobalFit(..)"
  //
  p_clust->nPhoton = 2;
  chiSq = GlobalFit(2, 1, p_clust);
    
  return chiSq;
};


Int_t Yiqun::FitEvent(Int_t nTows, Int_t &nClusts, Int_t &nRealClusts, Bool_t &junkyEvent)
{
  
	// possible alternative clusters for 1-photon fit: for catagory 0
	//

	HitCluster altClu ;
	

	// number of photons found by fit
	//
	Int_t nPh = nClusts = 0 ;

	TObjArray arrTow(nTows);

	// at the beginning of each event, first need to clear the TObjArray
	//

	arrTow.Clear();
	for(Int_t ic=0; ic<MAX_NUMER_CLUSTERS; ic++) {
	  if(clust[ic].tow)clust[ic].tow->Clear();
		clust[ic].Clear();
	}
	// re-add the TowerFPD objects to "arrTow" for each event!
	//
	for(Int_t itow=0; itow<nTows; itow++) {
	  arrTow.AddAt(&(towers->at(itow)), itow);
	}

	// find clusters
	//
	// "TObjArray::Sort()" sort array from small to large ( [0]<[1]<...<[48] )
	//
	// Need to first mark "TObjArray *arrTow" as unsorted! Then sort!
	//
	// calling "TObjArray::GetAbsLast()" will restore "fLast" to a value other than -2 (-2 is setho after "Sort()" is called) ?!
	//
	arrTow.UnSort();
	arrTow.Sort();
	nClusts = pTowerUtil->FindTowerCluster(&arrTow, clust);
			       
  std::list<HitCluster*> pClusters;
  for (int i(0); i < nClusts; ++i) {
    pClusters.push_back(&clust[i]);
  }  // for	
	std::list<HitCluster*> myClusters;
	copy_if(pClusters.begin(), pClusters.end(), std::back_inserter(myClusters),
	        IsGoodCluster(minRealClusterEne, maxHitsInRealCluster));

	// 2003-09-12
	//
	// Cluster energy should be at least 2 GeV (parameter "minRealClusterEne")
	//
	// 2003-09-16
	// New algorithm to move all clusters above "minRealClusterEne" to front, and
	//   clusters below "minRealClusterEne" to the end. This way, even if a low
	//   energy cluster were before a high energy one, we would not miscount!
	//
	// "ifront" movers from the 1st clusters, "jend" moves from the last cluster
	// 2008-03-19 SH generalizes to "move to end" clusters with too low energy OR to many hits
	Int_t ifront = 0;
	Int_t jend   = nClusts - 1 ;
	while(1) {

		// break out the loop when "ifront" meets "jend"
		//
		if( ifront > jend ) {
			break;
		}

		if( clust[ifront].energy >= minRealClusterEne 
		    && clust[ifront].tow->GetEntries()<= maxHitsInRealCluster 
		    ) {
		  ifront ++;
		}
		else {
			// search from "jend" and forward, until find a clusters that is higher than "minRealClusterEne"
			//
			while(jend>=ifront && ( clust[jend].energy < minRealClusterEne 
			       || clust[jend].tow->GetEntries()> maxHitsInRealCluster )) {
				jend --;
			}
			//
			// check that "jend" is in front of "iend", we have exhausted the array. Stop!
			//
			if( jend < ifront ) {
				break;
			}
			//
			// exchange the clusters: "ifront"<-->"jend", then increase "ifront" by 1 and decrease "jend" by 1;
			//
			else {
				altClu        = clust[ifront] ;
				clust[ifront] = clust[jend]   ;
				clust[jend]   = altClu        ;
				ifront ++;
				jend   --;
			}
		}
	}
	//
	// at the end of the above code, "ifront" counts the number of clusters above "minRealClusterEne"
	//
	nRealClusts = ifront ;

	if (nRealClusts != int(myClusters.size())) {
	  std::cerr << "Ya blew it!" << std::endl;
	  exit(-1);
	} else if (nRealClusts != nClusts) {
	  std::cerr << "huzzah! Steve's algorithm gives " << nRealClusts <<
	    ", mine gives " << myClusters.size() << " (out of " << nClusts <<
	    " total clusters)" << std::endl;
	}  // if
	// do mement analysis before catagozie!
	//
	//cout<<"test nRealClusts="<<nRealClusts<<endl;
	for(Int_t iic=0; iic<nRealClusts; iic++) {

	  clust[iic].FindClusterAxis(pTowerUtil->GetMomentEcutoff());
	}


	// loop over clusters, catagorize, guess the photon locations for cat 0 or 2 clusters
	// then fit, compare, and choose the best fit
	//
	Int_t clustCatag;

	// bad events?
	//
	junkyEvent = false;

	
	for(Int_t icc=0; icc<nRealClusts; icc++) {

		clustCatag = pTowerUtil->CatagBySigmXY( &clust[icc] );

		// point to the real TObjArray that contains the towers to be fitted
		// it is the same tower array for the cluster or all alternative clusters
		//
		fitter->tow2Fit = clust[icc].tow ;

		// store the goodness of fit
		//
		Double_t chiSq1, chiSq2;

		// Number of Degree of Freedom for the fit
		//
		Int_t ndf;
		ndf = clust[icc].tow->GetEntriesFast();

		// for catag-1 cluster, do 1-photon fit
		if( clustCatag == 1 ) {
			clust[icc].nPhoton   = 1 ;
			chiSq1 = FitOnePhoton(&clust[icc]);

			photons[nPh] = clust[icc].photon[0] ;

			ndf -= 3 ;
			if( ndf<=0 )
				ndf = 1;

			clust[icc].chiSquare = chiSq1 / ndf ;

			nPh ++ ;

		}
		//
		// for catagory-2 cluster, do 2-photon fit!
		// If the fit is good enough, it is 1-photon. Else also
		// try 2-photon fit, and find the best fit (including 1-photon fit).
		//
		else if( clustCatag == 2 ) {
			clust[icc].nPhoton   = 2 ;
			chiSq2 = Fit2PhotonClust(&clust[icc]);

			ndf -= 6 ;
			if( ndf<=0 ) 
				ndf = 1;

			clust[icc].chiSquare = chiSq2 / ndf ;

			if( clust[icc].chiSquare > MaxChi2Catag2 ) {
				junkyEvent = true;
			}

			// we need to fill in the information anyway
			//
			clust[icc].nPhoton   = 2 ;
			for(Int_t kp=0; kp<2; kp++) 
			  {
			    
			    photons[nPh] = clust[icc].photon[kp] ;

				nPh ++ ;
			}
		}
		//
		// for catagory-0 cluster, first try 1-photon fit!
		// If the fit is good enough, it is 1-photon. Else also
		// try 2-photon fit, and find the best fit (including 1-photon fit).
		//
		else if( clustCatag == 0 ) {
			// 2-photon or 1-photon
			//
			Bool_t is2Photon;
			is2Photon = true;

			altClu = clust[icc] ;
			altClu.nPhoton   = 1 ;
			chiSq1 = FitOnePhoton(&altClu);
			

			Int_t NDF1, NDF2;

			NDF1= ndf - 3 ;
			if( NDF1 <= 0 )
				NDF1 = 1;

			NDF2= ndf - 6 ;
			if( NDF2 <= 0 )
				NDF2 = 1;

			// decide if this 1-photon fit is good enough
			//
			if( chiSq1 < maxGood1PhChi2NDF * NDF1 ) {
				is2Photon = false ;
			}
			//
			// else try 2-photon fit
			//
			else {

				// 2003-09-15
				//
       // New logic for deciding if a 2-photon fit is good!
       // If one photon peak lies on top of low (compared to photon energy) or even zero tower,
       //    this photon is definitely bogus. This could happen if a nearby cluster took away
       //    towers that might make the fit having a large Chi-square (deviation).
       //
       // So for 2-photon fit, first check that the fitted photon of lower energy (higher energy
       //    photon should be fine) is over one of the non-zero tower of the cluster. First of all,
       //    this ensures that we don't have an "outside of cluster" bogus photon, i.e. a bogus
       //    photon could be the result of minimizing the chi-square over towers that do not include
       //    the supposed peak tower. 

			  chiSq2 = Fit2PhotonClust(&clust[icc]);
       // check that the fitted photon of lower energy is contained within one of the non-zero towers
       //    of the cluster
       //
			  Bool_t isBogus2ndPhoton;
			  isBogus2ndPhoton = true;
			  Int_t secondPhoton;
			  if( clust[icc].photon[0].energy < clust[icc].photon[1].energy ) 
			    {
			      secondPhoton = 0;
			    }
			  else 
			    {
			      secondPhoton = 1;
			    }
			  
			  // tower where the fitted photon of lower energy should hit
			  //
			  Int_t ix, iy;
			  ix = 1 + (Int_t) (clust[icc].photon[secondPhoton].xPos / widLG[0]) ;
			  iy = 1 + (Int_t) (clust[icc].photon[secondPhoton].yPos / widLG[1]) ;
			  TowerFPD tow2ndPhoton(clust[icc].photon[secondPhoton].energy, ix, iy, icc) ;
			  
       // now check whether this tower is one of the non-zero towers of the cluster
       //
			  for(Int_t itNZ=0; itNZ<clust[icc].numbTower; itNZ++) {
			    TowerFPD * nZTow = (TowerFPD *) clust[icc].tow->At(itNZ) ;
			    if( tow2ndPhoton.IsEqual( nZTow ) ) {
			      isBogus2ndPhoton = false;
			      tow2ndPhoton.energy = nZTow->energy ;
			      
			      break;
			    }
			  }
			  
			  // and the fitted energy is too large compared to the energy of the tower
			  //
			  if( !isBogus2ndPhoton ) {
			    if( tow2ndPhoton.energy < minHTEneOverPhoton * clust[icc].photon[secondPhoton].energy ) {
			      isBogus2ndPhoton = true;
			    }
			  }
			  
	// 2003-09-29
        // Check that the 2nd photon's "High-Tower" enery is too large compared to its fitted energy.
	//    If so, it is probably splitting one photon into two!
	//
			  if( !isBogus2ndPhoton ) 
			    {
			      Double_t eSS;

			      eSS = EnergyInTowerByPhoton(widLG[0], &tow2ndPhoton, 
							  &( clust[icc].photon[secondPhoton] ) );
			      if( tow2ndPhoton.energy > maxHTEneOverPhoton * eSS ) 
				{
				  isBogus2ndPhoton = true;
				}
			    }
			  
			  // 2003-09-29
			  // Check that the 2nd photon is not near the edge of another cluster.
			  // Namely, we check what would be the energy deposited in other clusters by this photon vs
			  //    energy deposited in its own cluster.
	// If the ratio is too high, this fitted photon is probably a bogus one!
	//

				if( !isBogus2ndPhoton ) {
					Double_t eneWInCluster = EnergyInClusterByPhoton(widLG[0], 
					       &clust[icc], &( clust[icc].photon[secondPhoton] ));

	// loop over all clusters, but skip its own cluster
	//
					for(Int_t kkc=0; kkc<nRealClusts; kkc++) 
					  {
					    if( kkc == icc )
					      continue;
					    
					    if( EnergyInClusterByPhoton(widLG[0], &clust[kkc], 
									&(clust[icc].photon[secondPhoton])) > 
						(maxRatioSpill * eneWInCluster) ) 
					      {
						isBogus2ndPhoton = true;
						break ;
					      }
					  }
				}

				
				if( isBogus2ndPhoton ) 
				  {
				    is2Photon = false ;
				  }
				else 
				  {
				    
				    // 2003-09-12
				    // Need to change!
				    //
				    // What is the criteria of choosing 2-photon fit over 1-photon fit?
				    // How about "chiSq2/NDF > chiSq1/NDF"?
				    //
				    // But there are cases when NDF<=0 !!!
				    //
				    // for catagory 0 cluster and when 1-photon fit is better
				    //
				    if( (chiSq1/NDF1) < (chiSq2/NDF2) ) {
				      is2Photon = false ;
				    }
				    else {
				      is2Photon = true ;
				    }
				    
				  } // isBogus2ndPhoton == true
				
			} // try 2-photon fit
			
			// now fill in the fit result
			//
			
			if( !is2Photon ) 
			  {
			    
			    // 1-photon fit
			    //			
			    clust[icc].nPhoton   = altClu.nPhoton   = 1 ;
			    clust[icc].chiSquare = altClu.chiSquare = chiSq1 / NDF1 ;
			    
			    photons[nPh] = clust[icc].photon[0] =  altClu.photon[0] ;
			    
			    nPh ++ ;
			  }
			else 
			  {
			    // 2-photon fit is better
			    //
			    clust[icc].nPhoton   = 2 ;
			    clust[icc].chiSquare = chiSq2 / NDF2 ;
			    
			    // check if this cluster has too much stuff in it!
			    //
			    if( clust[icc].chiSquare > MaxChi2Catag2 ) {
			      junkyEvent = true;
			    }
			    
			    
			    for(Int_t kp=0; kp<2; kp++) {
			      
			      photons[nPh] = clust[icc].photon[kp] ;
			      
			      nPh ++ ;
			    }
			} // else ( 2-photon fit is better)
			
		} // catag-0
		//
		// should not happen!
		//
		else {
		  std::cout << "Your logic of catagory is wrong! Something impossible happens! This a catagory-" << clustCatag;
		  std::cout << " clusters! Don't know how to fit it!\n" << "\n";
		}
		
	} // loop over all real clusters

	// 2003-10-04
	// Do a global fit for all events (even for 1 real cluster). Return
	//    Chi2/NDF instead of Chi2 in "FitEvent(...)".
	//
	// 2003-09-08
	// if there are more than 1 cluster, do a global fit!
	//
	if( nPh > FitTower::MAX_NUMB_PHOTONS ) 
	  {
	    //
	    // myFitter can only do up to "MAX_NUMB_PHOTONS"-photon fit
	    //
	    std::cout << "Can not fit " << nPh << " (more than " << FitTower::MAX_NUMB_PHOTONS << " photons!" << "\n";
	    return nPh;
	  }
	
	// 2003-09-08
	// for global fit, add all towers from all clusters
	//
	// first clear the TObjArray!
	//
	TObjArray allTow(nTows);
	allTow.Clear();
	Int_t ndfg = 0 ;
	for(Int_t jjc=0; jjc<nRealClusts; jjc++) 
	  {
	    allTow.AddAll(clust[jjc].tow);
	    ndfg += (clust[jjc].tow->GetEntriesFast() - 3 * clust[jjc].nPhoton);
	  }
	
	if( ndfg <= 0 ) 
	  {
	    ndfg = 1;
	  }
	
	fitter->tow2Fit = &allTow ;

	// 2003-10-13
	// only do global fit for 2 or more clusters (2-photon fit for one cluster already has global fit)
	//
	if(nRealClusts > 1 ) 
	  {
	    double chiSqG = GlobalFit(nPh, nRealClusts, clust) / ndfg ;
	    Int_t iph=0;
	    for(Int_t icl=0;icl<nRealClusts;icl++)
	      {
		
		for(Int_t iclph=0;iclph<clust[icl].nPhoton;iclph++)
		  {
		    
		    photons[iph] = clust[icl].photon[iclph];  
		    iph++;
		    if(iph>nPh){printf("ERROR total nPh=%d iph=%d \n",nPh,iph);break;};
		  }
	      };
	    if(iph!=nPh){printf("ERROR total nPh=%d iph=%d \n",nPh,iph);};

	  }
	else if( nRealClusts == 1 ) 
	  {
	    double chiSqG = clust[0].chiSquare ;
	  }
	else 
	  {
	    // temp revove comment
	    //std::cout << "Something wrong! Should NOT have " << nRealClusts << " real clusters! Error!" << "\n";
	    double chiSqG = -1 ;
	  };
	
	return nPh;
}



// Calculate the energy deposit in a cluster by a photon (from shower-shape function).
// In this case, we only need to consider non-zero towers.
//
Double_t Yiqun::EnergyInClusterByPhoton(Double_t widthLG, HitCluster *p_clust, PhotonHitFPD *p_photon)
{
	Double_t eSS = 0;
	for(Int_t it=0; it<p_clust->numbTower; it++) {
		eSS += EnergyInTowerByPhoton(widthLG, (TowerFPD *) p_clust->tow->At(it), p_photon);
	};

	return eSS ;
};


// Calculate the energy deposit in a tower by a photon (from shower-shape function).
//
Double_t Yiqun::EnergyInTowerByPhoton(Double_t widthLG, TowerFPD *p_tower, PhotonHitFPD* p_photon)
{
	Double_t xx, yy, eSS;
	//widthLG no longer used SH

	xx = ( (Double_t) p_tower->col - 0.5 ) * widLG[0] - p_photon->xPos ;
	yy = ( (Double_t) p_tower->row - 0.5 ) * widLG[1] - p_photon->yPos ;
	eSS = p_photon->energy * fitter->GetFunctShowShape()->Eval( xx, yy, 0 );

	return eSS;
};
Yiqun::Yiqun(TowerList* pEm,Geom* pgeom,Int_t iew,Int_t nstb)
{
  p_geom=pgeom;
  EW=iew;
  NSTB=nstb;
  Y(pEm);
};
void Yiqun::Y(TowerList* pEm)
{   
  towers = pEm;   
  NPh=0;
  NTower=0;
  NClusts=0;
  NRealClusts=0;

  pTowerUtil=new TowerUtil();
  pTowerUtil->SetMomentEcutoff(.5);  
  tow_Arr=0;
  widLG[0]=widLG[1]=(p_geom->FpdTowWid(EW,NSTB))[0];
  if(p_geom->FMSGeom)widLG[1]=(p_geom->FpdTowWid(EW,NSTB))[1];
  NTower=pEm->size();

  if(NTower>578)
    {
      printf("Too many towers for Fit\n");
      exit(-1);
    };

  Int_t cnt=-1;
  tow_Arr=new TObjArray(NTower);

  for (TowerList::iterator i = pEm->begin(); i != pEm->end(); ++i) {
    if (i->energy > 0.001) {
      tow_Arr->Add(&(*i));
    }  // if
  }  // for
    TIter next(tow_Arr);
    while(TowerFPD* tow=(TowerFPD*) next()){
      tow->SetContext(tow_Arr,EW,NSTB);
    }

  // choice of Chi-square function
  // 1:  Steve's
  // 2:  Larisa's
  //
  
  // clusters with just 1 photon
  //

  posDif_2PC=0.2;
  eneRat_2PC=0.05;
  dggPara[0]=18.0;
  dggPara[1]=2.2;
  dggPara[2]=0.5;
  dggPara[3]=60.0;
  dggPara[4]=0.085;
  dggPara[5]=3.5;
  
  thetaPara=2.8;
  
  // global fit (more than 1 cluster)
  //

  posDif_Gl=1.25;
  eneRat_Gl=0.3;
  maxGood1PhChi2NDF=5;
  minHTEneOverPhoton=0.25;
  maxHTEneOverPhoton=1.5;
  maxRatioSpill=0.2;
  MaxChi2Catag2=10.;
  minRealClusterEne=2.0;
  maxHitsInRealCluster=49;
  if(EW==2 && (NSTB==1 || NSTB==2))
    {
      maxHitsInRealCluster=10;
      minRealClusterEne=.75;
      //      maxHitsInRealCluster=4;
      maxHitsInRealCluster=25;
    };
  NClusts=0;
  NRealClusts=0;
  fitter=0;

  fitter = new FitTower(p_geom,EW,NSTB);

  fitter->SetTWidthCM(*(p_geom->FpdTowWid(EW,NSTB)));
  if(p_geom->FMSGeom)fitter->SetXYTWidthCM(p_geom->FpdTowWid(EW,NSTB));
  
  Bool_t badEvent(true);
  NPh = FitEvent(NTower, NClusts, NRealClusts, badEvent);
}

Yiqun::~Yiqun()
{
  fitter->GetFunctShowShape()->Draw("colz");
  gPad->SetLogz(true);
  gPad->Print("showerShape.eps");
  if(fitter!=0)delete fitter;

  if(pTowerUtil)delete pTowerUtil;
 
  if(tow_Arr)delete tow_Arr;
};

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

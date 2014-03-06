#include <iostream>
#include <stdlib.h>
#include <math.h>

#include "TMath.h"
#include "TVector2.h"

#include "StEvent/StFmsHit.h"

#include "TowerFPD.h"
#include "PhotonHitFPD.h"
#include "HitCluster.h"
#include "Yiqun.h"
using namespace std;
using namespace PSUGlobals;

ClassImp(HitCluster);



//HitCluster *clust2Fit;



HitCluster::HitCluster()
{

  tow=0;
  Ecutoff=.5;
  Clear();
  IsEUpdated=false;
}


HitCluster::~HitCluster()
{
  if(tow)delete tow;
};

void HitCluster::Clear(void)
{ 
  catag   = -1 ;
  numbTower = nPhoton =  0 ;
  energy  =  0 ;
  x0 = y0 = sigmaX = sigmaY = sigmaXY = chiSquare = sigmaMin = sigmaMax = -1 ;
  thetaAxis = -10 ;
  for(Int_t ip=0; ip<MAX_PHOTON_PER_CLUSTER; ip++) {
    photon[ip].Clear();
  }
  if(tow)delete tow;
  tow=new TObjArray(50);
}

void HitCluster::Print(void)
{
  if(tow)
    {
      
      for(Int_t i=0; i<tow->GetEntriesFast(); i++) 
	{
	  ( (TowerFPD *) tow->At(i) )->Print();
	}
    };
  
  cout << "catag=" << catag << "\t numbTower=" << numbTower << "\t nPhoton=" << nPhoton << "\t energy=" << energy ;
  cout << "\t (" << x0 << ", " << y0 << ")\t (" << sigmaX << ", " << sigmaY << ", " << sigmaXY << ")" << endl;
  
  cout << "thetaAxis=" << thetaAxis << ",\t sigmaMin=" << sigmaMin << ",\t sigmaMax=" << sigmaMax << endl;

  for(Int_t j=0; j<nPhoton; j++) {
    cout << "Photon #" << j << ":\n";
    photon[j].Print();
  }
  cout << "chiSquare = " << chiSquare << "\n" << endl;
  
}

void HitCluster::FindClusterAxis(void)
{
  Double_t dSigma2, aA, bB;
  dSigma2 = sigmaX * sigmaX - sigmaY * sigmaY ;
  aA = sqrt( dSigma2 * dSigma2 + 4.0 * sigmaXY * sigmaXY ) + dSigma2 ;
  bB = 2 * sigmaXY ;

	if( sigmaXY < 1e-10 ) {
		if( aA < 1e-10 ) {
			bB = sqrt( dSigma2 * dSigma2 + 4.0 * sigmaXY * sigmaXY ) - dSigma2 ;
			aA = 2 * sigmaXY ;
		}
	}


	thetaAxis = atan2(bB, aA);


	Double_t myPi = TMath::Pi() ; 
	while( thetaAxis > (myPi / 2.0) )
		thetaAxis -= myPi;

	while( thetaAxis < -(myPi / 2.0) )
		thetaAxis += myPi;

	sigmaMin = GetSigma(thetaAxis) ;
	sigmaMax = GetSigma(thetaAxis - TMath::Pi() / 2.0) ;
	//cout<<"sigmaMin="<<sigmaMin<<" sigmaMax="<<sigmaMax<<endl;
}



// calculate sigma w.r.t the axis going through the "center" and of an angle "theta" in x-y plane
//
Double_t HitCluster::GetSigma(Double_t theta){

	Double_t sigma = 0 ;

	// 2-d vector vaxis define the axis
	//
	TVector2 vaxis( cos(theta), sin(theta) );

	// loop over all towers pointer in cluster
	//
	TowerFPD * oneTower;
	float wnew =0;
	for(Int_t it=0; it<tow->GetEntriesFast(); it++) {

		oneTower = (TowerFPD *) tow->At(it);
		// the 2-d vector from the "center" of cluster to tower
		// "center" are at 0.5, 1.5, etc! Need shift of 0.5
		//
		TVector2 v1( oneTower->col - 0.5 - x0, oneTower->row - 0.5 - y0 );
		
		// perpendicular distance to the axis = length of the component of vector "v1" that is norm to "vaxis"
		//
		Double_t dis = (v1.Norm(vaxis)).Mod();
		
		// contribution to sigma
		//
		//sigma += oneTower->energy * dis * dis ;
		float wtmp = log(oneTower->hit->energy() + 1-Ecutoff)>0 ? log(oneTower->hit->energy() +1.-Ecutoff) : 0;
		wnew    += wtmp ;
		sigma += wtmp * dis * dis ;
	}

	//cout<<"here1 "<<(wnew>0? sqrt(sigma/wnew) : 0)<<endl;
	return wnew>0? sqrt(sigma/wnew) : 0;
}
void HitCluster::CalClusterMoment(Float_t Ecoff)
{
  Ecutoff=Ecoff;
  Float_t w0, w1, mtmp, mx, my, sigx, sigy, sigXY;
  w0 = w1 = mtmp = mx = my = sigx = sigy = sigXY = 0 ;
	
  TowerFPD * oneTow;
  TIter next(tow);
  while(oneTow=(TowerFPD*) next())
    {
      Float_t xxx, yyy;
      xxx = oneTow->col - 0.5 ;
      yyy = oneTow->row - 0.5 ;
      mtmp = log(oneTow->hit->energy()+1.-Ecoff)>0 ? log(oneTow->hit->energy()+1.-Ecoff) : 0;
      w1 += mtmp;
      w0    += oneTow->hit->energy() ;
      mx    += mtmp * xxx ;
      my    += mtmp * yyy ;
      sigx  += mtmp * xxx * xxx ;
      sigy  += mtmp * yyy * yyy ;
      sigXY += mtmp * xxx * yyy ;
    }
	
  energy  = w0 ;
  if(w1>0)
    {
      x0      = mx / w1 ;
      y0      = my / w1 ;
      sigmaX  = sqrt( fabs(sigx / w1 - x0 * x0 ) ) ;
      sigmaY  = sqrt( fabs(sigy / w1 - y0 * y0 ) ) ;
      sigmaXY = sigXY / w1 - x0 * y0 ;
    }
  else
    {
      x0      = 0 ;
      y0      = 0 ;
      sigmaX  = 0 ;
      sigmaY  = 0 ;
      sigmaXY = 0 ;
    }
  
};

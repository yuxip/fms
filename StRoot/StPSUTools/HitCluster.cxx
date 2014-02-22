#include <iostream>
#include <stdlib.h>
#include <math.h>

#include "TMath.h"
#include "TVector2.h"

#include "TowerFPD.h"
#include "PhotonHitFPD.h"
#include "HitCluster.h"
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

void HitCluster::EDepUpdate(WasExternal* pwe)
{
  if(IEW != 2 )return;
  if(pwe==0)return;
  Float_t Eratio=1;
  if(!IsEUpdated)
    {
	    if(NSTB<3)
	      {
		Eratio=pwe->EDepCorrection->Eval(20)/pwe->EDepCorrection->Eval(energy);
	      }
	    else
	      {
		Eratio=pwe->EDepCorrection->Eval(30.)/pwe->EDepCorrection->Eval(energy);
	      };
	energy=energy*Eratio;
	TIter next(tow);
	TowerFPD* hit;
	while(hit=(TowerFPD*) next())
	  {
	    hit->energy=hit->energy*Eratio;
	  };
    };
  IsEUpdated=true;
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


// 2003-10-13
// Use matrix algebra to find the axis. Basically, the 2nd moments form a symmetrical 2-D matrix. Diagonize it
//   and the two eigen-values are maximum moment and minimum moment, while two eigen-vectors are the vectors of
//   axis. To speed up the process, the solution is explicitly written in "C" (instead of using TMatrix).
// 
//
// do a 1-parameter fix (theta angle in x-y plane) to find the "least-2nd-sigma" axis of the cluster
//
void HitCluster::AddCluster(HitCluster* p_bclu)
{
  if(p_bclu)
    {
      TIter next(p_bclu->tow);
      TowerFPD* hit;
      while(hit=(TowerFPD*) next())
      {
	tow->Add(hit);
      };
      numbTower=tow->GetEntries();
      catag   = -1 ;
      numbTower = nPhoton =  0 ;
      energy  =  0 ;
      x0 = y0 = sigmaX = sigmaY = sigmaXY = chiSquare = sigmaMin = sigmaMax = -1 ;
      thetaAxis = -10 ;
      for(Int_t ip=0; ip<MAX_PHOTON_PER_CLUSTER; ip++) {
	photon[ip].Clear();
      }
      UpdateEnergy();
    }
};

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
		float wtmp = log(oneTower->energy + 1-Ecutoff)>0 ? log(oneTower->energy +1.-Ecutoff) : 0;
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
      mtmp = log(oneTow->energy+1.-Ecoff)>0 ? log(oneTow->energy+1.-Ecoff) : 0;
      w1 += mtmp;
      w0    += oneTow->energy ;
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

Float_t HitCluster::GetLastEnergySum(Int_t nclu)
{
  if(!EnergyUpdated)UpdateEnergy();
  if(tow==0)return 0;
  Int_t ntow=tow->GetEntries();
  if(ntow<=nclu)return energy;
  Float_t esum=0;
  for(int j=ntow;j>ntow-nclu;j--)
    {
      esum+=((TowerFPD*) tow->At(j-1))->energy;
    };
  return esum;
};
Bool_t HitCluster::UpdateEnergy()
{
  if(tow==0)return false;
  energy=0;
  //temporary change sign of energy for sorting
  TIter next(tow);
  while(TowerFPD* hit=(TowerFPD*) next())
    {
      hit->energy=-fabs(hit->energy);
    };
  tow->UnSort();
  tow->Sort();
  TIter next2(tow);
  while(TowerFPD* hit=(TowerFPD*) next2())
    {
      hit->energy=fabs(hit->energy);
      energy+=hit->energy;
    };
  EnergyUpdated=true;
  FindClusterAxis();
  
  return true;
};

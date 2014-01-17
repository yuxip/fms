////////////////////////////////////////////////////////////////////
//
//  definition of a HitCluster (cluster of FPD towers)
//
///////////////////////////////////////////////////////////////////

#ifndef HITCLUSTER_H
#define HITCLUSTER_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>

#include "TMath.h"
#include "TVector2.h"

#include "TObjArray.h"
#include "TMinuit.h"
#include "TVector3.h"



#include "PhotonHitFPD.h"
#include "WasExternal.h"
#define MAX_PHOTON_PER_CLUSTER 2

namespace PSUGlobals {//$NMSPC
class HitCluster : public TObject {

 public:

  PhotonHitFPD   photon[MAX_PHOTON_PER_CLUSTER];     // Photon-Hits in the cluster

  Int_t          catag  ;               // catagory of cluster (1: 1-photon,  2: 2-photon,  0: could be either 1- or 2-photon)
  Int_t          numbTower;             // number of non_zero towers in the cluster
  Int_t          nPhoton;               // number of photons contained in this cluster: ( nPhoton <= MAX_PHOTON_PER_CLUSTER ! )
  Float_t        energy ;               // total energy  contained in this cluster (0th moment)
  
  Float_t        x0     ;               // mean x ("center of gravity") in local grid coordinate (1st moment)
  Float_t        y0     ;               // mean y ("center of gravity") in local grid coordinate (1st moment)
  Float_t        sigmaX ;               // 2nd moment in x
  Float_t        sigmaY ;               // 2nd moment in y
  Float_t        sigmaXY  ;             // 2nd moment in x-y
  Float_t        thetaAxis;             // theta angle in x-y plane that define the direction of least-2nd-sigma axis
  Float_t        sigmaMin;              // 2nd sigma w.r.t to the least-2nd-sigma axis
  Float_t        sigmaMax;              // 2nd sigma w.r.t to the axis orthogonal to the least-2nd-sigma axis
  Float_t        chiSquare;             // Chi-square of the fitting
  
  TObjArray *    tow    ;               //!  TowerFPD objects that make the cluster
  void CalClusterMoment(Float_t Ecoff);
  Bool_t UpdateEnergy();
  Bool_t UpdateEnergy(Float_t Ecoff){Ecutoff=Ecoff;UpdateEnergy();};
  HitCluster() ;
  void AddCluster(HitCluster* p_bclu=0);
  ~HitCluster();
  Bool_t IsEUpdated;
  void EDepUpdate(WasExternal* pwe); 
  void Clear(void) ;
  
  HitCluster& operator=(const HitCluster& rhs) 
    {
      
      // HitCluster assignment operator
      //
      if (this != &rhs) 
	{
	  catag     = rhs.catag   ;
	  numbTower = rhs.numbTower;
	  nPhoton   = rhs.nPhoton ;
	  energy    = rhs.energy  ;
	  x0        = rhs.x0      ;
	  y0        = rhs.y0      ;
	  sigmaX    = rhs.sigmaX  ;
	  sigmaY    = rhs.sigmaY  ;
	  sigmaXY   = rhs.sigmaXY ;
	  thetaAxis = rhs.thetaAxis;
	  sigmaMin  = rhs.sigmaMin;
	  sigmaMax  = rhs.sigmaMax;
	  chiSquare = rhs.chiSquare;
	  if(tow)delete tow;
	  tow       = new TObjArray(*rhs.tow);
	  Ecutoff   = rhs.Ecutoff;
	  IEW       = rhs.IEW;
	  NSTB      = rhs.NSTB;
	  IsEUpdated= rhs.IsEUpdated;
	};
      
      for(Int_t i=0; i<nPhoton; i++) 
	{
	  photon[i] = rhs.photon[i] ;
	}
   
      return *this;
    };
  
  void Print(void);
  
  void FindClusterAxis(void);
  void FindClusterAxis(Float_t Ecoff){Ecutoff=Ecoff;FindClusterAxis();};
  
  Double_t GetSigma(Double_t theta);
  Float_t GetLastEnergySum(Int_t nclu);
  Int_t IEW;
  Int_t NSTB;
 private:
  Float_t Ecutoff;
  Bool_t EnergyUpdated;
  ClassDef (HitCluster, 7);
  
};
}
#endif

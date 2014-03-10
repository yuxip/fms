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
#define MAX_PHOTON_PER_CLUSTER 2

/*
 http://google-styleguide.googlecode.com/svn/trunk/cppguide.xml#Copy_Constructors
 A macro to disallow the copy constructor and operator= functions
 This should be used in the private: declarations for a class e.g.
  class Foo {
   public:
    Foo(int f);
    ~Foo();

   private:
    DISALLOW_COPY_AND_ASSIGN(Foo);
  };
*/
#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
  TypeName(const TypeName&);               \
  void operator=(const TypeName&)

namespace PSUGlobals {//$NMSPC
enum EClusterCategory {
  kAmbiguousCluster = 0,  // Could be 1- or 2-photon, needs to be fitted
  k1PhotonCluster = 1,
  k2PhotonCluster = 2,
  kInvalidClusterCategory
};
class HitCluster : public TObject {

 public:

  PhotonHitFPD   photon[MAX_PHOTON_PER_CLUSTER];     // Photon-Hits in the cluster

  Int_t          catag  ;               // catagory of cluster (1: 1-photon,  2: 2-photon,  0: could be either 1- or 2-photon)
  Int_t          numbTower;             // number of non_zero towers in the cluster
  Int_t          nPhoton;               // number of photons contained in this cluster: ( nPhoton <= MAX_PHOTON_PER_CLUSTER ! )
  Int_t          index;                 // cluster number in an event, counts from 0
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
  HitCluster() ;
  ~HitCluster();
  Bool_t IsEUpdated;
  void Clear(void) ;
  
  void Print(void);
  
  void FindClusterAxis(void);
  void FindClusterAxis(Float_t Ecoff){Ecutoff=Ecoff;FindClusterAxis();};
  
  Double_t GetSigma(Double_t theta);
  Float_t GetLastEnergySum(Int_t nclu);
  Int_t IEW;
  Int_t NSTB;
 private:
  DISALLOW_COPY_AND_ASSIGN(HitCluster);
  Float_t Ecutoff;
  Bool_t EnergyUpdated;
  ClassDef (HitCluster, 7);
  
};
}
#endif

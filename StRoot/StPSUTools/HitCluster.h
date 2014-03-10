#ifndef HITCLUSTER_H
#define HITCLUSTER_H

#include <TObject.h>

#include "StPSUTools/PhotonHitFPD.h"

class TObjArray;

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
  k1PhotonCluster = 1,  // A cluster created by 1 photon
  k2PhotonCluster = 2,  // A cluster created by 2 photons
  kInvalidClusterCategory
};  // enum EClusterCategory

/**
 A cluster of FMS towers
 */
class HitCluster : public TObject {
 public:
  HitCluster();
  ~HitCluster();
  void CalClusterMoment(Float_t Ecoff);
  void Clear();
  void FindClusterAxis(Float_t Ecoff);
  static const int mMaxPhotonsPerCluster = 2;
  PhotonHitFPD photon[mMaxPhotonsPerCluster];  ///< Photon-Hits in the cluster
  Bool_t IsEUpdated;
  Int_t catag;  ///< catagory of cluster (see EClusterCategory)
  Int_t numbTower;  ///< number of non_zero towers in the cluster
  Int_t nPhoton;  ///< number of photons contained in this cluster
  Int_t index;  ///< cluster number in an event, counts from 0
  Float_t energy;  ///< total energy  contained in this cluster (0th moment)
  Float_t x0;  ///< mean x in local grid coordinate (1st moment)
  Float_t y0;  ///< mean y in local grid coordinate (1st moment)
  Float_t sigmaX;  ///< 2nd moment in x
  Float_t sigmaY;  ///< 2nd moment in y
  Float_t sigmaXY;  ///< 2nd moment in x-y
  Float_t thetaAxis;  ///< theta angle in x-y plane that define the direction
                      ///< of least-2nd-sigma axis
  Float_t sigmaMin;  ///< 2nd sigma w.r.t to the least-2nd-sigma axis
  Float_t sigmaMax;  ///< 2nd sigma w.r.t to the axis orthogonal to the
                     ///< least-2nd-sigma axis
  Float_t chiSquare;  ///< Chi-square of the fitting
  TObjArray* tow;  //!<  TowerFPD objects that make the cluster

 private:
  DISALLOW_COPY_AND_ASSIGN(HitCluster);
  void FindClusterAxis();
  Double_t GetSigma(Double_t theta);
  Float_t Ecutoff;
  Bool_t EnergyUpdated;
  ClassDef(HitCluster, 7)
};  // class HitCluster

inline void HitCluster::FindClusterAxis(Float_t Ecoff) {
  Ecutoff = Ecoff;
  FindClusterAxis();
}
}  // namespace PSUGlobals
#endif

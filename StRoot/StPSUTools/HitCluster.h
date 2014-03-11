#ifndef HITCLUSTER_H
#define HITCLUSTER_H

#include <TObject.h>

#include "StEvent/StFmsCluster.h"
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
 
 This is an elaborated version of the simple StFmsCluster class, storing extra
 information needed during the clustering process.
 */
class HitCluster : public StFmsCluster {
 public:
  HitCluster();
  ~HitCluster();
  void CalClusterMoment(Float_t Ecoff);
  void Clear(const char* optionNotUsed = "");
  void FindClusterAxis(Float_t Ecoff);
  /** Return the index of this cluster in the event */
  Int_t index() const;
  /** Sets the index of this cluster in the event */
  void setIndex(Int_t index);
  /** Return the &chi;<sup>2</sup> of the photon fit for this cluster */
  Float_t chiSquare() const;
  /** Set the &chi;<sup>2</sup> of the photon fit for this cluster */
  void setChiSquare(Float_t chi2);
  Float_t thetaAxis() const;
  /** Return the list of towers in this cluster */
  TObjArray* towers();
  /** \overload */
  const TObjArray* towers() const;
  /** Return the array of photons creating this cluster */
  PhotonHitFPD* photons();
  /** \overload */
  const PhotonHitFPD* photons() const;
  /** Update an StFmsCluster with values from this cluster */
  void copyTo(StFmsCluster*) const;

 protected:
  Int_t mIndex;  ///< cluster number in an event, counts from 0
  Float_t mSigmaX;  ///< 2nd moment in x
  Float_t mSigmaY;  ///< 2nd moment in y
  Float_t mSigmaXY;  ///< 2nd moment in x-y
  Float_t mThetaAxis;  ///< theta angle in x-y plane that define the direction
                      ///< of least-2nd-sigma axis
  Float_t mChiSquare;  ///< Chi-square of the fitting
  TObjArray* mTowers;  //!<  TowerFPD objects that make the cluster
  static const int mMaxPhotonsPerCluster = 2;
  PhotonHitFPD mPhotons[mMaxPhotonsPerCluster];  ///< Photon-Hits in the cluster

 private:
  DISALLOW_COPY_AND_ASSIGN(HitCluster);
  void FindClusterAxis();
  Double_t GetSigma(Double_t theta);
  Float_t Ecutoff;
  ClassDef(HitCluster, 7)
};  // class HitCluster

inline Int_t HitCluster::index() const { return mIndex; }

inline void HitCluster::setIndex(Int_t index) { mIndex = index; }

inline Float_t HitCluster::chiSquare() const { return mChiSquare; }

inline void HitCluster::setChiSquare(Float_t chi2) {
  mChiSquare = chi2;
}

inline Float_t HitCluster::thetaAxis() const { return mThetaAxis; }

inline TObjArray* HitCluster::towers() { return mTowers; }

inline const TObjArray* HitCluster::towers() const { return mTowers; }

inline PhotonHitFPD* HitCluster::photons() { return mPhotons; }

inline const PhotonHitFPD* HitCluster::photons() const { return mPhotons; }

inline void HitCluster::FindClusterAxis(Float_t Ecoff) {
  Ecutoff = Ecoff;
  FindClusterAxis();
}
}  // namespace PSUGlobals
#endif

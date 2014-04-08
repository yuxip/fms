#ifndef STFMSTOWERCLUSTER_H
#define STFMSTOWERCLUSTER_H

#include <TObject.h>

#include "StPSUTools/StFmsFittedPhoton.h"

class TObjArray;
class StFmsCluster;

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
/**
 A cluster of FMS towers
 
 This is an elaborated version of the simple StFmsCluster class, storing extra
 information needed during the clustering process.
 */
class StFmsTowerCluster {
 public:
  explicit StFmsTowerCluster(StFmsCluster* cluster);
  ~StFmsTowerCluster();
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
  StFmsFittedPhoton* photons();
  /** \overload */
  const StFmsFittedPhoton* photons() const;
  /** Return the StEvent cluster structure */
  StFmsCluster* cluster();
  /** \overload */
  const StFmsCluster* cluster() const;

 protected:
  Int_t mIndex;  ///< cluster number in an event, counts from 0
  Float_t mSigmaX;  ///< 2nd moment in x
  Float_t mSigmaY;  ///< 2nd moment in y
  Float_t mSigmaXY;  ///< 2nd moment in x-y
  Float_t mThetaAxis;  ///< theta angle in x-y plane that define the direction
                      ///< of least-2nd-sigma axis
  Float_t mChiSquare;  ///< Chi-square of the fitting
  TObjArray* mTowers;  //!<  StFmsTower objects that make the cluster
  StFmsCluster* mCluster;  //!< Pointer to StEvent cluster structure
  static const int mMaxPhotonsPerCluster = 2;
  StFmsFittedPhoton mPhotons[mMaxPhotonsPerCluster];  ///< Photons in cluster

 private:
  DISALLOW_COPY_AND_ASSIGN(StFmsTowerCluster);
  void FindClusterAxis();
  Double_t GetSigma(Double_t theta);
  Float_t Ecutoff;
  ClassDef(StFmsTowerCluster, 7)
};  // class StFmsTowerCluster

inline Int_t StFmsTowerCluster::index() const { return mIndex; }

inline void StFmsTowerCluster::setIndex(Int_t index) { mIndex = index; }

inline Float_t StFmsTowerCluster::chiSquare() const { return mChiSquare; }

inline void StFmsTowerCluster::setChiSquare(Float_t chi2) {
  mChiSquare = chi2;
}

inline Float_t StFmsTowerCluster::thetaAxis() const { return mThetaAxis; }

inline TObjArray* StFmsTowerCluster::towers() { return mTowers; }

inline const TObjArray* StFmsTowerCluster::towers() const { return mTowers; }

inline StFmsCluster* StFmsTowerCluster::cluster() { return mCluster; }

inline const StFmsCluster* StFmsTowerCluster::cluster() const {
  return mCluster;
}

inline StFmsFittedPhoton* StFmsTowerCluster::photons() { return mPhotons; }

inline const StFmsFittedPhoton* StFmsTowerCluster::photons() const {
  return mPhotons;
}

inline void StFmsTowerCluster::FindClusterAxis(Float_t Ecoff) {
  Ecutoff = Ecoff;
  FindClusterAxis();
}
}  // namespace PSUGlobals
#endif

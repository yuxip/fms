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

namespace FMSCluster {  // $NMSPC
/**
 A cluster of FMS towers
 
 This is an elaborated version of the simple StFmsCluster class, storing extra
 information needed during the clustering process.
 */
class StFmsTowerCluster {
 public:
  /**
   Constructor
   
   Initialise with an StFmsCluster, which the StFmsTowerCluster continues to
   reference. However the StFmsTowerCluster does not own the StFmsCluster,
   so it must be deleted elsewhere if dynamically allocated.
   */
  explicit StFmsTowerCluster(StFmsCluster* cluster);
  /** Destructor */
  ~StFmsTowerCluster();
  /**
   Calculate cluster moments (mean and sigma of tower positions)
   
   Ignore towers below the energy cutoff
   */
  void calculateClusterMoments(Float_t energyCutoff);
  /* Clear photon and tower lists and reset other values to defaults */
  void Clear(const char* optionNotUsed = "");
  /* Determine cluster axis. Also sets energy cutoff for cluster moments */
  void findClusterAxis(Float_t Ecoff) {
    mEnergyCutoff = Ecoff;
    findClusterAxis();
  }
  /** Return the index of this cluster in the event */
  Int_t index() const { return mIndex; }
  /** Sets the index of this cluster in the event */
  void setIndex(Int_t index) { mIndex = index; }
  /** 2nd moment in x */
  float sigmaX() const { return mSigmaX; }
  /** 2nd moment in y */
  float sigmaY() const { return mSigmaY; }
  /** 2nd moment in x-y */
  float sigmaXY() const { return mSigmaXY; }
  /** angle in x-y plane that define the direction of least-2nd-sigma axis */
  Float_t thetaAxis() const { return mThetaAxis; }
  /** Return the &chi;<sup>2</sup> of the photon fit for this cluster */
  Float_t chiSquare() const { return mChiSquare; }
  /** Set the &chi;<sup>2</sup> of the photon fit for this cluster */
  void setChiSquare(Float_t chi2) { mChiSquare = chi2; }
  /** Cutoff on towers to use in moment calculations */
  float energyCutoff() const { return mEnergyCutoff; }
  /** Return the list of towers in this cluster */
  TObjArray* towers() { return mTowers; }
  /** \overload */
  const TObjArray* towers() const { return mTowers; }
  /** Return the array of photons creating this cluster */
  StFmsFittedPhoton* photons() { return mPhotons; }
  /** \overload */
  const StFmsFittedPhoton* photons() const { return mPhotons; }
  /** Return the StEvent cluster structure */
  StFmsCluster* cluster() { return mCluster; }
  /** \overload */
  const StFmsCluster* cluster() const{ return mCluster; }

 protected:
  static const int kMaxPhotonsPerCluster = 2;
  void findClusterAxis();
  Double_t getSigma(Double_t theta) const;
  Int_t mIndex;  ///< cluster number in an event, counts from 0
  Float_t mSigmaX;  ///< 2nd moment in x
  Float_t mSigmaY;  ///< 2nd moment in y
  Float_t mSigmaXY;  ///< 2nd moment in x-y
  Float_t mThetaAxis;  ///< theta angle in x-y plane that define the direction
                       ///< of least-2nd-sigma axis
  Float_t mChiSquare;  ///< Chi-square of the fitting
  Float_t mEnergyCutoff;  //!< Cutoff on towers to use in moment calculations
  TObjArray* mTowers;  //!<  StFmsTower objects that make the cluster
  StFmsCluster* mCluster;  //!< Pointer to StEvent cluster structure
  StFmsFittedPhoton mPhotons[kMaxPhotonsPerCluster];  ///< Photons in cluster

 private:
  DISALLOW_COPY_AND_ASSIGN(StFmsTowerCluster);
  ClassDef(StFmsTowerCluster, 7)
};  // class StFmsTowerCluster
}  // namespace FMSCluster
#endif

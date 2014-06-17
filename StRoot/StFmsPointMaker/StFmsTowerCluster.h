// $Id$
//
// $Log$
/**
 \file      StFmsTowerCluster.h
 \brief     Declaration of StFmsTowerCluster, a cluster of FMS towers
 \author    Steven Heppelmann <steveheppelmann@gmail.com>
 \author    Yuxi Pan <yuxipan@physics.ucla.edu>
 \author    Thomas Burton <tpb@bnl.gov>
 \date      2014
 \copyright Brookhaven National Lab
 */
#ifndef STROOT_STFMSPOINTMAKER_STFMSTOWERCLUSTER_H_
#define STROOT_STFMSPOINTMAKER_STFMSTOWERCLUSTER_H_

#include <TObject.h>

#include <list>
#include <memory>

#include "StFmsPointMaker/StFmsFittedPhoton.h"

class TObjArray;
class StFmsCluster;

namespace FMSCluster {  // $NMSPC
class StFmsTower;
/**
 A cluster of FMS towers.

 This is an elaborated version of the simple StFmsCluster class, storing extra
 information needed during the clustering process.

 The StFmsTower objects in its towers() list are not owned by the
 StFmsTowerCluster, so those towers must have a longer lifetime than it.
 The user must also take responsibility for deleting them if they are
 dynamically allocated.
 */
class StFmsTowerCluster {
 public:
  typedef std::list<StFmsTower*> Towers;  ///< Shorthand for tower collection
  /**
   Constructor.
   
   Initialise with a dynamically allocated StFmsCluster. The StFmsTowerCluster
   owns the StFmsCluster until the release() method is called, after which
   the StFmsTowerCluster no longer references a cluster and should not be used
   any longer.
   */
  explicit StFmsTowerCluster(StFmsCluster* cluster);
  // Use default copy constructor and assignment operator
  /** Destructor */
  ~StFmsTowerCluster();
  /* Clear photon and tower lists and reset other values to defaults */
  void Clear(const char* optionNotUsed = "");
  /**
   Calculate cluster moments (mean and sigma of tower (x, y) position).
   
   Ignore towers below the energy cutoff.
   */
  void calculateClusterMoments(Float_t energyCutoff);
  /**
   Determine cluster axis.

   Also sets energy cutoff for cluster moments.
   */
  void findClusterAxis(Float_t Ecoff) {
    mEnergyCutoff = Ecoff;
    findClusterAxis();
  }
  /** Return the index of this cluster in the event. */
  Int_t index() const { return mIndex; }
  /** Sets the index of this cluster in the event. */
  void setIndex(Int_t index) { mIndex = index; }
  /** 2nd moment in x. */
  float sigmaX() const { return mSigmaX; }
  /** 2nd moment in y. */
  float sigmaY() const { return mSigmaY; }
  /** 2nd moment in x-y. */
  float sigmaXY() const { return mSigmaXY; }
  /** Angle in x-y plane that defines the direction of least-2nd-sigma axis. */
  Float_t thetaAxis() const { return mThetaAxis; }
  /** Return the &chi;<sup>2</sup> of the photon fit for this cluster. */
  Float_t chiSquare() const { return mChiSquare; }
  /** Set the &chi;<sup>2</sup> of the photon fit for this cluster. */
  void setChiSquare(Float_t chi2) { mChiSquare = chi2; }
  /** Cutoff on towers to use in moment calculations. */
  float energyCutoff() const { return mEnergyCutoff; }
  /** Return the list of towers in this cluster. */
  Towers& towers() { return mTowers; }
  /** \overload */
  const Towers& towers() const { return mTowers; }
  /** Return the array of photons creating this cluster. */
  StFmsFittedPhoton* photons() { return mPhotons; }
  /** \overload */
  const StFmsFittedPhoton* photons() const { return mPhotons; }
  /** Return the StEvent cluster structure. */
  StFmsCluster* cluster() { return mCluster.get(); }
  /** \overload */
  const StFmsCluster* cluster() const { return mCluster.get(); }
  /** Return and give up ownership of the StEvent cluster structure. */
  StFmsCluster* release() { return mCluster.release(); }

 protected:
  static const int kMaxPhotonsPerCluster = 2;  ///< Support 2-photon clusters
  /** Determine cluster axis. */
  void findClusterAxis();
  /** Calculate sigma w.r.t the axis going through the "center" and of an angle
      "theta" in x-y plane. */
  Double_t getSigma(Double_t theta) const;
  Int_t mIndex;  ///< cluster number in an event, counts from 0
  Float_t mSigmaX;  ///< 2nd moment in x
  Float_t mSigmaY;  ///< 2nd moment in y
  Float_t mSigmaXY;  ///< 2nd moment in x-y
  Float_t mThetaAxis;  ///< theta angle in x-y plane that define the direction
                       ///< of least-2nd-sigma axis
  Float_t mChiSquare;  ///< Chi-square of the fitting
  Float_t mEnergyCutoff;  //!< Cutoff on towers to use in moment calculations
  Towers mTowers;  //!< Towers that make the cluster
  std::auto_ptr<StFmsCluster> mCluster;  //!< Pointer to StEvent cluster
  StFmsFittedPhoton mPhotons[kMaxPhotonsPerCluster];  ///< Photons in cluster

 private:
  ClassDef(StFmsTowerCluster, 0)
};  // class StFmsTowerCluster
}  // namespace FMSCluster
#endif  // STROOT_STFMSPOINTMAKER_STFMSTOWERCLUSTER_H_

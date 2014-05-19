//Yuxi Pan 03/31/2013
#ifndef STROOT_STEVENT_STFMSCLUSTER_H_
#define STROOT_STEVENT_STFMSCLUSTER_H_

#include <TLorentzVector.h>

#include "StObject.h"

#include "StEvent/StContainers.h"  // StPtrVecFmsHit, StPtrVecFmsPoint

enum EFmsClusterCategory {
  kAmbiguousCluster = 0,  // Could be 1- or 2-photon, needs to be fitted
  k1PhotonCluster = 1,  // A cluster created by 1 photon
  k2PhotonCluster = 2,  // A cluster created by 2 photons
  kInvalidClusterCategory
};  // enum EFmsClusterCategory

/** Represents a cluster of a group of towers */
class StFmsCluster : public StObject {
 public:
  /** Constructor */
  StFmsCluster();
  // Use default copy constructor and assignment operator
  /** Destructor */
  ~StFmsCluster();
  /** Print cluster information to LOG_INFO */
  void Print(Option_t* optionNotUsed = "") const;
  /** Sub-detector */
  UShort_t detectorId() const { return mDetectorId; }
  /** Cluster category (see EFmsClusterCategory) */
  Int_t category() const { return mCategory; }
  /** Number of towers (hits) in this cluster */
  Int_t nTowers() const { return mNTowers; }
  /** Number of points/photons forming this cluster */
  Int_t nPhotons() const { return mNPhotons; }
  /** Total cluster energy */
  Float_t energy() const { return mEnergy; }
  /** Mean x ("center of gravity") in local grid coordinate (1st moment) */
  Float_t x() const { return mX; }
  /** Mean y ("center of gravity") in local grid coordinate (1st moment) */
  Float_t y() const { return mY; }
  /** Maximum 2nd moment (along major axis) */
  Float_t sigmaMax() const { return mSigmaMax; }
  /** Minimum 2nd moment */
  Float_t sigmaMin() const { return mSigmaMin; }
  /** chi2 / ndf for 1-photon fit */
  Float_t chi2Ndf1Photon() const { return mChi2Ndf1Photon; }
  /** chi2 / ndf for 2-photon fit */
  Float_t chi2Ndf2Photon() const { return mChi2Ndf2Photon; }
  /** Eventwise cluster ID, also include nstb info. */
  Int_t id() const { return mId; }
  /** Cluster four-momentum (px, py, pz E) */
  TLorentzVector fourMomentum() const { return mFourMomentum; }
  /** Set sub-detector */
  void setDetectorId(UShort_t detector) { mDetectorId = detector; }
  /** Set cluster category (see EFmsClusterCategory) */
  void setCategory(Int_t catag) { mCategory = catag; }
  /** Set number of towers (hits) in this cluster */
  void setNTowers(Int_t numbTower) { mNTowers = numbTower; }
  /** Set number of points/photons forming this cluster */
  Bool_t setNPhotons(Int_t nPhoton);
  /** Set total cluster energy */
  void setEnergy(Float_t energy) { mEnergy = energy; }
  /** Set cluster mean x in local grid coordinates */
  void setX(Float_t x0) { mX = x0; }
  /** Set cluster mean y in local grid coordinates */
  void setY(Float_t y0) { mY = y0; }
  /** Set minimum 2nd moment */
  void setSigmaMin(Float_t sigmaMax) { mSigmaMin = sigmaMax; }
  /** Set maximum 2nd moment */
  void setSigmaMax(Float_t sigmaMax) { mSigmaMax = sigmaMax; }
  /** Set chi2 / ndf for 1-photon fit */
  void setChi2Ndf1Photon(Float_t chi2ndfph1) { mChi2Ndf1Photon = chi2ndfph1; }
  /** Set chi2 / ndf for 2-photon fit */
  void setChi2Ndf2Photon(Float_t chi2ndfph2) { mChi2Ndf2Photon = chi2ndfph2; }
  /** Set cluster ID */
  void setId(Float_t cluid) { mId = cluid; }
  /** Set cluster four-momentum (px, py, pz, E) */
  void setFourMomentum(TLorentzVector p4) { mFourMomentum = p4; }
  /** Towers/hits in this cluster */
  StPtrVecFmsHit& hits() { return mHits; }
  /** \overload */
  const StPtrVecFmsHit& hits() const { return mHits; }
  /** Points/photons forming this cluster */
  StPtrVecFmsPoint& points() { return mPhotons; }
  /** \overload */
  const StPtrVecFmsPoint& points() const { return mPhotons; }

 protected:
  UShort_t mDetectorId;  ///< Detector starts from 1
  Int_t mCategory;  ///< catagory of cluster (see EFmsClusterCategory)
  Int_t mNTowers;  ///< Number of non_zero towers in the cluster
  Int_t mNPhotons;  ///< Number of photons contained in this cluster
  Float_t mEnergy;  ///< Total energy  contained in this cluster (0th moment)
  Float_t mX;  ///< Mean x ("center of gravity") in local grid coordinate
                ///< (1st moment)
  Float_t mY;  ///< Mean y ("center of gravity") in local grid coordinate
                ///< (1st moment)
  Float_t mSigmaMin;  ///< Minimum 2nd moment
  Float_t mSigmaMax;  ///< Maximum 2nd moment (along major axis)
  Float_t mChi2Ndf1Photon;  ///< chi2 / ndf for 1-photon fit
  Float_t mChi2Ndf2Photon;  ///< chi2 / ndf for 2-photon fit
  Int_t mId;  ///< Eventwise cluster id, also include nstb info.
  TLorentzVector mFourMomentum;  ///< Cluster four momentum;
  StPtrVecFmsPoint mPhotons;  ///< Tower hits of the current cluster
  StPtrVecFmsHit mHits;  ///< Fitted points (photons) in the cluster
  ClassDef(StFmsCluster, 1)
};

#endif  // STROOT_STEVENT_STFMSCLUSTER_H_

//Yuxi Pan 03/31/2013
#ifndef StStFmsCluster_HH
#define StStFmsCluster_HH

#include <TLorentzVector.h>

#include "StObject.h"

#include "StEvent/StContainers.h"  // StPtrVecFmsHit, StPtrVecFmsPoint

/** Represents a cluster of a group of towers */
class StFmsCluster : public StObject {
 public:
  /** Constructor */
  StFmsCluster();
  /** Destructor */
  ~StFmsCluster();
  /** Print cluster information to LOG_INFO */
  void Print(Option_t* optionNotUsed = "") const;
  /** Sub-detector */
  Int_t GetNstb() const { return mNstb; }
  /** Cluster category (see EClusterCategory) */
  Int_t GetCatag() const { return mCatag; }
  /** Number of towers (hits) in this cluster */
  Int_t GetNTower() const { return mNumbTower; }
  /** Number of points/photons forming this cluster */
  Int_t GetNphoton() const { return mNphoton; }
  /** Total cluster energy */
  Float_t GetEnergy() const { return mEnergy; }
  /** Mean x ("center of gravity") in local grid coordinate (1st moment) */
  Float_t GetX0() const { return mX0; }
  /** Mean y ("center of gravity") in local grid coordinate (1st moment) */
  Float_t GetY0() const { return mY0; }
  /** Maximum 2nd moment (along major axis) */
  Float_t GetSigmaMax() const { return mSigmaMax; }
  /** Minimum 2nd moment */
  Float_t GetSigmaMin() const { return mSigmaMin; }
  /** chi2 / ndf for 1-photon fit */
  Float_t GetChi2NdfPh1() const { return mChi2NdfPh1; }
  /** chi2 / ndf for 2-photon fit */
  Float_t GetChi2NdfPh2() const { return mChi2NdfPh2; }
  /** Eventwise cluster ID, also include nstb info. */
  Int_t GetClusterId() const { return mCluId; }
  /** Cluster four-momentum (px, py, pz E) */
  TLorentzVector GetFourMomentum() const { return mFourMomentum; }
  /** Set sub-detector */
  void SetNstb(Int_t nstb) { mNstb = nstb; }
  /** Set cluster category (see EClusterCategory) */
  void SetCatag(Int_t catag) { mCatag = catag; }
  /** Set number of towers (hits) in this cluster */
  void SetNumbTower(Int_t numbTower) { mNumbTower = numbTower; }
  /** Set number of points/photons forming this cluster */
  Bool_t SetNphoton(Int_t nPhoton);
  /** Set total cluster energy */
  void SetClusterEnergy(Float_t energy) { mEnergy = energy; }
  /** Set cluster mean x in local grid coordinates */
  void SetX0(Float_t x0) { mX0 = x0; }
  /** Set cluster mean y in local grid coordinates */
  void SetY0(Float_t y0) { mY0 = y0; }
  /** Set maximum 2nd moment */
  void SetSigmaMax(Float_t sigmaMax) { mSigmaMax = sigmaMax; }
  /** Set minimum 2nd moment */
  void SetSigmaMin(Float_t sigmaMin) { mSigmaMin = sigmaMin; }
  /** Set chi2 / ndf for 1-photon fit */
  void SetChi2NdfPh1(Float_t chi2ndfph1) { mChi2NdfPh1 = chi2ndfph1; }
  /** Set chi2 / ndf for 2-photon fit */
  void SetChi2NdfPh2(Float_t chi2ndfph2) { mChi2NdfPh2 = chi2ndfph2; }
  /** Set cluster ID */
  void SetClusterId(Float_t cluid) { mCluId = cluid; }
  /** Set cluster four-momentum (px, py, pz, E) */
  void SetFourMomentum(TLorentzVector p4) { mFourMomentum = p4; }
  /** Towers/hits in this cluster */
  StPtrVecFmsHit& hits() { return mHits; }
  /** Points/photons forming this cluster */
  StPtrVecFmsPoint& points() { return mPhotons; }

 protected:
  Int_t mNstb;  ///< Nstb starts from 1
  Int_t mCatag;  ///< catagory of cluster (see EClusterCategory)
  Int_t mNumbTower;  ///< Number of non_zero towers in the cluster
  Int_t mNphoton;  ///< Number of photons contained in this cluster
  Float_t mEnergy;  ///< Total energy  contained in this cluster (0th moment)
  Float_t mX0;  ///< Mean x ("center of gravity") in local grid coordinate
                ///< (1st moment)
  Float_t mY0;  ///< Mean y ("center of gravity") in local grid coordinate
                ///< (1st moment)
  Float_t mSigmaMax;  ///< Maximum 2nd moment (along major axis)
  Float_t mSigmaMin;  ///< Minimum 2nd moment
  Float_t mChi2NdfPh1;  ///< chi2 / ndf for 1-photon fit
  Float_t mChi2NdfPh2;  ///< chi2 / ndf for 2-photon fit
  Int_t mCluId;  ///< Eventwise cluster id, also include nstb info.
  TLorentzVector mFourMomentum;  ///< Cluster four momentum;
  StPtrVecFmsPoint mPhotons;  ///< Tower hits of the current cluster
  StPtrVecFmsHit mHits;  ///< Fitted points (photons) in the cluster
  ClassDef(StFmsCluster, 1)
};

#endif

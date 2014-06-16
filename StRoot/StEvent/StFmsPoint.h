// $Id$
//
// $Log$
/**
 \file      StFmsPoint.h
 \brief     Declaration of StFmsPoint, the StEvent FMS photon structure
 \author    Yuxi Pan <yuxipan@physics.ucla.edu>
 \author    Thomas Burton <tpb@bnl.gov>
 \date      2014
 \copyright Brookhaven National Lab
 */
#ifndef STROOT_STEVENT_STFMSPOINT_H_
#define STROOT_STEVENT_STFMSPOINT_H_

#include <TLorentzVector.h>

#include "StRoot/St_base/StObject.h"

class StFmsCluster;

/** Represents a "point" (photon etc) fitted to a cluster of FMS towers. */
class StFmsPoint : public StObject {
 public:
  /** Constructor. */
  StFmsPoint();
  // Use default copy constructor and assignment operator
  /** Destructor. */
  ~StFmsPoint();
  /**
   Sub-detector in which the point was reconstructed.

   ID numbers are defined the STAR database.
   */
  UShort_t detectorId() const { return mDetectorId; }
  /** Energy of the point in GeV. */
  Float_t energy() const { return mEnergy; }
  /** x position in cm at which point intersects the sub-detector. */
  Float_t x() const { return mX; }
  /** y position in cm at which point intersects the sub-detector. */
  Float_t y() const { return mY; }
  /** ID of the point in the current event. */
  Int_t id() const { return mId; }
  /** Parent cluster of the photon. */
  StFmsCluster* cluster() { return mCluster; }
  /** Parent cluster of the photon. */
  const StFmsCluster* cluster() const { return mCluster; }
  /** ID of the parent cluster containing this point. */
  Int_t parentClusterId() const { return mParentClusterId; }
  /** Number of points in the parent cluster. */
  Int_t nParentClusterPhotons() const { return mNParentClusterPhotons; }
  /** 4-momentum of the point (px, py, pz, E). */
  TLorentzVector fourMomentum() const { return mFourMomentum; }
  /** Set the sub-detector in which the point was reconstructed. */
  void setDetectorId(UShort_t detector) { mDetectorId = detector; }
  /** Set the energy of the point in GeV. */
  void setEnergy(Float_t energy) { mEnergy = energy; }
  /** Set the x position in cm at which point intersects the sub-detector. */
  void setX(Float_t xpos) { mX = xpos; }
  /** Set the y position in cm at which point intersects the sub-detector. */
  void setY(Float_t ypos) { mY = ypos; }
  /** Set the ID of the point in the current event. */
  void setId(Int_t phid) { mId = phid; }
  /** Set the parent cluster of the photon. */
  void setCluster(StFmsCluster* cluster) { mCluster = cluster; }
  /** Set the ID of the parent cluster containing this point. */
  void setParentClusterId(Int_t cluid) { mParentClusterId = cluid; }
  /** Set the number of points in the parent cluster. */
  void setNParentClusterPhotons(Int_t nclph) { mNParentClusterPhotons = nclph; }
  /** Set the 4-momentum of the point (px, py, pz, E). */
  void setFourMomentum(const TLorentzVector& p4) { mFourMomentum = p4; }

 protected:
  UShort_t mDetectorId;  ///< Detector starts from 1
  Float_t mEnergy;  ///< Fitted energy
  Float_t mX;  ///< Fitted x-position
  Float_t mY;  ///< Fitted y-position
  Int_t mId;  ///< Photon ID within event
  Int_t mParentClusterId;  ///< ID of the parent cluster within event
  Int_t mNParentClusterPhotons;  ///< Number of photons in the parent cluster
  StFmsCluster* mCluster;  //!< Parent cluster of this photon
  TLorentzVector mFourMomentum;  ///< Photon 4-momentum
  ClassDef(StFmsPoint, 1)
};

#endif  // STROOT_STEVENT_STFMSPOINT_H_

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

/** Represents a point (photon etc) fitted to a cluster of FMS towers */
class StFmsPoint : public StObject {
 public:
  /** Constructor */
  StFmsPoint();
  // Use default copy constructor and assignment operator
  /** Destructor */
  ~StFmsPoint();
  /** Sub-detector in which the point was reconstructed */
  UShort_t detectorId() const { return mDetectorId; }
  /** Energy of the point in GeV */
  Float_t energy() const { return mEnergy; }
  /** x position in cm at which point intersects the sub-detector */
  Float_t x() const { return mX; }
  /** y position in cm at which point intersects the sub-detector */
  Float_t y() const { return mY; }
  /** ID of the point in the current event */
  Int_t id() const { return mId; }
  /** ID of the parent cluster containing this point */
  Int_t parentClusterId() const { return mParentClusterId; }
  /** Number of points in the parent cluster */
  Int_t nParentClusterPhotons() const { return mNParentClusterPhotons; }
  /** 4-momentum of the point */
  TLorentzVector fourMomentum() const { return mFourMomentum; }
  /** Set the sub-detector in which the point was reconstructed */
  void setDetectorId(UShort_t detector) { mDetectorId = detector; }
  /** Set the energy of the point in GeV */
  void setEnergy(Float_t energy) { mEnergy = energy; }
  /** Set the x position in cm at which point intersects the sub-detector */
  void setX(Float_t xpos) { mX = xpos; }
  /** Set the y position in cm at which point intersects the sub-detector */
  void setY(Float_t ypos) { mY = ypos; }
  /** Set the ID of the point in the current event */
  void setId(Int_t phid) { mId = phid; }
  /** Set the ID of the parent cluster containing this point */
  void setParentClusterId(Int_t cluid) { mParentClusterId = cluid; }
  /** Set the number of points in the parent cluster */
  void setNParentClusterPhotons(Int_t nclph) { mNParentClusterPhotons = nclph; }
  /** Set the 4-momentum of the point */
  void setFourMomentum(const TLorentzVector& p4) { mFourMomentum = p4; }

 protected:
  UShort_t mDetectorId;  ///< Detector starts from 1
  Float_t mEnergy;  ///< fitted energy
  Float_t mX;  ///< fitted(relative) x-position
  Float_t mY;  ///< fitted(relative) y-position
  Int_t mId;  ///< photon id within event, also include det12 info
  Int_t mParentClusterId;  ///< id of the parent cluster
  Int_t mNParentClusterPhotons;  ///< # of photons in by the parent cluster
  TLorentzVector mFourMomentum;  ///< Photon 4-momentum
  ClassDef(StFmsPoint, 1)
};

#endif  // STROOT_STEVENT_STFMSPOINT_H_

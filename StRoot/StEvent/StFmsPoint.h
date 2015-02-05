/*************************************************************************
 *
 * $Id$
 *
 * Author: Thomas Burton, Yuxi Pan, 2014
 *************************************************************************
 *
 * Description: Declaration of StFmsPoint, the StEvent FMS photon structure
 *
 **************************************************************************
 *
 * $Log$
 *
 **************************************************************************/
#ifndef STROOT_STEVENT_STFMSPOINT_H_
#define STROOT_STEVENT_STFMSPOINT_H_

#include "StLorentzVectorF.hh"

#include "StObject.h"

class StFmsCluster;

/** Represents a "point" (photon etc) fitted to a cluster of FMS towers. */
class StFmsPoint : public StObject {
 public:
  StFmsPoint();
  ~StFmsPoint();

  unsigned short detectorId() const;
  float energy() const;
  float x() const; // x position in cm at which point intersects the sub-detector.
  float y() const; // y position in cm at which point intersects the sub-detector.
  int id() const; // ID of the point in the current event.
  StFmsCluster* cluster(); //  Parent cluster of the photon.
  const StFmsCluster* cluster() const;
  int parentClusterId() const; // ID of the parent cluster containing this point.
  int nParentClusterPhotons() const; // Number of points in the parent cluster.
  StLorentzVectorF fourMomentum() const;
  void setDetectorId(unsigned short detector);
  void setEnergy(float energy);
  void setX(float xpos);
  void setY(float ypos);
  void setId(int phid);
  void setCluster(StFmsCluster* cluster);
  void setParentClusterId(int cluid);
  void setNParentClusterPhotons(int nclph);
  void setFourMomentum(const StLorentzVectorF& p4);

 private:
  UShort_t mDetectorId;  ///< Detector starts from 1
  Float_t mEnergy;  ///< Fitted energy
  Float_t mX;  ///< Fitted x-position
  Float_t mY;  ///< Fitted y-position
  Int_t mId;  ///< Photon ID within event
  Int_t mParentClusterId;  ///< ID of the parent cluster within event
  Int_t mNParentClusterPhotons;  ///< Number of photons in the parent cluster
  StFmsCluster* mCluster;  //!< Parent cluster of this photon
  StLorentzVectorF mFourMomentum;  ///< Photon 4-momentum
  ClassDef(StFmsPoint, 1)
};

  inline unsigned short StFmsPoint::detectorId() const { return mDetectorId; }
  inline float StFmsPoint::energy() const { return mEnergy; }
  inline float StFmsPoint::x() const { return mX; } // x position in cm at which point intersects the sub-detector.
  inline float StFmsPoint::y() const { return mY; } // y position in cm at which point intersects the sub-detector.
  inline int StFmsPoint::id() const { return mId; } // ID of the point in the current event.
  inline StFmsCluster* StFmsPoint::cluster() { return mCluster; } //  Parent cluster of the photon.
  inline const StFmsCluster* StFmsPoint::cluster() const { return mCluster; }
  inline int StFmsPoint::parentClusterId() const { return mParentClusterId; } // ID of the parent cluster containing this point.
  inline int StFmsPoint::nParentClusterPhotons() const { return mNParentClusterPhotons; } // Number of points in the parent cluster.
  inline StLorentzVectorF StFmsPoint::fourMomentum() const { return mFourMomentum; }
  inline void StFmsPoint::setDetectorId(unsigned short detector) { mDetectorId = detector; }
  inline void StFmsPoint::setEnergy(float energy) { mEnergy = energy; }
  inline void StFmsPoint::setX(float xpos) { mX = xpos; }
  inline void StFmsPoint::setY(float ypos) { mY = ypos; }
  inline void StFmsPoint::setId(int phid) { mId = phid; }
  inline void StFmsPoint::setCluster(StFmsCluster* cluster) { mCluster = cluster; }
  inline void StFmsPoint::setParentClusterId(int cluid) { mParentClusterId = cluid; }
  inline void StFmsPoint::setNParentClusterPhotons(int nclph) { mNParentClusterPhotons = nclph; }
  inline void StFmsPoint::setFourMomentum(const StLorentzVectorF& p4) { mFourMomentum = p4; }

#endif  // STROOT_STEVENT_STFMSPOINT_H_

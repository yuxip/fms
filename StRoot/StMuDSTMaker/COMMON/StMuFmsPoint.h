// $Id$
//
// $Log$
/**
 \file      StMuFmsPoint.h
 \brief     Declaration of StMuFmsPoint, the MuDST FMS "point" class
 \author    Thomas Burton <tpb@bnl.gov>
 \date      2014
 \copyright Brookhaven National Lab
 */
#ifndef STROOT_STMUDSTMAKER_COMMON_STMUFMSPOINT_H_
#define STROOT_STMUDSTMAKER_COMMON_STMUFMSPOINT_H_

#include <TLorentzVector.h>
#include <TObject.h>
#include <TRef.h>
#include <TVector3.h>

class StFmsPoint;  // The equivalent point structure in StEvent
class StMuFmsCluster;

/**
 Micro-DST FMS "point" class.

 Describes a "point" - the energy deposited by a single particle in a cluster.
 One or more points may be form a cluster of adjacent towers in the FMS.

 Maintains a persistent reference to the cluster formed by the point.
 The cluster is owned by the relevant TClonesArray in the micro-DST, not
 StMuFmsPoint, and should not be deleted.
 */
class StMuFmsPoint : public TObject {
 public:
  /** Constructor. */
  StMuFmsPoint(int detectorId = 0, float energy = 0.f,
               float x = 0.f, float y = 0.f, float z = 0.f);
  /** Construct from the equivalent StEvent point structure. */
  explicit StMuFmsPoint(const StFmsPoint&);
  /** Destructor. */
  virtual ~StMuFmsPoint();
  /** ID of the sub-detector with which the point is associated. */
  UShort_t detectorId() const { return mDetectorId; }
  /** Total point energy. */
  float energy() const { return mEnergy; }
  /** x "center of gravity" of the point (cm). */
  float x() const { return mX; }
  /** y "center of gravity" of the point (cm). */
  float y() const { return mY; }
  /** z position of front face of sub-detector (cm). */
  float z() const { return mZ; }
  /** (x, y, z) position of point at sub-detector face. */
  TVector3 xyz() const { return TVector3(mX, mY, mZ); }
  /**
   (px, py, pz) of point.

   Assumes some mass, which must be <= energy.
   */
  TVector3 momentum(float m = 0.f) const;
  /** (px, py, pz, E) of point. See also comments for momentum(). */
  TLorentzVector fourMomentum(float m = 0.f) const;
  /** Parent cluster of this photon (NULL if not known). */
  StMuFmsCluster* cluster();
  /** Parent cluster of this photon (NULL if not known). */
  const StMuFmsCluster* cluster() const;
  /** Set ID of the sub-detector with which the point is associated. */
  void setDetectorId(UShort_t detector) { mDetectorId = detector; }
  /** Set total point energy (sum over towers). */
  void setEnergy(float energy) { mEnergy = energy; }
  /** Set x "center of gravity" of the point. */
  void setX(float x) { mX = x; }
  /** Set y "center of gravity" of the point. */
  void setY(float y) { mY = y; }
  /** Set z position of front face of sub-detector (cm). */
  void setZ(float z) { mZ = z; }
  /** Set properties from an StFmsPoint. */
  void set(const StFmsPoint&);
  /** Set parent cluster of this photon. */
  void setCluster(StMuFmsCluster* cluster);

 protected:
  UShort_t mDetectorId;  ///< Detector ID as defined in database
  Float_t mEnergy;  ///< Total energy contained in the point
  Float_t mX;  ///< Mean x ("center of gravity")
  Float_t mY;  ///< Mean y ("center of gravity")
  Float_t mZ;  ///< z at front face of sub-detector
  TRef mCluster;  ///< Parent cluster of this photon

 private:
  /**
   Disallow copy construction.

   Duplication should only be done via Clone().
   */
  StMuFmsPoint(const StMuFmsPoint&);
  /**
   Disallow assignment.

   Duplication should only be done via Clone().
   */
  StMuFmsPoint& operator=(const StMuFmsPoint&);
  ClassDef(StMuFmsPoint, 1)
};
#endif  // STROOT_STMUDSTMAKER_COMMON_STMUFMSPOINT_H_

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

#include <TObject.h>

class StFmsPoint;  // The equivalent point structure in StEvent

/**
 Micro-DST FMS "point" class

 Describes a "point" - an energy deposition from a single particle.
 One or more points may be form a cluster of adjacent towers in the FMS.
 */
class StMuFmsPoint : public TObject {
 public:
  /** Constructor */
  StMuFmsPoint(int detectorId = 0, float energy = 0.f,
               float x = 0.f, float y = 0.f, float z = 0.f);
  /** Construct from the equivalent StEvent point structure */
  StMuFmsPoint(const StFmsPoint&);
  /** Destructor */
  virtual ~StMuFmsPoint();
  /** ID of the sub-detector with which the point is associated */
  UShort_t detectorId() const { return mDetectorId; }
  /** Total point energy */
  float energy() const { return mEnergy; }
  /** x "center of gravity" of the point (cm) */
  float x() const { return mX; }
  /** y "center of gravity" of the point (cm) */
  float y() const { return mY; }
  /** z position of front face of sub-detector (cm) */
  float z() const { return mZ; }
  /** Set ID of the sub-detector with which the point is associated */
  void setDetectorId(UShort_t detector) { mDetectorId = detector; }
  /** Set total point energy (sum over towers) */
  void setEnergy(float energy) { mEnergy = energy; }
  /** Set x "center of gravity" of the point */
  void setX(float x) { mX = x; }
  /** Set y "center of gravity" of the point */
  void setY(float y) { mY = y; }
  /** Set z position of front face of sub-detector (cm) */
  void setZ(float z) { mZ = z; }

 protected:
  UShort_t mDetectorId;  ///< Detector ID as defined in database
  Float_t mEnergy;  ///< Total energy contained in the point
  Float_t mX;  ///< Mean x ("center of gravity")
  Float_t mY;  ///< Mean y ("center of gravity")
  Float_t mZ;  ///< z at front face of sub-detector

 private:
  /**
   Disallow copy construction.

   Duplication should only be done via Clone()
   */
  StMuFmsPoint(const StMuFmsPoint&);
  /**
   Disallow assignment.

   Duplication should only be done via Clone()
   */
  StMuFmsPoint& operator=(const StMuFmsPoint&);
  ClassDef(StMuFmsPoint, 1)
};

#endif  // STROOT_STMUDSTMAKER_COMMON_STMUFMSPOINT_H_

// $Id$
//
// $Log$
/**
 \file      StMuFmsCluster.h
 \brief     Declaration of StMuFmsCluster, the MuDST FMS cluster class
 \author    Thomas Burton <tpb@bnl.gov>
 \date      2014
 \copyright Brookhaven National Lab
 */

#ifndef STROOT_STMUDSTMAKER_COMMON_STMUFMSCLUSTER_H_
#define STROOT_STMUDSTMAKER_COMMON_STMUFMSCLUSTER_H_

#include <TRefArray.h>

#include "StRoot/St_base/StObject.h"

class StFmsCluster;  // Equivalent class in StEvent

/**
 Micro-DST FMS cluster class

 Describes a cluster (collection of adjacent towers) in the FMS.
 May be created by one or more photons.
 */
class StMuFmsCluster : public StObject {
 public:
  /** Constructor */
  StMuFmsCluster(int detectorId = 0, int category = -1, float energy = 0.f,
                 float x = 0.f, float y = 0.f);
  /** Initialise from an equivalent StEvent cluster */
  StMuFmsCluster(const StFmsCluster&);
  /** Destructor */
  virtual ~StMuFmsCluster();
  /** ID of the sub-detector with which the cluster is associated */
  UShort_t detectorId() const { return mDetectorId; }
  /** Category of the cluster (see EFmsClusterCategory) */
  unsigned short category() const { return mCategory; }
  /** Total cluster energy (sum over towers) */
  float energy() const { return mEnergy; }
  /** x "center of gravity" of the cluster */
  float x() const { return mX; }
  /** y "center of gravity" of the cluster */
  float y() const { return mY; }
  /** The collection of hits in the cluster */
  TRefArray* hits() { return &mHits; }
  /** \overload */
  const TRefArray* hits() const { return &mHits; }
  /** Photons in this cluster */
  TRefArray* photons() { return &mPhotons; }
  /** \overload */
  const TRefArray* photons() const { return &mPhotons; }
  /** Set ID of the sub-detector with which the cluster is associated */
  void setDetectorId(UShort_t detector) { mDetectorId = detector; }
  /** Set category of the cluster (see EFmsClusterCategory) */
  void setCategory(unsigned short category) { mCategory = category; }
  /** Set total cluster energy (sum over towers) */
  void setEnergy(float energy) { mEnergy = energy; }
  /** Set x "center of gravity" of the cluster */
  void setX(float x) { mX = x; }
  /** Set y "center of gravity" of the cluster */
  void setY(float y) { mY = y; }

 protected:
  UShort_t mDetectorId;  ///< Detector ID as defined in database
  UShort_t mCategory;  ///< Category of cluster (see EFmsClusterCategory)
  Float_t mEnergy;  ///< Total energy contained in the cluster
  Float_t mX;  ///< Mean x ("center of gravity")
  Float_t mY;  ///< Mean y ("center of gravity")
  TRefArray mHits;  ///< StMuFmsHits in the current cluster
  TRefArray mPhotons;  ///< StMuFmsPoints in the cluster

 private:
  /**
   Disallow copy construction.

   Duplication should only be done via Clone()
   */
  StMuFmsCluster(const StMuFmsCluster&);
  /**
   Disallow assignment.

   Duplication should only be done via Clone()
   */
  StMuFmsCluster& operator=(const StMuFmsCluster&);
  ClassDef(StMuFmsCluster, 1)
};

#endif  // STROOT_STMUDSTMAKER_COMMON_STMUFMSCLUSTER_H_

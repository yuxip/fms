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
  /** Destructor */
  virtual ~StMuFmsCluster();

 protected:
  Int_t mDetector;  ///< Detector ID as defined in database
  Int_t mCategory;  ///< Category of cluster (see EFmsClusterCategory)
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

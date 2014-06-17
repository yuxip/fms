// $Id$
//
// $Log$
/**
 \file      StFmsFittedPhoton.h
 \brief     Declaration of StFmsFittedPhoton, a photon fitted to an FMS cluster
 \author    Steven Heppelmann <steveheppelmann@gmail.com>
 \author    Yuxi Pan <yuxipan@physics.ucla.edu>
 \author    Thomas Burton <tpb@bnl.gov>
 \date      2014
 \copyright Brookhaven National Lab
 */
#ifndef STROOT_STFMSPOINTMAKER_STFMSFITTEDPHOTON_H_
#define STROOT_STFMSPOINTMAKER_STFMSFITTEDPHOTON_H_

#include "Rtypes.h"  // For ClassDef macro

namespace FMSCluster {  // $NMSPC
/**
 Definition of a photon hit (SMD position info and reconstructed energy).

 The photon properties are determined by a shower-shape fit to a cluster of
 FMS towers.
 */
struct StFmsFittedPhoton {
  /** Constructor with optional position and energy. */
  StFmsFittedPhoton(Float_t x = -1., Float_t y = -1., Float_t e = 0.,
                    Float_t xerr = -1., Float_t yerr = -1., Float_t eerr = -1.);
  // Use default copy constructor and assignment operator
  /** Destructor */
  ~StFmsFittedPhoton() { }
  /** Reset all values to defaults. */
  void Clear();
  Float_t energy;  ///< Fitted energy
  Float_t errEne;  ///< Energy fit error
  Float_t xPos;  ///< Fitted (relative) x-position
  Float_t errXPos;  ///< x-position fit error
  Float_t yPos;  ///< Fitted (relative) y-position
  Float_t errYPos;  ///< y-position fit error
  ClassDef(StFmsFittedPhoton, 0)
};  // class StFmsFittedPhoton
}  // namespace FMSCluster

#endif  // STROOT_STFMSPOINTMAKER_STFMSFITTEDPHOTON_H_

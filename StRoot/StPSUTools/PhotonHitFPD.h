#ifndef PHOTONHITFPD_H
#define PHOTONHITFPD_H

namespace PSUGlobals {//$NMSPC
/**
 Definition of a photon hit (SMD position info and reconstructed energy)
 */
struct PhotonHitFPD {
  Float_t energy;  ///< Fitted energy
  Float_t errEne;  ///< Energy fit error
  Float_t xPos;  ///< Fitted (relative) x-position
  Float_t errXPos;  ///< x-position fit error
  Float_t yPos;  ///< Fitted (relative) y-position
  Float_t errYPos;  ///< y-position fit error
  /** Constructor with optional position and energy */
  PhotonHitFPD(Float_t x = -1., Float_t y = -1., Float_t e = 0.);
  /** Destructor */
  ~PhotonHitFPD() { }
  /** Reset all values to defaults */
  void Clear();
  ClassDef(PhotonHitFPD, 4)
};  // class PhotonHitFPD
}  // namespace PSUGlobals

#endif

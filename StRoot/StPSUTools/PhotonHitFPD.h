#ifndef PHOTONHITFPD_H
#define PHOTONHITFPD_H

#include <TObject.h>

namespace PSUGlobals {//$NMSPC
/**
 Definition of a photon hit (SMD position info and reconstructed energy)
 */
class PhotonHitFPD : public TObject {
 public:
  Float_t energy;  ///< Fitted energy
  Float_t errEne;  ///< Energy fit error
  Float_t xPos;  ///< Fitted (relative) x-position
  Float_t errXPos;  ///< x-position fit error
  Float_t yPos;  ///< Fitted (relative) y-position
  Float_t errYPos;  ///< y-position fit error
  /** Default constructor */
  PhotonHitFPD();
  /** Constructor with position and energy */
  PhotonHitFPD(const Float_t x, const Float_t y, const Float_t e);
  /** Destructor */
  ~PhotonHitFPD() { }
  /** Reset all values */
  void Clear();
  /** Set (x, y) position and energy */
  ClassDef(PhotonHitFPD, 4)
};  // class PhotonHitFPD
}  // namespace PSUGlobals

#endif

#ifndef STFMSFITTEDPHOTON_H
#define STFMSFITTEDPHOTON_H

#include <Rtypes.h>

namespace PSUGlobals {//$NMSPC
/**
 Definition of a photon hit (SMD position info and reconstructed energy)
 */
struct StFmsFittedPhoton {
  Float_t energy;  ///< Fitted energy
  Float_t errEne;  ///< Energy fit error
  Float_t xPos;  ///< Fitted (relative) x-position
  Float_t errXPos;  ///< x-position fit error
  Float_t yPos;  ///< Fitted (relative) y-position
  Float_t errYPos;  ///< y-position fit error
  /** Constructor with optional position and energy */
  StFmsFittedPhoton(Float_t x = -1., Float_t y = -1., Float_t e = 0.,
                    Float_t xerr = -1., Float_t yerr = -1., Float_t eerr = -1.);
  /** Destructor */
  ~StFmsFittedPhoton() { }
  /** Reset all values to defaults */
  void Clear();
  ClassDef(StFmsFittedPhoton, 4)
};  // class StFmsFittedPhoton
}  // namespace PSUGlobals

#endif

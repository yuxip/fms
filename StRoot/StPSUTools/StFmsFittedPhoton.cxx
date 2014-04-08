#include "StFmsFittedPhoton.h"

namespace PSUGlobals {
PhotonHitFPD::PhotonHitFPD(const Float_t x, const Float_t y,
                           const Float_t e, const Float_t xerr,
                           const Float_t yerr, const Float_t eerr)
    : energy(e), errEne(eerr), xPos(x), errXPos(xerr), yPos(y), errYPos(yerr) {
}

void PhotonHitFPD::Clear() {
  energy = 0;
  errEne = xPos = errXPos = yPos = errYPos = -1;
}
}  // namespace PSUGlobals
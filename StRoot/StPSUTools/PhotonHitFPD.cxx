#include "PhotonHitFPD.h"

namespace PSUGlobals {
PhotonHitFPD::PhotonHitFPD(const Float_t x, const Float_t y, const Float_t e)
    : energy(e), errEne(-1.), xPos(x), errXPos(-1.), yPos(y), errYPos(-1.) {
}

void PhotonHitFPD::Clear() {
  energy = 0;
  errEne = xPos = errXPos = yPos = errYPos = -1;
}
}  // namespace PSUGlobals

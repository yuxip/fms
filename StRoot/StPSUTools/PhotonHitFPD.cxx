#include "PhotonHitFPD.h"

#include <cmath>
#include <iostream>

using namespace PSUGlobals;

PhotonHitFPD::PhotonHitFPD() {
  Clear();
}

PhotonHitFPD::PhotonHitFPD(const Float_t x, const Float_t y, const Float_t e) {
  energy = e;
  xPos = x;
  yPos = y;
}

void PhotonHitFPD::Clear(void) {
  energy = 0;
  errEne = xPos = errXPos = yPos = errYPos = -1 ;
}

#include "PhotonHitFPD.h"

#include <cmath>
#include <iostream>

using namespace PSUGlobals;

ClassImp(PhotonHitFPD);

PhotonHitFPD::PhotonHitFPD() {
  Clear();
}

TVector3 PhotonHitFPD::xyz(Float_t zPos) {
  return TVector3(xPos,yPos,zPos);
}

PhotonHitFPD::PhotonHitFPD(const Float_t x, const Float_t y, const Float_t e) {
  energy = e ;
  xPos   = x ;
  yPos   = y ;
}

void PhotonHitFPD::Clear(void) {
  nSmdMatch = -1;
  smdXCen = smdXSum = smdXSig = smdYCen = smdYSum = smdYSig = -1;
  energy = 0;
  errEne = xPos = errXPos = yPos = errYPos = -1 ;
}

void PhotonHitFPD::Set(const Float_t x, const Float_t y, const Float_t e) {
  energy = e ;
  xPos   = x ;
  yPos   = y ;
}

void PhotonHitFPD::SetSMD(const Float_t smd[6]) {
  smdXCen   = smd[0]  ;
  smdXSum   = smd[1]  ;
  smdXSig   = smd[2]  ;
  smdYCen   = smd[3]  ;
  smdYSum   = smd[4]  ;
  smdYSig   = smd[5]  ;
}

// calculate the distance between "this" and another "hit"
Float_t PhotonHitFPD::Distance(const PhotonHitFPD& hit) {
  Float_t dx = xPos - hit.xPos;
  Float_t dy = yPos - hit.yPos;
  return sqrt(dx * dx + dy* dy) ;
}

Float_t PhotonHitFPD::WeighterDistance(const PhotonHitFPD& hit) {
  Float_t dx = xPos - hit.xPos;
  Float_t dy = yPos - hit.yPos;
  Float_t sumE = energy + hit.energy ;
  Float_t weight = 4 * energy * hit.energy / (sumE * sumE);
  return sqrt(dx * dx + dy* dy) * weight;
}

void PhotonHitFPD::Print() {
  std::cout << "xPos:" << xPos << "(" << errXPos << ")\t yPos:";
  std::cout << yPos << "(" << errYPos << ")\t " << "energy: " << energy <<"("
    << errEne << ")" << "\n";
}

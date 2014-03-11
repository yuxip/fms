#ifndef PHOTONHITFPD_H
#define PHOTONHITFPD_H

#include <TObject.h>
#include <TVector3.h>

namespace PSUGlobals {//$NMSPC
/**
 Definition of a photon hit (SMD position info and reconstructed energy)
 */
class PhotonHitFPD : public TObject {
 public:
  Int_t nSmdMatch;  ///< -1: no SMD data, 0: not matched, 1: matched  
  Float_t smdXCen;
  Float_t smdXSum;
  Float_t smdXSig;
  Float_t smdYCen;
  Float_t smdYSum;
  Float_t smdYSig;
  Float_t energy;  ///< Fitted energy
  Float_t errEne;  ///< Energy fit error
  Float_t xPos;  ///< Fitted (relative) x-position
  Float_t errXPos;  ///< x-position fit error
  Float_t yPos;  ///< Fitted (relative) y-position
  Float_t errYPos;  ///< y-position fit error
  TVector3 xyz(Float_t zPos = 0);
  /** Default constructor */
  PhotonHitFPD();
  /** Constructor with position and energy */
  PhotonHitFPD(const Float_t x, const Float_t y, const Float_t e);
  /** Destructor */
  ~PhotonHitFPD() {};
  /** Reset all values */
  void Clear();
  /** Print values to standard output */
  void Print();
  /** Set (x, y) position and energy */
  void Set(const Float_t x, const Float_t y, const Float_t e);
  void SetSMD(const Float_t smd[6]);
  Float_t Distance(const PhotonHitFPD& hit);
  Float_t WeighterDistance(const PhotonHitFPD& hit);
  PhotonHitFPD& operator=(const PhotonHitFPD& rhs) {
    if (this != &rhs) {
      nSmdMatch = rhs.nSmdMatch;
      smdXCen = rhs.smdXCen;
      smdXSum = rhs.smdXSum;
      smdXSig = rhs.smdXSig;
      smdYCen = rhs.smdYCen;
      smdYSum = rhs.smdYSum;
      smdYSig = rhs.smdYSig;
      energy = rhs.energy;
      errEne = rhs.errEne;
      xPos = rhs.xPos;
      errXPos = rhs.errXPos;
      yPos = rhs.yPos;
      errYPos = rhs.errYPos;
    }  // if
    return *this;
  }
 private:
  ClassDef(PhotonHitFPD, 4)
};  // class PhotonHitFPD
}  // namespace PSUGlobals

#endif

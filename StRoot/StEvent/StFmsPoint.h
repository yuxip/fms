#ifndef StStFmsPoint_HH
#define StStFmsPoint_HH

#include "StObject.h"

#include "TLorentzVector.h"

//class that represents a fitted photon
//Yuxi Pan 03/31/2013
class StFmsPoint : public StObject {
 public:
  /** Constructor */
  StFmsPoint();
  // Use default copy constructor and assignment operator
  /** Destructor */
  ~StFmsPoint();
  void Print(const Option_t* opt = "") const;
  /** Sub-detector */
  UShort_t detectorId() const { return mDetectorId; }
  Float_t energy() const { return mEnergy; }
  Float_t x() const { return mX; }  //in cm
  Float_t y() const { return mY; }  //in cm
  Int_t id() const { return mId; }
  Int_t parentClusterId() const { return mParentClusterId; }
  Int_t nParentClusterPhotons() const { return mNParentClusterPhotons; }
  TLorentzVector fourMomentum() const { return mFourMomentum; }
  /** Set sub-detector */
  void setDetectorId(UShort_t detector) { mDetectorId = detector; }
  void setEnergy(Float_t energy) { mEnergy = energy; }
  void setX(Float_t xpos) { mX = xpos; }
  void setY(Float_t ypos) { mY = ypos; }
  void setId(Int_t phid) { mId = phid; }
  void setParentClusterId(Int_t cluid) { mParentClusterId = cluid; }
  void setNParentClusterPhotons(Int_t nclph) { mNParentClusterPhotons = nclph; }
  void setFourMomentum(TLorentzVector p4) { mFourMomentum = p4; }

 protected:
  UShort_t mDetectorId;  ///< Detector starts from 1
  Float_t mEnergy;  ///< fitted energy
  Float_t mX;  ///< fitted(relative) x-position
  Float_t mY;  ///< fitted(relative) y-position
  Int_t mId;  ///< photon id within event, also include det12 info
  Int_t mParentClusterId;  ///< id of the parent cluster
  Int_t mNParentClusterPhotons;  ///< # of photons in by the parent cluster
  TLorentzVector mFourMomentum;
  ClassDef(StFmsPoint, 1)
};

ostream& operator<<(ostream&, const StFmsPoint&);

#endif

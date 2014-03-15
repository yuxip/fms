#ifndef Geom_
#define Geom_

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include <iostream>
#include <map>
#include <vector>

#include "TObject.h"
#include "TString.h"
#include "TObjArray.h"
#include "TVector3.h"

class fmsDetectorPosition_st;

namespace PSUGlobals {//$NMSPC
class Geom : public TObject {
 public:
  Geom();
  bool InitDBGeom();
  ~Geom();
  Float_t z(Int_t detectorId) const;
  Float_t xOffset(Int_t detectorId) const;
  Float_t yOffset(Int_t detectorId) const;
  /** Return [x, y] tower widths in cm */
  std::vector<Float_t> towerWidths(Int_t detectorId) const;
  /**
   Return the position information of a detector
   
   Return NULL if the detector ID is invalid, or database information is
   unavailable.
   */
  const fmsDetectorPosition_st* find(Int_t detectorId) const;

 private:
  typedef std::map<int, fmsDetectorPosition_st*> Table;
  Table mPositions;
  ClassDef(Geom,3);
};
}  // namespace PSUGlobals
#endif


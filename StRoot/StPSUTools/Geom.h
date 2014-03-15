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
class StFmsDbMaker;

namespace PSUGlobals {//$NMSPC
class Geom : public TObject {
 public:
  Geom();
  ~Geom();
  /**
   Initialise geometry from the FMS database
   
   If the argument is NULL, attempt to locate an StFmsDbMaker in the current
   chain and use that.
   Return true if the geometry is initialised, false if it is not.
   */
  Bool_t initialize(StFmsDbMaker* fmsDbMaker);
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


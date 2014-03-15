#ifndef Geom_
#define Geom_

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include <iostream>
#include <map>

#include "TObject.h"
#include "TString.h"
#include "TObjArray.h"
#include "TVector3.h"

class fmsDetectorPosition_st;

/*
  If the geom.txt file has 16 elements in the FpdTowWidth field it
  will be assumed that the file corresponds to the FMS data set. 
  FMSGeom will be set to true.
 */
namespace PSUGlobals {//$NMSPC
class Geom : public TObject {
 public:
  Geom();
  bool InitDBGeom();
  ~Geom();
  Bool_t FMSGeom;  // Only retained for backward compatibility. Can be deleted
                   // once other code's dependence on it is removed.
  const Float_t* ZFPD(Int_t detectorId) const;
  const Float_t* xOffset(Int_t detectorId) const;
  const Float_t* yOffset(Int_t detectorId) const;
  const Float_t* FpdTowWid(Int_t detectorId) const;
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


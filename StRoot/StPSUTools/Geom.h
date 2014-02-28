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
  Float_t* ZFPD(Int_t, Int_t );
  Float_t* xOffset(Int_t , Int_t );
  Float_t* yOffset(Int_t, Int_t );
  Float_t* FpdTowWid(Int_t , Int_t );

 private:
  /*
   This is a temporary function to get detector ID from east/west (1 or 2) and
   nstb (1-4), until we change the rest of the code to work with subdetectors
   via detector ID instead of ew/nstb. Returns an ID in the range [8, 11].
   */
  int ewNstbToDetectorId(int ew, int nstb);
  fmsDetectorPosition_st* find(int detectorId);
  typedef std::map<int, fmsDetectorPosition_st*> Table;
  Table mPositions;
  ClassDef(Geom,3);
};
}  // namespace PSUGlobals
#endif


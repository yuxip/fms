#ifndef Geom_
#define Geom_
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <iostream>
#include "TObject.h"
#include "TString.h"
#include "TObjArray.h"
#include "TVector3.h"
/*
  If the geom.txt file has 16 elements in the FpdTowWidth field it
  will be assumed that the file corresponds to the FMS data set. 
  FMSGeom will be set to true.
 */
namespace PSUGlobals {//$NMSPC
class Geom : public TObject
{
 public:
  Geom();
  void InitDBGeom();
  ~Geom();
  void  Print();
  void  printgeom(const char* file="geomDebug.txt");
  Bool_t Defined;
  Bool_t FMSGeom;
  TString datatype[14];
  Int_t dataSize[14];
  Float_t data[14][16];
  Int_t getNSTB(Float_t gx, Float_t gy);
  Float_t* ZFPD(Int_t, Int_t );
  Float_t* xOffset(Int_t , Int_t );
  Float_t* yOffset(Int_t, Int_t );
  Float_t* FpdTowWid(Int_t , Int_t );
  TVector3 GlobalXYZ(Int_t EW12,Int_t NSTB16,TVector3 Localxyz);
  TVector3 LocalXYZ(Int_t EW12,Int_t NSTB16,TVector3 Globalxyz,Bool_t BlockUnit=false);
 private:
    ClassDef(Geom,3);
    Int_t Date;
};
}
#endif


#ifndef CalibStr_
#define CalibStr_
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "TObject.h"
#include "TObjArray.h"
#include "TMatrix.h"
namespace PSUGlobals {//$NMSPC
class CalibStr :public TObject
{
 public:
  CalibStr();
    CalibStr(Int_t, const char*); 
    Int_t CalibStrBuild(Int_t,const  char*); //Modified to have logic value...
    Int_t MatrixDimensionRows[6][2];
    Int_t MatrixDimensionCols[6][2];
    Int_t SetValue(Int_t IEW,Int_t INSTB, Int_t Irow0, Int_t Icol0,Float_t val);
    TObjArray* M_array;
    TString filename;
    TString OutName;
    TMatrix* tm(Int_t ew,Int_t nstb)
      {return (TMatrix*)( M_array->At((ew-1)*6+nstb-1));};
    Int_t error;
    Int_t GetDBCalibration(); //!DataBase!
    Float_t ZERO;
    Float_t* value(Int_t index,Int_t i_ns,Int_t i_ew);
    Int_t SetValue(Int_t ,Int_t ,Int_t ,Float_t );
    Int_t UpdateFile(Int_t);
    Float_t GetValue(Int_t IEW,Int_t INSTB, Int_t Irow0, Int_t Icol0);
    Bool_t Compare(Int_t IEW,Int_t INSTB,CalibStr* other);
 private:
    Int_t CurrentRunNumber;
    ClassDef(CalibStr,2);
};
}
#endif

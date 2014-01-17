#ifndef WASEXTERNAL_H
#define WASEXTERNAL_H
#include "TObjArray.h"
#include "TF2.h"
#include "TMinuit.h"
#include "TRandom.h"
#include "TMatrix.h"
#define MAX_NUMB_PHOTONS 7

namespace PSUGlobals {//$NMSPC
class WasExternal: public TObject
{
 public:
  WasExternal()
    {
      Double_t step0[3*MAX_NUMB_PHOTONS+1]= {0, 0.1,0.1,0.2, 0.1,0.1,0.2, 0.1,
					     0.1,0.2, 0.1,0.1,0.2, 0.1,0.1,0.2, 
					     0.1,0.1,0.2, 0.1,0.1,0.2};
      for(int j=0;j<3*MAX_NUMB_PHOTONS+1;j++)step[j]=step0[j];
      fcnSS=0;
      UseThis_ab=false;
      UseThis_Err=false;
      NoCatag=false;
      DoGlobal=true;
      ForceCatag=0;
      widLG[0]=3.81;
      widLG[1]=3.81;
      dev=0;
      dchi2=0;
      Force2Mass=-1.;
      for(int j=0;j<10;j++)FreeGlobals[j]=-1;
      SubClu2=false;
      BlockFit=false;
      UseEDepCorrection=true;
      /*
      EDepCorrection=new TF1("EDepCorrection","(1-.11*exp(-(x)/[0])-.23*exp(-(x)/[1]))",1,250);
      EDepCorrection->SetParameter(0,5);
      EDepCorrection->SetParameter(1,70);
      */
      EDepCorrection=new TF1("EDepCorrection","(1.3-.15*exp(-(x)/[0])-.6*exp(-(x)/[1]))",1,250);
      EDepCorrection->SetParameter(0,10.);
      EDepCorrection->SetParameter(1,70);
      
    };
  
  ~WasExternal(){};
  //temporary SH
  TF1* EDepCorrection;
  Bool_t  UseEDepCorrection;
  Bool_t SubClu2;
  TMatrix* dev;
  TMatrix*  dchi2;
  //
  TObjArray* tow2Fit;
  TF2* fcnSS;
  UInt_t showerShapeFunc;
  UInt_t choiceChi2;
  Float_t showerWidthX;
  Float_t showerWidthY;
  Float_t widLG[2];//glass width X,Y
  UInt_t useLimitForPara;
  Bool_t DoGlobal;
  Bool_t NoCatag;
  //  Bool_t NoCatag2;
  Int_t ForceCatag;
  Double_t errFactor;
  Double_t errQ;
  Float_t Power1;
  Float_t Power2;
  TRandom *myRand;
  Double_t step[3*MAX_NUMB_PHOTONS+1];
  Bool_t UseThis_ab;
  Bool_t UseThis_Err;
  Float_t Force2Mass;
  Bool_t BlockFit;
  Float_t a1;
  Float_t a2;
  Float_t a3;
  Float_t b1;
  Float_t b2;
  Float_t b3;
  Double_t FreeGlobals[10];//there to be used;
  Float_t Energy_study;
 private:
 ClassDef(WasExternal,3);

};
}
#endif


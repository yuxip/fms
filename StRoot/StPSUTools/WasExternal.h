#ifndef WASEXTERNAL_H
#define WASEXTERNAL_H
#include "TObjArray.h"
#include "TF2.h"
#include "TMinuit.h"
#include "TRandom.h"
#include "TMatrix.h"
#define MAX_NUMB_PHOTONS 7

namespace PSUGlobals {//$NMSPC
/**
 This class collects a large number of generally unrelated flags and values.
 
 They presumably used to be globals in an earlier iteration of the code, hence
 the name "WasExternal"?
 \note
 All these values are in a sense really members of FitTower; that class has
 a static member WasExternal, which is the only instance of this class that is
 ever created. It would make more sense if these variables were moved to
 FitTower.
 */
class WasExternal: public TObject
{
 public:
  WasExternal()
    {
      a1 = 1.070804;
      a2 = 0.167773;
      a3 = -0.238578;
      b1 = 0.535845;
      b2 = 0.850233;
      b3 = 2.382637;
      Double_t step0[3*MAX_NUMB_PHOTONS+1]= {0, 0.1,0.1,0.2, 0.1,0.1,0.2, 0.1,
					     0.1,0.2, 0.1,0.1,0.2, 0.1,0.1,0.2, 
					     0.1,0.1,0.2, 0.1,0.1,0.2};
      for(int j=0;j<3*MAX_NUMB_PHOTONS+1;j++)step[j]=step0[j];
      fcnSS=0;
      widLG[0]=3.81;
      widLG[1]=3.81;
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
  /**
   \note
   Response function for nonlinear energy correction, based on cerenkov studies.
   It should be part of another class, but is used in both HitCluster and Yiqun.
   We therefore have to decide how to pass it to both HitCluster and Yiqun.
   Note, I'm not actually sure if the HitCluster usage is ever called… so
   could just be needed for Yiqun.
   */
  TF1* EDepCorrection;
  /**
   \note
   An array of TowerFPD objects to fit in the FitTower fitting functions. It
   makes more sense to me to just pass the tower array directly to FitTower as
   arguments to the relevant fitter functions, as this array isn't needed
   anywhere else.
   */
  TObjArray* tow2Fit;
  /**
   \note
   This function is allocated in FitTower, never deleted (maybe ROOT does that),
   then copied by FitTower (not sure what the need is for that). It is also used
   once by Yiqun, but as Yiqun has a member FitTower we could just make it a
   member of FitTower.
   */
  TF2* fcnSS;
  /**
   \note
   "Width lead-glass" i.e. cell (x, y) widths. It is set via Geom, from the
   database, which is good. There are functions to set it via FitTower, but I
   don't think these are needed, as it's only ever needed internally in Yiqun.
   As Yiqun already has a member Geom, we can just access this information from
   there.
   */
  Float_t widLG[2];//glass width X,Y
  /**
   \note
   Is only used in Yiqun, so it should located be there. Also, it is a TRandom
   object, which the ROOT developers themselves say is a bad generator. It may
   be that STAR specifies its own random generator, in which case we should use
   that, or else we just use the ROOT global gRandom, which is a much superior
   TRandom3; there doesn't seem a need to define a separate random generator
   for Yiqun.
   */
  TRandom *myRand;
  /**
   \note
   Is some kind of array used in the fitting procedure, used by both Yiqun and
   FitTower (FitTower copies the WasExternal array but doesn't modify it; Yiqun
   uses the WasExternal version directly; so both appear to use the same
   numbers). No need to have it external; probably best as a member of FitTower
   with an accessor method so Yiqun can get hold of it.
   */
  Double_t step[3*MAX_NUMB_PHOTONS+1];
  /**
   \note
   Only used in FitTower, and always a constant value, set in StFmsHitMaker via
   FitTower::Setwe_ab(). It is one of the parameters passed to fcnSS. If it is
   only ever going to be constant there doesn't seem a need to store it as a
   variable.
   */
  Float_t a1;
  /**
   \note
   See a1.
   */
  Float_t a2;
  /**
   \note
   See a1.
   */
  Float_t a3;
  /**
   \note
   See a1.
   */
  Float_t b1;
  /**
   \note
   See a1.
   */
  Float_t b2;
  /**
   \note
   See a1.
   */
  Float_t b3;
 private:
 ClassDef(WasExternal,3);

};
}
#endif


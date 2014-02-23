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
   "Width lead-glass" i.e. cell (x, y) widths. It is set via Geom, from the
   database, which is good. There are functions to set it via FitTower, but I
   don't think these are needed, as it's only ever needed internally in Yiqun.
   As Yiqun already has a member Geom, we can just access this information from
   there.
   */
  Float_t widLG[2];//glass width X,Y
 private:
 ClassDef(WasExternal,3);

};
}
#endif


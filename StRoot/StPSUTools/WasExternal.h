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
 private:
 ClassDef(WasExternal,3);

};
}
#endif


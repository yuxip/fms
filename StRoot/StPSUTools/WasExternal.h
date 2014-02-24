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
    };
  
  ~WasExternal(){};
  //temporary SH
 private:
 ClassDef(WasExternal,3);

};
}
#endif


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
  /**
   \note
   The function used when UseEDepCorrection is true.
   We could replace the flag with only allocating the function if it is
   required, and use the pointer as the flag. See also notes for
   UseEDepCorrection.
   */
  TF1* EDepCorrection;
  /**
   \note
   Defaults to true and is never changed, so maybe we scrap code that relies on
   it evaluating as false. It is used in both HitCluster and Yiqun, and is set
   via FitTower::Setwe_EDepCor(), but that function is never called. If it is
   retained we have to decide how to pass it to both HitCluster and Yiqun.
   Also, I'm not actually sure if the HitCluster usage is ever called… so
   could just be needed for Yiqun.
   */
  Bool_t  UseEDepCorrection;
  /**
   \note
   It's never used, so it can be removed.
   */
  Bool_t SubClu2;
  /**
   \note
   Is only used in FitTower, despite being allocated and deleted in Yiqun. It
   might as well just be a member of FitTower - in fact, it's only used in
   FitTower::Fcn1(), so we could just as well create it locally there.
   */
  TMatrix* dev;
  /**
   \note
   Basically the same deal as dev. Actually, even though I can see these
   matrices being filled, they never seem to be used. So probably just
   remove them.
   */
  TMatrix*  dchi2;
  //
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
   Is set to 0 at initialisation and never changed. There is code in FitTower
   that tests if showerShapeFunc==1, but this can never execute. We could make
   it a member of FitTower, but we could really just get rid of both it and any
   code that relies on it being non-zero.
   */
  UInt_t showerShapeFunc;
  /**
   \note
   Is initialised to 2. There is code in FitTower that depends on choiceChi2==2
   or choiceChi2==1, but only the "2" version will ever execute. It is also
   referenced in Yiqun, but never used for anything, so it can be made a member
   of FitTower, if it's needed at all (which I don't think it is).
   */
  UInt_t choiceChi2;
  /**
   \note
   Is initialised as 1, never changed, and is only used to divide another
   variable (obviously, to no effect). Does not seem to be needed.
   */
  Float_t showerWidthX;
  /**
   \note
   See showerWidthX
   */
  Float_t showerWidthY;
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
   Only used in Yiqun, initialised to 1. But, even though there is
   code in Yiqun::FitOnePhoton(), Yiqun::FitTwoPhoton(), Yiqun::GlobalFit()
   and Yiqun::Fit2PhotonClust() that have if(useLimitForPara ==<different
   values>), it is never changed from 1, so there doesn't seem much use for
   this flag at all.
   */
  UInt_t useLimitForPara;
  /**
   \note
   Is initialised to true and never changed. Code with if(DoGlobal) in Yiqun
   is always executed, so it seems pointless.
   */
  Bool_t DoGlobal;
  /**
   \note
   Is initialised to false and never changed. The only functions that use it
   are FitTower::SetNoCatag(), which sets it (and is never called), and
   Yiqun::FitEvent(), which has an if(NoCatag) clause, which is never executed,
   so it makes sense to get rid of it.
   */
  Bool_t NoCatag;
  //  Bool_t NoCatag2;
  /**
   \note
   Initialised to 0, never changes because FitTower::SetNoCatag() isn't called
   (that is the only place it is set), and is only used in case if(NoCatag)
   evaluates true, which it never does. Doesn't seem to serve a function.
   */
  Int_t ForceCatag;
  /**
   \note
   Always uses a default value because FitTower::Setwe_ErrFactors() is never
   called, so UseThis_Err always evaluates false. However it *is* used, so we
   need to keep it (at least the value, if not the actual variable). It is
   only used in FitTower, so it should be placed there.
   */
  Double_t errFactor;
  /**
   \note
   See errFactor.
   */
  Double_t errQ;
  /**
   \note
   Is used in FitTower::Fcn1(). It is set via FitTower::Setwe_ErrFactors(), but
   that function is never called in the code here (in fact, I can't find it
   being called in Steve's original package either). It's actually only ever
   used if we do Setwe_ErrFactors(), and we don't, so maybe it's not even
   needed…
   */
  Float_t Power1;
  /**
   \note
   See Power1.
   */
  Float_t Power2;
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
   This flag is only used by FitTower, and so can be part of it. It is always
   true in the current form of the code (because it is set so in StFmsHitMaker),
   so it's not really needed.
   */
  Bool_t UseThis_ab;
  /**
   \note
   Is only referenced by a function that is never used:
   FitTower::Setwe_ErrFactors(). It is always false in the current form of the
   code (this is the default value; it's not modified anywhere), so it's not
   really needed.
   */
  Bool_t UseThis_Err;
  /**
   \note
   Not clear what it is for. It can be set via FitTower::SetForceMass(), but
   that function is never called so it is never changed from the default value
   of -1. It is only used in FitTower::Fcn2(), so it can be made a member of
   FitTower, if indeed it is needed at all.
   */
  Float_t Force2Mass;
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
  /**
   \note
   Apparently a collection of global values available to the user as needed.
   They are never used, only set via FitTower::SetFreeGlobals(), which is never
   called. Might as well scrap them.
   */
  Double_t FreeGlobals[10];//there to be used;
 private:
 ClassDef(WasExternal,3);

};
}
#endif


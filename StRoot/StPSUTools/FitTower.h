#ifndef FIT_TOWER_H
#define FIT_TOWER_H

#include <iostream>

#include <TMinuit.h>
#include <TObject.h>

#include "StPSUTools/Geom.h"
#include "StPSUTools/TowerFPD.h"
#include "StPSUTools/TowerUtil.h"

class TF2;

namespace PSUGlobals {//$NMSPC
/**
 \todo
 It may be safer to make FitTower a singleton class, or something like that.
 The shower shape fit function is shared across all objects (by necessity, in
 order to interface with TMinuit), but each object updates the function itself
 with different parameters. Therefore bad things would happen if there were
 more than one object in existence at any time. There isn't ever more than one
 instance created in this code, but I think it would be good to enforce that.
 */
class FitTower : public TObject {
 public:
  FitTower(Geom* pgeom,Int_t iew,Int_t nstb);
  FitTower() { }
  ~FitTower();
  TF2* GetFunctShowShape();
  Int_t Fit(const Double_t *par, const Double_t *step, 
            const Double_t *low, const Double_t *up);
  Int_t Fit2Pin1Clust(const Double_t *para, const Double_t *step, 
                      const Double_t *low, const Double_t *up);
  static Double_t GGams(Double_t *x, Double_t *par);
  TowerUtil* pTowerUtil;
  TMinuit* fMn;  // Minuit fitter
  static const Int_t MAX_NUMB_PHOTONS = 7;
  Double_t step[3*MAX_NUMB_PHOTONS+1];
  static TObjArray* tow2Fit;

 private:  
  static Double_t FGams(Double_t *x, Double_t *par);
  static void  Fcn1(Int_t & npar, Double_t *grad,  
                    Double_t &fval, Double_t *par, Int_t iflag);
  static void Fcn2(Int_t & nparam, Double_t *grad, 
                   Double_t &fval, Double_t *param, Int_t iflag);
  void SetFCN(void (*fcn)(Int_t &, Double_t *, Double_t &, Double_t *, Int_t));
  void SetStep();
  Double_t fTWidthCM;  // width of one lead glass module
  Double_t fTXWidthCM;  // width of one lead glass module in X
  Double_t fTYWidthCM;  // width of one lead glass module in Y
  Double_t fNumbPhotons;  // number of photons to be fitted: Minuit wants a Double_t!
  Double_t fFitPara[3*MAX_NUMB_PHOTONS+1];  // Minuit fit parameter
  static Float_t widLG[2];  //glass width X,Y
  ClassDef(FitTower,4);
};  // class FitTower
}  // namespace PSUGlobals
#endif

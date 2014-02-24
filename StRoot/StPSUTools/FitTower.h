#ifndef FIT_TOWER_H
#define FIT_TOWER_H

#include <iostream>
#include  "TObject.h"
#include "Geom.h"
#include "TMatrix.h"
#include "TMinuit.h"
#include "WasExternal.h"
#include "TF2.h"
#include "TowerUtil.h"
#include "TMath.h"
#include "TowerFPD.h"
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
class FitTower : public TObject
{
 public:
  
  TMinuit*    fMn;    // Minuit fitter ma
  FitTower(TMatrix* ptm, Geom* pgeom,Int_t iew,Int_t nstb);
  /*
    constructor removed by SH  8/2009
  FitTower(const Int_t dim[2], const Double_t wd, TF2 *ssFunct);
  */
  FitTower(){};

  ~FitTower();
  static TObjArray* tow2Fit;
  static  WasExternal we;
  Double_t step[3*MAX_NUMB_PHOTONS+1];
  TowerUtil* pTowerUtil;
  static Double_t GGams(Double_t *x, Double_t *par);
  static Double_t FGams(Double_t *x, Double_t *par);
  
  static void  Fcn1(Int_t & npar, Double_t *grad,  
		    Double_t &fval, Double_t *par, Int_t iflag);
  static void Fcn2(Int_t & nparam, Double_t *grad, 
		   Double_t &fval, Double_t *param, Int_t iflag);
  void SetRow(const Int_t row) {fRow=row;};
  void SetCol(const Int_t col) {fCol=col;};
  void SetTWidthCM(Float_t tw) 
    {
      fTWidthCM = tw;
      fTXWidthCM=fTYWidthCM=tw;
      FitTower::widLG[0] = FitTower:: widLG[1] = tw;
    };
  void SetXYTWidthCM(Float_t*  ptw)
    {
      fTXWidthCM = FitTower::widLG[0] = ptw[0];
      fTYWidthCM = FitTower::widLG[1] = ptw[1];
    };
  void SetFCN(void (*fcn)(Int_t &, Double_t *, Double_t &, Double_t *, Int_t));
  void SetNumberPhoton(const Int_t nP);

  void GetRow(Int_t &row) {row=fRow;};
  void GetCol(Int_t &col) {col=fCol;};
  void GetParameter(Double_t *ap){/* not implemeted*/};
  void GetTWidthCM(Double_t &tw) {tw=fTWidthCM;};
  TF2* GetFunctShowShape();
  Int_t Fit(const Double_t *par, const Double_t *step, 
	    const Double_t *low, const Double_t *up);
  Int_t Fit2Pin1Clust(const Double_t *para, const Double_t *step, 
		      const Double_t *low, const Double_t *up);
 private:	
  
  ClassDef(FitTower,4);
  void SetStep();
  
  Int_t       fCol;           // number of columns for lead glass array
  Int_t       fRow;           // number of rows for lead glass array
  Double_t    fTWidthCM;      // width of one lead glass module
  Double_t    fTXWidthCM;      // width of one lead glass module in X
  Double_t    fTYWidthCM;      // width of one lead glass module in Y
  Double_t    fNumbPhotons;   // number of photons to be fitted: Minuit wants a Double_t!
  Double_t    fFitPara[3*MAX_NUMB_PHOTONS+1];  // Minuit fit parameter
  static Float_t widLG[2];//glass width X,Y
};
}
#endif

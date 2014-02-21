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
  TF2* fcnSS;
  TObjArray* tow2Fit;
  static  WasExternal we;
  static Bool_t Setwe_EDepCor(Bool_t useEdepCor=true);
  static Bool_t Setwe_ab(Float_t a1=.8,Float_t a2=.3,Float_t a3=-.1,Float_t b1=.8,Float_t b2=.2,Float_t b3=7.6);
  static Bool_t Setwe_ErrFactors(float errQ,float errFactor,Float_t p1= 1.95,Float_t p2=2.72,float energy_study=100.);
  static Bool_t SetDoGlobal(Bool_t sdg=true);
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
      we.widLG[0]=we.widLG[1]=tw;
    };
  void SetXYTWidthCM(Float_t*  ptw)
    {
      fTXWidthCM=we.widLG[0]=ptw[0];
      fTYWidthCM=we.widLG[1]=ptw[1];
    };
  void SetFCN(void (*fcn)(Int_t &, Double_t *, Double_t &, Double_t *, Int_t));
  void SetFunctShowShape(TF2 *fSS) {fFunctSS = fSS;};
  void SetNumberPhoton(const Int_t nP);

  void GetRow(Int_t &row) {row=fRow;};
  void GetCol(Int_t &col) {col=fCol;};
  void GetParameter(Double_t *ap){/* not implemeted*/};
  void GetTWidthCM(Double_t &tw) {tw=fTWidthCM;};
  TF2* GetFunctShowShape(void) {return fFunctSS;};
  Int_t Fit(const Double_t *par, const Double_t *step, 
	    const Double_t *low, const Double_t *up);
  Int_t Fit2Pin1Clust(const Double_t *para, const Double_t *step, 
		      const Double_t *low, const Double_t *up);
  static Bool_t SetForceMass(Float_t fmass=-1.);
 private:	
  
  ClassDef(FitTower,4);
  void SetStep();
  
  Int_t       fCol;           // number of columns for lead glass array
  Int_t       fRow;           // number of rows for lead glass array
  Double_t    fTWidthCM;      // width of one lead glass module
  Double_t    fTXWidthCM;      // width of one lead glass module in X
  Double_t    fTYWidthCM;      // width of one lead glass module in Y
  TF2 *       fFunctSS;       // shower-shape function
  Double_t    fNumbPhotons;   // number of photons to be fitted: Minuit wants a Double_t!
  Double_t    fFitPara[3*MAX_NUMB_PHOTONS+1];  // Minuit fit parameter

};
}
#endif

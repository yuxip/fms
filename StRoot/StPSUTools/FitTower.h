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
  /** Constructor using detector geometry */
  FitTower(Geom* pgeom,Int_t iew,Int_t nstb);
  /** Default constructor */
  FitTower() { }
  /** Destructor */
  ~FitTower();
  /**
   Return the shower shape function
   
   The shower shape gives the fractional energy deposition by a photon in a
   tower as a function of the distance of the photon from the tower center.
   */
  TF2* GetFunctShowShape();
  /**
   Fit photons to the list of towers.
   
   par, step, low and up are all arrays of size 3N+1 for an N-photon fit
    - par: start values for the fit variables
    - step: step size when fitting
    - low: lower bound on fit variable
    - up: upper bound on fit varaible
   The first element in each array refers to the number of photons to fit.
   The next 3 are the x position, y position and energy of a photon.
   e.g. for 1 photon [1, x0, y0, E0]; for 2 photons [2, x0, y0, E0, x1, y1, E1]
   */
  Int_t Fit(const Double_t *par, const Double_t *step, 
            const Double_t *low, const Double_t *up);
  /**
   Specialized fit function for exactly 2-photon fit
   
   Arguments are as for Fit(). However the input arrays should always have 7
   (3 * 2 photons + 1) elements, and have a special meaning here:
    - 0: still a constant parameter, should be set to 2 for 2-photon fitting
    - 1: xPi      (x-position of pi^0)
    - 2: yPi      (y-position of pi^0)
    - 3: d_gg     (distance between 2 photons)
    - 4: theta    (theta angle of displacement vector from photon 2 to photon 1)
    - 5: z_gg     (this z_gg can go from -1 to +1, so we do not set E1>E2)
    - 6: E_gg     (total energy of two photons)
   */
  Int_t Fit2Pin1Clust(const Double_t *para, const Double_t *step, 
                      const Double_t *low, const Double_t *up);
  /**
   Shower shape function for use with TF2
   
   Yields the fraction of energy deposited in a tower by a photon as a function
   of the distance of the photon from the tower center.
   Arguments:
    - x: array with x[0] = x distance from tower center, x[1] = y distance
    - par: array of fit parameters (fixed, derived from FMS studies)
   Integrates F(x,y) over a tower, with F(x,y) defined as here:
   https://drupal.star.bnl.gov/STAR/blog/leun/2010/aug/02/fms-meeting-20100802
   \todo Provide LaTeX math function in documentation
   */
  static Double_t GGams(Double_t *x, Double_t *par);
  TowerUtil* pTowerUtil;
  TMinuit* fMn;  // Minuit fitter
  static const Int_t MAX_NUMB_PHOTONS = 7;
  Double_t step[3*MAX_NUMB_PHOTONS+1];
  static TObjArray* tow2Fit;

 private:  
  /**
   Shower shape helper function
   
   Used by GGams to evaluate one component of the shower shape function.
   Evaluates F(x,y) as defined here:
   https://drupal.star.bnl.gov/STAR/blog/leun/2010/aug/02/fms-meeting-20100802
   \todo Provide LaTeX math function in documentation
   */
  static Double_t FGams(Double_t *x, Double_t *par);
  /**
   Minuit minimization function for Fit() routine
   
   For the purpose of this function and a description of its arguments see
   https://wwwasdoc.web.cern.ch/wwwasdoc/minuit/node14.html
   For its use in ROOT via TMinuit see
   http://root.cern.ch/root/htmldoc/TMinuit.html#TMinuit:SetFCN
   */
  static void Fcn1(Int_t & npar, Double_t *grad,  
                    Double_t &fval, Double_t *par, Int_t iflag);
  /**
   Minuit minimization function for Fit2Pin1Clust() routine
   
   See also Fcn1().
   */
  static void Fcn2(Int_t & nparam, Double_t *grad, 
                   Double_t &fval, Double_t *param, Int_t iflag);
  void SetFCN(void (*fcn)(Int_t &, Double_t *, Double_t &, Double_t *, Int_t));
  /**
   Set default steps for fitting functions Fit() and Fit2Pin1Clust()
   */
  void SetStep();
  Double_t fTWidthCM;  ///< width of one lead glass module
  Double_t fTXWidthCM;  ///< width of one lead glass module in X
  Double_t fTYWidthCM;  ///< width of one lead glass module in Y
  Double_t fNumbPhotons;  ///< number of photons to be fitted: Minuit wants a Double_t!
  Double_t fFitPara[3*MAX_NUMB_PHOTONS+1];  ///< Minuit fit parameter
  static Float_t widLG[2];  ///< glass width X,Y
  ClassDef(FitTower,4);
};  // class FitTower
}  // namespace PSUGlobals
#endif

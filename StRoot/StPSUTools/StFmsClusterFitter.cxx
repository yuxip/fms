#include "StFmsClusterFitter.h"

#include <vector>

#include <TF2.h>
#include <TMath.h>

#include "StRoot/St_base/StMessMgr.h"
#include "StEvent/StFmsHit.h"

#include "StPSUTools/StFmsGeometry.h"
#include "StPSUTools/StFmsTower.h"

namespace {
Int_t numbPara = 10;
TF2 showerShapeFitFunction("showerShapeFitFunction",
                       &PSUGlobals::StFmsClusterFitter::energyDepositionInTower,
                      -25.0, 25.0, -25.0, 25.0, numbPara);
}  // unnamed namespace

namespace PSUGlobals {
// Instantiate static members
Float_t StFmsClusterFitter::mTowerWidthXY[2];
TObjArray* StFmsClusterFitter::mTowers(NULL);

TF2* StFmsClusterFitter::showerShapeFunction() {
  return &showerShapeFitFunction;
}

StFmsClusterFitter::StFmsClusterFitter(StFmsGeometry* pgeom, Int_t detectorId)
    : mMinuit(3 * MAX_NUMB_PHOTONS + 1) {
  // Set steps for Minuit fitting
  const Double_t step[3 * MAX_NUMB_PHOTONS + 1]= {
    0.0, 0.1, 0.1, 0.2, 0.1, 0.1, 0.2, 0.1, 0.1, 0.2, 0.1,
    0.1, 0.2, 0.1, 0.1, 0.2, 0.1, 0.1, 0.2, 0.1, 0.1, 0.2
  };
  for(int j = 0; j < 3 * MAX_NUMB_PHOTONS + 1; j++) {
    mSteps[j] = step[j];
  }  // for
  std::vector<Float_t> towerWidth = pgeom->towerWidths(detectorId);
  mTowerWidth = towerWidth[0];
  StFmsClusterFitter::mTowerWidthXY[0] = towerWidth[0];
  StFmsClusterFitter::mTowerWidthXY[1] = towerWidth[1];
  Double_t para[numbPara];
  para[0] = mTowerWidth;
  para[1] = 1.070804;
  para[2] = 0.167773;
  para[3] = -0.238578;
  para[4] = 0.535845;
  para[5] = 0.850233;
  para[6] = 2.382637;
  para[7] = 0.0;
  para[8] = 0.0;
  para[9] = 1.0;
  showerShapeFitFunction.SetParameters(para); 
  // create a Minuit instance
  mMinuit.SetPrintLevel(-1);  // Quiet, including suppression of warnings
}

StFmsClusterFitter::~StFmsClusterFitter() { }

// Calculate fractional photon energy deposition in a tower based on its (x, y)
// position relative to the tower center
Double_t StFmsClusterFitter::energyDepositionDistribution(
    Double_t* xy,
    Double_t* parameters) {
  Double_t f = 0;
  Double_t x = xy[0];
  Double_t y = xy[1];
  for (Int_t i = 1; i <= 3; i++) {
    // The parameter array has 10 elements, but we only use 6
    // Parameters 1 to 6 are a1, a2, a3, b1, b2, b3 as defined in
    // https://drupal.star.bnl.gov/STAR/blog/leun/2010/aug/02/fms-meeting-20100802
    Double_t a = parameters[i];
    Double_t b = parameters[i + 3];
    f += a * atan(x * y / (b * sqrt(b * b + x * x + y * y)));
  }  // for
  return f / TMath::TwoPi();
}

// xy array contains (x, y) position of the photon relative to the tower center
Double_t StFmsClusterFitter::energyDepositionInTower(Double_t* xy,
                                                     Double_t* para) {
  Double_t gg(0);
  // Calculate the energy deposited in a tower
  // Evaluate energyDepositionDistribution at x+/-d/2 and y+/-d/2, for tower
  // width d. The double-loop below is equivalent to
  // F(x+d/2, y+d/2) + F(x-d/2, y-d/2) - F(x-d/2, y+d/2) - F(x+d/2, y-d/2)
  for (Int_t ix = 0; ix < 2; ix++) {
    for (Int_t iy = 0; iy < 2; iy++) {
      Double_t signX = pow(-1.0, ix);  // 1 or -1
      Double_t signY = pow(-1.0, iy);  // 1 or -1
      // para[0] is the cell width
      // para[7] and para[8] are offsets that are normally zero
      Double_t s[2];
      s[0] = xy[0] - para[7] + signX * para[0] / 2.0;  // x +/- d/2
      s[1] = xy[1] - para[8] + signY * para[0] / 2.0;  // y +/- d/2
      gg += signX * signY * energyDepositionDistribution(s, para);
    }  // for
  }  // for
  return gg * para[9];
}

// Uses the signature needed for TMinuit interface:
// http://root.cern.ch/root/htmldoc/TMinuit.html#TMinuit:SetFCN
void StFmsClusterFitter::Fcn1(Int_t& npara, Double_t* grad, Double_t& fval,
                              Double_t* para, Int_t iflag) {
  // Number of expected photons should ALWAYS be the first parameter "para[0]"
  Int_t numbPh = (Int_t)para[0];
  StFmsTower* oneTow;
  // Sum energy of all towers being studied
  Double_t sumCl = 0;
  TIter next(StFmsClusterFitter::mTowers);
  while(oneTow=(StFmsTower*)next()) {
    sumCl += oneTow->hit()->energy();
  }  // while
  // Loop over all towers that are involved in the fit
  fval = 0;  // Stores sum of chi2 over each tower
  TIter nextTower(StFmsClusterFitter::mTowers);
  while(oneTow=(StFmsTower*) nextTower()) {
    // The shower shape function expects the centers of towers in units of cm
    // Tower centers are stored in row/column i.e. local coordinates
    // Therefore convert to cm, remembering to subtract 0.5 from row/column to
    // get centres not edges
    const Double_t x = (oneTow->column() - 0.5) *
                       StFmsClusterFitter::mTowerWidthXY[0];
    const Double_t y = (oneTow->row() - 0.5) *
                       StFmsClusterFitter::mTowerWidthXY[1];
    // Measured energy
    const Double_t eMeas = oneTow->hit()->energy();
    // Expected energy from Shower-Shape
    Double_t eSS = 0;
    for (Int_t iph = 0; iph < numbPh; iph++) {
      Int_t j = 3 * iph;
      Double_t Eshape = para[j + 3] *
                        showerShapeFitFunction.Eval(x - para[j + 1],
                                                    y - para[j + 2], 0);
      eSS += Eshape;
    }  // for
    const Double_t dev = eMeas - eSS;
    const Double_t errFactor = 0.03;
    const Double_t errQ = 0.01;
    // Larisa's Chi2 function
    const Double_t err = (errFactor * pow(eMeas / sumCl, 1. - 0.001 * sumCl) *
                         pow(1 - eMeas / sumCl, 1. - 0.007 * sumCl)) * sumCl
                         + errQ;
    const Double_t dchi2 = dev * dev / err;  // Calculate chi^2 of this tower
    fval += dchi2;  // Add chi2 for this tower to the sum
  }  // while
  // require that the fraction be positive!
  if (fval < 0) {
    fval = 0;
  }  // if
}

Double_t StFmsClusterFitter::fit(const Double_t* para, const Double_t* step,
                                 const Double_t* low, const Double_t* up,
                                 PhotonList* photons) {
  if (!step) {
    step = mSteps;
  }  // if
  Double_t chiSq(-1.);  // Return value
  // Check that there is a pointer to TObjArray of towers
  if(!StFmsClusterFitter::mTowers) {
    LOG_ERROR << "no tower data available! return -1!" << endm;
    return chiSq;
  }  // if
  mMinuit.SetFCN(Fcn1);  // Must set the function for Minuit to use
  Int_t nPh = (Int_t)para[0];  // Get the number of photons from parameters
  if (nPh < 1 || nPh > MAX_NUMB_PHOTONS) {
    LOG_ERROR << "nPh = " << nPh << "! Number of photons must be between 1 and "
      << MAX_NUMB_PHOTONS << "! Set it to be 1!" << endm;
    nPh = 1;
  }  // if
  mMinuit.mncler();  // Clear old parameters, so we can define the new parameters
  // The first parameter tells Minuit how many photons to fit!
  // It should be a fixed parameter, and between 1 and the max number of photons
  Int_t ierflg = 0;
  Double_t nPhotons(nPh);  // Minuit needs a double argument
  mMinuit.mnparm(0, "nph", nPhotons, 0, 0.5, 4.5, ierflg);
  // Set the rest of parameters: 3 parameters per photon
  for (Int_t i = 0; i < nPh; i++) {
    Int_t j = 3 * i + 1;  // Need to set 3 parameters per photon
    mMinuit.mnparm(j, Form("x%d", i+1), para[j], step[j], low[j], up[j], ierflg);
    j++;
    mMinuit.mnparm(j, Form("y%d", i+1), para[j], step[j], low[j], up[j], ierflg);
    j++;
    mMinuit.mnparm(j, Form("E%d", i+1), para[j], step[j], low[j], up[j], ierflg);
  }  // if
  Double_t arglist[10];
  arglist[0] = 1000;
  arglist[1] = 1.;
  ierflg = 0;
  mMinuit.mnexcm("MIGRAD", arglist ,2,ierflg);
  // Populate the list of photons
  if (0 == mMinuit.GetStatus() && photons) {
    // Get the fit results for starting positions and errors
    Double_t param[1 + 3 * MAX_NUMB_PHOTONS];
    Double_t error[1 + 3 * MAX_NUMB_PHOTONS];
    mMinuit.GetParameter(0, param[0], error[0]);
    // There are 3 parameters per photon, plus the 1st parameter
    nPh = (Int_t)param[0];  // Shouldn't have changed, but just to be safe...
    const Int_t nPar = 3 * nPh + 1;
    for (Int_t i(1); i < nPar; ++i) {  // Get remaining fit parameters (x, y, E)
      mMinuit.GetParameter(i, param[i], error[i]);
    }  // for
    for (Int_t par(1); par < nPar; par += 3) {  // Fill photons from parameters
      photons->push_back(  // x, y, E, error x, error y, error E
        StFmsFittedPhoton(param[par], param[par + 1], param[par + 2],
                          error[par], error[par + 1], error[par + 2]));
    }  // for
    // Evaluate chi-square (*not* chi-square per degree of freedom)
    Int_t iflag = 1;  // Don't calculate 1st derivatives, 2nd argument unneeded
    mMinuit.Eval(photons->size(), NULL, chiSq, param, iflag);
  }  // for
  return chiSq;
}

// a different set of parameters for 2-photon clusters only:
//    param[0]:  still a constant parameter, should be set to 2 for 2-photon fitting
//    param[1]:  xPi      (x-position of pi^0)
//    param[2]:  yPi      (y-position of pi^0)
//    param[3]:  d_gg     (distance between 2 photons)
//    param[4]:  theta    (theta angle of displacement vector from photon 2 to photon 1)
//    param[5]:  z_gg     (this z_gg can go from -1 to +1, so we do not set E1>E2)
//    param[6]:  E_gg     (total energy of two photons)
// Thus, in more conventional parameterization: x1, y1, E1, x2, y2, E2:
//     E1 = E_gg * (1 + z_gg)/2
//     E2 = E_gg * (1 - z_gg)/2
//     x1 = xPi + cos(theta) * d_gg * (1 - z_gg)/2
//     y1 = yPi + sin(theta) * d_gg * (1 - z_gg)/2
//     x2 = xPi - cos(theta) * d_gg * (1 + z_gg)/2
//     y2 = yPi - sin(theta) * d_gg * (1 + z_gg)/2
// The advantage of the new parameterization is that for 2-photon cluster fitting, we can
//    ensure that the two photons never get to close. The old parameterization certainly
//    suffers from this shortcoming if we let the parameters vary freely.
// What we already know about the limits of the new parameters:
//    xPi and yPi:   rarely do they go beyond 0.3 unit of lgd
//    theta:         have a narrow theta range (for r=sigmaMin/sigmaMax, |theta|<0.5*r/0.65
//                      when r<0.65, and linear increase from 0.5 to Pi/2 for 0.65<r<1)
//    E_gg:          given by Ec (+/- 20% or less)
//    z_gg:          should just let it vary from -1 to 1.
//    d_gg:          a lower bound is given by r=sqrt(sigmaX^2+sigmaY^2). 
//                      d_gg > Max( 2.5*(r-0.6), 0.5 )
Int_t StFmsClusterFitter::fit2PhotonCluster(const Double_t* para,
                                            const Double_t* step,
                                            const Double_t* low,
                                            const Double_t* up,
                                            PhotonList* photons) {
  if (!step) {
    step = mSteps;
  }  // if
  Double_t chiSq(-1.);  // Return value
  // Check that there is a pointer to TObjArray of towers
  if (!StFmsClusterFitter::mTowers) {
    LOG_ERROR << "no tower data available! return -1!" << endm;
    return chiSq;
  }  // if
  mMinuit.SetFCN(Fcn2);  // Must set the function for Minuit to use
  Int_t nPh = (Int_t)para[0];
  if (nPh != 2) {
    LOG_ERROR << "number of photons must be 2 for special 2-photon cluster "
      << "fitter \"Int_t StFmsClusterFitter::Fit2Pin1Clust(...)\"!"
      << " Set it to be 2!" << endm;
    nPh = 2;
  }  // if
  mMinuit.mncler();  // Clear old parameters, so we can define the new parameters
  // The first parameter tells Minuit how many photons to fit!
  // It should be a fixed parameter, in this case 2
  Int_t ierflg = 0;
  Double_t nPhotons(nPh);  // Minuit needs a double argument
  mMinuit.mnparm(0, "nph", nPhotons, 0, 1.5, 2.5, ierflg);
  mMinuit.mnparm(1, "xPi"  , para[1], step[1], low[1], up[1], ierflg);
  mMinuit.mnparm(2, "yPi"  , para[2], step[2], low[2], up[2], ierflg);
  mMinuit.mnparm(3, "d_gg" , para[3], step[3], low[3], up[3], ierflg);
  mMinuit.mnparm(4, "theta", para[4], step[4], low[4], up[4], ierflg);
  mMinuit.mnparm(5, "z_gg" , para[5], step[5], low[5], up[5], ierflg);
  mMinuit.mnparm(6, "E_gg" , para[6], step[6], low[6], up[6], ierflg);
  // Fix E_total and theta, we don't want these to be free parameters
  mMinuit.FixParameter(6);
  mMinuit.FixParameter(4);
  Double_t arglist[10];
  arglist[0] = 1000;
  arglist[1] = 1.;
  ierflg = 0;
  mMinuit.mnexcm("MIGRAD", arglist ,2,ierflg);
  mMinuit.mnfree(0);  // Free fixed parameters before next use of mMinuit
  if (0 == mMinuit.GetStatus() && photons) {
    // Get the fit results
    Double_t param[7];  // 3 * nPhotons + 1 parameters = 7 for 2 photons
    Double_t error[7];
    mMinuit.GetParameter(0, param[0], error[0]);
    Int_t nPar = 3 * (Int_t)param[0] + 1;  // Should be 7
    for (Int_t ipar = 1; ipar < nPar; ipar++) {  // Get the remaining parameters
      mMinuit.GetParameter(ipar, param[ipar], error[ipar]);
    }  // for
    // Put the fit result back in "clust". Need to translate the special
    // parameters for 2-photon fit into x, y, E, which looks a bit complicated!
    // First photons
    double x = param[1] + cos(param[4]) * param[3] * (1 - param[5]) / 2.0;
    double xErr = error[1] + (cos(param[4])*error[3]-error[4]*sin(param[4])*param[3])*(1-param[5])/2 - cos(param[4])*param[3]*error[5]/2.0;
    double y = param[2] + sin(param[4]) * param[3] * (1 - param[5]) / 2.0;
    double yErr = error[2] + (sin(param[4])*error[3]+error[4]*cos(param[4])*param[3])*(1-param[5])/2 - sin(param[4])*param[3]*error[5]/2.0;
    double E = param[6] * (1 + param[5]) / 2.0;
    double EErr = error[6]*(1+param[5])/2.0 + param[6]*error[5]/2.0;
    photons->push_back(StFmsFittedPhoton(x, y, E, xErr, yErr, EErr));
    // Second photon
    x = param[1] - cos(param[4]) * param[3] * (1 + param[5]) / 2.0;
    xErr = error[1] + (-cos(param[4])*error[3]+error[4]*sin(param[4])*param[3])*(1+param[5])/2 - cos(param[4])*param[3]*error[5]/2.0;
    y = param[2] - sin(param[4]) * param[3] * (1 + param[5]) / 2.0;
    yErr = error[2] + (sin(param[4])*error[3]-error[4]*cos(param[4])*param[3])*(1+param[5])/2 - sin(param[4])*param[3]*error[5]/2.0;
    E = param[6] * (1 - param[5]) / 2.0;
    EErr = error[6]*(1-param[5])/2.0 - param[6]*error[5]/2.0;
    photons->push_back(StFmsFittedPhoton(x, y, E, xErr, yErr, EErr));
    // Evaluate the Chi-square function
    Int_t iflag = 1;  // Don't calculate 1st derivatives...
    mMinuit.Eval(7, NULL, chiSq, param, iflag);  // ... so 2nd argument unneeded
  }  // if
  return chiSq;
}

void StFmsClusterFitter::Fcn2(Int_t& nparam, Double_t* grad, Double_t& fval,
                              Double_t* param, Int_t iflag) {
  // Only need to translate into the old parameterization
  Double_t oldParam[7];
  float dd = param[3];
  oldParam[0] = param[0];  // Number of photons, unchanged
  oldParam[1] = param[1] + cos(param[4]) * dd * (1 - param[5]) / 2.0;  // x 1
  oldParam[2] = param[2] + sin(param[4]) * dd * (1 - param[5]) / 2.0;  // y 1
  oldParam[3] = param[6] * (1 + param[5]) / 2.0;  // Energy 1
  oldParam[4] = param[1] - cos(param[4]) * dd * (1 + param[5]) / 2.0;  // x 2
  oldParam[5] = param[2] - sin(param[4]) * dd * (1 + param[5]) / 2.0;  // y 2
  oldParam[6] = param[6] * (1 - param[5]) / 2.0;  // Energy 2
  // Now we can call the regular Fcn1 with the translated parameters
  Fcn1(nparam, grad, fval, oldParam, iflag);
}
}  // namespace PSUGlobals

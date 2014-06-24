// $Id$
//
// $Log$
/**
 \file      StFmsClusterFitter.cxx
 \brief     Implementation of StFmsClusterFitter, shower-shape fitting routine
 \author    Steven Heppelmann <steveheppelmann@gmail.com>
 \author    Yuxi Pan <yuxipan@physics.ucla.edu>
 \author    Thomas Burton <tpb@bnl.gov>
 \date      2014
 \copyright Brookhaven National Lab
 */
#include "StFmsPointMaker/StFmsClusterFitter.h"

#include <cmath>
#include <functional>
#include <numeric>
#include <vector>

#include "TF2.h"
#include "TMath.h"
#include "TString.h"

#include "StRoot/St_base/StMessMgr.h"
#include "StEvent/StFmsHit.h"

#include "StFmsPointMaker/StFmsGeometry.h"
#include "StFmsPointMaker/StFmsTower.h"

namespace {
const Int_t kMaxNPhotons = 7;  // Maximum number of photons that can be fitted
const Int_t kNFitParameters = 10;  // Parameters for shower-shape function
TF2 showerShapeFitFunction("showerShapeFitFunction",
                       &FMSCluster::StFmsClusterFitter::energyDepositionInTower,
                      -25.0, 25.0, -25.0, 25.0, kNFitParameters);
/*
 Compose Minuit step size in each fit variable.
 The first value is for the number of photons.
 Each subsequent triplet is for the (x, y, E) of a photon, up to kMaxNPhotons.
 */
std::vector<double> defaultMinuitStepSizes() {
  std::vector<double> steps(1, 0.);  // Initialise with nPhoton step
  // Append default (x, y, E) steps for each photon
  for (int i(0); i < kMaxNPhotons; ++i) {
    steps.insert(steps.end(), {0.1, 0.1, 0.2});
  }  // for
  return steps;
}
std::vector<float> towerWidths;  // Tower (x, y) width in cm

// Helper function for accumulating tower energies
double addTowerEnergy(double energy, const FMSCluster::StFmsTower* tower) {
  return energy + tower->hit()->energy();
}

// Returns a * f(x,y,b) as defined here:
// https://drupal.star.bnl.gov/STAR/blog/leun/2010/aug/02/fms-meeting-20100802
double showerShapeComponent(double x, double y, double a, double b) {
  return a * atan(x * y / (b * sqrt(b * b + x * x + y * y)));
}
}  // unnamed namespace

namespace FMSCluster {
// Instantiate static members
StFmsTowerCluster::Towers* StFmsClusterFitter::mTowers(nullptr);

StFmsClusterFitter::StFmsClusterFitter(const StFmsGeometry* geometry,
                                       Int_t detectorId)
    : mMinuit(3 * kMaxNPhotons + 1) {
  // Set tower (x, y) widths for this detector
  towerWidths = geometry->towerWidths(detectorId);
  double parameters[kNFitParameters] = {towerWidths.at(0), 1.070804, 0.167773,
                                        -0.238578, 0.535845, 0.850233, 2.382637,
                                        0.0, 0.0, 1.0};
  showerShapeFitFunction.SetParameters(parameters);
  mMinuit.SetPrintLevel(-1);  // Quiet, including suppression of warnings
}

StFmsClusterFitter::~StFmsClusterFitter() { }

TF2* StFmsClusterFitter::showerShapeFunction() {
  return &showerShapeFitFunction;
}

Double_t StFmsClusterFitter::fit(const std::vector<double>& para,
                                 const std::vector<double>& step,
                                 const std::vector<double>& low,
                                 const std::vector<double>& up,
                                 PhotonList* photons) {
  Double_t chiSq(-1.);  // Return value
  // Check that there is a pointer to TObjArray of towers
  if (!StFmsClusterFitter::mTowers) {
    LOG_ERROR << "no tower data available! return -1!" << endm;
    return chiSq;
  }  // if
  mMinuit.SetFCN(minimizationFunctionNPhoton);
  int nPh = para.size() / 3;
  if (nPh < 1 || nPh > kMaxNPhotons) {
    LOG_ERROR << "nPh = " << nPh << "! Number of photons must be between 1 and "
      << kMaxNPhotons << "! Set it to be 1!" << endm;
    nPh = 1;
  }  // if
  mMinuit.mncler();  // Clear old parameters so we can define the new parameters
  // The first parameter tells Minuit how many photons to fit!
  // It should be a fixed parameter, and between 1 and the max number of photons
  setMinuitParameter(0, "nph", para, step, low, up);
  // Set the rest of parameters: 3 parameters per photon
  for (Int_t i = 0; i < nPh; i++) {
    Int_t j = 3 * i + 1;  // Need to set 3 parameters per photon
    setMinuitParameter(j++, Form("x%d", i + 1), para, step, low, up);
    setMinuitParameter(j++, Form("y%d", i + 1), para, step, low, up);
    setMinuitParameter(j++, Form("E%d", i + 1), para, step, low, up);
  }  // if
  std::vector<double> arglist = {1000., 1.};
  // Perform the minimization
  int ierflg = 0;
  mMinuit.mnexcm("MIGRAD", arglist.data(), arglist.size(), ierflg);
  // Populate the list of photons from the fit results
  if (0 == mMinuit.GetStatus() && photons) {
    // Get the fit results and errors
    std::vector<double> param(para.size(), 0.), error(para.size(), 0.);
    for (unsigned i(0); i < para.size(); ++i) {
      mMinuit.GetParameter(i, param.at(i), error.at(i));
    }  // for
    for (unsigned i(1); i < para.size(); i += 3) {
      photons->push_back(  // x, y, E, error x, error y, error E
        StFmsFittedPhoton(param.at(i), param.at(i + 1), param.at(i + 2),
                          error.at(i), error.at(i + 1), error.at(i + 2)));
    }  // for
    // Evaluate chi-square (*not* chi-square per degree of freedom)
    Int_t iflag = 1;  // Don't calculate 1st derivatives, 2nd argument unneeded
    mMinuit.Eval(photons->size(), nullptr, chiSq, param.data(), iflag);
  }  // for
  return chiSq;
}

Double_t StFmsClusterFitter::fit(const std::vector<double>& parameters,
                                 const std::vector<double>& lower,
                                 const std::vector<double>& upper,
                                 PhotonList* photons) {
  return fit(parameters, defaultMinuitStepSizes(), lower, upper, photons);
}

/*
 A different set of parameters for 2-photon clusters only:
  0: still a constant parameter, should be set to 2 for 2-photon fitting
  1: xPi, x-position of pi^0
  2: yPi, y-position of pi^0
  3: d_gg, distance between 2 photons
  4: theta, angle of displacement vector from photon 2 to photon 1
  5: z_gg, can go from -1 to +1, so we do not set E1 > E2
  6: E_gg, total energy of two photons
 Thus, in the more conventional fit() parameterization: x1, y1, E1, x2, y2, E2:
  E1 = E_gg * (1 + z_gg) / 2
  E2 = E_gg * (1 - z_gg) / 2
  x1 = xPi + cos(theta) * d_gg * (1 - z_gg) / 2
  y1 = yPi + sin(theta) * d_gg * (1 - z_gg) / 2
  x2 = xPi - cos(theta) * d_gg * (1 + z_gg) / 2
  y2 = yPi - sin(theta) * d_gg * (1 + z_gg) / 2

 The advantage of this parameterization is that for 2-photon cluster fitting
 we can ensure that the two photons never get to close. The old parameterization
 suffers from this shortcoming if we let the parameters vary freely.

 What we already know about the limits of the new parameters:
  xPi and yPi: rarely do they go beyond 0.3 unit of lgd
  theta:       have a narrow theta range (for r = sigmaMax / sigmaMax,
               |theta| < 0.5 * r / 0.65 when r < 0.65, and linear increase
               from 0.5 to pi/2 for 0.65 < r < 1)
  E_gg:        given by Ec (+/- 20% or less)
  z_gg:        should just let it vary from -1 to 1
  d_gg:        a lower bound is given by r = sqrt(sigmaX^2 + sigmaY^2). 
               d_gg > Max(2.5 * (r - 0.6), 0.5)
 */
Int_t StFmsClusterFitter::fit2PhotonCluster(const std::vector<double>& para,
                                            std::vector<double> step,
                                            const std::vector<double>& low,
                                            const std::vector<double>& up,
                                            PhotonList* photons) {
  if (step.empty()) {
    step = defaultMinuitStepSizes();
  }  // if
  Double_t chiSq(-1.);  // Return value
  if (!StFmsClusterFitter::mTowers) {
    LOG_ERROR << "no tower data available! return -1!" << endm;
    return chiSq;
  }  // if
  mMinuit.SetFCN(minimizationFunction2Photon);
  int nPh = para.size() / 3;
  if (nPh != 2) {
    LOG_ERROR << "number of photons must be 2 for special 2-photon cluster "
      << "fitter \"Int_t StFmsClusterFitter::fit2PhotonCluster(...)\"!"
      << " Set it to be 2!" << endm;
    nPh = 2;
  }  // if
  mMinuit.mncler();  // Clear old parameters so we can define the new parameters
  // The first parameter tells Minuit how many photons to fit, in this case 2
  const std::vector<TString> names = {
    "nph", "xPi", "yPi", "d_gg", "theta", "z_gg", "E_gg"
  };
  for (unsigned i = 0; i < names.size(); ++i) {
    setMinuitParameter(i, names.at(i), para, step, low, up);
  }  // for
  // Fix E_total and theta, we don't want these to be free parameters
  mMinuit.FixParameter(6);
  mMinuit.FixParameter(4);
  std::vector<double> arglist = {1000., 1.};
  // Perform the minimization
  int ierflg = 0;
  mMinuit.mnexcm("MIGRAD", arglist.data(), arglist.size(), ierflg);
  mMinuit.mnfree(0);  // Free fixed parameters before next use of mMinuit
  if (0 == mMinuit.GetStatus() && photons) {
    // Get the fit results for starting positions and errors
    // 3 * nPhotons + 1 parameters = 7 for 2 photons
    std::vector<double> param(7, 0.), error(7, 0.);
    for (unsigned i = 0; i < param.size(); ++i) {
      mMinuit.GetParameter(i, param.at(i), error.at(i));
    }  // for
    // Put the fit result back in "clust". Need to translate the special
    // parameters for 2-photon fit into x, y, E, which looks a bit complicated!
    // First photon
    double x = param[1] + cos(param[4]) * param[3] * (1 - param[5]) / 2.0;
    double xErr = error[1] +
      (cos(param[4]) * error[3] - error[4] * sin(param[4]) * param[3]) *
      (1 - param[5]) / 2 - cos(param[4]) * param[3] * error[5] / 2.0;
    double y = param[2] + sin(param[4]) * param[3] * (1 - param[5]) / 2.0;
    double yErr = error[2] +
      (sin(param[4]) * error[3] + error[4] * cos(param[4]) * param[3]) *
      (1 - param[5]) / 2 - sin(param[4]) * param[3] * error[5] / 2.0;
    double E = param[6] * (1 + param[5]) / 2.0;
    double EErr = error[6] * (1 + param[5]) / 2.0 + param[6] * error[5] / 2.0;
    photons->push_back(StFmsFittedPhoton(x, y, E, xErr, yErr, EErr));
    // Second photon
    x = param[1] - cos(param[4]) * param[3] * (1 + param[5]) / 2.0;
    xErr = error[1] +
           (-cos(param[4]) * error[3] + error[4] * sin(param[4]) * param[3]) *
           (1 + param[5]) / 2 - cos(param[4]) * param[3] * error[5] / 2.0;
    y = param[2] - sin(param[4]) * param[3] * (1 + param[5]) / 2.0;
    yErr = error[2] +
           (sin(param[4]) * error[3] - error[4] * cos(param[4]) * param[3]) *
           (1 + param[5]) / 2 - sin(param[4]) * param[3] * error[5] / 2.0;
    E = param[6] * (1 - param[5]) / 2.0;
    EErr = error[6] * (1 - param[5]) / 2.0 - param[6] * error[5] / 2.0;
    photons->push_back(StFmsFittedPhoton(x, y, E, xErr, yErr, EErr));
    // Evaluate the chi-square function
    Int_t iflag = 1;  // Don't calculate 1st derivatives, 2nd argument unneeded
    mMinuit.Eval(7, nullptr, chiSq, param.data(), iflag);
  }  // if
  return chiSq;
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
      Double_t s[2] = {xy[0] - para[7] + signX * para[0] / 2.0,   // x +/- d/2
                       xy[1] - para[8] + signY * para[0] / 2.0};  // y +/- d/2
      gg += signX * signY * energyDepositionDistribution(s, para);
    }  // for
  }  // for
  return gg * para[9];
}

int StFmsClusterFitter::maxNFittedPhotons() {
  return kMaxNPhotons;
}

// Calculate fractional photon energy deposition in a tower based on its (x, y)
// position relative to the tower center
Double_t StFmsClusterFitter::energyDepositionDistribution(
    Double_t* xy,
    Double_t* parameters) {
  double f = 0;
  // The parameter array has 10 elements, but we only use 6 here
  // 1 to 6 are a1, a2, a3, b1, b2, b3 as defined in
  // https://drupal.star.bnl.gov/STAR/blog/leun/2010/aug/02/fms-meeting-20100802
  for (int i = 1; i < 4; i++) {  // 1, 2, 3
    f += showerShapeComponent(xy[0], xy[1], parameters[i], parameters[i + 3]);
  }  // for
  return f / TMath::TwoPi();
}

// Uses the signature needed for TMinuit interface:
// http://root.cern.ch/root/htmldoc/TMinuit.html#TMinuit:SetFCN
void StFmsClusterFitter::minimizationFunctionNPhoton(Int_t& npara,
                                                     Double_t* grad,
                                                     Double_t& fval,
                                                     Double_t* para,
                                                     Int_t iflag) {
  const double energySum = std::accumulate(mTowers->begin(), mTowers->end(),
                                           0., addTowerEnergy);
  fval = 0;  // Stores sum of chi2 over each tower
  const int nPhotons = static_cast<int>(para[0]);
  for (auto i = mTowers->begin(); i != mTowers->end(); ++i) {
    const StFmsTower* tower = *i;
    // The shower shape function expects the centers of towers in units of cm
    // Tower centers are stored in row/column i.e. local coordinates
    // Therefore convert to cm, remembering to subtract 0.5 from row/column to
    // get centres not edges
    const double x = towerWidths.at(0) * (tower->column() - 0.5);
    const double y = towerWidths.at(1) * (tower->row() - 0.5);
    // Add expected energy in tower from each photon, according to shower-shape
    double expected = 0;
    for (int j = 0; j < nPhotons; ++j) {  // Recall there are 3 paras per photon
      int k = 3 * j;
      expected += para[k + 3] *  // total energy
                  showerShapeFitFunction.Eval(x - para[k + 1], y - para[k + 2]);
    }  // for
    const double measured = tower->hit()->energy();
    const double deviation = measured - expected;
    // Larisa's chi2 function definition
    const Double_t err = 0.03 *
                         pow(measured / energySum, 1. - 0.001 * energySum) *
                         pow(1 - measured / energySum, 1. - 0.007 * energySum) *
                         energySum + 0.01;
    fval += pow(deviation, 2.) / err;
  }  // for
  fval = std::max(fval, 0.);  // require that the fraction be positive
}

// Uses the signature needed for TMinuit interface:
// http://root.cern.ch/root/htmldoc/TMinuit.html#TMinuit:SetFCN
void StFmsClusterFitter::minimizationFunction2Photon(Int_t& nparam,
                                                     Double_t* grad,
                                                     Double_t& fval,
                                                     Double_t* param,
                                                     Int_t iflag) {
  // Only need to translate into the old parameterization
  const float dd = param[3];
  Double_t oldParam[7] = {
    param[0],  // Number of photons, unchanged
    param[1] + cos(param[4]) * dd * (1 - param[5]) / 2.0,  // x 1
    param[2] + sin(param[4]) * dd * (1 - param[5]) / 2.0,  // y 1
    param[6] * (1 + param[5]) / 2.0,  // Energy 1
    param[1] - cos(param[4]) * dd * (1 + param[5]) / 2.0,  // x 2
    param[2] - sin(param[4]) * dd * (1 + param[5]) / 2.0,  // y 2
    param[6] * (1 - param[5]) / 2.0  // Energy 2
  };
  // Now call the regular minimization function with the translated parameters
  minimizationFunctionNPhoton(nparam, grad, fval, oldParam, iflag);
}

int StFmsClusterFitter::setMinuitParameter(int index, const TString& name,
                                           const std::vector<double>& par,
                                           const std::vector<double>& step,
                                           const std::vector<double>& low,
                                           const std::vector<double>& up) {
  int error = 0;
  mMinuit.mnparm(index, name, par.at(index), step.at(index),
                 low.at(index), up.at(index), error);
  return error;
}
}  // namespace FMSCluster

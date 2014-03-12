#include "FitTower.h"

#include <TF2.h>
#include <TMath.h>

#include "StEvent/StFmsHit.h"

using namespace PSUGlobals;

ClassImp(FitTower)

namespace {
Int_t numbPara = 10;
TF2 showerShapeFitFunction("showerShapeFitFunction", &FitTower::GGams,
                           -25.0, 25.0, -25.0, 25.0, numbPara);
}  // unnamed namespace

// Instantiate static members
Float_t FitTower::widLG[2];
TObjArray* FitTower::tow2Fit(NULL);

TF2* FitTower::GetFunctShowShape() {
  return &showerShapeFitFunction;
}

FitTower::FitTower(Geom* pgeom,Int_t iew,Int_t nstb) {
  SetStep();
  pTowerUtil=new TowerUtil();
  Float_t* towerWidth = pgeom->FpdTowWid(iew, nstb);
  fTWidthCM = towerWidth[0];
  fTXWidthCM = FitTower::widLG[0] = towerWidth[0];
  fTYWidthCM = FitTower::widLG[1] = towerWidth[1];
  Double_t para[numbPara];
  para[0] = fTWidthCM ;
  para[1] =  1.070804;
  para[2] =  0.167773;
  para[3] =  -0.238578;
  para[4] =  0.535845;
  para[5] =  0.850233;
  para[6] =  2.382637;
  para[7] =  0.0 ;
  para[8] =  0.0 ;
  para[9] =  1.0 ;
  showerShapeFitFunction.SetParameters(para); 
  // create a Minuit instance
  fMn = new TMinuit(3*MAX_NUMB_PHOTONS+1);
  fMn->SetPrintLevel(-1);  // Quiet, including suppression of warnings
}

FitTower::~FitTower() {
  if (fMn) {
    delete fMn;
  }  // if
  delete pTowerUtil;
}

void FitTower::SetStep() {
  fMn=0;
  const Double_t step0[3 * MAX_NUMB_PHOTONS + 1]= {
    0.0, 0.1, 0.1, 0.2, 0.1, 0.1, 0.2, 0.1, 0.1, 0.2, 0.1,
    0.1, 0.2, 0.1, 0.1, 0.2, 0.1, 0.1, 0.2, 0.1, 0.1, 0.2
  };
  for(int j=0;j<3*MAX_NUMB_PHOTONS+1;j++)step[j]= step0[j];
}

Double_t FitTower::FGams(Double_t *x, Double_t *para) {
  Double_t f=0;
  Double_t xx=x[0];
  Double_t yy=x[1];
  for(Int_t i=1;i<=3;i++) {
    Int_t j;
    j = i + 3 ;
    f += para[i]*( atan( xx*yy/ (para[j]* sqrt(para[j]*para[j]+xx*xx+yy*yy) ) ) );
  };
  return f/(2 * TMath::Pi() ) ;
}

Double_t FitTower::GGams(Double_t *x, Double_t *para) {
  Double_t gg, s[2];
  gg = 0 ;  
  for (Int_t ix=0; ix<2; ix++) {
    for (Int_t iy=0; iy<2; iy++) {
      Double_t ax, ay;
      ax = pow(-1.0, ix);
      ay = pow(-1.0, iy);
      s[0] = x[0] - para[7] + ax * para[0] / 2.0 ;
      s[1] = x[1] - para[8] + ay * para[0] / 2.0 ;
      gg += ax * ay * FGams(s, para) ;
    }  // for
  }  // for
  return gg*para[9];
}

void FitTower::Fcn1(Int_t& npara, Double_t* grad, Double_t& fval,
                    Double_t* para, Int_t iflag) {
  // number of expected photons
  // should ALWAYS be the first parameter "para[0]"
  Int_t numbPh;
  numbPh = (Int_t) para[0] ;
  fval = 0 ;
  Double_t err ;
  Double_t xx, yy;
  Double_t eSS, eMeas;
  Double_t dev;
  // we are in GeV, not ADC count
  TowerFPD * oneTow;
  // first get cluster sum
  Double_t sumCl = 0;
  TIter next(FitTower::tow2Fit);
  while(oneTow=(TowerFPD*) next()) {
    sumCl+=oneTow->hit()->energy();
  }  // while
  // loop over all towers that are involved in the fit
  TIter nextTower(FitTower::tow2Fit);
  while(oneTow=(TowerFPD*) nextTower()) {
    // center of tower in unit of "cm"
    // my towers are center at 0.5 to 6.5, as Steve Heppelmann
    //Note from SFH need more clever position for FMS
    xx = ( (Double_t) oneTow->column() - 0.5 ) * FitTower::widLG[0] ;
    yy = ( (Double_t) oneTow->row() - 0.5 ) * FitTower::widLG[1] ;
    // measured energy
    eMeas = oneTow->hit()->energy();
    // expected energy from Shower-Shape
    eSS = 0 ;
    for(Int_t iph=0; iph<numbPh; iph++) {
      Int_t j;
      j = 3 * iph ;
      // shower-shape function calculate the fraction of energy
      // in coords of center of tower relative to photon
      Double_t Eshape = para[j+3] * showerShapeFitFunction.Eval(xx-para[j+1], yy-para[j+2], 0);
      eSS+=Eshape;
    }  // for
    dev = eMeas - eSS ;
    Double_t dchi2;
    const double errFactor = 0.03;
    const double errQ = 0.01;
    // Larisa'e Chi2 function
    err = (errFactor * pow(eMeas/sumCl,1.-.001*sumCl) *
           pow(1 - eMeas/sumCl,1.-.007*sumCl))*sumCl
           +errQ;
    dchi2 = dev * dev / err;
    float dsign=1.;
    if(dev<0)dsign=-1.;
    fval += dchi2;
  }  // while
  // require that the fraction be positive!
  if (fval < 0) {
    fval = 0;
  }  // if
}

Double_t FitTower::Fit(const Double_t *para, const Double_t *step,
                       const Double_t *low, const Double_t *up,
                       PhotonList* photons) {
  Double_t chiSq(-1.);
  // check that there is a pointer to TObjArray of towers
  if( !(FitTower::tow2Fit) ) {
    std::cerr << "no tower data available! return -1!" << "\n";
    return chiSq;
  }  // if
  // must set the function to "Fcn1"!
  fMn->SetFCN(Fcn1);
  Int_t nPh = (Int_t) para[0];
  if( nPh < 1 || nPh > MAX_NUMB_PHOTONS ) {
    fNumbPhotons = 1;
    std::cerr << "nPh = " << nPh << "! Number of photons must be between 1 and " << MAX_NUMB_PHOTONS << "! Set it to be 1!" << "\n";
  } else {
    fNumbPhotons = nPh ;
  }  // if
  // clear old parameters, so we can define the new parameters
  fMn->mncler();
  // The first parameter tells Minuit how many photons to fit!
  // It should be a fixed parameter!
  Int_t ierflg = 0;
  fMn->mnparm(0, "nph", fNumbPhotons, 0, 0.5, 4.5, ierflg);
  // set the rest of parameters, 3 parameters for 1 photon
  for (Int_t i=0; i<fNumbPhotons; i++) {
    Int_t j ;
    j = 3*i+1 ;
    fMn->mnparm(j, Form("x%d", i+1), para[j], step[j], low[j], up[j], ierflg);
    j++ ;
    fMn->mnparm(j, Form("y%d", i+1), para[j], step[j], low[j], up[j], ierflg);
    j++ ;
    fMn->mnparm(j, Form("E%d", i+1), para[j], step[j], low[j], up[j], ierflg);
  }  // if
  Double_t arglist[10];
  arglist[0] = 1000;
  arglist[1] = 1.;
  ierflg = 0;
  fMn->mnexcm("MIGRAD", arglist ,2,ierflg);
  // Populate the list of photons
  if (0 == fMn->GetStatus() && photons) {
    // Get the fit results for starting positions and errors
    Double_t param[1 + 3 * MAX_NUMB_PHOTONS];
    Double_t error[1 + 3 * MAX_NUMB_PHOTONS];
    fMn->GetParameter(0, param[0], error[0]);
    // There are 3 parameters per photon, plus the 1st parameter,
    // which gives the number of photons
    nPh = (Int_t)param[0];
    Int_t nPar = 3 * nPh + 1;
    for (Int_t i(1); i < nPar; ++i) {
      fMn->GetParameter(i, param[i], error[i]);
    }  // for
    // There are 3 parameters per photon, starting at parameter 1
    for (Int_t par(1); par < nPar; par += 3) {
      photons->push_back(
        PhotonHitFPD(param[par], param[par + 1], param[par + 2],
                     error[par], error[par + 1], error[par + 2]));
    }  // for
    // Evaluate chi-square (*not* chi-square per degree of freedom)
    Double_t gradient[1 + 3 * MAX_NUMB_PHOTONS];
    Int_t iflag = 1;
    fMn->Eval(photons->size(), gradient, chiSq, param, iflag);
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
Int_t FitTower::Fit2Pin1Clust(const Double_t *para, const Double_t *step,
                              const Double_t *low, const Double_t *up,
                              PhotonList* photons) {
  Double_t chiSq(-1.);
  // check that there is a pointer to TObjArray of towers
  if (!(FitTower::tow2Fit)) {
    std::cerr << "no tower data available! return -1!" << "\n";
    return chiSq;
  }  // if
  // must set the function to "Fcn2"!
  fMn->SetFCN(Fcn2);
  Int_t nPh = (Int_t) para[0];
  fNumbPhotons = 2;
  if (nPh != 2) {
    std::cerr << "number of photons must be 2 for special 2-photon cluster fitter \"Int_t FitTower::Fit2Pin1Clust(...)\"!";
    std::cerr << " Set it to be 2!" << "\n";
    fNumbPhotons = 2;
  }  // if
  // clear old parameters, so we can define the new parameters
  fMn->mncler();
  // The first parameter tells Minuit how many photons to fit!
  // It should be a fixed parameter!
  Int_t ierflg = 0;
  fMn->mnparm(0, "nph", fNumbPhotons, 0, 1.5, 2.5, ierflg);
  fMn->mnparm(1, "xPi"  , para[1], step[1], low[1], up[1], ierflg);
  fMn->mnparm(2, "yPi"  , para[2], step[2], low[2], up[2], ierflg);
  fMn->mnparm(3, "d_gg" , para[3], step[3], low[3], up[3], ierflg);
  fMn->mnparm(4, "theta", para[4], step[4], low[4], up[4], ierflg);
  fMn->mnparm(5, "z_gg" , para[5], step[5], low[5], up[5], ierflg);
  fMn->mnparm(6, "E_gg" , para[6], step[6], low[6], up[6], ierflg);
  Double_t arglist[10];
  arglist[0] = 1000;
  arglist[1] = 1.;
  ierflg = 0;
  // fix E_total and theta
  fMn->FixParameter(6);
  fMn->FixParameter(4);
  fMn->mnexcm("MIGRAD", arglist ,2,ierflg);
  fMn->mnfree(0);
  if (0 == fMn->GetStatus() && photons) {
    // Get the fit result
    Double_t param[7];  // 3 * nPhotons + 1 parameters
    Double_t error[7];
    fMn->GetParameter(0, param[0], error[0]);
    Int_t nPar = 3*((Int_t) param[0])+1;
    for (Int_t ipar = 1; ipar < nPar; ipar++) {
      fMn->GetParameter(ipar, param[ipar], error[ipar]);
    }  // for
    // Put the fit result back in "clust"
    double x = param[1] + cos(param[4]) * param[3] * (1 - param[5]) / 2.0 ;
    double xErr = error[1] + (cos(param[4])*error[3]-error[4]*sin(param[4])*param[3])*(1-param[5])/2 - cos(param[4])*param[3]*error[5]/2.0;
    double y = param[2] + sin(param[4]) * param[3] * (1 - param[5]) / 2.0 ;
    double yErr = error[2] + (sin(param[4])*error[3]+error[4]*cos(param[4])*param[3])*(1-param[5])/2 - sin(param[4])*param[3]*error[5]/2.0;
    double E = param[6] * (1 + param[5]) / 2.0 ;
    double EErr = error[6]*(1+param[5])/2.0 + param[6]*error[5]/2.0;
    photons->push_back(PhotonHitFPD(x, y, E, xErr, yErr, EErr));
    // Second photon
    x = param[1] - cos(param[4]) * param[3] * (1 + param[5]) / 2.0 ;
    xErr = error[1] + (-cos(param[4])*error[3]+error[4]*sin(param[4])*param[3])*(1+param[5])/2 - cos(param[4])*param[3]*error[5]/2.0;
    y = param[2] - sin(param[4]) * param[3] * (1 + param[5]) / 2.0 ;
    yErr = error[2] + (sin(param[4])*error[3]-error[4]*cos(param[4])*param[3])*(1+param[5])/2 - sin(param[4])*param[3]*error[5]/2.0;
    E = param[6] * (1 - param[5]) / 2.0 ;
    EErr = error[6]*(1-param[5])/2.0 - param[6]*error[5]/2.0;
    photons->push_back(PhotonHitFPD(x, y, E, xErr, yErr, EErr));
    // Evaluate the Chi-square function
    Double_t gradient[7];
    Int_t iflag = 1;
    fMn->Eval(7, gradient, chiSq, param, iflag);
  }  // if
  return chiSq;
}

void FitTower::Fcn2(Int_t & nparam, Double_t *grad, Double_t &fval,
                    Double_t *param, Int_t iflag) {
  // only need to translate into the old parameterization
  Double_t oldParam[7];
  float dd=param[3];
  oldParam[0] = param[0] ;
  oldParam[1] = param[1] + cos(param[4]) * dd * (1 - param[5]) / 2.0 ;
  oldParam[2] = param[2] + sin(param[4]) * dd * (1 - param[5]) / 2.0 ;
  oldParam[3] = param[6] * (1 + param[5]) / 2.0 ;
  oldParam[4] = param[1] - cos(param[4]) * dd * (1 + param[5]) / 2.0 ;
  oldParam[5] = param[2] - sin(param[4]) * dd * (1 + param[5]) / 2.0 ;
  oldParam[6] = param[6] * (1 - param[5]) / 2.0 ;
  // then just call "Fcn1(...)"
  Fcn1(nparam, grad, fval, oldParam, iflag);
}

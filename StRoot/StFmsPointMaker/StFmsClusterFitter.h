// $Id$
//
// $Log$
/**
 \file      StFmsClusterFitter.h
 \brief     Declaration of StFmsClusterFitter, shower-shape fitting routine
 \author    Steven Heppelmann <steveheppelmann@gmail.com>
 \author    Yuxi Pan <yuxipan@physics.ucla.edu>
 \author    Thomas Burton <tpb@bnl.gov>
 \date      2014
 \copyright Brookhaven National Lab
 */
#ifndef STROOT_STFMSPOINTMAKER_STFMSCLUSTERFITTER_H_
#define STROOT_STFMSPOINTMAKER_STFMSCLUSTERFITTER_H_

#include <list>

#include "TMinuit.h"
#include "TObject.h"

#include "StFmsPointMaker/StFmsFittedPhoton.h"
#include "StFmsPointMaker/StFmsTowerCluster.h"

class TF2;

namespace FMSCluster {  // $NMSPC
typedef std::list<StFmsFittedPhoton> PhotonList;
class StFmsGeometry;
class StFmsTower;
/**
 Photon shower-shape fitting routine for FMS clusters.

 Fits tower clusters with the photon energy-deposition shower-shape
 distribution function i.e. how a photon hitting the FMS distributes its energy
 over towers.
 This serves to:
  - find the best-fit photon properties for a cluster.
  - distinguish clusters made by a single photon from those made by two
    nearby photons merging into a single cluster, which can frequently happen
    at the far-forward rapidity of the FMS.

 \todo
 It may be safer to make StFmsClusterFitter a singleton class, or something like that.
 The shower shape fit function is shared across all objects (by necessity, in
 order to interface with TMinuit), but each object updates the function itself
 with different parameters. Therefore bad things would happen if there were
 more than one object in existence at any time. There isn't ever more than one
 instance created in this code, but I think it would be good to enforce that.
 */
class StFmsClusterFitter : public TObject {
 public:
  /** Constructor using detector geometry for a single sub-detector */
  StFmsClusterFitter(const StFmsGeometry* geometry, Int_t detectorId);
  /**
   Default constructor.

   \todo Actually initialize things!
   */
  StFmsClusterFitter() { }
  /** Destructor */
  ~StFmsClusterFitter();
  /**
   Return the shower shape function.
   
   The shower shape gives the fractional energy deposition by a photon in a
   tower as a function of the distance of the photon from the tower center.
   */
  TF2* showerShapeFunction();
  /** Set the tower list to fit when calling fit() or fit2PhotonCluster() */
  void setTowers(StFmsTowerCluster::Towers* towers) { mTowers = towers; }
  /**
   Fit photons to the list of towers.

   All array arguments are of size 3N+1 for an N-photon fit:
    - par: initial guess values for each fit variable.
    - step: step size when fitting (if NULL use default values).
    - low: lower bound on each fit variable.
    - up: upper bound on each fit variable.

   In each 3N+1D array, the first element is the number of photons to fit.
   Each subsequent triplet is the x position, y position and energy of a photon.
   e.g.
    - for a 1-photon fit: [1, x0, y0, E0]
    - for a 2-photon fit: [2, x0, y0, E0, x1, y1, E1]

   Returns the &chi;<sup>2</sup> of the fit.
   */
  Double_t fit(const Double_t* par, const Double_t* step,
               const Double_t* low, const Double_t* up, PhotonList* photons);
  /**
   Specialized fit function for exactly 2-photon fit.
   
   Argument arrays are as for fit().
   However as this is for 2-photon fits only, the input arrays should always
   have 7 (3 * 2 photons + 1) elements. Additionally each element has a
   different meaning here:
    - 0: Still a constant parameter, should be set to 2 for 2-photon fitting.
    - 1: x-position of pi<sup>0</sup>.
    - 2: y-position of pi<sup>0</sup>.
    - 3: Distance between 2 photons.
    - 4: Theta angle of displacement vector from photon 2 to photon 1.
    - 5: z_gg, from -1 to +1, so we do not set E1 > E2.
    - 6: Total energy of two photons.

   Returns the &chi;<sup>2</sup> of the fit.
   */
  Int_t fit2PhotonCluster(const Double_t* para, const Double_t* step,
                          const Double_t* low, const Double_t* up,
                          PhotonList* photons);
  /**
   Energy-deposition shower-shape function, for use with a TF2.

   Yields the fraction of energy deposited in a tower by a photon, as a function
   of the distance from the photon to the tower center.
   Arguments:
    - x: array with x[0] = x distance from tower center, x[1] = y distance.
    - par: array of fit parameters (fixed, derived from FMS studies).

   Integrates F(x, y) over a tower, with F(x, y) defined as here:
   https://drupal.star.bnl.gov/STAR/blog/leun/2010/aug/02/fms-meeting-20100802

   \todo Provide LaTeX math function in documentation
   */
  static Double_t energyDepositionInTower(Double_t* x, Double_t* par);
  /** Maximum number of photons that can be fit at once. */
  static int maxNFittedPhotons() { return kMaxNPhotons; }

 private:
  static const Int_t kMaxNPhotons = 7;  ///< Maximum number that can be fitted
  /** Disallow copy construction. */
  StFmsClusterFitter(const StFmsClusterFitter&);
  /** Disallow assignment. */
  StFmsClusterFitter& operator=(const StFmsClusterFitter&);
  /**
   Shower shape helper function.

   Evaluates F(x,y) as defined here:
   https://drupal.star.bnl.gov/STAR/blog/leun/2010/aug/02/fms-meeting-20100802

   Used by energyDepositionInTower() to integrate F(x,y) over a tower.

   \todo Provide LaTeX math function in documentation
   */
  static Double_t energyDepositionDistribution(Double_t* x, Double_t* par);
  /**
   Minuit minimization function for fit() routine.
   
   For the purpose of this function and a description of its arguments see
   https://wwwasdoc.web.cern.ch/wwwasdoc/minuit/node14.html

   For its use in ROOT via TMinuit see
   http://root.cern.ch/root/htmldoc/TMinuit.html#TMinuit:SetFCN
   */
  static void minimizationFunctionNPhoton(Int_t& npar, Double_t* grad,
                                          Double_t& fval, Double_t* par,
                                          Int_t iflag);
  /**
   Minuit minimization function for fit2PhotonCluster() routine
   
   Also see comments for minimizationFunction().
   */
  static void minimizationFunction2Photon(Int_t& nparam, Double_t* grad,
                                          Double_t& fval, Double_t* param,
                                          Int_t iflag);
  Double_t mSteps[3 * kMaxNPhotons + 1];  ///< Step size in each fit variable
  Double_t mTowerWidth;  ///< width of one lead glass module
  TMinuit mMinuit;  ///< Minuit fitting interface
  static StFmsTowerCluster::Towers* mTowers;  ///< List of towers to fit
  static Float_t mTowerWidthXY[2];  ///< glass width X, Y in cm
  ClassDef(StFmsClusterFitter, 0)
};  // class StFmsClusterFitter
}  // namespace FMSCluster
#endif  // STROOT_STFMSPOINTMAKER_STFMSCLUSTERFITTER_H_

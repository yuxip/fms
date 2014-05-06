// $Id$
//
// $Log$
/**
 \file      StFmsPointMaker.h
 \brief     Declaration of StFmsPointMaker, the FMS cluster/photon maker
 \author    Yuxi Pan <yuxipan@physics.ucla.edu>
 \author    Thomas Burton <tpb@bnl.gov>
 \date      2014
 \copyright Brookhaven National Lab
 */
#ifndef STROOT_STFMSPOINTMAKER_STFMSPOINTMAKER_H_
#define STROOT_STFMSPOINTMAKER_STFMSPOINTMAKER_H_

#include <map>
#include <vector>

#include "StRoot/StChain/StMaker.h"
#include "StRoot/StPSUTools/StFmsTower.h"
#include "StRoot/StPSUTools/StFmsGeometry.h"

class StFmsCollection;
class StFmsDbMaker;
class StFmsPoint;

namespace FMSCluster {
class StFmsTowerCluster;
class StFmsFittedPhoton;
}  // namespace FMSCluster

/**
 Find FMS clusters and fit clusters with photon hypothesis (shower fit)

 A cluster is a collection of adjacent FMS towers with energy depositions.
 "Point" is a generic term for the energy deposited by individual particle.
 This will typically be a photon, but may be e.g. an electron or energy left
 by a hadron.
 A single cluster may be formed by more than one point, so clusters are fitted
 with a shower shape function to disentangle depositions from different
 particles.
 */
class StFmsPointMaker : public StMaker {
 public:
  /** Constructor */
  StFmsPointMaker(const char* name = "StFmsPointMaker");
  /** Destructor */
  ~StFmsPointMaker();
  /** Called by StMaker when switch to a new run number */
  Int_t InitRun(Int_t runNumber);
  /** Called once per event to process the event */
  Int_t Make();
  /** Called after each event to reset values */
  void Clear(Option_t* option = "");

 private:
  // Define a collection of towers
  typedef std::vector<FMSCluster::StFmsTower> TowerList;
  // Define a collection of tower lists, a TowerList per sub-detector
  // keyed by detector ID
  typedef std::map<int, TowerList> TowerMap;
  /** Disallow copy construction */
  StFmsPointMaker(const StFmsPointMaker&);
  /** Disallow assignment */
  StFmsPointMaker& operator=(const StFmsPointMaker&);
  /**
   Get the StFmsCollection from the current StEvent

   Print messages to LOG_ERROR if StEvent/StFmsCollection cannot be found
   */
  StFmsCollection* getFmsCollection();
  /**
   Perform photon reconstruction in all sub-detectors for a single event

   Populate StFmsCollection with the generated clusters and photons.

   Return kStOk upon success, kStErr in case of an error.
   */
  int clusterEvent();
  /**
   Perform photon reconstruction on a single sub-detector

   Cluster all towers for a sub-detector with ID "detector".
   Update the cluster and photon lists in the provided StFmsCollection with
   the generated clusters and photons.

   Return standard STAR error codes (kStOk, kStWarn, kStErr).
   */
  int clusterDetector(TowerList* towers, int detectorId,
                      StFmsCollection* fmsCollection);
  /**
   Verify that the sum of tower energies is sensible

   Return true if the sum is non-negative and does not exceed the
   center-of-mass energy. Return false otherwise.
   */
  bool validateTowerEnergySum(const TowerList& towers) const;
  /**
   Process a single StFmsTowerCluster and store it as an StFmsCluster

   Update StFmsCollection with the cluster and any photons therein.
   Return true if the cluster is processed, or false if it is skipped due to
   bad values (e.g. unphysical coordinates).
   */
  bool processTowerCluster(FMSCluster::StFmsTowerCluster& towerCluster,
                           int detectorId, StFmsCollection* fmsCollection);
  /** Create a new StFmsPoint from an StFmsFittedPhoton */
  StFmsPoint* makeFmsPoint(const FMSCluster::StFmsFittedPhoton& photon,
                           int detectorId);
  /** Read hits from StEvent and prepare them for clustering */
  bool populateTowerLists();
  /**
   Test channel validity

   Return true if a detector/row/column number physically exists
   Detector values should be as defined as in the database and row and column
   numbers are in the range [1, N].
  */
  bool isValidChannel(int detector, int row, int col);
  StFmsDbMaker* mFmsDbMaker;  //!< Access to FMS database information
  FMSCluster::StFmsGeometry mGeometry;  //!< Access to current FMS geometry
  TowerMap mTowers;  //!< One for each sub-detector, keyed by detector ID
  ClassDef(StFmsPointMaker, 0)
};
#endif  // STROOT_STFMSPOINTMAKER_STFMSPOINTMAKER_H_

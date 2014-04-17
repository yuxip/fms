#ifndef StFmsPointMaker_HH
#define StFmsPointMaker_HH

#include <map>
#include <vector>

#include "StChain/StMaker.h"

#include "StPSUTools/StFmsTower.h"

class StFmsCollection;
class StFmsDbMaker;

namespace FMSCluster { class StFmsGeometry; }
/**
 Find FMS clusters and fit clusters with photon hypothesis (shower fit)
 adapted from PSU code by Yuxi Pan --03/31/2013
 */
class StFmsPointMaker : public StMaker {
 public:
  /** Constructor */
  StFmsPointMaker(const char* name = "StFmsPointMaker");
  /** Destructor */
  ~StFmsPointMaker();
  /** Called at the start to perform one-time initialization steps */
  Int_t Init();
  /** Called by StMaker when switch to a new run number */
  Int_t InitRun(Int_t runNumber);
  /** Called once per event to process the event */
  Int_t Make();
  /** Called after each event to reset values */
  void Clear(Option_t* option = "");
  /** Called at the end to perform one-time cleanup steps */
  Int_t Finish();

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
  /** Perform the actual photon reconstruction */
  int doClustering();
  /**
   Test channel validity

   Return true if a detector/row/column number physically exists
   Detector values should be as defined as in the database and row and column
   numbers are in the range [1, N].
  */
  bool isValidChannel(int detector, int row, int col);
  /** Read hits from StEvent and prepare them for clustering */
  bool populateTowerLists();
  /**
   Verify that the sum of tower energies is sensible

   Return true if the sum is non-negative and does not exceed the
   center-of-mass energy. Return false otherwise.
   */
  bool validateTowerEnergySum(const TowerList& towers) const;
  StFmsDbMaker* mFmsDbMaker;  //!< Access to FMS database information
  FMSCluster::StFmsGeometry* mGeometry;  //!< Access to current FMS geometry
  TowerMap mTowers;  //!< One for each sub-detector, keyed by detector ID
  ClassDef(StFmsPointMaker, 0)
};
#endif

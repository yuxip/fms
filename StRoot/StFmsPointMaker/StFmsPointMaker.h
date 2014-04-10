#ifndef StFmsPointMaker_HH
#define StFmsPointMaker_HH

#include <vector>

#include "StChain/StMaker.h"

#include "StPSUTools/StFmsTower.h"

class StFmsDbMaker;
class StFmsPointCollection;
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
  /** Interface to the actual photon reconstruction */
  Int_t FindPoint();
  /** Return true if a detector/row/column number physically exists */
  Bool_t Legal(Int_t iew, Int_t nstb, Int_t row0, Int_t col0);
  StFmsDbMaker* mFmsDbMaker;  //!< Access to FMS database information
  typedef std::vector<FMSCluster::StFmsTower> TowerList;
  std::vector<TowerList> mTowers; ///< One for each of four FMS sub-detectors
  /** Read hits from StEvent and prepare them for clustering */
  Bool_t populateTowerLists();
  /** all the points (photons) extracted from clusters */
  //StFmsPointCollection*   mFmsPtsColl;  //!
  FMSCluster::StFmsGeometry* fmsgeom;  //!< Access to current FMS geometry
  ClassDef(StFmsPointMaker, 0)
};
#endif

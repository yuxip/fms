// $Id$
//
// $Log$
/**
 \file      StFmsPointMaker.cxx
 \brief     Implementation of StFmsPointMaker, the FMS cluster/photon maker
 \author    Yuxi Pan <yuxipan@physics.ucla.edu>
 \author    Thomas Burton <tpb@bnl.gov>
 \date      2014
 \copyright Brookhaven National Lab
 */
#include "StRoot/StFmsPointMaker/StFmsPointMaker.h"

#include <TLorentzVector.h>
#include <TProcessID.h>

#include "StRoot/St_base/StMessMgr.h"
#include "StRoot/StEvent/StEvent.h"
#include "StRoot/StEvent/StFmsCluster.h"
#include "StRoot/StEvent/StFmsCollection.h"
#include "StRoot/StEvent/StFmsHit.h"
#include "StRoot/StEvent/StFmsPoint.h"
#include "StRoot/StEvent/StRunInfo.h"
#include "StRoot/StFmsDbMaker/StFmsDbMaker.h"

#include "StRoot/StFmsPointMaker/StFmsEventClusterer.h"
#include "StRoot/StFmsPointMaker/StFmsFittedPhoton.h"
#include "StRoot/StFmsPointMaker/StFmsTower.h"
#include "StRoot/StFmsPointMaker/StFmsTowerCluster.h"

namespace {
// Calculate a 4 momentum from a direction/momentum vector and energy
// assuming zero mass i.e. E = p
TLorentzVector compute4Momentum(const TVector3& xyz, Double_t energy) {
  TVector3 mom3 = xyz.Unit() * energy;  // Momentum vector with m = 0
  return TLorentzVector(mom3, energy);
}
}  // unnamed namespace

StFmsPointMaker::StFmsPointMaker(const char* name)
    : StMaker(name), mFmsDbMaker(nullptr), mObjectCount(0) { }

StFmsPointMaker::~StFmsPointMaker() { }

Int_t StFmsPointMaker::InitRun(Int_t runNumber) {
  // Ensure we can access database information
  mFmsDbMaker = static_cast<StFmsDbMaker*>(GetMaker("fmsDb"));
  if (!mFmsDbMaker) {
    return kStErr;
  }  // if
  // Set up geometry, which stays constant for each run
  if (!mGeometry.initialize(mFmsDbMaker)) {
    // Return an error if geometry initialization fails
    return kStErr;
  }  // if
  return StMaker::InitRun(runNumber);
}

Int_t StFmsPointMaker::Make() {
  // Cache the current count of referenced objects, as discussed here
  // http://root.cern.ch/root/htmldoc/TRef.html
  /** \note Object count caching should probably be done in e.g. StMuDstMaker
      however we do it in StFmsPointMaker for now until we can verify changes
      with the MuDST coordinator. It should work OK as I don't think any other
      STAR makers use TRef. */
  mObjectCount = TProcessID::GetObjectCount();
  if (!populateTowerLists()) {
    LOG_ERROR << "StFmsPointMaker::Make() - failed to initialise tower " <<
      "lists for the event" << endm;
  }  // if
  if (clusterEvent() == kStOk) {
     return StMaker::Make();
  }  // if
  return kStErr;
}

void StFmsPointMaker::Clear(Option_t* option) {
  mTowers.clear();
  // Reset the count of referenced objects to the value before the previous
  // call to Make(), in order prevent ever-growing table of objects
  TProcessID::SetObjectCount(mObjectCount);
  StMaker::Clear(option);
}

StFmsCollection* StFmsPointMaker::getFmsCollection() {
  StEvent* event = static_cast<StEvent*>(GetInputDS("StEvent"));
  StFmsCollection* fms = nullptr;
  if (event) {
    fms = event->fmsCollection();
  }  // if
  if (!fms) {
    LOG_ERROR << "StFmsPointMaker did not find "
              << "an StFmsCollection in StEvent" << endm;
  }  // if
  return fms;
}

int StFmsPointMaker::clusterEvent() {
  StFmsCollection* fmsCollection = getFmsCollection();
  if (!fmsCollection) {
    return kStErr;
  }  // if
  for (auto i = mTowers.begin(); i != mTowers.end(); ++i) {
    if (!validateTowerEnergySum(i->second)) {
      continue;  // To remove LED trails
    }  // if
    clusterDetector(&i->second, i->first, fmsCollection);
  }  // for
  return kStOk;
}

/* Perform photon reconstruction on a single sub-detector */
int StFmsPointMaker::clusterDetector(TowerList* towers, const int detectorId,
                                     StFmsCollection* fmsCollection) {
  FMSCluster::StFmsEventClusterer clustering(&mGeometry, detectorId);
  // Perform tower clustering, skip this subdetector if an error occurs
  if (!clustering.cluster(towers)) {  // Cluster tower list
    return kStWarn;
  }  // if
  // Saved cluster info into StFmsCluster
  auto& clusters = clustering.clusters();
  for (auto cluster = clusters.begin(); cluster != clusters.end(); ++cluster) {
    processTowerCluster(cluster->get(), detectorId, fmsCollection);
  }  // for
  return kStOk;
}

bool StFmsPointMaker::validateTowerEnergySum(const TowerList& towers) const {
  // Attempt to get center-of-mass energy from StRunInfo.
  // If it can't be accessed assume 500 GeV running.
  float centerOfMassEnergy(500.);
  const StEvent* event = static_cast<const StEvent*>(GetInputDS("StEvent"));
  if (event) {
    if (event->runInfo()) {
      centerOfMassEnergy = event->runInfo()->centerOfMassEnergy();
    }  // if
  }  // if
  // Sum tower energies and test validity of the sum
  float Esum = 0.f;
  typedef TowerList::const_iterator TowerIter;
  for (TowerIter i = towers.begin(); i != towers.end(); ++i) {
    Esum += i->hit()->energy();
  }  // for
  return Esum >= 0.f && Esum <= centerOfMassEnergy;
}

bool StFmsPointMaker::processTowerCluster(
    FMSCluster::StFmsTowerCluster* towerCluster,
    const int detectorId,
    StFmsCollection* fmsCollection) {
  // Update the StFmsCluster object we want to store in StEvent with information
  // not automatically propagated via StFmsTowerCluster
  StFmsCluster* cluster = towerCluster->cluster();
  // Skip clusters that don't have physically sensible coordinates
  if (!(cluster->x() > 0. && cluster->y() > 0.)) {
    return false;
  }  // if
  cluster->setDetectorId(detectorId);
  // Cluster id is id of the 1st photon, not necessarily the highest-E photon
  cluster->setId(305 + 20 * detectorId + fmsCollection->numberOfPoints());
  // Cluster locations are in column-row coordinates so convert to cm
  TVector3 xyz = mGeometry.columnRowToGlobalCoordinates(
    cluster->x(), cluster->y(), detectorId);
  cluster->setFourMomentum(compute4Momentum(xyz, cluster->energy()));
  // Save photons reconstructed from this cluster
  for (Int_t np = 0; np < towerCluster->photons().size(); np++) {
    StFmsPoint* point = makeFmsPoint(towerCluster->photons()[np], detectorId);
    point->setDetectorId(detectorId);
    point->setId(305 + 20 * detectorId + fmsCollection->numberOfPoints());
    point->setParentClusterId(cluster->id());
    point->setNParentClusterPhotons(towerCluster->photons().size());
    point->setCluster(cluster);
    // Add it to both the StFmsCollection and StFmsCluster
    // StFmsCollection owns the pointer, the cluster merely references it
    fmsCollection->points().push_back(point);
    cluster->points().push_back(point);
  }  // for
  // Save the tower hit info.
  auto& towers = towerCluster->towers();
  for (auto i = towers.begin(); i != towers.end(); ++i) {
    if ((*i)->hit()->adc() >= 1) {  // Min ADC = 1
      cluster->hits().push_back((*i)->hit());
    }  // if
  }  // for
  // Release StFmsCluster held by towerCluster to pass ownership to
  // StFmsCollection (and hence StEvent).
  fmsCollection->addCluster(towerCluster->release());
  return true;
}

StFmsPoint* StFmsPointMaker::makeFmsPoint(
    const FMSCluster::StFmsFittedPhoton& photon, const int detectorId) {
  StFmsPoint* point = new StFmsPoint;
  point->setEnergy(photon.energy);
  // Calculate photon 4 momentum
  // StFmsFittedPhoton position is in detector-local (x, y) cm coordinates
  // Convert to global STAR coordinates for StFmsPoint
  TVector3 xyz = mGeometry.localToGlobalCoordinates(
    photon.x, photon.y, detectorId);
  point->setX(xyz.x());
  point->setY(xyz.y());
  point->setFourMomentum(compute4Momentum(xyz, point->energy()));
  return point;
}

bool StFmsPointMaker::populateTowerLists() {
  StFmsCollection* fmsCollection = getFmsCollection();
  if (!fmsCollection) {
      return false;
  }  // if
  auto& hits = fmsCollection->hits();
  for (auto i = hits.begin(); i != hits.end(); ++i) {
    StFmsHit* hit = *i;
    const int detector = hit->detectorId();
    const int row = mFmsDbMaker->getRowNumber(detector, hit->channel());
    const int column = mFmsDbMaker->getColumnNumber(detector, hit->channel());
    if (!isValidChannel(detector, row, column)) {
      continue;
    }  // if
    if (hit->adc() > 0) {
      // Insert a tower list for this detector ID if there isn't one already
      // This method is faster than using find() followed by insert()
      // http://stackoverflow.com/questions/97050/stdmap-insert-or-stdmap-find
      auto low = mTowers.lower_bound(detector);
      if (low == mTowers.end() || mTowers.key_comp()(detector, low->first)) {
        mTowers.insert(TowerMap::value_type(detector, TowerList()));
      }  // if
      FMSCluster::StFmsTower tower(hit);
      // Ensure tower information is valid before adding
      if (tower.initialize(mFmsDbMaker)) {
        mTowers[detector].push_back(tower);
      }  // if
    }  // if
  }  // for
  return true;
}

/* Test channel validity by detector and row, column in the range [1, N] */
bool StFmsPointMaker::isValidChannel(int detector, int row, int column) {
  // Simplest check first, test lower bounds are valid
  if (row < 1 || column < 1) {
    return false;
  }  // if
  // Omit gaps in the detector
  switch (detector) {
    case FMSCluster::kFmsNorthLarge:  // Deliberate fall-through
    case FMSCluster::kFmsSouthLarge:  // Large-cell FMS sub-detector
      if (fabs(row - 17.5) < 8 && column < 9) {  // Central hole
        return false;
      }  // if
      // This cuts off a 7x7 triangle from the corners
      if (fabs(17.5 - row) + column > 27.) {
        return false;
      }  // if
      break;
    case FMSCluster::kFmsNorthSmall:  // Deliberate fall-through
    case FMSCluster::kFmsSouthSmall:  // Small-cell FMS sub-detector
      if (fabs(row - 12.5) < 5 && column < 6) {  // Central hole
        return false;
      }  // if
      break;
    default:  // Don't currently support non-FMS sub-detectors
      return false;
  }  // switch (detector)
  // Test row and column number against the numbers stored in the database for
  // this detector. Leave this to last to avoid database calls when possible.
  // Also serves as a double-check on detector, as the database will
  // return -1 for both numbers in case of an invalid detector number.
  if (mFmsDbMaker) {
    const int nRows = mFmsDbMaker->nRow(detector);
    if (nRows < 0 || row > nRows) {
      return false;
    }  // if
    const int nColumns = mFmsDbMaker->nColumn(detector);
    if (nColumns < 0 || column > nColumns) {
      return false;
    }  // if
  }  // if
  return true;
}

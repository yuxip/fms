#include "StFmsPointMaker.h"

#include <algorithm>

#include <boost/foreach.hpp>

#include <TLorentzVector.h>

#include "St_base/StMessMgr.h"
#include "StEvent/StEvent.h"
#include "StEvent/StFmsCluster.h"
#include "StEvent/StFmsCollection.h"
#include "StEvent/StFmsHit.h"
#include "StEvent/StFmsPoint.h"
#include "StEvent/StTriggerData.h"
#include "StFmsDbMaker/StFmsDbMaker.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"

#include "StPSUTools/StFmsGeometry.h"
#include "StPSUTools/StFmsTowerCluster.h"
#include "StPSUTools/StFmsTower.h"
#include "StPSUTools/StFmsClusterFinder.h"  // Defines ClusterList
#include "StPSUTools/StFmsEventClusterer.h"

namespace {
// Calculate a 4 momentum from a direction/momentum vector and energy
// assuming zero mass i.e. E = p
TLorentzVector compute4Momentum(const TVector3& xyz, Double_t energy) {
  TVector3 mom3 = xyz.Unit() * energy;  // Momentum vector with m = 0
  return TLorentzVector(mom3, energy);
}
}  // unnamed namespace

StFmsPointMaker::StFmsPointMaker(const char* name)
    : StMaker(name), mFmsDbMaker(NULL), mGeometry(NULL) { }

StFmsPointMaker::~StFmsPointMaker() {
  LOG_DEBUG << "StFmsPointMaker:: destructor " << endm;
}

void StFmsPointMaker::Clear(Option_t* option) {
  LOG_DEBUG << "StFmsPointMaker::Clear() " << endm;
  LOG_DEBUG << "after StFmsPointMaker::Clear()" << endm;
  mTowers.clear();
  StMaker::Clear(option);
}

Int_t StFmsPointMaker::Init() {
  LOG_DEBUG << "StFmsPointMaker::Init() " << endm;
  return StMaker::Init();
}

Int_t StFmsPointMaker::InitRun(Int_t runNumber) {
  // Ensure we can access database information
  mFmsDbMaker = static_cast<StFmsDbMaker*>(GetMaker("fmsDb"));
  if (!mFmsDbMaker) {
    return kStErr;
  }  // if
  // Set up geometry, which stays constant for each run
  // Only allocate new space in the beginning, not in between runs
  if (!mGeometry) {
    mGeometry = new FMSCluster::StFmsGeometry;
    if (!mGeometry->initialize(mFmsDbMaker)) {
      // Return an error if geometry initialization fails
      return kStErr;
    }  // if
  }  // if
  return StMaker::InitRun(runNumber);
}

Int_t StFmsPointMaker::Finish() {
  LOG_DEBUG << "StFmsPointMaker::Finish() " << endm;
  return StMaker::Finish();
}

Int_t StFmsPointMaker::Make() {
  LOG_DEBUG << "StFmsPointMaker::Make() " << endm;
  if (!populateTowerLists()) {
    LOG_ERROR << "StFmsPointMaker::Make() - failed to initialise tower " <<
      "lists for the event" << endm;
  }  // if
  if (doClustering() == kStOk){
     LOG_DEBUG << "Cluster finder returns successfully" <<endm;
     return StMaker::Make();
  }  // if
  LOG_INFO << " cluster finder returns error!!!" << endm;
  return kStErr;
}

StFmsCollection* StFmsPointMaker::getFmsCollection() {
  StEvent* event = static_cast<StEvent*>(GetInputDS("StEvent"));
  StFmsCollection* fms(NULL);
  if (event) {
    fms = event->fmsCollection();
  } else {
    LOG_ERROR << "StFmsPointMaker did not find StEvent" << endm;
  }  // if
  if (!fms) {
    LOG_ERROR << "StFmsPointMaker did not find "
              << "an StFmsCollection in StEvent" << endm;
  }  // if
  return fms;
}

int StFmsPointMaker::doClustering() {
  LOG_DEBUG << " StFmsPointMaker::FindPoint() " << endm;
  StFmsCollection* fmsCollection = getFmsCollection();
  if (!fmsCollection) {
    return kStErr;
  }  // if
  BOOST_FOREACH(TowerMap::value_type& subdetector, mTowers) {
    if (!validateTowerEnergySum(subdetector.second)) {
      continue;  // To remove LED trails
    }  // if
    const int detectorId = subdetector.first;
    using namespace FMSCluster;
    StFmsEventClusterer clustering(mGeometry, detectorId);
    // Perform tower clustering, skip this subdetector if an error occurs
    if (!clustering.cluster(&subdetector.second)) {  // Cluster tower list
      continue;
    }  // if
    // Saved cluser info into StFmsCluster
    Int_t iPh = 0;  // Sequence # in StFmsEventClusterer::photons[]
    BOOST_FOREACH(StFmsTowerCluster& towerCluster, clustering.clusters()) {
      StFmsCluster* cluster = towerCluster.cluster();
      // Cluster id = id of the 1st photon, not necessarily the highE photon
      cluster->setDetector(detectorId);
      cluster->setId(305 + 20 * detectorId + iPh);
      // Skip clusters that don't have physically sensible coordinates
      if (!(cluster->x() > 0. && cluster->y() > 0.)) {
        continue;
      }  // if
      // Cluster locations are in column-row coordinates (not cm)
      TVector3 xyz = mGeometry->columnRowToGlobalCoordinates(
        cluster->x(), cluster->y(), clustering.detector());
      cluster->setFourMomentum(compute4Momentum(xyz, cluster->energy()));
      // Save photons reconstructed from this cluster
      for (Int_t np = 0; np < cluster->nPhotons(); np++) {
        StFmsPoint* clpoint = new StFmsPoint;
        clpoint->setEnergy(towerCluster.photons()[np].energy);
        clpoint->setId(305 + 20 * detectorId + iPh);
        iPh++;
        // Calculate photon 4 momentum
        // Photon position is in local (x, y) cm coordinates
        TVector3 xyzph = mGeometry->localToGlobalCoordinates(
          towerCluster.photons()[np].xPos, towerCluster.photons()[np].yPos,
          clustering.detector());
        clpoint->setXYZLab(xyzph);
        clpoint->setFourMomentum(compute4Momentum(xyzph, clpoint->energy()));
        clpoint->setParentClusterId(cluster->id());
        clpoint->setNParentClusterPhotons(cluster->nPhotons());
        // Add it to both the StFmsCollection and StFmsCluster
        // StFmsCollection owns the pointer, the cluster merely references it
        fmsCollection->points().push_back(clpoint);
        cluster->points().push_back(clpoint);
      }  // for
      // Save the tower hit info.
      BOOST_FOREACH(const StFmsTower* tow, towerCluster.towers()) {
        if (tow->hit()->adc() >= 1) {  // Min ADC = 1
          cluster->hits().push_back(tow->hit());
        }  // if
      }  // BOOST_FOREACH(towers)
      fmsCollection->addCluster(towerCluster.release());
    }  // BOOST_FOREACH(clusters)
  }  // BOOST_FOREACH(subdetectors)
  LOG_DEBUG << "StFmsPointMaker::FindPoint() --StFmsCluster collections filled "
    << endm;
  return kStOk;
}

bool StFmsPointMaker::populateTowerLists() {
  StFmsCollection* fmsCollection = getFmsCollection();
  if (!fmsCollection) {
      return false;
  }  // if
  StSPtrVecFmsHit& hits = fmsCollection->hits();
  for (StSPtrVecFmsHitIterator i = hits.begin(); i != hits.end(); ++i) {
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
      TowerMap::iterator low = mTowers.lower_bound(detector);
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
  LOG_INFO << "Tower summary:" << endm;
  BOOST_FOREACH(TowerMap::value_type& subdetector, mTowers) {
    LOG_INFO << "detector " << subdetector.first << " nTowers = " <<
      subdetector.second.size() << endm;
  }  // BOOST_FOREACH
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
    case 8:  // Deliberate fall-through
    case 9:  // Large-cell FMS sub-detector
      if (fabs(row - 17.5) < 8 && column < 9) {  // Central hole
        return false;
      }  // if
      // This cuts off a 7x7 triangle from the corners
      if (fabs(17.5 - row) + column > 27.) {
        return false;
      }  // if
      break;
    case 10:  // Deliberate fall-through
    case 11:  // Small-cell FMS sub-detector
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
  BOOST_FOREACH(const FMSCluster::StFmsTower& tower, towers) {
    Esum += tower.hit()->energy();
  }  // BOOST_FOREACH
  return Esum >= 0.f && Esum <= centerOfMassEnergy;
}

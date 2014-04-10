#include "StFmsPointMaker.h"

#include <algorithm>
#include <cassert>

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

#ifndef __CINT__
typedef FMSCluster::ClusterList::iterator ClusterIter;
typedef FMSCluster::ClusterList::const_iterator ClusterConstIter;
#endif  // __CINT__

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
  return kStOk;
}

Int_t StFmsPointMaker::Make() {
  LOG_DEBUG << "StFmsPointMaker::Make() " << endm;
  if (!populateTowerLists()) {
    LOG_ERROR << "StFmsPointMaker::Make() - failed to initialise tower " <<
      "lists for the event" << endm;
  }  // if
  if (doClustering() == kStOk){
     LOG_DEBUG << "Cluster finder returns successfully" <<endm;
     return kStOk;
  }  // if
  LOG_INFO << " cluster finder returns error!!!" << endm;
  return kStErr;
}
  
int StFmsPointMaker::doClustering() {
  LOG_DEBUG << " StFmsPointMaker::FindPoint() " << endm;
  StEvent* event = static_cast<StEvent*>(GetDataSet("StEvent"));
  if (!event) {
    LOG_ERROR << "StFmsPointMaker::populateTowerLists() did not find "
      << "an StEvent" << endm;
      return false;
  }  // if
  StFmsCollection* fmsCollection = event->fmsCollection();
  if (!fmsCollection) {
    LOG_ERROR << "StFmsPointMaker::populateTowerLists() did not find "
      << "an StFmsCollection in StEvent" << endm;
      return false;
  }  // if
  for (Int_t instb = 0; instb < 4; instb++) {
    TowerList& towers = mTowers.at(instb);
    Float_t Esum = 0.f;
    for (TowerList::const_iterator i = towers.begin(); i != towers.end(); ++i) {
      Esum += i->hit()->energy();
    }  // for
    if (Esum == 0 || Esum > 500) {
      continue;  // To remove LED trails, for pp500 GeV
    }  // if
    Int_t detectorId = instb + 8;  // FMS IDs from 8 to 11
    FMSCluster::StFmsEventClusterer clustering(mGeometry, detectorId);
    // Perform tower clustering, skip this subdetector if an error occurs
    if (!clustering.cluster(&towers)) {
      continue;
    }  // if
    // Saved cluser info into StFmsCluster
    Int_t iPh = 0;  // Sequence # in StFmsEventClusterer::photons[]
    FMSCluster::ClusterList& clusters = clustering.clusters();
    for (ClusterIter ci = clusters.begin(); ci != clusters.end(); ++ci) {
      StFmsCluster* cluster = ci->cluster();
      // Cluster id = id of the 1st photon, not necessarily the highE photon
      cluster->SetNstb(instb + 1);
      cluster->SetClusterId(305 + 20 * instb + iPh);
      // Skip clusters that don't have physically sensible coordinates
      if (!(cluster->GetX0() > 0. && cluster->GetY0() > 0.)) {
        continue;
      }  // if
      // Cluster locations are in column-row coordinates (not cm)
      TVector3 xyz = mGeometry->columnRowToGlobalCoordinates(
        cluster->GetX0(), cluster->GetY0(), clustering.detector());
      cluster->SetFourMomentum(compute4Momentum(xyz, cluster->GetEnergy()));
      // Save photons reconstructed from this cluster
      for (Int_t np = 0; np < cluster->GetNphoton(); np++) {
        StFmsPoint* clpoint = new StFmsPoint;
        clpoint->SetEnergy(ci->photons()[np].energy);
        clpoint->SetPhotonId(305 + 20 * instb + iPh);
        iPh++;
        // Calculate photon 4 momentum
        // Photon position is in local (x, y) cm coordinates
        TVector3 xyzph = mGeometry->localToGlobalCoordinates(
          ci->photons()[np].xPos, ci->photons()[np].yPos,
          clustering.detector());
        clpoint->SetPointXYZLab(xyzph);
        clpoint->SetFourMomentum(compute4Momentum(xyzph, clpoint->GetEnergy()));
        clpoint->SetParentCluId(cluster->GetClusterId());
        clpoint->SetParentNclPh(cluster->GetNphoton());
        // Add it to both the StFmsCollection and StFmsCluster
        // StFmsCollection owns the pointer, the cluster merely references it
        fmsCollection->points().push_back(clpoint);
        cluster->points().push_back(clpoint);
      }  // for
      // Save the tower hit info.
      BOOST_FOREACH(const FMSCluster::StFmsTower* tow, ci->towers()) {
        if (tow->hit()->adc() >= 1) {  // Min ADC = 1
          cluster->hits().push_back(tow->hit());
          // Make sure the hit is in the original collection
          assert(std::find(fmsCollection->hits().begin(),
                           fmsCollection->hits().end(),
                           cluster->hits().back()) !=
                 fmsCollection->hits().end());
        }  // if
      }  // BOOST_FOREACH
      fmsCollection->addCluster(cluster);
    }  // for loop over clusters
  }  // for loop over NSTB
  LOG_DEBUG << "StFmsPointMaker::FindPoint() --StFmsCluster collections filled "
    << endm;
  return kStOk;
}

bool StFmsPointMaker::populateTowerLists() {
  StEvent* event = static_cast<StEvent*>(GetDataSet("StEvent"));
  if (!event) {
    LOG_ERROR << "StFmsPointMaker::populateTowerLists() did not find "
      << "an StEvent" << endm;
      return false;
  }  // if
  StFmsCollection* fmsCollection = event->fmsCollection();
  if (!fmsCollection) {
    LOG_ERROR << "StFmsPointMaker::populateTowerLists() did not find "
      << "an StFmsCollection in StEvent" << endm;
      return false;
  }  // if
  mTowers.assign(4, TowerList());
  StSPtrVecFmsHit& hits = fmsCollection->hits();
  for (StSPtrVecFmsHitIterator i = hits.begin(); i != hits.end(); ++i) {
    StFmsHit* hit = *i;
    Int_t row = mFmsDbMaker->getRowNumber(hit->detectorId(), hit->channel());
    Int_t column = mFmsDbMaker->getColumnNumber(hit->detectorId(),
                                                hit->channel());
    Int_t nstb = 0;
    if (hit->detectorId() > 7 && hit->detectorId() < 12) {
      nstb = hit->detectorId() - 7;
    }  // if
    if (nstb == 1 || nstb == 2) {
      // because channel geometry in the database assigns row1 as the bottom row
      row = 35 - row;
    } else if (nstb == 3 || nstb == 4) {
      row = 25 - row;
    }  // if
    Int_t ew  = 2;  // east=1, west=2
    if (!isValidChannel(ew, nstb, row - 1, column - 1)) {
      continue;
    }  // if
    unsigned index = hit->detectorId() - 8;  // FMS IDs range from 8 to 11
    if (hit->adc() > 0 && index >= 0 && index < mTowers.size()) {
      FMSCluster::StFmsTower tower(hit);
      // Ensure tower information is valid before adding
      if (tower.initialize(mFmsDbMaker)) {
        mTowers.at(index).push_back(tower);
      }  // if
    }  // if
  }  // for
  return true;
}

bool StFmsPointMaker::isValidChannel(int iew, int nstb, int row0, int col0) {
  // nstb starts from 1, row0, col0 starts from 0  
  if (iew > 0 && iew < 2) {
    return false;
  }  // if
  if (nstb < 1 || nstb > 4) {
    return false;
  }  // if
  if (nstb > 2) {
    if (row0 < 0 || row0 > 23) {
      return false;
    }  // if
    if (col0 < 0 || col0 > 11) {
      return false;
    }  // if
    if (fabs(1. * row0 - 11.5) < 5 && col0 < 5) {
      return false;
    }  // if
  } else {
    if (row0 < 0 || row0 > 33) {
      return false;
    }  // if
    if (col0 < 0 || col0 > 16) {
      return false;
    }  // if
    if (fabs(1. * row0 - 16.5) < 8 && col0 < 8) {
      return false;
    }  // if
    if (row0 < col0 - 9.5) {
      return false;
    }  // if
    if (33 - row0 < col0 - 9.5) {
      return false;
    }  // if
  } // if (nstb > 2)
  return true;
}

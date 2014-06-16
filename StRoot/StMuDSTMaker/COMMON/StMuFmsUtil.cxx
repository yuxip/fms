/***************************************************************************
 *
 * $Id: StMuFmsUtil.cxx,v 1.1 2010/01/25 03:57:39 tone421 Exp $
 *
 * Author: Jingguo Ma, Jan 2010
 ***************************************************************************
 *
 * Description: FMS Util to convert between StEvent and MuDst
 *
 ***************************************************************************
 *
 * $Log: StMuFmsUtil.cxx,v $
 * Revision 1.1  2010/01/25 03:57:39  tone421
 * Added FMS and Roman pot arrays
 *
 **************************************************************************/
#include "StFmsCluster.h"
#include "StMuFmsCluster.h"
#include "StMuFmsHit.h"
#include "StMuFmsPoint.h"
#include "StMuFmsUtil.h"
#include "StMuFmsCollection.h"
#include "StEvent.h"
#include "StMessMgr.h"
#include "StEventTypes.h"

#include <algorithm>  // For std::find
#include <iterator>  // For std::distance

#include <TCollection.h>  // For TIter
#include <TRefArray.h>
#include <TVector3.h>

ClassImp(StMuFmsUtil)

namespace {
/*
 Return the index of an element in an StPtrVec-like array.

 The element should be of the pointer type in the StPtrVec.
 Return -1 if the element cannot be located in the array.
 See StEvent/StArray.h for the definition of St[S]PtrVec arrays.
 */
template<class StPtrVec, class Element>
int findElementIndex(const StPtrVec& array, const Element* element) {
  // This typedef is equivalent to an iterator to an Element in the StPtrVec,
  // as defined in StArray.h. The StPtrVec definition does not include a
  // template iterator typedef itself, but this the same as what
  // StPtrVec::const_iterator would do.
  typedef Element* const * ElementConstIter;
  // Find where the element is in the array
  ElementConstIter location = std::find(array.begin(), array.end(), element);
  // Ensure the desired element actually is in the array, as the behaviour
  // of std::distance is undefined otherwise.
  // Return -1 if the element isn't in the array, otherwise return its index
  if (location == array.end()) {  // Not in the array
    return -1;
  } else {
    return std::distance(array.begin(), location);
  }  // if
}
}  // namespace

StMuFmsUtil::StMuFmsUtil()
{
}
StMuFmsUtil::~StMuFmsUtil()
{
}

StMuFmsCollection* StMuFmsUtil::getMuFms(StFmsCollection *fmscol)
{
  if(!fmscol) return NULL;
  StMuFmsCollection* muFms=new StMuFmsCollection();
  fillMuFms(muFms,fmscol);
  return muFms;
}  

StFmsCollection* StMuFmsUtil::getFms(StMuFmsCollection* muFms)
{
  if(!muFms) return NULL;
  
  StFmsCollection *fms=new StFmsCollection();
  fillFms(fms,muFms);
  return fms;
}

void StMuFmsUtil::fillMuFms(StMuFmsCollection *muFms,StFmsCollection *fmscol)
{
  if(!fmscol) return;
  if(!muFms) return;
  // Do hits and points before clusters, so that the hit and point lists are
  // populated before we try to set hit- and photon-in-cluster information
  // during the cluster loop
  fillMuFmsHits(muFms, fmscol);
  fillMuFmsPoints(muFms, fmscol);
  fillMuFmsClusters(muFms, fmscol);
  // Now we need to go back and set parent cluster of each point, now that the
  // cluster list in StMuFmsCollection is populated (as these are the clusters
  // we reference).
  setMuFmsPointParentClusters(muFms, fmscol);
}

void StMuFmsUtil::fillFms(StFmsCollection* fmscol,StMuFmsCollection* muFms)
{
  if(!muFms) return;
  if(!fmscol) return;
  fillFmsHits(fmscol, muFms);
  fillFmsPoints(fmscol, muFms);
  fillFmsClusters(fmscol, muFms);
  // Now we need to go back and set parent cluster of each point, now that the
  // cluster list in StFmsCollection is populated (as these are the clusters
  // we reference).
  setFmsPointParentClusters(fmscol, muFms);
}

void StMuFmsUtil::fillMuFmsHits(StMuFmsCollection* muFms,
                                StFmsCollection* fmscol) {
  StSPtrVecFmsHit vecHit = fmscol->hits();
  for(unsigned int i=0; i<fmscol->numberOfHits(); i++){
    unsigned short detId = vecHit[i]->detectorId();
    unsigned short ch    = vecHit[i]->channel();
    unsigned short crate = vecHit[i]->qtCrate();
    unsigned short slot  = vecHit[i]->qtSlot();
    unsigned short qtch  = vecHit[i]->qtChannel();
    unsigned short adc   = vecHit[i]->adc();
    unsigned short tdc   = vecHit[i]->tdc();
    float          ene   = vecHit[i]->energy();
    muFms->addHit();
    StMuFmsHit* muFmsHit = muFms->getHit(i);
    muFmsHit->setDetectorId(detId);
    muFmsHit->setChannel(ch);
    muFmsHit->setQtCrate(crate);
    muFmsHit->setQtSlot(slot);
    muFmsHit->setQtChannel(qtch);
    muFmsHit->setAdc(adc);
    muFmsHit->setTdc(tdc);
    muFmsHit->setEnergy(ene);
  }
}

void StMuFmsUtil::fillMuFmsClusters(StMuFmsCollection* muFms,
                                    StFmsCollection* fmscol) {
  // Fill clusters
  for (unsigned i(0); i < fmscol->numberOfClusters(); ++i) {
    const StFmsCluster* cluster = fmscol->clusters()[i];
    muFms->addCluster();  // Expand StMuFmsCollection cluster array by 1
    StMuFmsCluster* muCluster = muFms->getCluster(i);
    muCluster->setDetectorId(cluster->detectorId());
    muCluster->setCategory(cluster->category());
    muCluster->setEnergy(cluster->energy());
    muCluster->setX(cluster->x());
    muCluster->setY(cluster->y());
    // Propagate hits-in-cluster information
    // Remember, clusters don't *own* hits, they just reference them.
    // For each StFmsHit in the cluster, find the index of that hit in the main
    // StFmsCollection hit array. Then add the StMuFmsHit (from the main
    // StMuFmsCollection hit array) at the same index to the StMuFmsCluster.
    StPtrVecFmsHitConstIterator hit;  // Iterate over StFmsHits
    for (hit = cluster->hits().begin(); hit != cluster->hits().end(); ++hit) {
      const int index = findElementIndex(fmscol->hits(), *hit);
      if (index != -1) {
        muCluster->hits()->Add(muFms->getHit(index));
      }  // if
    }  // for
    // Do the same procedure for photon-in-cluster information
    StPtrVecFmsPointConstIterator p;
    for (p = cluster->points().begin(); p != cluster->points().end(); ++p) {
      const int index = findElementIndex(fmscol->points(), *p);
      if (index != -1) {
        muCluster->photons()->Add(muFms->getPoint(index));
      } // if
    }  // for
  }  // for
}

void StMuFmsUtil::fillMuFmsPoints(StMuFmsCollection* muFms,
                                  StFmsCollection* fmscol) {
  for (unsigned i(0); i < fmscol->numberOfPoints(); ++i) {
    const StFmsPoint* point = fmscol->points()[i];
    StMuFmsPoint* muPoint = muFms->addPoint();
    if (point && muPoint) {
      muPoint->set(*point);
    }  // if
  }  // for
}

void StMuFmsUtil::setMuFmsPointParentClusters(StMuFmsCollection* muFms,
                                              StFmsCollection* fmscol) {
  for (unsigned i(0); i < muFms->numberOfPoints(); ++i) {
    // Points and clusters in the StMuFmsCollection and StFmsCollection are in
    // the same order, so we get the corresponding objects just by index
    const StFmsPoint* point = fmscol->points().at(i);
    if (!point) {
      continue;
    }  // if
    // Find the index of the point's parent cluster in the main cluster list
    const int index = findElementIndex(fmscol->clusters(), point->cluster());
    // If we found it, set the StMuFmsPoint's parent cluster to be the
    // corresponding cluster in the StMuFmsCollection
    if (index != -1) {
      StMuFmsPoint* muPoint = muFms->getPoint(i);
      if (muPoint) {
        muPoint->setCluster(muFms->getCluster(index));
      }   // if
    }  // if
  }  // for
}

void StMuFmsUtil::fillFmsHits(StFmsCollection* fmscol,
                              StMuFmsCollection* muFms) {
  // Using TIter to iterate is safe in the case of hits being NULL
  TIter next(muFms->getHitArray());
  StMuFmsHit* muHit(NULL);
  while ((muHit = static_cast<StMuFmsHit*>(next()))) {
    fmscol->addHit(new StFmsHit);
    StFmsHit* hit = fmscol->hits().back();
    hit->setDetectorId(muHit->detectorId());
    hit->setChannel(muHit->channel());
    hit->setQtCrate(muHit->qtCrate());
    hit->setQtSlot(muHit->qtSlot());
    hit->setQtChannel(muHit->qtChannel());
    hit->setAdc(muHit->adc());
    hit->setTdc(muHit->tdc());
    hit->setEnergy(muHit->energy());
  }  // while
}

void StMuFmsUtil::fillFmsClusters(StFmsCollection* fmscol,
                                  StMuFmsCollection* muFms) {
  // Using TIter to iterate is safe in the case of clusters being NULL
  TIter next(muFms->getClusterArray());
  StMuFmsCluster* muCluster(NULL);
  while ((muCluster = static_cast<StMuFmsCluster*>(next()))) {
    // Create an StFmsCluster from this StMuFmsCluster
    fmscol->addCluster(new StFmsCluster);
    StFmsCluster* cluster = fmscol->clusters().back();
    cluster->setDetectorId(muCluster->detectorId());
    cluster->setCategory(muCluster->category());
    cluster->setNTowers(muCluster->hits()->GetEntries());
    cluster->setNPhotons(muCluster->photons()->GetEntries());
    cluster->setEnergy(muCluster->energy());
    cluster->setX(muCluster->x());
    cluster->setY(muCluster->y());
    // StMuFmsCluster does not store all the information in StFmsCluster, so
    // sigmaMin, sigmaMax, chi2Ndf1Photon, chi2Ndf2Photon and id will not be
    // filled
    /** \todo fill 4-momentum. Requires adding z field to StMuFmsPoint */
    /** \todo propagate hit- and photon-in-cluster information */
  }  // while
}

void StMuFmsUtil::fillFmsPoints(StFmsCollection* fmscol,
                                StMuFmsCollection* muFms) {
  // Using TIter to iterate is safe in the case of points being NULL
  TIter next(muFms->getPointArray());
  StMuFmsPoint* muPoint(NULL);
  while ((muPoint = static_cast<StMuFmsPoint*>(next()))) {
    // Create an StFmsPoint from this StMuFmsPoint
    fmscol->addPoint(new StFmsPoint);
    StFmsPoint* point = fmscol->points().back();
    point->setDetectorId(muPoint->detectorId());
    point->setEnergy(muPoint->energy());
    point->setX(muPoint->x());
    point->setY(muPoint->y());
  }  // while
}

void StMuFmsUtil::setFmsPointParentClusters(StFmsCollection* fmscol,
                                            StMuFmsCollection* muFms) {
  // Points and clusters in the StMuFmsCollection and StFmsCollection are in
  // the same order, so we get the corresponding objects just by index
  for (unsigned i(0); i < muFms->numberOfClusters(); ++i) {
    const StMuFmsPoint* muPoint = muFms->getPoint(i);
    if (!muPoint) {
      continue;
    }  // if
    const int index = muFms->getClusterArray()->IndexOf(muPoint->cluster());
    if (index != -1) {
      StFmsPoint* point = fmscol->points().at(i);
      if (point) {
        point->setCluster(fmscol->clusters().at(index));
      }  // if
    }  // if
  }  // for
}

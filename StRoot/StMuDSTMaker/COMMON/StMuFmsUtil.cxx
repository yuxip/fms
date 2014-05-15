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
}

void StMuFmsUtil::fillFms(StFmsCollection* fmscol,StMuFmsCollection* muFms)
{
  if(!muFms) return;
  if(!fmscol) return;
  fillFmsHits(fmscol, muFms);
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
  for (int i(0); i < fmscol->numberOfClusters(); ++i) {
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
  for (int i(0); i < fmscol->numberOfPoints(); ++i) {
    const StFmsPoint* point = fmscol->points()[i];
    muFms->addPoint();
    StMuFmsPoint* muPoint = muFms->getPoint(i);
    muPoint->setDetectorId(point->detectorId());
    muPoint->setEnergy(point->energy());
    muPoint->setX(point->x());
    muPoint->setY(point->y());
  }  // for
}

void StMuFmsUtil::fillFmsHits(StFmsCollection* fmscol,
                              StMuFmsCollection* muFms) {
  TClonesArray* arrHit = muFms->getHitArray();
  for(unsigned int i=0; i<muFms->numberOfHits(); i++){
    unsigned short detId = ((StMuFmsHit*)(arrHit->At(i)))->detectorId();
    unsigned short ch    = ((StMuFmsHit*)(arrHit->At(i)))->channel();
    unsigned short crate = ((StMuFmsHit*)(arrHit->At(i)))->qtCrate();
    unsigned short slot  = ((StMuFmsHit*)(arrHit->At(i)))->qtSlot();
    unsigned short qtch  = ((StMuFmsHit*)(arrHit->At(i)))->qtChannel();
    unsigned short adc   = ((StMuFmsHit*)(arrHit->At(i)))->adc();
    unsigned short tdc   = ((StMuFmsHit*)(arrHit->At(i)))->tdc();
    float          ene   = ((StMuFmsHit*)(arrHit->At(i)))->energy();
    StFmsHit* hit = new StFmsHit();
    hit->setDetectorId(detId);
    hit->setChannel(ch);
    hit->setQtCrate(crate);
    hit->setQtSlot(slot);
    hit->setQtChannel(qtch);
    hit->setAdc(adc);
    hit->setTdc(tdc);
    hit->setEnergy(ene);
    fmscol->addHit(hit);
  }  // for
}

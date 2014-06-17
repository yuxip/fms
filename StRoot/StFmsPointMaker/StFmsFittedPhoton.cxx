// $Id$
//
// $Log$
/**
 \file      StFmsFittedPhoton.cxx
 \brief     Implementation of StFmsFittedPhoton,
            a photon fitted to an FMS cluster
 \author    Steven Heppelmann <steveheppelmann@gmail.com>
 \author    Yuxi Pan <yuxipan@physics.ucla.edu>
 \author    Thomas Burton <tpb@bnl.gov>
 \date      2014
 \copyright Brookhaven National Lab
 */
#include "StFmsPointMaker/StFmsFittedPhoton.h"

namespace FMSCluster {
StFmsFittedPhoton::StFmsFittedPhoton(const Float_t x, const Float_t y,
                                     const Float_t e, const Float_t xerr,
                                     const Float_t yerr, const Float_t eerr)
    : energy(e), errEne(eerr), xPos(x), errXPos(xerr), yPos(y), errYPos(yerr) {
}

void StFmsFittedPhoton::Clear() {
  energy = 0;
  errEne = xPos = errXPos = yPos = errYPos = -1;
}
}  // namespace FMSCluster

/****************************************************************************
 *
 * $Id$
 *
 * Author: Thomas Burton, Yuxi Pan, 2014
 ****************************************************************************
 *
 * Description: Implementation of StFmsPoint, the StEvent FMS photon structure
 *
 ****************************************************************************
 *
 * $Log$
 *
 ***************************************************************************/
#include "StFmsPoint.h"

StFmsPoint::StFmsPoint()
    : mDetectorId(0), mEnergy(-1.0), mX(-99.0), mY(-99.0), mId(-1) { }

StFmsPoint::~StFmsPoint() { }

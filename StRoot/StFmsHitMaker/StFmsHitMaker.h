#ifndef ST_FMS_HIT_MAKER_H
#define ST_FMS_HIT_MAKER_H

/*!
 *                                                                     
 * \class  StFmsHitMaker
 * \author Jingguo Ma
 * \date   2009/12/12
 * \brief  StFmsHitMaker
 *
 *
 * This maker makes Fms Hit, Cluster and photon lists and fill them in to
 * StFmsCollection in StEvent
 *
 *
 */                                                                      

//modified by Yuxi Pan 03/28/2013
class StFmsDbMaker;
class StFmsCollection;
class StMuFmsCollection;

/*class FilesSet;
class Qt;
class RunDepMgr;
class RunDepCor;
class CellTDep;
*/

#include "TMatrix.h"
#include "StMaker.h"
using namespace std;
namespace PSUGlobals { } using namespace PSUGlobals;

class StFmsHitMaker : public StMaker {
public:
  StFmsHitMaker(const char* name = "StFmsHitMaker");
  ~StFmsHitMaker();

  void   Clear(Option_t* option = "");
  Int_t  Init();
  Int_t  InitRun(Int_t runNumber);    //called by StMaker when switch to a new run#
  Int_t  Make();
  Int_t  Finish();

  TMatrix**	GetEnergyMatrices();
  Bool_t	Legal(Int_t iew,Int_t nstb,Int_t row0,Int_t col0);
  virtual const char *GetCVS() const
  {static const char cvs[]="Tag $Name:  $ $Id: StFmsHitMaker.h,v 1.1 2010/02/02 21:29:13 jgma Exp $ built "__DATE__" "__TIME__ ; return cvs;}

private:
  StFmsDbMaker*      mFmsDbMaker;    //! DB maker provides FMS geometry and calibration data
  StFmsCollection*   mFmsCollection; //! FMS data structure for StEvent
  StMuFmsCollection* mMuFmsColl;     //! FMS data structure for StMuEvent

  Int_t mCurrentRunNumber;
  ClassDef(StFmsHitMaker,1);
};

#endif

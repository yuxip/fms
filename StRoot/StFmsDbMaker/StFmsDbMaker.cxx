/***************************************************************************
 * $Id: StFmsDbMaker.cxx,v 1.3 2011/01/13 02:56:34 jgma Exp $
 * \author: akio ogawa
 ***************************************************************************
 *
 * Description: This maker is the interface between FMS and the STAR database
 *
 ***************************************************************************
 *
 * $Log: StFmsDbMaker.cxx,v $
 * Revision 1.3  2011/01/13 02:56:34  jgma
 * Fixed bug in function nRow and nColumn
 *
 * Revision 1.2  2010/01/11 20:35:30  jgma
 * Added reversed map and some other minor updates
 *
 * Revision 1.1  2009/10/28 16:11:15  jgma
 * This is the first check in of the code.
 *
 **************************************************************************/


#include "StFmsDbMaker.h"
#include "St_db_Maker/St_db_Maker.h"
#include "StMessMgr.h"
#include "tables/St_fmsDetectorPosition_Table.h"
#include "tables/St_fmsChannelGeometry_Table.h"
#include "tables/St_fmsMap_Table.h"
#include "tables/St_fmsPatchPanelMap_Table.h"
#include "tables/St_fmsQTMap_Table.h"
#include "tables/St_fmsGain_Table.h"
#include "tables/St_fmsGainCorrection_Table.h"
#include "tables/St_fmsRec_Table.h"
StFmsDbMaker* gStFmsDbMaker=NULL; 

ClassImp(StFmsDbMaker)

StFmsDbMaker::StFmsDbMaker(const Char_t *name) : StMaker(name), mDebug(0),mChannelGeometry(0),mDetectorPosition(0),mMap(0),mmMap(0),mPatchPanelMap(0),mQTMap(0),mGain(0),mmGain(0),mGainCorrection(0),mmGainCorrection(0),mRecPar(0),mRecConfig(StFmsDbConfig::Instance()){
	gStFmsDbMaker = this;
//	mRecConfig = StFmsDbConfig::Instance();
}
StFmsDbMaker::~StFmsDbMaker() {deleteArrays(); gStFmsDbMaker = 0;}
Int_t StFmsDbMaker::Init(){
	LOG_DEBUG<<"StFmsDbMaker Init Start"<<endm; 
	return StMaker::Init();
}
Int_t StFmsDbMaker::Make(){LOG_DEBUG<<"StFmsDbMaker Make"<<endm; return kStOK;}
void StFmsDbMaker::Clear(const Char_t*){LOG_DEBUG<<"StFmsDbMaker Clear"<<endm; StMaker::Clear();}
Int_t StFmsDbMaker::Finish(){LOG_DEBUG<<"StFmsDbMaker Finish"<<endm; return kStOK;}

Int_t StFmsDbMaker::InitRun(Int_t runNumber) {
  LOG_DEBUG << "StFmsDbMaker::InitRun - run = " << runNumber << endm;
  deleteArrays();

  //! Accessing DBs
//  if(mDebug>0) {
    St_db_Maker* dbmaker = (St_db_Maker*)GetMaker("db");
    LOG_INFO << "StFmsDbMaker::InitRun - Date&time from St_db_Maker="<<dbmaker->GetDate()<<","<< dbmaker->GetTime() << endm;
//    }

  TDataSet *DBgeom = 0;
  TDataSet *DBmapping = 0;
  TDataSet *DBcalibration = 0;
  DBgeom  = GetInputDB("Geometry/fms");
  DBmapping = GetInputDB("Calibrations/fms/mapping");
  DBcalibration= GetInputDB("Calibrations/fms");
  if(!DBgeom)    {LOG_ERROR << "StFmsDbMaker::InitRun - No Geometry/fms"<<endm;            return kStFatal;}
  if(!DBmapping) {LOG_ERROR << "StFmsDbMaker::InitRun - No Calibration/fms/mapping"<<endm; return kStFatal;} 
  if(!DBcalibration) {LOG_ERROR << "StFmsDbMaker::InitRun - No Calibration/fms"<<endm; return kStFatal;}

  //!Getting DB tables
  St_fmsChannelGeometry *dbChannelGeometry = 0;
  St_fmsDetectorPosition *dbDetectorPosition = 0;
  St_fmsMap             *dbMap             = 0;
  St_fmsPatchPanelMap   *dbPatchPanelMap   = 0;
  St_fmsQTMap           *dbQTMap           = 0;
  St_fmsGain            *dbGain            = 0;
  St_fmsGainCorrection  *dbGainCorrection  = 0;
  St_fmsRec		*dbRec		   = 0;
  dbChannelGeometry = (St_fmsChannelGeometry*) DBgeom->Find("fmsChannelGeometry");
  dbDetectorPosition = (St_fmsDetectorPosition*) DBgeom->Find("fmsDetectorPosition");
  dbMap             = (St_fmsMap*)             DBmapping->Find("fmsMap");
  dbPatchPanelMap   = (St_fmsPatchPanelMap*)   DBmapping->Find("fmsPatchPanelMap");
  dbQTMap           = (St_fmsQTMap*)           DBmapping->Find("fmsQTMap");
  dbGain            = (St_fmsGain*)            DBcalibration->Find("fmsGain");
  dbGainCorrection  = (St_fmsGainCorrection*)  DBcalibration->Find("fmsGainCorrection");
  dbRec		    = (St_fmsRec*)	       DBcalibration->Find("fmsRec");
  if(!dbChannelGeometry){LOG_ERROR << "StFmsDbMaker::InitRun - No Geometry/fms/fmsChannelGeometry"         <<endm; return kStFatal;}
  if(!dbDetectorPosition){LOG_ERROR << "StFmsDbMaker::InitRun - No Geometry/fms/fmsDetectorPosition"         <<endm; return kStFatal;}
  if(!dbMap)            {LOG_ERROR << "StFmsDbMaker::InitRun - No Calibration/fms/mapping/fmsMap"          <<endm; return kStFatal;}
  if(!dbPatchPanelMap)  {LOG_ERROR << "StFmsDbMaker::InitRun - No Calibration/fms/mapping/fmsPatchPanelMap"<<endm; return kStFatal;}
  if(!dbQTMap)          {LOG_ERROR << "StFmsDbMaker::InitRun - No Calibration/fms/mapping/fmsQTMap"        <<endm; return kStFatal;}
  if(!dbGain)           {LOG_ERROR << "StFmsDbMaker::InitRun - No Calibration/fms/fmsGain"                 <<endm; return kStFatal;}
  if(!dbGainCorrection) {LOG_ERROR << "StFmsDbMaker::InitRun - No Calibration/fms/fmsGainCorrection"       <<endm; return kStFatal;}
  if(!dbRec) 		{LOG_ERROR << "StFmsDbMaker::InitRun - No Calibration/fms/fmsRec"	           <<endm; return kStFatal;}

  //!fmsChannelGeometry
  fmsChannelGeometry_st *tChannelGeometry = 0;
  tChannelGeometry = (fmsChannelGeometry_st*) dbChannelGeometry->GetTable();
  Int_t max = dbChannelGeometry->GetNRows();
  mMaxDetectorId = 0;
  for(Int_t i=0; i<max; i++){
    if(mMaxDetectorId < tChannelGeometry[i].detectorId) mMaxDetectorId = tChannelGeometry[i].detectorId;     
  }
  mChannelGeometry = new fmsChannelGeometry_st[mMaxDetectorId+1];
  memset(mChannelGeometry,0,sizeof(fmsChannelGeometry_st)*(mMaxDetectorId+1));
  for(Int_t i=0; i<max; i++){ 
    memcpy(&mChannelGeometry[tChannelGeometry[i].detectorId], &tChannelGeometry[i], sizeof(fmsChannelGeometry_st));
  }
  LOG_DEBUG << "StFmsDbMaker::InitRun - Got Geometry/fms/fmsChannelGeometry with maxDetectorId = "<<mMaxDetectorId<< endm;

  //!fmsDetectorPosition
  fmsDetectorPosition_st *tDetectorPosition = 0;
  tDetectorPosition = (fmsDetectorPosition_st*) dbDetectorPosition->GetTable();
  mDetectorPosition = new fmsDetectorPosition_st[mMaxDetectorId+1];
  memset(mDetectorPosition,0,sizeof(fmsDetectorPosition_st)*(mMaxDetectorId+1));
  max = dbDetectorPosition->GetNRows();
  for(Int_t i=0; i<max; i++){
    memcpy(&mDetectorPosition[tDetectorPosition[i].detectorId], &tDetectorPosition[i], sizeof(fmsDetectorPosition_st));
  }
  LOG_DEBUG << "StFmsDbMaker::InitRun - Got Geometry/fms/fmsDetectorPosition with  "<<max<<" detectors"<< endm;

  //!fmsPatchPanelMap
  mPatchPanelMap = (fmsPatchPanelMap_st*) dbPatchPanelMap->GetTable();
  mMaxModule = dbPatchPanelMap->GetNRows();
  //!LOG_INFO << mMaxModule << "***********" << endm;
  LOG_DEBUG << "StFmsDbMaker::InitRun - Got Calibration/fms/mapping/fmsPatchPanelMap with mMaxModule = "<<mMaxModule<< endm;

  //!fmsQTMap
  mQTMap = (fmsQTMap_st*) dbQTMap->GetTable();
  mMaxNS = dbQTMap->GetNRows();
  LOG_DEBUG << "StFmsDbMaker::InitRun - Got Calibration/fms/mapping/fmsQTMap with mMaxNS = "<<mMaxNS<< endm;
  
  //!fmsMap
  mMap = (fmsMap_st*) dbMap->GetTable();
  mMaxMap = dbMap->GetNRows();
  mmMap = new fmsMap_st* [mMaxDetectorId+1];
  memset(mmMap,0,sizeof(fmsMap_st*)*(mMaxDetectorId+1));
  memset(mReverseMapDetectorId,0,sizeof(mReverseMapDetectorId));
  memset(mReverseMapChannel,0,sizeof(mReverseMapChannel));
  for(Int_t i=0; i<mMaxMap; i++){
    Int_t d=mMap[i].detectorId;
    Int_t c=mMap[i].ch;
    if(d<0 || d>mMaxDetectorId){
      LOG_ERROR << "StFmsDbMaker::InitRun - Calibration/fms/mapping/fmsMap detectorId="<<d<<" exceed max="<<mMaxDetectorId<<endm; 
      return kStFatal;
    }
    if(c<1 || c>maxChannel(d)){
      LOG_ERROR << "StFmsDbMaker::InitRun - Calibration/fms/mapping/fmsMap ch="<<c<<" exceed max="<<maxChannel(d)<<endm; 
      return kStFatal;
    }
    if(mmMap[d]==0){
      mmMap[d] = new fmsMap_st [maxChannel(d)];
      memset(mmMap[d],0,sizeof(fmsMap_st)*maxChannel(d));
    }
    memcpy(&mmMap[d][c-1],&mMap[i],sizeof(fmsMap_st));
    //creating reverse mapping
    Int_t crt,slot,ch;
    getMap(d,c,&crt,&slot,&ch);
    mReverseMapDetectorId[crt][slot][ch]=d;
    mReverseMapChannel[crt][slot][ch]=c;
  }
  
  LOG_DEBUG << "StFmsDbMaker::InitRun - Got Geometry/fms/mapping/fmsMap with mMaxMap = "<<mMaxMap<< endm;
  
  //!fmsGain
  mGain = (fmsGain_st*) dbGain->GetTable();
  mMaxGain = dbGain->GetNRows();
  mmGain = new fmsGain_st* [mMaxDetectorId+1];
  memset(mmGain,0,sizeof(fmsGain_st*)*(mMaxDetectorId+1));
  for(Int_t i=0; i<mMaxGain; i++){
    Int_t d=mGain[i].detectorId;
    Int_t c=mGain[i].ch;
    if(maxChannel(d)<1){
      LOG_ERROR << "StFmsDbMaker::InitRun - Calibration/fms/fmsGain invalid max number of channel = "<<maxChannel(d)<<endm; 
      continue;
    }
    if(d<0 || d>mMaxDetectorId){
      LOG_ERROR << "StFmsDbMaker::InitRun - Calibration/fms/fmsGain detectorId="<<d<<" exceed max = "<<mMaxDetectorId<<endm; 
      continue;
    }
    if(c<1 || c>maxChannel(d)){
      LOG_ERROR << "StFmsDbMaker::InitRun - Calibration/fms/fmsGain detectorId="<<d<<" ch="<<c<<" exceed max = "<<maxChannel(d)<<endm; 
      continue;
    }
    if(mmGain[d]==0){
      mmGain[d] = new fmsGain_st [maxChannel(d)];
      memset(mmGain[d],0,sizeof(fmsGain_st)*maxChannel(d));
    }
    memcpy(&mmGain[d][c-1],&mGain[i],sizeof(fmsGain_st));
  }
  LOG_DEBUG << "StFmsDbMaker::InitRun - Got Calibration/fms/fmsGain with mMaxGain = "<<mMaxGain<< endm;
  
  //!fmsGainCorrection
  mGainCorrection = (fmsGainCorrection_st*) dbGainCorrection->GetTable();
  mMaxGainCorrection = dbGainCorrection->GetNRows();
  mmGainCorrection = new fmsGainCorrection_st* [mMaxDetectorId+1];
  memset(mmGainCorrection,0,sizeof(fmsGainCorrection_st*)*(mMaxDetectorId+1));
  for(Int_t i=0; i<mMaxGainCorrection; i++){
    Int_t d=mGainCorrection[i].detectorId;
    Int_t c=mGainCorrection[i].ch;
    if(maxChannel(d)<1){
      LOG_ERROR << "StFmsDbMaker::InitRun - Calibration/fms/fmsGainCorrection invalid max number of channel = "<<maxChannel(d)<<endm; 
      continue;
    }
    if(d<0 || d>mMaxDetectorId){
      LOG_ERROR << "StFmsDbMaker::InitRun - Calibration/fms/fmsGainCorrection detectorId="<<d<<" exceed max="<<mMaxDetectorId<<endm; 
      continue;
    }
    if(c<1 || c>maxChannel(d)){
      LOG_ERROR << "StFmsDbMaker::InitRun - Calibration/fms/fmsGainCorrection ch="<<c<<" exceed max="<<maxChannel(d)<<endm; 
      continue;
    }
    if(mmGainCorrection[d]==0){
      mmGainCorrection[d] = new fmsGainCorrection_st [maxChannel(d)];
      memset(mmGainCorrection[d],0,sizeof(fmsGainCorrection_st)*maxChannel(d));
    }
    memcpy(&mmGainCorrection[d][c-1],&mGainCorrection[i],sizeof(fmsGainCorrection_st));
  }
  LOG_DEBUG << "StFmsDbMaker::InitRun - Got Geometry/fms/fmsGainCorrection with mMaxGainCorrection = "<<mMaxGainCorrection<< endm;
  
  //!fmsRec
  mMaxRecPar = 80; //dummy
  mRecPar = (fmsRec_st*)dbRec->GetTable();
  mRecConfig.readMap(*mRecPar); //read recPar into internal memory
  LOG_DEBUG << "StFmsDbMaker::InitRun - Got Calibration/fms/fmsRec "<< endm;
  //!Debug
  if(mDebug>0){
    dumpFmsChannelGeometry();
    dumpFmsDetectorPosition();
    dumpFmsMap();
    dumpFmsPatchPanelMap();
    dumpFmsQTMap();
    dumpFmsGain();
    dumpFmsGainCorrection();
    dumpFmsRec();
  }
  return kStOK;
}

StFmsDbConfig& StFmsDbMaker::getRecConfig(){ return mRecConfig; }

void StFmsDbMaker::deleteArrays(){
  if(mChannelGeometry) delete [] mChannelGeometry;  
  if(mDetectorPosition) delete [] mDetectorPosition;  
  if(mmMap){
    for(Int_t d=0; d<=mMaxDetectorId; d++){
      if(mmMap[d]) delete [] mmMap[d];
    }
    delete [] mmMap;
  }
    
  if(mmGain){
    for(Int_t d=0; d<=mMaxDetectorId; d++){ 
      if(mmGain[d]) delete [] mmGain[d];
    }
    delete [] mmGain;
  }
  if(mmGainCorrection){
    for(Int_t d=0; d<=mMaxDetectorId; d++){ 
      if(mmGainCorrection[d]) delete [] mmGainCorrection[d];
    }
    delete [] mmGainCorrection;
  }
}

//! get coordinates in STAR frame, will be implemented when the detector position table is fill into the database and available
StThreeVectorF StFmsDbMaker::getStarXYZ(Int_t detectorId,Float_t FmsX, Float_t FmsY)
{
  Float_t x = 0;
  Float_t y = 0;
  Float_t z = 0;
  y = mDetectorPosition[detectorId].yoffset - FmsY*mDetectorPosition[detectorId].ywidth;
  z = mDetectorPosition[detectorId].zoffset;
  if(northSouth(detectorId) == 0) //! north side
    x = mDetectorPosition[detectorId].xoffset - FmsX*mDetectorPosition[detectorId].xwidth;
  else  //! south side
    x = mDetectorPosition[detectorId].xoffset + FmsX*mDetectorPosition[detectorId].xwidth;
  return StThreeVectorF(x,y,z);
}
Float_t StFmsDbMaker::getPhi(Int_t detectorId,Float_t FmsX, Float_t FmsY){ return (getStarXYZ(detectorId,FmsX,FmsY)).phi();}
Float_t StFmsDbMaker::getEta(Int_t detectorId,Float_t FmsX, Float_t FmsY, Float_t Vertex) { return (getStarXYZ(detectorId,FmsX,FmsY)).pseudoRapidity();}

void StFmsDbMaker::setDebug(Int_t debug) {mDebug = debug;}

//!getting the whole table
fmsDetectorPosition_st* StFmsDbMaker::DetectorPosition() {if(mDetectorPosition) return mDetectorPosition; return 0;}
fmsChannelGeometry_st* StFmsDbMaker::ChannelGeometry() {if(mChannelGeometry) return mChannelGeometry; return 0;}
fmsMap_st*             StFmsDbMaker::Map()             {if(mMap)             return mMap;             return 0;}
fmsPatchPanelMap_st*   StFmsDbMaker::PatchPanelMap()   {if(mPatchPanelMap)   return mPatchPanelMap;   return 0;}
fmsQTMap_st*           StFmsDbMaker::QTMap()           {if(mQTMap)           return mQTMap;           return 0;}
fmsGain_st*            StFmsDbMaker::Gain()            {if(mGain)            return mGain;            return 0;}
fmsGainCorrection_st*  StFmsDbMaker::GainCorrection()  {if(mGainCorrection)  return mGainCorrection;  return 0;}
fmsRec_st*	       StFmsDbMaker::RecPar()	       {if(mRecPar)	     return mRecPar;	      return 0;}

//!ChannelGeometry
Int_t StFmsDbMaker::maxDetectorId()             {return mMaxDetectorId;}
Int_t StFmsDbMaker::eastWest(Int_t detectorId){
  if(detectorId>=0 && detectorId<=mMaxDetectorId && maxChannel(detectorId)>0) return mChannelGeometry[detectorId].ew;
  else{
    LOG_WARN<<"StFmsDbMaker::eastWest: Corresponding channel geometry not found."<<endm;
    return -1;
  }
}

Int_t StFmsDbMaker::northSouth(Int_t detectorId){
  if(detectorId>=0 && detectorId<=mMaxDetectorId && maxChannel(detectorId)>0) return mChannelGeometry[detectorId].ns;
  else{
    LOG_WARN<<"StFmsDbMaker::northSouth: Corresponding channel geometry not found."<<endm;
    return -1;
  }
}

Int_t StFmsDbMaker::type(Int_t detectorId){
  if(detectorId>=0 && detectorId<=mMaxDetectorId && maxChannel(detectorId)>0) return mChannelGeometry[detectorId].type;
  else{
    LOG_WARN<<"StFmsDbMaker::type: Corresponding channel geometry not found."<<endm;
    return -1;
  }
}

Int_t StFmsDbMaker::nRow(Int_t detectorId){
  if(detectorId>=0 && detectorId<=mMaxDetectorId && maxChannel(detectorId)>0) return mChannelGeometry[detectorId].nY;
  else{
    LOG_WARN<<"StFmsDbMaker::nRow: Corresponding channel geometry not found for detectorId = "<<detectorId<<endm;
    return -1;
  }
}

Int_t StFmsDbMaker::nColumn(Int_t detectorId){
  if(detectorId>=0 && detectorId<=mMaxDetectorId && maxChannel(detectorId)>0)
    return mChannelGeometry[detectorId].nX;
  else{
    LOG_WARN<<"StFmsDbMaker::nColumn: Corresponding channel geometry not found."<<endm;
    return -1;
  }
}

Int_t StFmsDbMaker::maxChannel(Int_t detectorId){
  if(detectorId>=0 && detectorId<=mMaxDetectorId && mChannelGeometry[detectorId].nX>0)
    return mChannelGeometry[detectorId].nX*mChannelGeometry[detectorId].nY;
  else{
    LOG_WARN<<"StFmsDbMaker::maxChannel: Corresponding channel geometry not found."<<endm;
    return -1;
  }
}

Int_t StFmsDbMaker::detectorId(Int_t ew, Int_t ns, Int_t type){
  for(Int_t i=0; i<=mMaxDetectorId; i++)
    if((mChannelGeometry+i)){
      if(mChannelGeometry[i].ew   == ew && mChannelGeometry[i].ns   == ns && mChannelGeometry[i].type == type)
	return mChannelGeometry[i].detectorId;
    }
  LOG_WARN<<"StFmsDbMaker::detectorId: Corresponding channel geometry not found."<<endm;
  return -1;
}

Int_t StFmsDbMaker::getRowNumber(Int_t detectorId, Int_t ch){
  if(maxChannel(detectorId)>0) return mChannelGeometry[detectorId].nY - (ch-1)/mChannelGeometry[detectorId].nX;
  return -1;
}

Int_t StFmsDbMaker::getColumnNumber(Int_t detectorId, Int_t ch){
  if(maxChannel(detectorId)>0) return (ch-1)%mChannelGeometry[detectorId].nX + 1;
  return -1;
}

Int_t StFmsDbMaker::getChannelNumber(Int_t detectorId, Int_t row, Int_t column){
  if(maxChannel(detectorId)>0) return column + mChannelGeometry[detectorId].nX * (mChannelGeometry[detectorId].nY - row);
  return -1;
}

StThreeVectorF StFmsDbMaker::getDetectorOffset(Int_t detectorId){
  if(detectorId>=0 && detectorId<=mMaxDetectorId && maxChannel(detectorId)>0)
    return StThreeVectorF(mDetectorPosition[detectorId].xoffset, mDetectorPosition[detectorId].yoffset, mDetectorPosition[detectorId].zoffset);
  return StThreeVectorF(0, 0, 0);
}

Float_t StFmsDbMaker::getXWidth(Int_t detectorId){
  if(detectorId>=0 && detectorId<=mMaxDetectorId)
  return mDetectorPosition[detectorId].xwidth; return -1;
}

Float_t StFmsDbMaker::getYWidth(Int_t detectorId){
  if(detectorId>=0 && detectorId<=mMaxDetectorId)
  return mDetectorPosition[detectorId].ywidth; return -1;
}

//!fmsMap
Int_t StFmsDbMaker::maxMap() {return mMaxMap;}
void StFmsDbMaker::getMap(Int_t detectorId, Int_t ch, Int_t* qtCrate, Int_t* qtSlot, Int_t* qtChannel){
  if(detectorId<0 || detectorId>mMaxDetectorId || ch<1 || ch>maxChannel(detectorId) || mmMap[detectorId]==0 ){
    *qtCrate=0; *qtSlot=0; *qtChannel=0;
    return;
  }
  *qtCrate   = mmMap[detectorId][ch-1].qtCrate;
  *qtSlot    = mmMap[detectorId][ch-1].qtSlot;
  *qtChannel = mmMap[detectorId][ch-1].qtChannel;
}
void StFmsDbMaker::getReverseMap(Int_t qtCrate, Int_t qtSlot, Int_t qtChannel, Int_t* detectorId, Int_t* ch){
  if(qtCrate==0 && qtSlot==0 && qtChannel==0) {
    *detectorId = 0;
    *ch         = 0;
  }else{
    *detectorId = mReverseMapDetectorId[qtCrate][qtSlot][qtChannel];
    *ch         = mReverseMapChannel[qtCrate][qtSlot][qtChannel];
  }
}

//!fmsPatchPanelMap
Int_t StFmsDbMaker::maxModule() {return mMaxModule;}

//!fmsQTMap()
Int_t StFmsDbMaker::maxNS() {return mMaxNS;}

//!fmsGain/GainCorrection
Int_t StFmsDbMaker::maxGain() {return mMaxGain;}
Int_t StFmsDbMaker::maxGainCorrection() {return mMaxGainCorrection;}
Float_t StFmsDbMaker::getGain(Int_t detectorId, Int_t ch){
  if(detectorId<0 || detectorId>mMaxDetectorId || ch<1 || ch>maxChannel(detectorId) || mmGain[detectorId]==0) return 0;
  return mmGain[detectorId][ch-1].gain;
}
Float_t StFmsDbMaker::getGainCorrection(Int_t detectorId, Int_t ch){
  if(detectorId<0 || detectorId>mMaxDetectorId || ch<1 || ch>maxChannel(detectorId) || mmGainCorrection[detectorId]==0) return 0;
  return mmGainCorrection[detectorId][ch-1].corr;
}

//!text dump for debugging
void StFmsDbMaker::dumpFmsChannelGeometry(const Char_t* filename) {
  FILE* fp;
  LOG_INFO << "Writing "<<filename<<endm;
  if((fp=fopen(filename,"w"))){
    fprintf(fp,"maxDetectorId = %d\n",maxDetectorId());
    fprintf(fp,"    i detiid  ew   ns type nRow nCol maxCh\n");
    for(Int_t i=0; i<mMaxDetectorId+1; i++){    
      fprintf(fp,"%5d%7d%4d%5d%5d%5d%5d%6d\n",
	      i,detectorId(eastWest(i),northSouth(i),type(i)),eastWest(i),northSouth(i),type(i),
	      nRow(i),nColumn(i),maxChannel(i));
    }
    for(Int_t i=0; i<mMaxDetectorId+1; i++){    
      fprintf(fp,"DetectorId=%d\n",i);
      fprintf(fp,"detiid  ch   getCh  getRow getCol\n");
      for(Int_t j=1; j<=maxChannel(i); j++){
	fprintf(fp,"%6d%4d%8d%8d%7d\n",
		i,j,getChannelNumber(i,getRowNumber(i,j),getColumnNumber(i,j)),getRowNumber(i,j),getColumnNumber(i,j));
      }	    
    }
    fclose(fp);
  }
}

void StFmsDbMaker::dumpFmsDetectorPosition(const Char_t* filename) {
  FILE* fp;
  LOG_INFO << "Writing "<<filename<<endm;
  if((fp=fopen(filename,"w"))){
    fprintf(fp,"maxDetectorId = %d\n",maxDetectorId());
    fprintf(fp,"  detiid   zoffset   xoffset yoffset    xwidth    ywidth\n");
    for(Int_t i=0; i<mMaxDetectorId+1; i++)
      if((mDetectorPosition+i))
	fprintf(fp,"%8d%10.1f%10.2f%8.1f%10.3f%10.3f\n", i,getDetectorOffset(i).z(),getDetectorOffset(i).x(),getDetectorOffset(i).y(),getXWidth(i),getYWidth(i));
    fclose(fp);
  }
}

void StFmsDbMaker::dumpFmsMap(const Char_t* filename) {
  FILE* fp;
  LOG_INFO << "Writing "<<filename<<endm;
  if((fp=fopen(filename,"w"))){
    fprintf(fp,"maxMap = %d\n",maxMap());
    fprintf(fp,"    i DetId    ch  crt  slt qtch    getmap()   getReverseMap\n");
    for(Int_t i=0; i<mMaxMap; i++){
      Int_t d=mMap[i].detectorId;
      Int_t c=mMap[i].ch;
      Int_t crt,slot,ch,dd,cc;
      getMap(d,c,&crt,&slot,&ch);
      getReverseMap(crt,slot,ch,&dd,&cc);
      fprintf(fp,"%5d%6d%6d%5d%5d%5d%5d%5d%5d%5d%5d\n",
	      i,d,c,mMap[i].qtCrate,mMap[i].qtSlot,mMap[i].qtChannel,crt,slot,ch,dd,cc);
      if(mMap[i].qtCrate>0 && (d-dd!=0 || c-cc!=0)) fprintf(fp,"Problem in reverse map!\n");
    }
    fclose(fp);
  }      
}

void StFmsDbMaker::dumpFmsPatchPanelMap(const Char_t* filename) {
  FILE* fp;
  LOG_INFO << "Writing "<<filename<<endm;
  if((fp=fopen(filename,"w"))){
    fprintf(fp,"  mod channel ppPanel ppRow ppColumn\n");
    for(Int_t i=0; i<mMaxModule; i++)
      for(Int_t j=0; j<maxChannel(i+8); j++)
	fprintf(fp,"%5d%8d%8d%6d%9d\n",i+1,j+1,mPatchPanelMap[i].ppPanel[j],mPatchPanelMap[i].ppRow[j],mPatchPanelMap[i].ppColumn[j]);
    fclose(fp);
  }
}

void StFmsDbMaker::dumpFmsQTMap(const Char_t* filename) {
  FILE* fp;
  LOG_INFO << "Writing "<<filename<<endm;
  if((fp=fopen(filename,"w"))){
    fprintf(fp,"ns ppPanel   row column crate slot channel\n");
    for(Int_t ns=0; ns<2; ns++)
      for(Int_t pp=0; pp<2; pp++)
	for(Int_t row=0; row<20; row++)
	  for(Int_t col=0; col<16; col++){
	    if(mQTMap[ns].qtCrate[pp][row][col]==0 && mQTMap[ns].qtSlot[pp][row][col]==0 && mQTMap[ns].qtChannel[pp][row][col]==0)
	      fprintf(fp,"-1      -1    -1     -1    -1   -1      -1\n");
	    else
	      fprintf(fp,"%2d%8d%6d%7d%6d%5d%8d\n",ns+1, pp+1, row+1, col+1,mQTMap[ns].qtCrate[pp][row][col],mQTMap[ns].qtSlot[pp][row][col],
		      mQTMap[ns].qtChannel[pp][row][col]);
	  }
    fclose(fp);
  }
}

void StFmsDbMaker::dumpFmsGain(const Char_t* filename) {
  FILE* fp;
  LOG_INFO << "Writing "<<filename<<endm;
  if((fp=fopen(filename,"w"))){
    fprintf(fp,"maxGain = %d\n",maxGain());
    fprintf(fp,"    i DetId    ch    gain  getGain()\n");
    for(Int_t i=0; i<mMaxGain; i++){
      Int_t d=mGain[i].detectorId;
      Int_t c=mGain[i].ch;
      fprintf(fp,"%5d%6d%6d%8.3f%11.3f\n",
	      i,d,c,mGain[i].gain,getGain(d,c));
    }
    fclose(fp);
  }      
}

void StFmsDbMaker::dumpFmsGainCorrection(const Char_t* filename) {
  FILE* fp;
  LOG_INFO << "Writing "<<filename<<endm;
  if((fp=fopen(filename,"w"))){
    fprintf(fp,"maxGainCorrection = %d\n",maxGainCorrection());
    fprintf(fp,"    i DetId    ch    gain  getGainCorrection()\n");
    for(Int_t i=0; i<mMaxGainCorrection; i++){
      Int_t d=mGainCorrection[i].detectorId;
      Int_t c=mGainCorrection[i].ch;
      fprintf(fp,"%5d%6d%6d%8.3f%21.3f\n",
	      i,d,c,mGainCorrection[i].corr,getGainCorrection(d,c));
    }
    fclose(fp);
  }      
}

void StFmsDbMaker::dumpFmsRec(const Char_t* filename) {
	
  LOG_INFO << "writing "<<filename<<endm;
  mRecConfig.writeMap(filename);

}

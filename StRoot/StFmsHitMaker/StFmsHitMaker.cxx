// ROOT
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"


#include "StMessMgr.h"
#include "StEventTypes.h"
#include "StMuDSTMaker/COMMON/StMuTypes.hh"

#include "StFmsDbMaker/StFmsDbMaker.h"
#include "StFmsCollection.h"
#include "StFmsHit.h"
#include "StFmsHitMaker.h"
#include "StTriggerData2009.h"

// PSU-FMS package
#include "StPSUTools/Geom.h"
#include "StPSUTools/FitTower.h"

#include <iostream>
#include <fstream>
#include <stdlib.h>
using namespace std;
using namespace PSUGlobals;

ClassImp(StFmsHitMaker);

StFmsHitMaker::StFmsHitMaker(const char* name) : StMaker(name){
	mFmsDbMaker = NULL;
	mFmsCollection = NULL;
	mMuFmsColl = NULL;
	LOG_INFO << "StFmsHitMaker initializing...."<<endm;
	LOG_DEBUG << "StFmsHitMaker::constructor." << endm;

	//asy. shower shape not available in StPSUTools
	Float_t a1[3]={1.070804, 0.167773, -0.238578};//top frac 85
	Float_t b1[3]={0.535845, 0.850233, 2.382637};
	FitTower::Setwe_ab(a1[0],a1[1],a1[2],b1[0],b1[1],b1[2]);
	for(Int_t i = 0; i < 4; i++){
		mAdc[i] = mEnergy[i] = 0;
	}
}

StFmsHitMaker::~StFmsHitMaker(){
  
  LOG_INFO << "StFmsHitMaker::destructor." << endm;
  for(Int_t i = 0; i < 4; i++){
	if(mAdc[i]) { delete mAdc[i]; mAdc[i] = 0; }
	if(mEnergy[i]){ delete mEnergy[i]; mEnergy[i] = 0;}
  }
}

void StFmsHitMaker::Clear(Option_t* option){
  LOG_DEBUG << "StFmsHitMaker::Clear()" << endm;
  for(Int_t i = 0; i < 4; i++){
	if(mAdc[i]) *(mAdc[i]) = 0.;
	if(mEnergy[i]) *(mEnergy[i]) = 0.;
  }
  StMaker::Clear(option);
}

int StFmsHitMaker::Init() {
  
	LOG_INFO<<"StFmsHitMaker::Init() "<<endm;

	return StMaker::Init();

}

int StFmsHitMaker::InitRun(Int_t runNumber) {
  
	LOG_INFO << "StFmsHitMaker::InitRun --run# changed to " << runNumber << endm;
	//initial energy/adc matrices
	for(Int_t j = 0; j < 4; j++){
                Int_t detectorId = j+8;
                LOG_INFO<<"nstb = "<<j+1<<", NROWS = "<<gStFmsDbMaker->nRow(detectorId)<<", NCOLS = "<<gStFmsDbMaker->nColumn(detectorId)<<endm;
		//only allocate new space in the begining, not in between runs
		if(!(mAdc[j]))mAdc[j] = new TMatrix(gStFmsDbMaker->nRow(detectorId),gStFmsDbMaker->nColumn(detectorId));
		if(!(mEnergy[j]))mEnergy[j] = new TMatrix(gStFmsDbMaker->nRow(detectorId),gStFmsDbMaker->nColumn(detectorId));
                *(mAdc[j]) = 0.;
                *(mEnergy[j]) = 0.;

        }
	mFmsDbMaker = gStFmsDbMaker;
	mCurrentRunNumber = runNumber; //called by maker's method :InitRun(run); when the run# changes

	return kStOK;
}

//// This is StFmsHitMaker Make, it reads in data and make hit and fill StFmsCollection
int StFmsHitMaker::Make(){

	LOG_DEBUG<<"StFmsHitMaker::Make start"<<endm;
	int flag = 0;
	StTriggerData* triggerData = 0;

	//first try to get StTriggerData from StTriggerDataMaker (works for proudction) and create StFmsCollection
	TObjectSet *os = (TObjectSet*)GetDataSet("StTriggerData");
	if (os) {
		triggerData = (StTriggerData*)os->GetObject();
		if(triggerData){
			flag=1;
			// mFmsCollection = new StFmsCollection();
			LOG_DEBUG<<"StFmsHitMaker::Make Found StTriggerData from StTriggerDataMaker"<<endm;
		}
	}
  
	//2nd try to get StTriggerData from StEvent
	//but once FMS data is killed in StEvent, this will not work and all you see is empty data
	StEvent* stEvent = (StEvent*) GetInputDS("StEvent");
	//StEvent* stEvent = (StEvent*) GetInputDS("StEvent");
	if(flag==0){
		if(stEvent){
			triggerData = stEvent->triggerData();
			if(triggerData) {
				flag=2;
				LOG_DEBUG<<"StFmsHitMaker::Make Found StTriggerData from StEvent"<<endm;
			}
			else{
				mFmsCollection = stEvent->fmsCollection();
				if(mFmsCollection){
	  				flag=3;
	  				LOG_DEBUG<<"StFmsHitMaker::Make Found StFmsCollection from StEvent"<<endm;
				}
			}
		} //found StEvent
	}
	LOG_DEBUG<<"before checking MuDst, flag is: "<<flag<<endl;
  
	//3rd try to get StTriggerData from StMuEvent, works for produced data (.MuDst.root) --Yuxi
	StMuDst* muDst = (StMuDst*)GetInputDS("MuDst");
	mCurrentRunNumber = muDst->event()->runNumber();
	if(flag==0){
		if(muDst){
			triggerData = (StTriggerData*)StMuDst::event()->triggerData();
			if(triggerData){
				flag = 4; //Yuxi
				LOG_INFO<<"StFmsHitMaker::Make Found StFmsTriggerData in MuDst"<<endm;
			}
			else LOG_ERROR << "Finally, no StFmsTriggerData in MuDst " <<endm;
		}
	}
   
 
	LOG_DEBUG<<"after checking MuDst, flag is: "<<flag<<endm;
	//after this step triggerData is pointing to StTriggerData block of StEvent

	if(flag>0){
	
		mFmsCollection = new StFmsCollection();
		//create StFmsHit and add it to StFmsCollection
		for(unsigned short crt=1; crt<=4; crt++){
			for(unsigned short slot=1; slot<=16; slot++){
				for(unsigned short ch=0; ch<32; ch++){
	  				unsigned short adc=0;
	  				unsigned short tdc=0;
	  				if(flag<=4){ //wont work when flag=3
	    					adc=triggerData->fmsADC(crt,slot-1,ch);
	    					tdc=triggerData->fmsTDC(crt,slot-1,ch);
	  				}
				
					if(adc>0 || tdc>0){
					//	LOG_INFO<<"adc of crt "<<crt<<", slot "<<slot<<", channel "<<ch<<" is: "<<adc<<endm;
					//	LOG_INFO<<"tdc=====================================================is: "<<tdc<<endm;
	    					StFmsHit* hit = new StFmsHit();
	    					if(!hit){
	      						LOG_ERROR <<"Failed to create FMS hit, skip this hit."<<endm;
	      						continue;
	    					}
						hit->setDetectorId(0);
						hit->setChannel(0);
						hit->setQtCrate(crt);
						hit->setQtSlot(slot);
						hit->setQtChannel(ch);
						hit->setAdc(adc);
						hit->setTdc(tdc);
						hit->setEnergy(0.0);
						mFmsCollection->addHit(hit);
						if(GetDebug()>0) hit->print();
	  				}
				}
      			}	
    		}
    

  		// Read DB and put DetectorId, channeel and apply Calibration to get Energy
  		LOG_DEBUG<<"2nd pass"<<endm;
    		for(unsigned int i=0; i<mFmsCollection->numberOfHits(); i++){
      			int d,c; //detector Id, channel# (starts from 1)
      			StFmsHit* fmsHit = (mFmsCollection->hits())[i];
      			int crt   =fmsHit->qtCrate();
      			int slot  =fmsHit->qtSlot();
      			int ch    =fmsHit->qtChannel();
      			float adc =fmsHit->adc();
      			mFmsDbMaker->getReverseMap(crt,slot,ch,&d,&c);
		//	int ns = mFmsDbMaker->northSouth(d);
		//	int type = mFmsDbMaker->type(d);
			int row = mFmsDbMaker->getRowNumber(d,c);
			int col = mFmsDbMaker->getColumnNumber(d,c);
			int nstb = 0;
			int ew  = 2; //east=1, west=2
			switch(d){
				case 9: //south large
					nstb = 2;
					break;
				case 8: //north large
					nstb = 1;
					break;
				case 11: //south small
					nstb = 4;
					break;
				case 10: //north small
					nstb = 3;
					break;
				default:
					nstb = 0;
					break;
			}			

      			float e=0.0;
      			if(d>0 || c>0){
				float g1 = mFmsDbMaker->getGain(d,c);
				float g2 = mFmsDbMaker->getGainCorrection(d,c);
				e = adc*g1*g2;
			}
      		
      			fmsHit->setDetectorId(d);
      			fmsHit->setChannel(c);
      			fmsHit->setEnergy(e);
		//	fmsHit->setRow(row); //seems impossible to change StFmsHit definition in StEvent
		//	fmsHit->setCol(col);
		//	fmsHit->setNs(ns);
		//	fmsHit->setType(type);
			if(nstb==1||nstb==2){
				row = 35 - row; //because channel geometry in the database assigns row1 as the bottom row
			}
			if(nstb==3||nstb==4){
				row = 25 - row;
			}			
			if(!Legal(ew,nstb,row-1,col-1)){
		//		LOG_INFO<<" illegal (nstb,row,col), dId = "<<dId<<", crt = "<<crt<<", slot = "<<slot<<", ch = "<<ch<<", nstb1 = "<<nstb<<", row0 = "<<row-1<<", col0 = "<<col-1<<endm;
				continue;
			}
		//	populate adc and energy matrices
			if(adc>0){
				(*(mAdc[nstb-1]))[row-1][col-1] = adc;
				(*(mEnergy[nstb-1]))[row-1][col-1] = e;
		//	LOG_INFO<<"adc of dId "<<dId<<", crt "<<crt<<", slot "<<slot<<", ch "<<ch<<", nstb1 "<<nstb<<", row0 "<<row-1<<", col0 "<<col-1<<" is: "<<adc<<endm;
			}
			

			if(GetDebug()>0) fmsHit->print();
		}

		/*
		cout<<"debug adc, energy arrays"<<endl;
		for(Int_t j = 0; j < 4; j++){
			cout<<"adc_"<<j<<endl;
			mAdc[j]->Print();
			cout<<"energy_"<<j<<endl;
			mEnergy[j]->Print();
		}*/
  
  		LOG_INFO<<"StFmsHitMaker::Make got "<<mFmsCollection->numberOfHits()<<" hits in StFmsCollection"<<endm;
	

						

	}//flag, received fms data from mudst, triggerdata, etc.

 
/*   
	//flag = 0; //yuxi debug
	if(flag==0) { //read hits created during production
		if(muDst){
			mMuFmsColl = StMuDst::muFmsCollection();
			if(mMuFmsColl){
				int nhits = mMuFmsColl->numberOfHits();
				LOG_DEBUG<<"StFmsHitMaker::Make Found "<<nhits<<" hits in MuDst"<<endm;
				cout<<"StFmsHitMaker::Make Found "<<nhits<<" hits in MuDst"<<endl;
				for(int i=0; i<nhits; i++){
		  			StMuFmsHit* mHit = mMuFmsColl->getHit(i);
		  			if(mHit && GetDebug()>0) mHit->print();
		  			mHit->print();
				}
				flag=5;
			}
		}
	}
*/
	if(flag>0 && flag<=2) {
		//Adding StFmsCollection to StEvent
		LOG_DEBUG<<"StFmsHitMaker::Make Adding StFmsCollection as FMSCOLLECTION"<<endm;
		stEvent->setFmsCollection(mFmsCollection);
	}

	return kStOk;
}

int StFmsHitMaker::Finish(){
  
  LOG_INFO << "StFmsHitMaker::Finish() " << endm;
  return kStOk;
}

TMatrix** StFmsHitMaker::GetEnergyMatrices() {
	
	return mEnergy;
}

Bool_t StFmsHitMaker::Legal(Int_t iew,Int_t nstb,Int_t row0,Int_t col0){
	//nstb starts from 1
	//row0,col0 starts from 0	

	if(iew>0 && iew<2)return false;
	if(nstb<1 || nstb>4)return false;
	if(nstb>2){
      		if(row0<0 || row0>23)return false;
		if(col0<0 || col0>11)return false;
		if(fabs(1.*row0-11.5)<5 && col0<5)return false;
	}
	else{
		if(row0<0 || row0>33)return false;
		if(col0<0 || col0>16)return false;
		if(fabs(1.*row0-16.5)<8 && col0<8)return false;
		if(row0<col0-9.5)return false;
		if(33-row0<col0-9.5)return false;
	}	
  return true;
}


#include "StFmsPointMaker.h"

#include "StEvent/StEvent.h"
#include "StEvent/StFmsCollection.h"
#include "StEvent/StFmsHit.h"
#include "StEvent/StTriggerData.h"
#include "StFmsDbMaker/StFmsDbMaker.h"
#include "StFmsHitMaker/StFmsHitMaker.h"
#include "StFmsClusterCollection.h"
#include "StFmsPointCollection.h"
#include "StFmsClHitCollection.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"

#include "StMessMgr.h"

#include "TFile.h"
#include "TMatrix.h"
#include "TTree.h"
#include "TMinuit.h"

#include "StPSUTools/Yiqun.h"
#include "StPSUTools/TowerFPD.h"
#include "StPSUTools/Geom.h"
#include "StPSUTools/HitCluster.h"
#include "StPSUTools/TowerUtil.h"  // Defines ClusterList

#ifndef __CINT__
typedef ClusterList::const_iterator ClusterCIter;
#endif  // __CINT__

using namespace std;
using namespace PSUGlobals;

ClassImp(StFmsPointMaker)

StFmsPointMaker::StFmsPointMaker(const char* name = "StFmsPointMaker")
: StMaker(name) { mFmsClColl = 0; }

StFmsPointMaker::~StFmsPointMaker() {
	
	LOG_INFO << "StFmsPointMaker:: destructor " << endm;

//	if(mFmsClColl){
	//	mFmsClColl->Clear();
	//	delete mFmsClColl;
	//	mFmsClColl = 0;
//	}

/*	if(mFmsPtsColl){
		delete mFmsPtsColl;
		mFmsPtsColl = 0;
	}
*/

}

void StFmsPointMaker::Clear( const char* opt ) {

	LOG_INFO << "StFmsPointMaker::Clear() " << endm;
	if(mFmsClColl){
                mFmsClColl->Clear();
		delete mFmsClColl;
		mFmsClColl = 0;
        }

/*	if(mFmsPtsColl){
		mFmsPtsColl->Clear();
		delete mFmsPtsColl;
		mFmsPtsColl = 0;
        }
*/
	LOG_INFO <<"after StFmsPointMaker::Clear()" <<endm;
	StMaker::Clear();	

}

Int_t StFmsPointMaker::Init() {
	
	LOG_INFO << "StFmsPointMaker::Init() " << endm;
	//mFmsClColl  = new StFmsClusterCollection();		
	//Get run# from StFmsHitMaker and access status files
	fmsgeom = 0;
	return StMaker::Init();
}

Int_t StFmsPointMaker::InitRun(Int_t runNumber){ //gStFmsDbMaker is filled after InitRun
	//setting up Geom, CalibStr, there are all internal functions that help the default cluster finder origanize detector
	//geometries and calibration tables, which would stay constant for each Run
	//only allocate new space in the begining, not in between runs
	if(!fmsgeom)fmsgeom = new Geom();
	// Ensure we can access database information
	mFmsDbMaker = static_cast<StFmsDbMaker*>(GetMaker("fmsDb"));
	if (!mFmsDbMaker) {
	  return kStErr;
	}  // if
}

Int_t StFmsPointMaker::Finish() {
	
	LOG_INFO << "StFmsPointMaker::Finish() " << endm;
	return kStOk;
}

Int_t StFmsPointMaker::Make() {
	
	LOG_INFO << "StFmsPointMaker::Make() " << endm;
	mFmsClColl  = new StFmsClusterCollection();
//	mFmsPtsColl = new StFmsPointCollection();
  if (!populateTowerLists()) {
    LOG_ERROR << "StFmsPointMaker::Make() - failed to initialise tower " <<
      "lists for the event" << endm;
  }  // if
  if(FindPoint()==kStOk){
		 LOG_INFO << "Cluster finder returns successfully" <<endm;
		 return kStOk;
	}
	LOG_INFO << " cluster finder returns error!!!" << endm;
	return kStErr;
}
	
Int_t StFmsPointMaker::FindPoint() {
	
	LOG_INFO << " StFmsPointMaker::FindPoint() " << endm;
	
	for(Int_t instb = 0; instb < 4; instb++){
		TowerList& towers = mTowers.at(instb);
		Float_t Esum = 0.f;
		for (TowerList::const_iterator i = towers.begin(); i != towers.end(); ++i) {
		  Esum += i->hit()->energy();
		}  // for
		if(Esum==0||Esum>500) continue; //to remove LED trails, for pp500 GeV
		Int_t detectorId = instb + 8;  // FMS IDs from 8 to 11
    Yiqun clustering(fmsgeom, detectorId);
    // Perform tower clustering, skip this subdetector if an error occurs
    if (!clustering.cluster(&towers)) {
      continue;
    }  // if
		//Saved cluser info into StFmsCluster
		Int_t iPh = 0;	//sequence # in Yiqun::photons[];
		const ClusterList& clusters = clustering.clusters();
		for (ClusterCIter ci = clusters.begin(); ci != clusters.end(); ++ci) {
			StFmsCluster* cluster = new StFmsCluster();
			
			Int_t cluid = 305 + 20*(instb) + iPh; //cluster id = id of the 1st photon, not necessarily the highE photon
			
			cluster->SetNstb(instb + 1);
			
      ci->copyTo(cluster);
			//calculate cluster four momentum
			Geom* p_geom = fmsgeom;
			if(!p_geom){
				LOG_ERROR << " StFmsPointMaker::Make() Geom missing! " << endm;
				return kStErr;
			}
			Float_t widLG[2];
			widLG[0] = (p_geom->FpdTowWid(clustering.mDetectorId))[0];//lead glass x width
			widLG[1] = (p_geom->FpdTowWid(clustering.mDetectorId))[1];//lead glass y width
			TVector3 xyz;
			xyz[2] = p_geom->z(clustering.mDetectorId);
			xyz[0] = cluster->GetX0()*widLG[0];
			xyz[1] = cluster->GetY0()*widLG[1];
			if (clustering.mDetectorId == 8 ||
			    clustering.mDetectorId == 10) { //north, negative x axis
        xyz[0] = p_geom->xOffset(clustering.mDetectorId) - xyz[0];
			}
			else{
        xyz[0] = p_geom->xOffset(clustering.mDetectorId) + xyz[0]; //south, positive x axis
			}
      xyz[1] = p_geom->yOffset(clustering.mDetectorId) - xyz[1];
			Double_t dist = xyz.Mag();
			TVector3 uvec(0.,0.,0.);
			if(dist!=0)uvec=(1./dist)*xyz;
			Float_t clEnergy = cluster->GetEnergy();
			TVector3 mom3;
			mom3 = clEnergy*uvec;
			TLorentzVector mom4(mom3,clEnergy);
			cluster->SetFourMomentum(mom4);
			
			//save photons reconstructed from this cluster
			for(Int_t np = 0; np < cluster->GetNphoton(); np++){
				
				StFmsPoint* clpoint = new StFmsPoint();
				clpoint->SetEnergy(ci->photons()[np].energy);
				clpoint->SetXpos(ci->photons()[np].xPos);//in cm
				clpoint->SetYpos(ci->photons()[np].yPos);//in cm
				
				Int_t phid = 305 + 20*(instb) + iPh;
				clpoint->SetPhotonId(phid);
				iPh++;
			
				//calculate photon 4 momentum;
				TVector3 xyzph;
				xyzph[2] = p_geom->z(clustering.mDetectorId);
				xyzph[0] = clpoint->GetXpos();	//in cm, towWidth*x0(fit)
				xyzph[1] = clpoint->GetYpos();  
				
        if (clustering.mDetectorId == 8 ||
            clustering.mDetectorId == 10) { //north, negative x axis
          xyzph[0] = p_geom->xOffset(clustering.mDetectorId) - xyzph[0];
        }
				else{
					xyzph[0] = p_geom->xOffset(clustering.mDetectorId) + xyzph[0];
				}
				xyzph[1] = p_geom->yOffset(clustering.mDetectorId) - xyzph[1];
					
				clpoint->SetPointXYZLab(xyzph);
				Double_t distph = xyzph.Mag();
				TVector3 uvecph(0.,0.,0.);
				if(distph!=0)uvecph=(1./distph)*xyzph;
				Float_t phEnergy = clpoint->GetEnergy();
				TVector3 phmom3;
				phmom3 = phEnergy*uvecph;
				TLorentzVector phmom4(phmom3,phEnergy);
				clpoint->SetFourMomentum(phmom4);
				Int_t cluid = cluster->GetClusterId();
				clpoint->SetParentCluId(cluid);
				Int_t nclph = cluster->GetNphoton();
				clpoint->SetParentNclPh(nclph);
				//Add it to StFmsCluster mPhotons
				cluster->GetPointCollection()->AddPoint(clpoint);

			/*	//make a copy of clpoint (use the default copy constructor), add it to mFmsPtsColl;
				StFmsPoint* ptspoint = new StFmsPoint(*clpoint);
				this->mFmsPtsColl->AddPoint(ptspoint);
			*/
			}
			
			//save the tower hit info.
			TIter next(ci->towers());
			while(TowerFPD* tow = (TowerFPD*)next()){
				if (tow->hit()->adc() >= 1) {			//minADC=1
					Int_t snstb = instb + 1; 	//starts from 1
					Int_t srow = (tow->row()) - 1;		//srow starts from 0
					Int_t scol = (tow->column()) - 1;		//scol starts from 0
					UInt_t adc = tow->hit()->adc();
					Float_t tenergy = tow->hit()->energy();
					UChar_t status = 0; 				//no status code for now --04/01/2013					
					StFmsClHit* hit = new StFmsClHit(snstb,srow,scol,adc,tenergy,status);	
					cluster->GetClHitCollection()->AddHit(hit);
				}
			}
		//LOG_INFO << "StFmsPointMaker::FindPoint() --StFmsCluster created: " << endm;
		//cluster->Print();
		this->mFmsClColl->AddCluster(cluster);
		
		}//loop over clusters
		
	}//loop over NSTB
	LOG_INFO << "StFmsPointMaker::FindPoint() --StFmsCluster collections filled " << endm;
	LOG_INFO << "nClusters = " << this->mFmsClColl->NumberOfClusters() << endm;
	//this->mFmsClColl->Print();
	AddData(mFmsClColl);
	LOG_INFO << "StFmsPointMaker::FindPoint() --StFmsClusterCollection added to TDataSet " << endm;
	
//	AddData(mFmsPtsColl);
//	LOG_INFO << "StFmsPointMaker::FindPoint() --StFmsPointCollection added to TDataSet " << endm;

	return kStOk;
}

Bool_t StFmsPointMaker::populateTowerLists() {
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
  const StSPtrVecFmsHit& hits = fmsCollection->hits();
  for (StSPtrVecFmsHitConstIterator i = hits.begin(); i != hits.end(); ++i) {
    const StFmsHit* hit = *i;
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
    if (!Legal(ew, nstb, row - 1, column - 1)) {
      continue;
    }  // if
    unsigned index = hit->detectorId() - 8;  // FMS IDs range from 8 to 11
    if (hit->adc() > 0 && index >= 0 && index < mTowers.size()) {
      PSUGlobals::TowerFPD tower(hit);
      // Ensure tower information is valid before adding
      if (tower.initialize(mFmsDbMaker)) {
        mTowers.at(index).push_back(tower);
      }  // if
    }  // if
  }  // for
  return true;
}

Bool_t StFmsPointMaker::Legal(Int_t iew,Int_t nstb,Int_t row0,Int_t col0){
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

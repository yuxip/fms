
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
using namespace std;
using namespace PSUGlobals;

ClassImp(StFmsPointMaker)

StFmsPointMaker::StFmsPointMaker(const char* name = "StFmsPointMaker", StFmsHitMaker* fmsHitMaker = 0 )
: StMaker(name), mFmsHitMaker(fmsHitMaker) { mFmsClColl = 0; }

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
  std::vector<TMatrix>::iterator iter;
  for (iter = mEnergyMatrices.begin(); iter != mEnergyMatrices.end(); ++iter) {
    *iter = 0.f;
  }  // for
	LOG_INFO <<"after StFmsPointMaker::Clear()" <<endm;
	StMaker::Clear();	

}

Int_t StFmsPointMaker::Init() {
	
	LOG_INFO << "StFmsPointMaker::Init() " << endm;
	if(!mFmsHitMaker){
		LOG_ERROR << " StFmsHitMaker is null " <<endm;
		return kStErr;
	}
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
	// Create energy matrices of the correct (row, column) dimensions for each
	// FMS subdetector, and initialise all elements to zero.
	mEnergyMatrices.assign(4, TMatrix());
	for(Int_t i = 0; i < 4; i++){
    Int_t detectorId = i + 8;  // FMS detector ID in range [8, 11]
    Int_t nRows = mFmsDbMaker->nRow(detectorId);
    Int_t nCols = mFmsDbMaker->nColumn(detectorId);
    mEnergyMatrices.at(i).ResizeTo(nRows, nCols);
    mEnergyMatrices.at(i) = 0.f;
  }  // for
}

Int_t StFmsPointMaker::Finish() {
	
	LOG_INFO << "StFmsPointMaker::Finish() " << endm;
	return kStOk;
}

Int_t StFmsPointMaker::Make() {
	
	LOG_INFO << "StFmsPointMaker::Make() " << endm;
	mFmsClColl  = new StFmsClusterCollection();
//	mFmsPtsColl = new StFmsPointCollection();
	
  initialiseEnergyMatrices();
  if(FindPoint()==kStOk){
		 LOG_INFO << "Cluster finder returns successfully" <<endm;
		 return kStOk;
	}
	LOG_INFO << " cluster finder returns error!!!" << endm;
	return kStErr;
}
	
Int_t StFmsPointMaker::FindPoint() {
	
	LOG_INFO << " StFmsPointMaker::FindPoint() " << endm;
	
	Yiqun* p_rec[4];
	for(Int_t instb = 0; instb < 4; instb++){
		
		Float_t Esum = mEnergyMatrices.at(instb).Sum();
		if(Esum==0||Esum>500) continue; //to remove LED trails, for pp500 GeV

		//call the cluster finder for each nstb
		p_rec[instb] = new Yiqun(&mEnergyMatrices.at(instb),fmsgeom,2,instb+1);
		
		//Saved cluser info into StFmsCluster
		Int_t iPh = 0;	//sequence # in Yiqun::photons[];
		for(Int_t ncl = 0; ncl<p_rec[instb]->NRealClusts; ncl++){

			StFmsCluster* cluster = new StFmsCluster();
			
			Int_t cluid = 305 + 20*(p_rec[instb]->NSTB-1) + iPh; //cluster id = id of the 1st photon, not necessarily the highE photon
			
			cluster->SetNstb(p_rec[instb]->NSTB);
			cluster->SetClusterId(cluid);
			cluster->SetCatag(p_rec[instb]->clust[ncl].catag);
			cluster->SetNumbTower(p_rec[instb]->clust[ncl].numbTower);
			cluster->SetNphoton(p_rec[instb]->clust[ncl].nPhoton);
			cluster->SetClusterEnergy(p_rec[instb]->clust[ncl].energy);
			cluster->SetX0(p_rec[instb]->clust[ncl].x0);	//in units of tower width
			cluster->SetY0(p_rec[instb]->clust[ncl].y0);	//in units of tower width
			cluster->SetSigmaMax(p_rec[instb]->clust[ncl].sigmaMax);
			cluster->SetSigmaMin(p_rec[instb]->clust[ncl].sigmaMin);
		//	cluster->SetChi2NdfPh1(p_rec[instb]->clust[ncl].Chi2NdfPh1);
		//	--no such funtion in SH's package --Yuxi
		//	cluster->SetChi2NdfPh2(p_rec[instb]->clust[ncl].Chi2NdfPh2);

			//calculate cluster four momentum
			Geom* p_geom = fmsgeom;
			if(!p_geom){
				LOG_ERROR << " StFmsPointMaker::Make() Geom missing! " << endm;
				return kStErr;
			}
			Float_t widLG[2];
			widLG[0] = (p_geom->FpdTowWid(2,p_rec[instb]->NSTB))[0];//lead glass x width
			widLG[1] = (p_geom->FpdTowWid(2,p_rec[instb]->NSTB))[1];//lead glass y width
			TVector3 xyz;
			xyz[2] = *(p_geom->ZFPD(2,p_rec[instb]->NSTB));
			xyz[0] = cluster->GetX0()*widLG[0];
			xyz[1] = cluster->GetY0()*widLG[1];
			if((p_rec[instb]->NSTB)==1||(p_rec[instb]->NSTB)==3){ //north, negative x axis
				xyz[0] = (*(p_geom->xOffset(2,p_rec[instb]->NSTB))) - xyz[0];
			}
			else{
				xyz[0] = (*(p_geom->xOffset(2,p_rec[instb]->NSTB))) + xyz[0]; //south, positive x axis
			}
			xyz[1] = (*(p_geom->yOffset(2,p_rec[instb]->NSTB))) - xyz[1];
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
				clpoint->SetEnergy(p_rec[instb]->clust[ncl].photon[np].energy);
				clpoint->SetXpos(p_rec[instb]->clust[ncl].photon[np].xPos);//in cm
				clpoint->SetYpos(p_rec[instb]->clust[ncl].photon[np].yPos);//in cm
				
				Int_t phid = 305 + 20*(p_rec[instb]->NSTB-1) + iPh;
				clpoint->SetPhotonId(phid);
				iPh++;
			
				//calculate photon 4 momentum;
				TVector3 xyzph;
				xyzph[2] = *(p_geom->ZFPD(2,p_rec[instb]->NSTB));
				xyzph[0] = clpoint->GetXpos();	//in cm, towWidth*x0(fit)
				xyzph[1] = clpoint->GetYpos();  
				
				if((p_rec[instb]->NSTB)==1||(p_rec[instb]->NSTB)==3){ //north, negative x axis
                                	xyzph[0] = (*(p_geom->xOffset(2,p_rec[instb]->NSTB))) - xyzph[0];
                        	}
				else{
					xyzph[0] = (*(p_geom->xOffset(2,p_rec[instb]->NSTB))) + xyzph[0];
				}
				xyzph[1] = (*(p_geom->yOffset(2,p_rec[instb]->NSTB))) - xyzph[1];
					
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
			TIter next(p_rec[instb]->clust[ncl].tow);
			while(TowerFPD* tow = (TowerFPD*)next()){
				if(tow->adc_over_ped>=1){			//minADC=1
					Int_t snstb = p_rec[instb]->NSTB; 	//starts from 1
					Int_t srow = (tow->row) - 1;		//srow starts from 0
					Int_t scol = (tow->col) - 1;		//scol starts from 0
					UInt_t adc = (UInt_t)(tow->adc_over_ped);
					Float_t tenergy = tow->energy;
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

Bool_t StFmsPointMaker::initialiseEnergyMatrices() {
  const StTriggerData* triggerData = 
      static_cast<const StTriggerData*>(StMuDst::event()->triggerData());
  // Loop over QT crate, slot and channel in the trigger data
  // and fill FMS energy matrices for this event.
  for (unsigned short crt = 1; crt <= 4; crt++){
    for (unsigned short slot = 1; slot <= 16; slot++){
      for (unsigned short ch = 0; ch < 32; ch++){
        unsigned short adc=triggerData->fmsADC(crt, slot - 1, ch);
        unsigned short tdc=triggerData->fmsTDC(crt, slot - 1, ch);
        // Determine the detector ID, row and channel
        Int_t detector(0);
        Int_t channel(0);
        mFmsDbMaker->getReverseMap(crt, slot, ch, &detector, &channel);
        Int_t row = mFmsDbMaker->getRowNumber(detector, channel);
        Int_t column = mFmsDbMaker->getColumnNumber(detector, channel);
        Int_t nstb = 0;
        switch (detector) {
          case 9:  // south large
            nstb = 2;
            break;
          case 8:  // north large
            nstb = 1;
            break;
          case 11:  // south small
            nstb = 4;
            break;
          case 10:  // north small
            nstb = 3;
            break;
          default:
            nstb = 0;
            break;
        }  // switch	
        Float_t energy = 0.f;
        if (detector > 0 || channel > 0) {
          Float_t g1 = mFmsDbMaker->getGain(detector, channel);
          Float_t g2 = mFmsDbMaker->getGainCorrection(detector, channel);
          energy = adc * g1 * g2;
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
        if (adc > 0) {
          mEnergyMatrices.at(nstb - 1)[row - 1][column - 1] = energy;
        }  // if
      }  // for
    }  // for
  }  // for
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

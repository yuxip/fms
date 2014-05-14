#include "StFmsQAHistoMaker.h"

#include <sstream>
#include <string>

#include "StEvent/StFmsCluster.h"
#include "StEvent/StFmsPoint.h"

#include "TObjArray.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include <TMath.h>
#include "TTree.h"
#include "TLorentzVector.h"

// STAR
#include "StMessMgr.h"
#include "StEventTypes.h"

//StMuDstMaker
#include "StMuDSTMaker/COMMON/StMuTypes.hh"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuEmcCollection.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuEmcUtil.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuPrimaryVertex.h"
#include "StMuDSTMaker/COMMON/StMuTriggerIdCollection.h"

//StEmc
#include "StEmcClusterCollection.h"
#include "StEmcPoint.h"
#include "StEmcUtil/geometry/StEmcGeom.h"
#include "StEmcUtil/others/emcDetectorName.h"
#include "StEmcADCtoEMaker/StBemcData.h"
#include "StEmcADCtoEMaker/StEmcADCtoEMaker.h"
#include "StEmcRawMaker/defines.h"
#include "StEmcRawMaker/StBemcRaw.h"
#include "StEmcRawMaker/StBemcTables.h"
#include "StEmcRawMaker/StEmcRawMaker.h"
#include "StEmcRawMaker/defines.h"
#include "StEmcUtil/database/StBemcTables.h"

//EEMC
#include "StEEmcUtil/database/StEEmcDb.h"
#include "StEEmcUtil/database/EEmcDbItem.h"
#include "StEEmcUtil/EEmcGeom/EEmcGeomSimple.h"

#include "StRoot/StFmsDbMaker/StFmsDbMaker.h"
#include "StRoot/StMuDSTMaker/COMMON/StMuFmsCluster.h"
#include "StRoot/StMuDSTMaker/COMMON/StMuFmsPoint.h"

#include <assert.h>

ClassImp(StFmsQAHistoMaker)

StFmsQAHistoMaker::StFmsQAHistoMaker(const char* name):StMaker(name) {
	
	SetFmsQA();
	SetEmcQA();
	mEmcEt = 0.2; //Et threshold for EMC qa
	SetTrackQA();
}

StFmsQAHistoMaker::~StFmsQAHistoMaker() {
	
	LOG_INFO << " StFmsQAHistoMaker:: destructor " << endm;
	
/*	if(hbbcratevsevt) { delete hbbcratevsevt; hbbcratevsevt = 0;}
	if(hbunchid) { delete hbunchid; hbunchid = 0;}
	
	if(hfmsNhitvsevt) { delete hfmsNhitvsevt; hfmsNhitvsevt = 0;}
	if(hfmshitEvsevt) { delete hfmshitEvsevt; hfmshitEvsevt = 0;}
	if(hfmsNcluvsevt) { delete hfmsNcluvsevt; hfmsNcluvsevt = 0;}
	if(hfmscluEvsevt) { delete hfmscluEvsevt; hfmscluEvsevt = 0;}
	if(hfmsNphovsevt) { delete hfmsNphovsevt; hfmsNphovsevt = 0;}
	if(hfmsphoEvsevt) { delete hfmsphoEvsevt; hfmsphoEvsevt = 0;}
//	if(hfmshitEvseta) { delete hfmshitEvseta; hfmshitEvseta = 0;}
//	if(hfmshitEvsphi) { delete hfmshitEvsphi; hfmshitEvsphi = 0;}
	if(hfmscluEvseta) { delete hfmscluEvseta; hfmscluEvseta = 0;}
	if(hfmscluEvsphi) { delete hfmscluEvsphi; hfmscluEvsphi = 0;}
	if(hfmsphoEvseta) { delete hfmsphoEvseta; hfmsphoEvseta = 0;}
	if(hfmsphoEvsphi) { delete hfmsphoEvsphi; hfmsphoEvsphi = 0;}

	if(hemcNhitvsevt) { delete hemcNhitvsevt; hemcNhitvsevt = 0;}
	if(hemchitEvsevt) { delete hemchitEvsevt; hemchitEvsevt = 0;}
	if(hemchitEvseta) { delete hemchitEvseta; hemchitEvseta = 0;}
	if(hemchitEvsphi) { delete hemchitEvsphi; hemchitEvsphi = 0;}
	if(hemchitEtvsevt) { delete hemchitEtvsevt; hemchitEtvsevt = 0;}
	if(hemchitEtvseta) { delete hemchitEtvseta; hemchitEtvseta = 0;}
	if(hemchitEtvsphi) { delete hemchitEtvsphi; hemchitEtvsphi = 0;}
	if(hemchitEtavsevt) { delete hemchitEtavsevt; hemchitEtavsevt = 0;}
	if(hemchitPhivsevt) { delete hemchitPhivsevt; hemchitPhivsevt = 0;}
	if(hemchitEtavsphi) { delete hemchitEtavsphi; hemchitEtavsphi = 0;}
	
	if(htpcverzvsevt) { delete htpcverzvsevt; htpcverzvsevt = 0;}
	if(htpcntrkvsevt) { delete htpcntrkvsevt; htpcntrkvsevt = 0;}
	if(htpcptvsevt) { delete htpcptvsevt; htpcptvsevt = 0;}
	if(htpcetavsevt) { delete htpcetavsevt; htpcetavsevt = 0;}
	if(htpcphivsevt) { delete htpcphivsevt; htpcphivsevt = 0;}
	if(htpcetavsphi) { delete htpcetavsphi; htpcetavsphi = 0;}
*/	
	
}

void StFmsQAHistoMaker::Clear( const char* opt ) {
	
}

void StFmsQAHistoMaker::SetOutputFile( Char_t* filename = "default.root" ) {
	
	mFilename = filename;
}

Int_t StFmsQAHistoMaker::InitRun(Int_t runNumber) {
  // Ensure we can access database information
  StFmsDbMaker* fmsDbMaker = static_cast<StFmsDbMaker*>(GetMaker("fmsDb"));
  if (!fmsDbMaker) {
    return kStErr;
  }  // if
  // Set up geometry, which stays constant for each run
  if (!mGeometry.initialize(fmsDbMaker)) {
    // Return an error if geometry initialization fails
    return kStErr;
  }  // if
  return StMaker::InitRun(runNumber);
}

Int_t StFmsQAHistoMaker::Init() {
  LOG_INFO << "StFmsQAHistoMaker::Init() " << endm;
  if (mEmcQA) {
    // If data, StEmcADCtoEMaker will be in the chain
    if (StEmcADCtoEMaker* adc2e = (StEmcADCtoEMaker*)StMaker::GetChain()->GetMakerInheritsFrom("StEmcADCtoEMaker")) {
      mBemcTables = adc2e->getBemcData()->getTables();
    }  // if
    assert(mBemcTables);
  }  // if
	mFile = new TFile(mFilename,"recreate");
	assert(mFile);
	
	hbbcratevsevt = new TH2F("hbbcratevsevt","bbc coincidence rate vs event number",2e3,0,2e6,500,1e6,4e6);
	hbunchid = new TH1F("hbunchid","bunchXing id",120,0,120);
	if(mFmsQA){
		hfmsNhitvsevt = new TH2F("hfmsNhitvsevt","FMS #hits vs event number",2e3,0,2e6,1264,0,1264);
		hfmshitEvsevt = new TH2F("hfmshitEvsevt","FMS hit energy vs event number",2e3,0,2e6,250,0,250);
		hfmsNcluvsevt = new TH2F("hfmsNcluvsevt","FMS #clusters vs event number",2e3,0,2e6,50,0,50);
		hfmscluEvsevt = new TH2F("hfmscluEvsevt","FMS cluster energy vs event number",2e3,0,2e6,250,0,250);
		hfmsNphovsevt = new TH2F("hfmsNphovsevt","FMS #photons vs event number",2e3,0,2e6,50,0,50);
		hfmsphoEvsevt = new TH2F("hfmsphoEvsevt","FMS photon energy vs event number",2e3,0,2e6,250,0,250);
	//	hfmshitEvseta = new TH2F("hfmshitEvseta","FMS hit energy vs eta",100,2.5,4.5,250,0,250);
	//	hfmshitEvsphi = new TH2F("hfmshitEvsphi","FMS hit energy vs phi",100,-TMath::Pi(),TMath::Pi(),250,0,250);	
    hfmshitEvsChannel = new TH2F("hfmshitEvsChannel", "StEvent FMS hit energy vs channel", 300, 0, 600, 250, 0.1, 250);
    // Create via clone to ensure identical binning
    hmufmshitEvsChannel = static_cast<TH2F*>(hfmshitEvsChannel->Clone("hmufmshitEvsChannel"));
    hmufmshitEvsChannel->SetTitle("StMuDst FMS hit energy vs channel");
		hfmscluEvseta = new TH2F("hfmscluEvseta","FMS cluster energy vs eta",100,2.5,4.5,250,0,250);
		hmufmscluEvseta = static_cast<TH2F*>(hfmscluEvseta->Clone("hmufmscluEvseta"));
		hmufmscluEvseta->SetTitle("StMuDst FMS cluster energy vs eta");
		hfmscluEvsphi = new TH2F("hfmscluEvsphi","FMS cluster energy vs phi",100,-TMath::Pi(),TMath::Pi(),250,0,250);	
		hmufmscluEvsphi = static_cast<TH2F*>(hfmscluEvsphi->Clone("hmufmscluEvsphi"));
		hmufmscluEvsphi->SetTitle("StMuDst FMS cluster energy vs phi");
		hfmsphoEvseta = new TH2F("hfmsphoEvseta","FMS photon energy vs eta",100,2.5,4.5,250,0,250);
		hmufmsphoEvseta = static_cast<TH2F*>(hfmsphoEvseta->Clone("hmufmsphoEvseta"));
		hmufmsphoEvseta->SetTitle("StMuDst FMS photon energy vs eta");
		hfmsphoEvsphi = new TH2F("hfmsphoEvsphi","FMS photon energy vs phi",100,-TMath::Pi(),TMath::Pi(),250,0,250);	
		hmufmsphoEvsphi = static_cast<TH2F*>(hfmsphoEvsphi->Clone("hmufmsphoEvsphi"));
		hmufmsphoEvsphi->SetTitle("StMuDst FMS photon energy vs phi");
    for (int nstb(1); nstb < 5; ++nstb) {
      std::ostringstream oss;
      oss << "hfmsy0x0_" << nstb;
      TH2F* hist = new TH2F(oss.str().c_str(), ";local row;local column",
                            18, 0, 18, 36, 0, 36);
      hfmsy0x0.insert(std::make_pair(nstb, hist));
      oss.str("");
      oss << "Row vs. column, NSTB = " << nstb;
      hfmsy0x0[nstb]->SetTitle(oss.str().c_str());
    }  // for
	}
	if(mEmcQA){
	//	hemcNhitvsevt = new TH2F("hemcNhitvsevt","EMC #towers vs event number",2e3,0,2e6,5520,0,5520);
	//	hemchitEvsevt = new TH2F("hemchitEvsevt","EMC tower energy vs event number",2e3,0,2e6,250,0,250);
	//	hemchitEvseta = new TH2F("hemchitEvseta","EMC tower energy vs eta",300,-1,2,250,0,250);
	//	hemchitEvsphi = new TH2F("hemchitEvsphi","EMC tower energy vs phi",100,-TMath::Pi(),TMath::Pi(),250,0,250);
		hbemchitEtvsevt = new TH2F("hbemchitEtvsevt","BEMC tower Et vs event number",2e3,0,2e6,250,0,250);
		heemchitEtvsevt = new TH2F("heemchitEtvsevt","EEMC tower Et vs event number",2e3,0,2e6,250,0,250);
		hbemchitAdcvsevt = new TH2F("hbemchitAdcvsevt","BEMC tower Adc vs event number",2e3,0,2e6,4096,0,4096);
		heemchitAdcvsevt = new TH2F("heemchitAdcvsevt","EEMC tower Adc vs event number",2e3,0,2e6,4096,0,4096);
		hbemchitIdvsevt = new TH2F("hbemchitIdvsevt","BEMC tower Id vs event number",2e3,0,2e6,4800,0,4800);
		heemchitIdvsevt = new TH2F("heemchitIdvsevt","EEMC tower Id vs event number",2e3,0,2e6,720,0,720);
	//	hemchitEtvseta = new TH2F("hemchitEtvseta","EMC tower Et vs eta",300,-1,2,250,0,250);
	//	hemchitEtvsphi = new TH2F("hemchitEtvsphi","EMC tower Et vs phi",100,-TMath::Pi(),TMath::Pi(),250,0,250);
		hbemchitAdcvsid = new TH2F("hbemchitAdcvsid","BEMC tower Adc vs Id",4800,0,4800,4096,0,4096);
		heemchitAdcvsid = new TH2F("heemchitAdcvsid","EEMC tower Adc vs Id",720,0,720,4096,0,4096);
		hbemchitEtvsid = new TH2F("hbemchitEtvsid","BEMC tower Et vs Id",4800,0,4800,500,0,250);
		heemchitEtvsid = new TH2F("heemchitEtvsid","EEMC tower Et vs Id",720,0,720,500,0,250);
		hbemchitEtvseta = new TH2F("hbemchitEtvseta","BEMC tower Et vs eta",300,-1,2,250,0,250);
		hbemchitEtvsphi = new TH2F("hbemchitEtvsphi","BEMC tower Et vs phi",120,-TMath::Pi(),TMath::Pi(),250,0,250);
		heemchitEtvseta = new TH2F("heemchitEtvseta","EEMC tower Et vs eta",300,-1,2,250,0,250);
		heemchitEtvsphi = new TH2F("heemchitEtvsphi","EEMC tower Et vs phi",120,-TMath::Pi(),TMath::Pi(),250,0,250);
	//	hemchitEtavsevt = new TH2F("hemchitEtavsevt","EMC hit eta vs event number",2e3,0,2e6,300,-1,2);
	//	hemchitPhivsevt = new TH2F("hemchitPhivsevt","EMC hit phi vs event number",2e3,0,2e6,100,-TMath::Pi(),TMath::Pi());
	//	hemchitEtavsphi = new TH2F("hemchitEtavsphi","EMC hit eta vs phi",100,-TMath::Pi(),TMath::Pi(),300,-1,2);
		hbemchitEtavsphi = new TH2F("hbemchitEtavsphi","BEMC hit eta vs phi",120,-TMath::Pi(),TMath::Pi(),300,-1,2);
		heemchitEtavsphi = new TH2F("heemchitEtavsphi","EEMC hit eta vs phi",120,-TMath::Pi(),TMath::Pi(),300,-1,2);
	}
	if(mTrackQA){
		htpcverzvsevt = new TH2F("htpcverzvsevt","TPC vertex z vs event number",2e3,0,2e6,200,-100,100);
		htpcntrkvsevt = new TH2F("htpcntrkvsevt","TPC #tracks vs event numer",2e3,0,2e6,100,0,100);
		htpcptvsevt = new TH2F("htpcptvsevt","TPC track pT vs event number",2e3,0,2e6,100,0,100);
		htpcetavsevt = new TH2F("htpcetavsevt","TPC track eta vs event number",2e3,0,2e6,100,-2,2);
		htpcphivsevt = new TH2F("htpcphivsevt","TPC track phi vs event number",2e3,0,2e6,100,-TMath::Pi(),TMath::Pi());
		htpcetavsphi = new TH2F("htpcetavsphi","TPC track eta vs phi",100,-TMath::Pi(),TMath::Pi(),100,-2,2);
	}
	mEeDb = (StEEmcDb*)StMaker::GetChain()->GetDataSet("StEEmcDb");
	if(!mEeDb){
		LOG_ERROR << "StEEmcDb not found " <<endm;
		return kStErr;
	}
	
	if(mEeDb) mEeDb->setThreshold(3);
//	LOG_INFO << "StFmsQAHistoMaker Init() finished.." << endm;	
	return kStOk;
}
	
	
Int_t StFmsQAHistoMaker::Make() {

//	LOG_INFO << "StFmsQAHistoMaker::Make(): mFmsQA = "<<mFmsQA<<", mEmcQA = "<<mEmcQA<<", mTrackQA = "<<mTrackQA << endm;
	StMuDst* muDst = (StMuDst*)GetInputDS("MuDst");
	if(!muDst){
		LOG_WARN << "StFmsQAHistoMaker::Make() no MuDst" << endm;
		return kStWarn;
	}

	Rnum = muDst->event()->runNumber();
	ievt = muDst->event()->eventNumber();

	//remove FMS LED events;
	Bool_t checkLed = 0;
	if (muDst->event()->triggerIdCollection().nominal().isTrigger(320225))checkLed = true;
	
	Bool_t LedEvent = false;
	StL0Trigger& l0 = muDst->event()->l0Trigger();
	UInt_t mLedBit = l0.lastDsmArray(4);
	if(mLedBit & 0x0001) LedEvent = true;
	
	if(checkLed||LedEvent){
		LOG_INFO << " FMS LED fired, abort " <<endm;
		return kStOk;
	}
	
	StRunInfo* runInfo = &(muDst->event()->runInfo());
	Int_t fillnumber = runInfo->beamFillNumber(blue);
	Int_t bunchid = l0.bunchCrossingId7bit(Rnum);
	Int_t bbccoincrate = runInfo->bbcCoincidenceRate();
	
	hbbcratevsevt->Fill(ievt,bbccoincrate);
	hbunchid->Fill(bunchid);

	if(mFmsQA){
	  fmsEventQa();
    fmsMuDstQa();
	}//mFmsQA
//	LOG_INFO << "begin reading EMC" << endm;
	if(mEmcQA){

		StEvent *mEvent = static_cast<StEvent*>(this->GetDataSet("StEvent"));
		if(!mEvent){
			LOG_ERROR << " no StEvent " << endm;
			return kStErr;
		}		

		StEmcGeom* geom = StEmcGeom::instance("bemc"); // for towers
		if (!geom) LOG_WARN << " No StEmcGeom!" << endm;
		assert(geom);
		StEmcCollection *emc = mEvent->emcCollection();
		if (!emc) LOG_WARN << " No StEmcCollection" << endm;
		assert(emc);
		Int_t emcnhits = 0;
		Float_t Etcut = mEmcEt; //GeV, only look at hits with Et > 0.2 GeV
//		LOG_INFO <<" Et cut = "<<Etcut<<endm;
		//barrel
		StEmcDetector* detector = emc->detector(kBarrelEmcTowerId);
		if (!detector) LOG_WARN << " No StEmcDetector!" << endm;
		if (detector) {
			
			for(Int_t m = 1; m <= 120; m++){
				StEmcModule* module = detector->module(m);
				if(module){
					StSPtrVecEmcRawHit& rawHit=module->hits();
			//		emcnhits += rawHit.size();
					for(UInt_t k = 0; k < rawHit.size(); ++k){
						if(rawHit[k]){
			//				LOG_INFO <<" bemc k = "<<k<<endm;
							Int_t did;
							Int_t m=rawHit[k]->module();
							Int_t e=rawHit[k]->eta();
							Int_t s=abs(rawHit[k]->sub());
							Int_t adc=rawHit[k]->adc();
							Float_t energy=rawHit[k]->energy();
							
							Float_t tower_eta, tower_phi;
							geom->getId(m,e,s,did);
							geom->getEtaPhi(did,tower_eta,tower_phi); // to convert software id into eta/phi					
						/*	LOG_INFO <<" did = "<<did<<", tower_eta = "<<tower_eta<<", tower_phi = "<<tower_phi<<endm;
							//2nd way
							float towerX, towerY, towerZ;
							StEmcGeom::instance("bemc")->getXYZ(did, towerX, towerY, towerZ);
							TVector3 tower(towerX, towerY, towerZ);
							Float_t tmpeta = tower.Eta();
							Float_t tmpphi = tower.Phi();
							LOG_INFO <<" tmpeta = "<<tmpeta<<", tmpphi = "<<tmpphi<<endm;					*/		

							Float_t Et = (energy/cosh(tower_eta));
	
							Float_t pedestal;
							Float_t rms;
							Int_t status;
							Int_t CAP(0);
							mBemcTables->getPedestal(BTOW,did,CAP,pedestal,rms);
							mBemcTables->getStatus(BTOW,did,status);
							//same condition as jet finder
							if(status!=1)continue;
							if(Et<Etcut)continue;
							if(adc-pedestal <= 4 )continue;
							if(adc-pedestal <= 3*rms)continue;
					//		hemchitEvsevt->Fill(ievt,energy);
					//		hemchitEvseta->Fill(tower_eta,energy);
					//		hemchitEvsphi->Fill(tower_phi,energy);
							hbemchitAdcvsevt->Fill(ievt,(adc-pedestal));
							hbemchitEtvsevt->Fill(ievt,Et);
							hbemchitIdvsevt->Fill(ievt,did);
					//		hemchitEtvseta->Fill(tower_eta,Et);
					//		hemchitEtvsphi->Fill(tower_phi,Et);
							hbemchitEtvseta->Fill(tower_eta,Et);
							hbemchitEtvsphi->Fill(tower_phi,Et);
					//		hemchitEtavsevt->Fill(ievt,tower_eta);
					//		hemchitPhivsevt->Fill(ievt,tower_phi);
					//		hemchitEtavsphi->Fill(tower_phi,tower_eta);
							hbemchitEtavsphi->Fill(tower_phi,tower_eta);
							hbemchitAdcvsid->Fill(did,adc);
					//		hbemchitPedvsid->Fill(did,pedestal);
					//		hbemchitPedRmsvsid->Fill(did,rms);
							hbemchitEtvsid->Fill(did,Et);
							emcnhits++;
					//		LOG_INFO <<"pedestal = "<<pedestal<<", rms = "<<rms<<", status = "<<status<<endm;
							
						}
					}
				}
			}
		}
		
		//endcap --only include towers which passed basic EEMC QA ( accepted by the jet maker )
		StMuEmcCollection* muEmc = StMuDst::muEmcCollection();
		for (Int_t id = 0; id < muEmc->getNEndcapTowerADC(); ++id) {
		//	LOG_INFO<<"eemc id = "<<id<<endm;
			Int_t rawadc, sec, sub, etabin;
			muEmc->getEndcapTowerADC(id, rawadc, sec, sub, etabin);

			// Sanity check
			if (rawadc < 0 || rawadc >= 4095) continue;

			assert(1 <= sec && sec <= 12);
			assert(1 <= sub && sub <= 5);
			assert(1 <= etabin && etabin <= 12);

			const EEmcDbItem *dbItem = mEeDb->getT(sec,sub-1+'A',etabin);
			if(dbItem->fail) continue; //drop broken channels
			if(dbItem->stat) continue; // drop not working channels and jumpy pedestal channels
			if(dbItem->gain<=0.) continue; // drop it, unless you work with ADC spectra
			if(rawadc<dbItem->thr) continue; // drop raw ADC < ped+N*sigPed, N==3 in init
			Double_t adc = rawadc - (dbItem->ped);
			//same as jet finder
			if(adc <= 4)continue;
			if(adc <= 3*dbItem->sigPed)continue;
			Double_t energy = adc/(dbItem->gain);
			const EEmcGeomSimple& geom = EEmcGeomSimple::Instance();
			TVector3 towerLocation = geom.getTowerCenter(sec-1,sub-1,etabin-1);
			Float_t eemc_eta = towerLocation.Eta();		
			Float_t eemc_phi = towerLocation.Phi();		
			Float_t Et_eemc = (energy/cosh(towerLocation.Eta()));
			
			/*
			//to make eemc_id = 175;
			Int_t tmpsec = 3;
			Int_t tmpsub = 5;
			Int_t tmpetabin = 8;
			*/
		/*	//to make eemc_id = 548;
			Int_t tmpsec = 10;
			Int_t tmpsub = 1;
			Int_t tmpetabin = 9;

			TVector3 tmploc = geom.getTowerCenter(tmpsec-1,tmpsub-1,tmpetabin-1);
			Float_t tmpeta = tmploc.Eta();
			Float_t tmpphi = tmploc.Phi();
		//	LOG_INFO <<"========== EEMC_id = 175, eta = "<<tmpeta<<", phi = "<<tmpphi << endm;
		//	LOG_INFO <<"========== EEMC_id = 548, eta = "<<tmpeta<<", phi = "<<tmpphi << endm;
		*/	

			if(Et_eemc < Etcut)continue;
			Int_t eemc_id = (sec-1)*60+(sub-1)*12+(etabin-1);
			
	//		hemchitEvsevt->Fill(ievt,energy);
	//		hemchitEvseta->Fill(eemc_eta,energy);
	//		hemchitEvsphi->Fill(eemc_phi,energy);
			heemchitEtvsevt->Fill(ievt,Et_eemc);
			heemchitAdcvsevt->Fill(ievt,adc);
			heemchitIdvsevt->Fill(ievt,eemc_id);
			heemchitEtvseta->Fill(eemc_eta,Et_eemc);
			heemchitEtvsphi->Fill(eemc_phi,Et_eemc);
	//		hemchitEtavsevt->Fill(ievt,eemc_eta);
	//		hemchitPhivsevt->Fill(ievt,eemc_phi);
			heemchitEtavsphi->Fill(eemc_phi,eemc_eta);
			heemchitEtvsid->Fill(eemc_id,Et_eemc);
			heemchitAdcvsid->Fill(eemc_id,adc);

			emcnhits++;
		}
			
	//	hemcNhitvsevt->Fill(ievt,emcnhits);
	//	LOG_INFO << "debug ievt = "<<ievt<<", emcnhits = "<<emcnhits<<endm;
		

	}//mEmcQA
	
	if(mTrackQA){
		
		Float_t vertex_z_cut = 60.0; //cm
		
		Int_t nVertices = muDst->numberOfPrimaryVertices();
		for(Int_t j = 0; j < nVertices; j++){
			
			StMuPrimaryVertex* vertex = muDst->primaryVertex(j);
			assert(vertex);
			muDst->setVertexIndex(j);
			StThreeVectorF r = vertex->position();
			Float_t verz = r.z();
			if(vertex->ranking() < 0 || j > 1)continue;
			htpcverzvsevt->Fill(ievt,verz); //after ranking cut
			if(fabs(r.z()) > vertex_z_cut) continue; //skip bad vertices
			Int_t nTracks = muDst->GetNPrimaryTrack();
			htpcntrkvsevt->Fill(ievt,nTracks);
			for(Int_t i = 0; i < nTracks; i++){
				
				StMuTrack* track = muDst->primaryTracks(i);
				assert(track);
				StThreeVectorF momentum = track->momentum();
				Bool_t flagTrack = isUsableTrack(*track);
          			if(flagTrack == false)continue;
				Float_t trackpT = track->pt();
				Float_t trackEta = track->eta();
				Float_t trackPhi = track->phi();
				htpcptvsevt->Fill(ievt,trackpT);
				htpcetavsevt->Fill(ievt,trackEta);
				htpcphivsevt->Fill(ievt,trackPhi);
				htpcetavsphi->Fill(trackPhi,trackEta);
				
			}
		}	
	}//mTrackQA
	
	return kStOk;
}

Int_t StFmsQAHistoMaker::Finish(){
	mFile->cd();
	TH2F* diffEvsChannel = static_cast<TH2F*>(hfmshitEvsChannel->Clone("diffEvsChannel"));
	diffEvsChannel->Add(hmufmshitEvsChannel, -1.);
	diffEvsChannel->SetOption();
	mFile->Write();
	mFile->Close();
	
	return kStOk;
}

Bool_t StFmsQAHistoMaker::isUsableTrack(const StMuTrack& track)
{
  if(track.flag() < 0)
    return false;

  if(track.dcaGlobal().mag()> 3.)
    return false;

  int dcaFlag = 1;
  if(true){
    Double_t limit = 3.-2.*track.pt();
    if(!((track.pt()<0.5&&track.dcaGlobal().mag()<=2.)||
         ((track.pt()>=0.5&&track.pt()<1.0)&&track.dcaGlobal().mag()<=limit)||(track.pt()>=1.0&&track.dcaGlobal().mag()<=1.0))) dcaFlag =0;
  }
  if(dcaFlag == 0)
    return false;
  if (track.topologyMap().trackFtpcEast() || track.topologyMap().trackFtpcWest()) {
            return false;
        }
  //    if(track->eta() < GetEtaLow()) { //GetEtaLow defined in StFourPMaker.cxx
  //        return false;
  //    }
  //    if(track->eta() > GetEtaHigh()) {
  //      return false;
  //    }
        if(static_cast<double>(track.nHits())/static_cast<double>(track.nHitsPoss()) < .51)
          return false;


   return true;
}

void StFmsQAHistoMaker::fmsEventQa() {
  StEvent *event = static_cast<StEvent*>(GetDataSet("StEvent"));
  if(!event){
    LOG_ERROR << " no StEvent " << endm;
    return;
  }
  StFmsCollection* fmsCollection = event->fmsCollection();
  if (!fmsCollection) {
    LOG_ERROR << "StFmsQAHistoMaker did not find "
      << "an StFmsCollection in StEvent" << endm;
      return;
  }  // if
  const StSPtrVecFmsHit& fmshits = fmsCollection->hits();
  const StSPtrVecFmsCluster& fmsclusters = fmsCollection->clusters();
  hfmsNcluvsevt->Fill(ievt,fmsclusters.size());
  const StSPtrVecFmsPoint& fmspoints = fmsCollection->points();
  Int_t nphotons = fmspoints.size();
  Int_t nhits = 0;
  for (StSPtrVecFmsHitConstIterator i = fmshits.begin(); i != fmshits.end(); ++i) {
    StFmsHit* hit = *i;
    if (hit->detectorId() >= 8 && hit->detectorId() <= 11 && hit->adc() > 0) {
      hfmshitEvsChannel->Fill(hit->channel(), hit->energy());
    }  // if
  }  // for
  for(StSPtrVecFmsClusterConstIterator iclu = fmsclusters.begin(); iclu != fmsclusters.end(); iclu++){
    nhits += (*iclu)->nTowers();
    Float_t clusterE = (*iclu)->energy();
    hfmscluEvsevt->Fill(ievt,clusterE);
    Float_t clusterEta = ((*iclu)->fourMomentum()).Eta();
    Float_t clusterPhi = ((*iclu)->fourMomentum()).Phi();
    hfmscluEvseta->Fill(clusterEta,clusterE);
    hfmscluEvsphi->Fill(clusterPhi,clusterE);
    const int nstb = (*iclu)->detectorId();
    if (hfmsy0x0.find(nstb) != hfmsy0x0.end()) {
      if ((*iclu)->x() > 0. && (*iclu)->y() > 0.) {
        hfmsy0x0[nstb]->Fill((*iclu)->x(), (*iclu)->y());
      }  // if
    }  // if
    //loop over hits
    for(StPtrVecFmsHitConstIterator ihit = (*iclu)->hits().begin(); ihit != (*iclu)->hits().end(); ihit++){
      Float_t hitE = (*ihit)->energy();
      hfmshitEvsevt->Fill(ievt,hitE);
    }
  }
  for(StSPtrVecFmsPointConstIterator ipts = fmspoints.begin(); ipts != fmspoints.end(); ipts++){
    //(*ipts)->Print();
    Float_t photonE = (*ipts)->energy();
    hfmsphoEvsevt->Fill(ievt,photonE);
    Float_t photonEta = ((*ipts)->fourMomentum()).Eta();
    Float_t photonPhi = ((*ipts)->fourMomentum()).Phi();
    hfmsphoEvseta->Fill(photonEta,photonE);
    hfmsphoEvsphi->Fill(photonPhi,photonE);
  }
  hfmsNhitvsevt->Fill(ievt,nhits);
  hfmsNphovsevt->Fill(ievt,nphotons);
}

void StFmsQAHistoMaker::fmsMuDstQa() {
	StMuDst* muDst = (StMuDst*)GetInputDS("MuDst");
	if(!muDst){
		return;
	}  // if
  // Fill histograms using hits from StMuDst
  StMuFmsCollection* mufmsCollection = muDst->muFmsCollection();
  for (int i(0); i < mufmsCollection->numberOfHits(); ++i) {
    StMuFmsHit* hit = mufmsCollection->getHit(i);
    if (hit->detectorId() >= 8 && hit->detectorId() <= 11 && hit->adc() > 0) {
      hmufmshitEvsChannel->Fill(hit->channel(), hit->energy());
    }  // if
  }  // for
  for (int i(0); i < mufmsCollection->numberOfClusters(); ++i) {
    StMuFmsCluster* cluster = mufmsCollection->getCluster(i);
    if (cluster) {
      TVector3 xyz = mGeometry.columnRowToGlobalCoordinates(
        cluster->x(), cluster->y(), cluster->detectorId());
      hmufmscluEvseta->Fill(xyz.Eta(), cluster->energy());
      hmufmscluEvsphi->Fill(xyz.Phi(), cluster->energy());
    }  // if
  }  // for
  for (int i(0); i < mufmsCollection->numberOfPoints(); ++i) {
    StMuFmsPoint* point = mufmsCollection->getPoint(i);
    TVector3 xyz = TVector3(point->x(), point->y(),
                            mGeometry.z(point->detectorId()));
    hmufmsphoEvseta->Fill(xyz.Eta(), point->energy());
    hmufmsphoEvsphi->Fill(xyz.Phi(), point->energy());
  }  // for
}

/*
void StFmsQAHistoMaker::SetFmsPhotonEcut ( Float_t ecutsmall, Float_t ecutlarge ) {
	
	mFmsEcutSmall = ecutsmall;
	mFmsEcutLarge = ecutlarge;
}
*/


#include "TH1.h"
#include "TChain.h"
#include "TSystem.h"
#include "TFile.h"
#include <iostream>

class StMuDstMaker;
StMuDstMaker* maker;

void testFms_stfms( Int_t ibegin = 1, Int_t iend = 100, const char* file = "fms7.list", const char* outfile = "stfmsAnal_run12098008.root", const char* qafile = "stfmsQAhisto_run12098008.root", Bool_t qa = true ){

	gROOT->Macro("loadMuDst.C");
	gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
	loadSharedLibraries();
	gROOT->Macro("LoadLogger.C");
	gSystem->Load("StTpcDb");
	gSystem->Load("StDetectorDbMaker");
	gSystem->Load("StDbUtilities");
	gSystem->Load("StDaqLib");
	gSystem->Load("StEmcRawMaker");
	gSystem->Load("StEmcADCtoEMaker");
	gSystem->Load("StPreEclMaker");
	gSystem->Load("StEpcMaker");
	gSystem->Load("StDbBroker");
	gSystem->Load("St_db_Maker");
	gSystem->Load("StEEmcUtil");
	gSystem->Load("StEEmcDbMaker");
	gSystem->Load("StSpinDbMaker");
	gSystem->Load("StEmcTriggerMaker");
	gSystem->Load("StTriggerUtilities");
	gSystem->Load("StMCAsymMaker");
	gSystem->Load("StRandomSelector");
	

	// Load your shared libraries here
	gSystem->Load("libMinuit.so");
	gSystem->Load("StFmsDbMaker");
	//gSystem->Load("StTriggerDataMaker");
	gSystem->Load("StTriggerFilterMaker");
	gSystem->Load("StPSUTools");
	using namespace PSUGlobals;
	gSystem->Load("StFmsHitMaker");
	gSystem->Load("StFmsPointMaker");
	gSystem->Load("StFmsQAHistoMaker");
	//gSystem->Load("StFmsAnalysisMaker");

	StChain* chain = new StChain("StChain");
	chain->SetDEBUG(0);
	StMuDstMaker* muDstMaker = new StMuDstMaker(0,0,"",file,".",1000,"MuDst");
	StMuDbReader* muDB = StMuDbReader::instance();
	gMessMgr->SwitchOff("D");
	gMessMgr->SwitchOff("I");	

	// Trigger filter
	StTriggerFilterMaker* filterMaker = new StTriggerFilterMaker;
	//filterMaker->addTrigger(89);//FMS*BEMC*JP0
	filterMaker->addTrigger(320220); //   FMSJP1
	filterMaker->addTrigger(320231); //   FMSJP2
	filterMaker->addTrigger(320227); //   FMSLgBS2
	filterMaker->addTrigger(320226); //   FMSLgBS1
	filterMaker->addTrigger(320222); //   FMSSmBS1
	filterMaker->addTrigger(320223); //   FMSSmBS2
	St_db_Maker *dbMk = new St_db_Maker("db","MySQL:StarDb","$STAR/StarDb");
	dbMk->SetDEBUG(0);
	dbMk->SetDateTime(20110601, 0);

	StFmsDbMaker* fmsdb = new StFmsDbMaker("fmsDb");
	fmsdb->setDebug(0);

	// Endcap database
	StEEmcDbMaker* eemcDb = new StEEmcDbMaker;
	// Barrel ADC to energy maker
        StEmcADCtoEMaker* adc = new StEmcADCtoEMaker;
        adc->saveAllStEvent(true);
	
	
	StFmsHitMaker* fmshitMk = new StFmsHitMaker();
	//fmshitMk->SetDEBUG();
	
	
	StFmsPointMaker* fmsptMk = new StFmsPointMaker("StFmsPointMaker",fmshitMk);
	if(qa){
		StFmsQAHistoMaker* fmsQa = new StFmsQAHistoMaker();
		fmsQa->SetOutputFile(qafile);
	}
	
	/*StFmsAnalysisMaker* fms = new StFmsAnalysisMaker("StFmsAnalysisMaker",jetmaker);
	fms->EnableEdepCor = true;
	fms->SetOutputFile(outfile);
	*/
	chain->Init();
	chain->EventLoop(ibegin,iend);
	chain->Finish();
	delete chain;
}


#ifndef ST_FMS_QAHISTO_MAKER_H
#define ST_FMS_QAHISTO_MAKER_H
#include <stdio.h>
#include <stdlib.h>

#include <map>

// ROOT
class TFile;
class TH1F;
class TH2F;
class TTree;
class TObjArray;
class TLorentzVector;

// STAR
#include "StMaker.h"
class StMuEvent;
class StMuTrack;
class StMuDstMaker;
class StBemcTables;
class StEEmcDb;

#include "StRoot/StFmsPointMaker/StFmsGeometry.h"

//only save event QA histograms (run-by-run). Lighter version of StFmsQAMaker
class StFmsQAHistoMaker : public StMaker {

public:
	StFmsQAHistoMaker( const char* name = "StFmsQAMaker" );

	~StFmsQAHistoMaker();
	
	void Clear( const char* opt = "" );
	void SetOutputFile( Char_t* filename );
	void SetFmsQA ( Bool_t fmsqa = true ) { mFmsQA = fmsqa; }
	void SetEmcQA ( Bool_t emcqa = true ) { mEmcQA = emcqa; }
	void SetEmcEt ( Float_t emcet = 0.2 ) { mEmcEt = emcet; }
	void SetTrackQA ( Bool_t trackqa = true ) { mTrackQA = trackqa; }
	Int_t InitRun(Int_t runNumber);
	Int_t Init();
        Int_t Make();
        Int_t Finish();
	
private:
	
	TFile* mFile;
	Bool_t mFmsQA;
	Bool_t mEmcQA;
	Float_t mEmcEt;
	Bool_t mTrackQA;

	UInt_t Rnum;			//run number
	UInt_t ievt;			//event number
	Char_t* mFilename;
	
	StBemcTables* mBemcTables;
	StEEmcDb* mEeDb;
	
	//BbcQA
	TH2F* hbbcratevsevt;	//bbc coincidence rate vs event#
	TH1F* hbunchid;

	//FmsQA
	TH2F* hfmsNhitvsevt;
	TH2F* hfmshitEvsevt;
	TH2F* hfmsNcluvsevt;
	TH2F* hfmscluEvsevt;
	TH2F* hfmsNphovsevt;
	TH2F* hfmsphoEvsevt;
//	TH2F* hfmshitEvseta;
//	TH2F* hfmshitEvsphi;
  TH2F* hfmshitEvsChannel;
  TH2F* hmufmshitEvsChannel;  // Like hfmshitEvsChannel but for StMuDst hits
	TH2F* hfmscluEvseta;
	TH2F* hmufmscluEvseta;  // Like hfmscluEvseta but for StMuDst clusters
	TH2F* hfmscluEvsphi;
	TH2F* hmufmscluEvsphi;  // Like hfmscluEvsphi but for StMuDst clusters
	TH2F* hfmsphoEvseta;
	TH2F* hmufmsphoEvseta;  // Like hfmsphoEvseta but for StMuDst photons
	TH2F* hfmsphoEvsphi;
	TH2F* hmufmsphoEvsphi;  // Like hfmsphoEvsphi but for StMuDst photons
  std::map<int, TH2F*> hfmsy0x0;
	
	//EmcQA
	//TH2F* hemchitEtavsevt = new TH2F("hemchitEtavsevt","emc hit #eta vs evt",2e4,0,2e6,100,-1,1);
//	TH2F* hemcNhitvsevt; //hit with good status
//	TH2F* hemchitEvsevt;
//	TH2F* hemchitEvseta;
//	TH2F* hemchitEvsphi;
	TH2F* hbemchitEtvsevt;
	TH2F* heemchitEtvsevt;
	TH2F* hbemchitAdcvsevt;
	TH2F* heemchitAdcvsevt;
	TH2F* hbemchitIdvsevt;
	TH2F* heemchitIdvsevt;
	
//	TH2F* hemchitEtvseta;
//	TH2F* hemchitEtvsphi;
//	TH2F* hemchitEtavsevt;	//for Et > 0.2 GeV with good tower status
//	TH2F* hemchitPhivsevt;
//	TH2F* hemchitEtavsphi;
	
	TH2F* hbemchitAdcvsid;
	TH2F* heemchitAdcvsid;
	TH2F* hbemchitEtvsid;
	TH2F* heemchitEtvsid;
	TH2F* hbemchitEtvseta;
	TH2F* hbemchitEtvsphi;
	TH2F* heemchitEtvseta;
	TH2F* heemchitEtvsphi;
	TH2F* hbemchitEtavsphi;
	TH2F* heemchitEtavsphi;	

	//TPCQA
	TH2F* htpcverzvsevt;
	TH2F* htpcntrkvsevt;
	TH2F* htpcptvsevt;
	TH2F* htpcetavsevt;
	TH2F* htpcphivsevt;
	TH2F* htpcetavsphi;

	Bool_t isUsableTrack(const StMuTrack& track);
  FMSCluster::StFmsGeometry mGeometry;
	ClassDef(StFmsQAHistoMaker,0)
};

#endif

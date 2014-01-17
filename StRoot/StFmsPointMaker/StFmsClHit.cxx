#include "StMessMgr.h"
#include "TMath.h"
#include "StFmsClHit.h"

#include <stdlib.h>

ClassImp(StFmsClHit)

StFmsClHit::StFmsClHit():StObject(){
	
	mNstb	= 0;
	mRow	= 0;
	mCol	= 0;
	mAdc	= 0;
	mEnergy	= 0;
	mStatus	= 0;

}

//___________________________________________________________________
StFmsClHit::StFmsClHit(Int_t nstb, Float_t row0, Float_t col0, UInt_t adc, Float_t energy, UChar_t status){

	if(!Legal(2,nstb,row0,col0)){
		LOG_ERROR << "illegal (nstb, row ,col) combo: "<<" ( "<<nstb<<", "<<row0<<", "<<col0<<" )"<<endm;
		exit(0);
	}

	mNstb	= nstb;
	mRow	= row0;
	mCol	= col0;

	SetEnergy(energy);
	SetAdc(adc);
	SetStatus(status);
}

//___________________________________________________________________
StFmsClHit::~StFmsClHit(){
;
}

//___________________________________________________________________
Bool_t StFmsClHit::Legal(Int_t iew, Int_t nstb, Float_t row0, Float_t col0 ){

	if(iew>0 && iew<2)                              return false;
        if(nstb<1 || nstb>4)                            return false;
        if(nstb>2){
                if(row0<0 || row0>23)                   return false;
                if(col0<0 || col0>11)                   return false;
                if(fabs(1.*row0-11.5)<5 && col0<5)      return false;
        }
        else{
                if(row0<0 || row0>33)                   return false;
                if(col0<0 || col0>16)                   return false;
                if(fabs(1.*row0-16.5)<8 && col0<8)      return false;
                if(row0<col0-9.5)                       return false;
                if(33-row0<col0-9.5)                    return false;
        }

        return true;
}

//___________________________________________________________________
void StFmsClHit::Print(Option_t *option)	const { LOG_INFO << *this << endm; }

//___________________________________________________________________
ostream& operator<<(ostream& os, const StFmsClHit& v){

	return os <<"StFmsClHit:\n\tnstb:\t"<<v.GetNstb()
		  <<"\n\trow:\t"<<v.GetRow()
		  <<"\n\tcol:\t"<<v.GetCol()
		  <<"\n\teff.adc:\t"<<v.GetAdc()
		  <<"\n\tenergy:\t"<<v.GetEnergy()
		  <<"\n\tstatus:\t0x"<<hex<<(Int_t)v.GetStatus()<<dec<<"\n";
}

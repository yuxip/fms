#include <iostream>
#include <stdlib.h>
#include <math.h>


#include "TowerFPD.h"

#include "StEvent/StFmsHit.h"

using namespace PSUGlobals;
ClassImp(TowerFPD);

TowerFPD::TowerFPD() 
{
  hit = NULL;
  col = row =  cluster = -1 ;
  adc_over_ped=9999;
  IEW=0;// East/West undefined
  NSTB=0;// Detector Undefined
  Lnk_LRUD=0;
}

TowerFPD::TowerFPD(const StFmsHit* fmsHit, Int_t towX, Int_t towY, Int_t clu)
{
  hit = fmsHit;
  col = towX;
  row = towY;
  cluster = clu;
  adc_over_ped=9999;
 
  IEW=0;// East/West undefined
  NSTB=0;// Detector Undefined
  Lnk_LRUD=0;
};

Int_t TowerFPD::Compare(const TObject *obj) const
{
  if( hit->energy() < ((TowerFPD *) obj)->hit->energy() )
    return -1;
  else if( hit->energy() > ((TowerFPD *) obj)->hit->energy() )
    return 1;
  else
    return 0;
};

Bool_t TowerFPD::IsNeighbor(TowerFPD *a) 
{
  if(!a)return false;
  return ( abs(this->col - a->col) + abs(this->row - a->row) == 1 )  ;
};

Bool_t TowerFPD::SetContext(TObjArray* towers,Int_t iEW,Int_t nSTB)
{
  IEW=iEW;
  NSTB=nSTB;
  TIter next(towers);
  if(Lnk_LRUD)delete Lnk_LRUD;
  Lnk_LRUD=new TObjArray(4);  
  Lnk_LRUD->SetOwner(0);
  while(TowerFPD* tower=(TowerFPD*) next())
    {
      if(IsNeighbor(tower))
	{
	  if((tower->col)<col)        {Lnk_LRUD->AddAt(tower,0);}
	  else if ((tower->col)>col)  {Lnk_LRUD->AddAt(tower,1);}
	  else if ((tower->row)<row)  {Lnk_LRUD->AddAt(tower,2);}
	  else if ((tower->row)>row)  {Lnk_LRUD->AddAt(tower,3);}
	};
    };
  return true;
};

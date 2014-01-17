#include <iostream>
#include <stdlib.h>
#include <math.h>


#include "TowerFPD.h"
using namespace PSUGlobals;
ClassImp(TowerFPD);

TowerFPD::TowerFPD() 
{
  energy = 0;
  col = row =  cluster = -1 ;
  adc_over_ped=9999;
  LinkDefined=false;
  CalibSet=false;  
  Gain=GainCorr=0;
  IEW=0;// East/West undefined
  NSTB=0;// Detector Undefined
  ClaimedID=-1;// Mark as unclaimed
  Lnk_LRUD=0;
}

TowerFPD::TowerFPD(Float_t ene, Int_t towX, Int_t towY, Int_t clu)
{
  energy = ene;
  col = towX;
  row = towY;
  cluster = clu;
  adc_over_ped=9999;
 
  LinkDefined=false;
  CalibSet=false;  
  Gain=GainCorr=0;
  IEW=0;// East/West undefined
  NSTB=0;// Detector Undefined
  ClaimedID=-1;// Mark as unclaimed
  Lnk_LRUD=0;
};

Int_t TowerFPD::Compare(const TObject *obj) const
{
  if( energy < ((TowerFPD *) obj)->energy )
    return -1;
  else if( energy > ((TowerFPD *) obj)->energy )
    return 1;
  else
    return 0;
};

Bool_t TowerFPD::IsNeighbor(TowerFPD *a) 
{
  if(!a)return false;
  return ( abs(this->col - a->col) + abs(this->row - a->row) == 1 )  ;
};

void TowerFPD::Print(void) 
{
  
  TString HiTow="";
  printf("energy=%5.2f col=%3d row=%3d cluster=%3d cluster2=%3d claimed=%3d adc=%5d",energy,col,row,cluster,cluster2, ClaimedID,adc_over_ped);
  Bool_t hi=LocalHighTower();
  if(hi)HiTow="Hi";
  printf(" IEW=%2d NSTB=%2d %s\n",IEW,NSTB,(const char*) HiTow);
}


Bool_t TowerFPD::SetContext(TObjArray* towers,CalibStr* gain,CalibStr* gaincorr,Int_t iEW,Int_t nSTB)
{
  IEW=iEW;
  NSTB=nSTB;
  CalibSet=true;
  TowerSet=towers;
  if(gain)
    {
      Gain=gain;
    }
  else CalibSet=false;
  if(gaincorr)
    {
      GainCorr=gaincorr;
    }
  else CalibSet=false;
  Live=false;
  Float_t gn=Gain->GetValue(IEW,NSTB,row-1,col-1);
  Float_t gnc=GainCorr->GetValue(IEW,NSTB,row-1,col-1);
  if(gn*gnc>0)Live=true;
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
  return  LinkDefined=true;
};

Bool_t TowerFPD::LocalHighTower()
{
 if(!LinkDefined || !CalibSet)return false;
 if(!Live)return false;
 if(!TowerFilter())return false;
  Bool_t retval=true;
  for(int j=0;j<4;j++)
    {
      if( (TowerFPD*) (Lnk_LRUD->At(j)) )
	{
	  if( ((TowerFPD*) (Lnk_LRUD->At(j)))->energy > energy)retval=false;
	};
    };
  //high tower may be a tie with neighbor
  return retval;
};

Bool_t TowerFPD::IsolatedDeadCell()
{
  if(Live)return false;
  if(!LinkDefined || !CalibSet)return false;
  Bool_t retval=true;
  for(int j=0;j<4;j++)
    {
      if((TowerFPD*) (Lnk_LRUD->At(j)))
	{
	  if(!(((TowerFPD*) (Lnk_LRUD->At(j)))->Live))retval=false;
	}
      else {retval=false;};
    };
  return retval;
};

void TowerFPD::ClearClusterInfo(TObjArray* TowerS)
{
  if(!TowerS)return;
  TIter next(TowerS);
  while(TowerFPD* tow = (TowerFPD*) next())
    {
      tow->ClaimedID=-1;	
    };
};

Bool_t TowerFPD::TowerFilter()
{
  if(energy<.002)return false;
  return true;
};
 
Int_t TowerFPD::AddToContiguous(Int_t set,Int_t hitcounter)
{
  ClaimedID=set;
  hitcounter++;
  if(!LinkDefined)return 0;
  if(hitcounter>50)return hitcounter; //break out after 50 hits in cluster

  for(int j=0;j<4;j++)
    {
      if(TowerFPD* tow=((TowerFPD*) Lnk_LRUD->At(j)))
	{

	  if(tow->Live && tow->TowerFilter() && ((tow->ClaimedID)<0))
	    {	
	      hitcounter=tow->AddToContiguous(set,hitcounter);

	    };
	};
    };
  return hitcounter;
};


Int_t TowerFPD::CreateContiguous()
{
  if(!LinkDefined)return 0;
  ClearClusterInfo(TowerSet);
  TIter next(TowerSet);
  Int_t set=0;
  while(TowerFPD* tow=(TowerFPD*) next())
    {
      if(tow->Live && tow->TowerFilter() && tow->ClaimedID==-1)
	{
	  tow->ClaimedID=set;
	  
	  tow->AddToContiguous(set,0);
	  set++;
	};
    }
  return set; // return number of Clusters
}
Bool_t TowerFPD::NeighborsAlive(Float_t dcell)
{
  if(!CalibSet)return false; // info not available so answer false
  for(int ir=-(int) dcell;ir<=(int) dcell;ir++)
    {
      for(int ic=-(int) dcell;ic<= (int) dcell;ic++)
	{
	  Float_t dist=sqrt(ic*ic+ir*ir);
	  if(dist<=dcell)
	    {
	      Int_t Row=row+ir-1;
	      Int_t Col=col+ic-1;
	      Int_t NRows=Gain->tm(IEW,NSTB)->GetNrows();
	      Int_t NCols=Gain->tm(IEW,NSTB)->GetNcols();
	      
	      if( (Row>=0)&&(Row<NRows))
	      {
		if( (Col>=0) && (Col<NCols))
		  {
		    Float_t gn=Gain->GetValue(IEW,NSTB,Row,Col);
		    Float_t gncor=GainCorr->GetValue(IEW,NSTB,Row,Col);
		    if(gn*gncor<=0)return false;
		  };
	      };
	    };
	};
    };
  return true;
};

Int_t TowerFPD::ClusterCategory(TObjArray* ta)
{
  if(!LinkDefined)return 0;
  TIter next(ta);
  Int_t cat=0;
  Float_t highEnergy=0;
  TowerFPD* hightow=0;
  while(TowerFPD* tow=(TowerFPD*) next())
    {
      if(tow->LocalHighTower())
	{
	  if(!(tow->energy==highEnergy && tow->IsNeighbor(hightow) ) )cat++;
	  
	  if(tow->energy>highEnergy)
	    {
	      highEnergy=tow->energy;
	      hightow=tow;
	    };
	};
    };
  ;
  return cat; // return category
}


TObjArray* TowerFPD::CreateClusterList(Int_t NCl)
{
  TObjArray* to=new TObjArray();

  //to is created here for use outside of this class
  //must be deleted outside too

  TIter next(TowerSet);
  Int_t count;
  Float_t esum=0;
  while(TowerFPD* tow=(TowerFPD*) next())
    {
      if(tow->ClaimedID==NCl)
	{
	  count++;
	  to->Add(tow);
	  esum+=tow->energy;
	};
      
    };
  if(((count<1) && (esum<2)) || (esum<.01))
    {
      delete to;
      return 0;
    };
  return to;
};


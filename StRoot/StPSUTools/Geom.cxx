#include "Geom.h"
#include "StFmsDbMaker/StFmsDbMaker.h"

#include "tables/St_fmsDetectorPosition_Table.h"
#include "tables/St_fmsChannelGeometry_Table.h"

using namespace std;
using namespace PSUGlobals;
ClassImp(Geom)

Geom::Geom(){
  InitDBGeom();
}
Geom::~Geom()
{

}
void Geom::printgeom(const char* file)
{
FILE* fp=fopen(file,"w");
if(fp){
	for(int i=0; i<13;i++){
	fprintf(fp,"* %s\n", (const char*)datatype[i]);
		for(int j=0; j<dataSize[i]; j++){
		(i==0) ? fprintf(fp, "%.1f", data[i][j]) :
		fprintf(fp, "%.3f", data[i][j]);
		int next = j+1;
		if(next<dataSize[i])fprintf(fp,"  ");
		}
	fprintf(fp,"\n*\n");
	}
}
else{
printf("FAILED TO OPEN FILE: %s", file);
    }
fclose(fp);
};

TVector3 Geom::LocalXYZ(Int_t EW12,Int_t NSTB16,TVector3 Globalxyz,Bool_t BlockUnit)
{
  TVector3 DirX(1.,0,0);
  TVector3 DirY(0,-1.,0);
  TVector3 DirZ(0,0,1);
  if(NSTB16==1 || (EW12==2 && NSTB16==3))DirX=-1*DirX;
  TVector3 voffset(*xOffset(EW12,NSTB16),*yOffset(EW12,NSTB16),
		  *ZFPD(EW12,NSTB16));
  TVector3 shiftxyz=Globalxyz-voffset;
  TVector3 locxyz=shiftxyz;
  locxyz(0)=shiftxyz.Dot(DirX);
  locxyz(1)=shiftxyz.Dot(DirY);
  locxyz(2)=shiftxyz.Dot(DirZ);
  Float_t scalefactorx=1./(*FpdTowWid(EW12,NSTB16));
  if(FMSGeom)
    {
      if(BlockUnit)
	{
	  Float_t scalefactory=1./(FpdTowWid(EW12,NSTB16)[1]);
	  locxyz(0)*=scalefactorx;
	  locxyz(1)*=scalefactory;
	};
    }
  else
    {
      if(BlockUnit)
	{
	  locxyz*=scalefactorx;
	};
    };
  return locxyz; 
};

TVector3 Geom::GlobalXYZ(Int_t EW12,Int_t NSTB16,TVector3 Localxyz)
{
  TVector3 DirX(1.,0,0);
  TVector3 DirY(0,-1.,0);
  TVector3 DirZ(0,0,1);
  if(NSTB16==1 || (EW12==2 && NSTB16==3))DirX=-1*DirX;
  TVector3 voffset(*xOffset(EW12,NSTB16),*yOffset(EW12,NSTB16),
		  *ZFPD(EW12,NSTB16));
  TVector3 locvec=Localxyz.X()*DirX+Localxyz.Y()*DirY+Localxyz.Z()*DirZ;
  return locvec+voffset;
};

Int_t Geom::getNSTB(Float_t gx, Float_t gy)
{
  Float_t xoff[4],yoff[4];
  for(int i=0; i<4; i++)
    {
      xoff[i] = *xOffset(2,i+1);
      yoff[i] = *yOffset(2,i+1);
    }
  Float_t xywidth12 = *FpdTowWid(2,1);
  Float_t xwidth34 = *FpdTowWid(2,3);
  Float_t ywidth34 = FpdTowWid(2,3)[1];
  if(gx>xoff[1]+xywidth12*17 || gx<xoff[0]-xywidth12*17 || gy>yoff[0] || gy<yoff[0]-xywidth12*34)
    return 0;
  else{
    if(gx<xoff[0] && (gx<xoff[0]-xywidth12*8 || gy>yoff[0]-xywidth12*9 || gy<yoff[0]-xywidth12*25))
      return 1;
    else{
      if(gx>xoff[1] && (gx>xoff[1]+xywidth12*8 || gy>yoff[1]-xywidth12*9 || gy<yoff[1]-xywidth12*25))
	return 2;
      else{
	if(gx>xoff[3]+xwidth34*12 || gx<xoff[2]-xwidth34*12 || gy>yoff[2] || gy<yoff[2]-ywidth34*24)
	  return 0;
	else{
	  if(gx<xoff[2] && (gx<xoff[2]-xwidth34*5 || gy>yoff[2]-ywidth34*5 || gy<yoff[2]-ywidth34*17))
	    return 3;
	  else{
	    if(gx>xoff[3] && (gx>xoff[3]+xwidth34*5 || gy>yoff[3]-ywidth34*5 || gy<yoff[3]-ywidth34*17))
	      return 4;
	    else return 0;
	  }
	}
      }
    }
  }
}

Float_t* Geom::ZFPD(Int_t ew, Int_t nstb)
{
  Int_t wcnt;
  if(ew==1 || ew==2)
    {
      if(nstb>0 && nstb<5)
	{
	  wcnt=(ew-1)*4+nstb-1;
	  return &(data[0][wcnt]);
	};
    };
  printf("Undefined reference in Geom\n");
  return 0;
};

Float_t* Geom::xOffset(Int_t ew, Int_t nstb)
{
  Int_t wcnt;
  if(ew==1 || ew==2)
    {
      if(nstb>0 && nstb<5)
	{
	  wcnt=(ew-1)*4+nstb-1;
	  return &(data[1][wcnt]);
	};
    };
  printf("Undefined reference in Geom\n");
  return 0;
};
Float_t* Geom::yOffset(Int_t ew, Int_t nstb)
{
  Int_t wcnt;
  if(ew==1 || ew==2)
    {
      if(nstb>0 && nstb<5)
	{
	  wcnt=(ew-1)*4+nstb-1;
	  
	  return &(data[2][wcnt]);
	};
    };
  printf("Undefined reference in Geom\n");
  return 0;
};
void Geom::Print()
{
  for(Int_t i=0;i<14;i++)
    {
      if(dataSize[i]>0)
	{
	  printf("%s : ",(const char*) datatype[i]);
	  for(Int_t j=0;j<dataSize[i];j++)
	    {
	      printf("%f ",(data[i][j]));
	    };
	  printf("\n");
	};
    };
};
Float_t* Geom::FpdTowWid(Int_t ew, Int_t nstb)
{
  Int_t wcnt;
  Int_t Nxy=1;
  if(FMSGeom)Nxy=2;
  // for FMS geometry file this routine returns an array of Float_t's 
  // first: horizontal width
  // second: vertical width
  if(ew==1 || ew==2)
    {
      if(nstb>0 && nstb<5)
	{
	  wcnt=(ew-1)*4+nstb-1;//ew==1, wcnt = {0-3} : ew==2, wcnt = {4-8}
	  return &(data[3][Nxy*wcnt]);
	};
    };
  printf("Undefined reference in Geom\n");
  return 0;
};
void Geom::InitDBGeom()
{

/*
//We don't use most of the geom information so the DB 
//call will only grab the stuff we care about.
//(ZFPD, xOffSet, yOffset,FpdTowWid)
*/

  fmsDetectorPosition_st* dbgeom;
  for(Int_t i=0;i<14;i++)
  	{
  		dataSize[i]=0;
      		datatype[i]="none";
    	};
  datatype[0]="ZFPD";
  dataSize[0]=8;
  datatype[1]="xOffset";
  dataSize[1]=8;
  datatype[2]="yOffset";
  dataSize[2]=8;
  datatype[3]="FpdTowWid";
  dataSize[3]=16;
  datatype[4]="SmdSep";
  dataSize[4]=0;
  datatype[5]="SmdOffset";
  dataSize[5]=0;
  datatype[6]="digitBits";
  dataSize[6]=0;
  datatype[7]="numebin";
  dataSize[7]=0;
  datatype[8]="ebin";
  dataSize[8]=0;
  datatype[9]="numetabin";
  dataSize[9]=0;
  datatype[10]="etabin";
  dataSize[10]=0;
  datatype[11]="absphibin";
  dataSize[11]=0;
  datatype[12]="end";
  dataSize[12]=0;

  Defined=false;
  FMSGeom=false;

  /*
  //Notes for the next for loop:
  //The idea for the for loop is that the file reads: data of [type] for detectorid (data[counter][0-16])
  //The new idea is that we have ~11 detectorids, however only 1-2 and 8-11 have anything in them...
  //this corresponds to the fpd and the fms detectors. The DB struct is index such that type becomes
  //invarient and id must shift (Reverse of the file). i.e. data[0-3][counter] where 0-3 correspond to the 4 data types.
  //Also note that TowWid is strange as it is a variation of xwidth and ywidth so we must do a 
  //data[3][0,2,4,6,8,10,12,14]=xwidth and data[3][1,3,5,7,9,11,13,15] = ywidth. 
  //This strangeness makes for some funky looking code that involves a second index updated in the for loop;
  //Also there are 11 detectorids in the DB, the code assumes we have 8. So I will make that fit by trimming
  // 2 of the useless detectorids via 2 for loops.
  */
  dbgeom=gStFmsDbMaker->DetectorPosition();
  if(dbgeom)
  	{
  		int j;
  		for(int i=0; i <= 3; i++)
  			{
  				j = 2*i;
  				data[0][i]=dbgeom[i].zoffset;//ZFPD of detector i
  				data[1][i]=dbgeom[i].xoffset;//xOffset of detector i 
  				data[2][i]=dbgeom[i].yoffset;//yOffset of detector i
  				data[3][j]=dbgeom[i].xwidth;//x tow width
  				data[3][j+1]=dbgeom[i].ywidth;//y tow width
  			}
  		for(int i=8; i <=11; i++)
  			{
  				j= (i-4)*2;
  				data[0][i-4]=dbgeom[i].zoffset;//ZFPD of detector i
  				data[1][i-4]=dbgeom[i].xoffset;//xOffset of detector i 
  				data[2][i-4]=dbgeom[i].yoffset;//yOffset of detector i
  				data[3][j]=dbgeom[i].xwidth;//x tow width
  				data[3][j+1]=dbgeom[i].ywidth;//y tow width
  			}
  		


  		FMSGeom=true;//Got the fms stuff, so we are good to go
  	}
  else{
  		LOG_ERROR << "Geom::InitDBGeom() - Failed to find Geom Info" <<endm;
  	}
}

#include <iostream>
#include <stdlib.h>
#include <math.h>


#include "TowerFPD.h"

#include "StEvent/StFmsHit.h"
#include "StFmsDbMaker/StFmsDbMaker.h"

using namespace PSUGlobals;
ClassImp(TowerFPD);

TowerFPD::TowerFPD()
    : hit(NULL), col(-1), row(-1), cluster(-1), Lnk_LRUD(NULL) { }

TowerFPD::TowerFPD(const StFmsHit* fmsHit)
    : hit(fmsHit), col(-1), row(-1), cluster(-1), Lnk_LRUD(NULL) { }

TowerFPD& TowerFPD::operator=(const TowerFPD& rhs) {
  if (this != &rhs) {
    hit = rhs.hit;
    col = rhs.col;
    row = rhs.row;
    cluster = rhs.cluster;
  }  // if
  return *this;
}

TowerFPD::~TowerFPD() {
  if (Lnk_LRUD) {
    delete Lnk_LRUD;
  }  // if
}

Bool_t TowerFPD::initialize(StFmsDbMaker* database) {
  if (!hit || !database) {  // Check for invalid input
    return false;
  }  // if
  // Get row and column from the database
  row = database->getRowNumber(hit->detectorId(), hit->channel());
  col = database->getColumnNumber(hit->detectorId(), hit->channel());
  // The database counts row number starting at the bottom, but the internals
  // of this code are set up to count from the top-down (for historical reasons)
  // so recalculate the row number. Valid FMS detector IDs are [8, 11]
  if (hit->detectorId() > 7 && hit->detectorId() < 12) {
    // Add 1 to maintain [1, N] row range rather than [0, N-1]
    row = database->nRow(hit->detectorId()) - row + 1;
  } else {  // Invalid detector ID, reset to invalid values
    row = col = -1;
  }  // if
  return row > -1 && col > -1;
}

Int_t TowerFPD::Compare(const TObject *obj) const {
  const TowerFPD* other = static_cast<const TowerFPD*>(obj);
  if (hit->energy() < other->hit->energy()) {
    return -1;
  } else if (hit->energy() > other->hit->energy()) {
    return 1;
  } else {
    return 0;
  }  // if
}

Bool_t TowerFPD::IsNeighbor(TowerFPD *a) {
  if (!a) {
    return false;
  }  // if
  return abs(this->col - a->col) + abs(this->row - a->row) == 1;
}

Bool_t TowerFPD::SetContext(TObjArray* towers) {
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

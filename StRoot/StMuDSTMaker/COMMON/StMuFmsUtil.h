/***************************************************************************
 *
 * $Id: StMuFmsUtil.h,v 1.1 2010/01/25 03:57:39 tone421 Exp $
 *
 * Author: Jingguo Ma, Jan 2010
 ***************************************************************************
 *
 * Description: FMS Util to convert between StEvent and MuDst
 *
 ***************************************************************************
 *
 * $Log: StMuFmsUtil.h,v $
 * Revision 1.1  2010/01/25 03:57:39  tone421
 * Added FMS and Roman pot arrays
 *
 **************************************************************************/
#ifndef StMuFmsUtil_h
#define StMuFmsUtil_h
#include "TObject.h"

class StMuFmsCollection;
class StFmsCollection;

class StMuFmsUtil : public TObject
{
  protected:
    
  public:
                       StMuFmsUtil();
                       ~StMuFmsUtil();
    StMuFmsCollection* getMuFms(StFmsCollection*);
    StFmsCollection*   getFms(StMuFmsCollection*);
    void               fillMuFms(StMuFmsCollection*,StFmsCollection*);
    void               fillFms(StFmsCollection*,StMuFmsCollection*);

 protected:
  /** Create StMuFmsHits from StFmsHits and fill StMuFmsCollection */
  void fillMuFmsHits(StMuFmsCollection*, StFmsCollection*);
  /** Create StMuFmsClusters from StFmsClusters and fill StMuFmsCollection */
  void fillMuFmsClusters(StMuFmsCollection*, StFmsCollection*);
  /** Create StMuFmsPoints from StFmsPoints and fill StMuFmsCollection */
  void fillMuFmsPoints(StMuFmsCollection*, StFmsCollection*);
  /** Create StFmsHits from StMuFmsHits and fill StFmsCollection */
  void fillFmsHits(StFmsCollection*, StMuFmsCollection*);
  /** Create StFmsClusters from StMuFmsClusters and fill StFmsCollection */
  void fillFmsClusters(StFmsCollection*, StMuFmsCollection*);
  /** Create StFmsPoints from StMuFmsPoints and fill StFmsCollection */
  void fillFmsPoints(StFmsCollection*, StMuFmsCollection*);
  ClassDef(StMuFmsUtil,1)
};

#endif

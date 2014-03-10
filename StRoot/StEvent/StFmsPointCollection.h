//container of StFmsPoint
//Yuxi Pan --04/01/2013
#ifndef StFmsPointCollection_HH
#define StFmsPointCollection_HH

#include "StObject.h"
#include "TDataSet.h"

#include "StFmsPoint.h"

class StFmsPointCollection : public TDataSet {

public:
        StFmsPointCollection();
        ~StFmsPointCollection();
        void Clear(const char* opt="");
        void Print(Option_t *option="") const;

        void AddPoint(StFmsPoint*);
        unsigned int NumberOfPoints() const;

        StPtrVecFmsPoint&             points();
        const StPtrVecFmsPoint&       points() const;

private:
        StPtrVecFmsPoint mPoints; //StFmsPointCollection owns the points

        ClassDef(StFmsPointCollection,1)
};
#endif

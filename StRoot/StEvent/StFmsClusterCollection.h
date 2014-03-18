//container of StFmsCluster
//Yuxi Pan --03/31/2013
#ifndef StFmsClusterCollection_HH
#define StFmsClusterCollection_HH

#include "StObject.h"
#include "TDataSet.h"

#include "StEvent/StContainers.h"
#include "StFmsCluster.h"

class StFmsClusterCollection : public TDataSet {

public:
	StFmsClusterCollection();
	~StFmsClusterCollection();
	void Clear(const char* opt="");
	void Print(Option_t *option="") const;
	
	void AddCluster(StFmsCluster*);
	unsigned int NumberOfClusters() const;
	
	StPtrVecFmsCluster&		clusters();
	const StPtrVecFmsCluster&	clusters() const;

private:
	StPtrVecFmsCluster mClusters; //StFmsClusterCollection owns the clusters
	
	ClassDef(StFmsClusterCollection,1)
};
#endif	

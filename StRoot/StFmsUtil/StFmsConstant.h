//stores all the detector/reconstruction related constants

//the following numbers are used by StFmsPointMaker::isValidChannel()
//to define detector boundaries
//row/column lower limit
#define ROW_LOW_LIMIT 1
#define COL_LOW_LIMIT 1
//central hole in large cells
#define CEN_ROW_LRG 17.5
#define CEN_ROW_WIDTH_LRG 8
#define CEN_UPPER_COL_LRG 9
//central hole in small cells
#define CEN_ROW_SML 12.5
#define CEN_ROW_WIDTH_SML 5
#define CEN_UPPER_COL_SML 6
//cuts off 7x7 triangle from the corners
#define CORNER_ROW 17.5
#define CORNER_LOW_COL 27.0

//the following are used for cluster/point id encodings
//in case the parent cluster of a point cannot be accessed
//by TRef
#define CLUSTER_BASE 305
#define CLUSTER_ID_FACTOR_DET 20

//the following are used by StFmsEventClusterer
#define TOTAL_TOWERS 578
#define PEAK_TOWER_FACTOR 1.6
#define TOWER_E_THRESHOLD 0.01
#define BAD_2PH_CHI2 10.0
#define BAD_MIN_E_LRG  0.75
#define BAD_MAX_TOW_LRG 25
#define BAD_MIN_E_SML  2.0
#define BAD_MAX_TOW_SML 49
#define VALID_FT 0.25
#define VALID_2ND_FT 1.5
#define VALID_E_OWN 0.2


//these are parameters of the shower shape
#define SS_C 0.0
#define SS_A1 1.070804
#define SS_A2 0.167773
#define SS_A3 -0.238578
#define SS_B1 0.535845
#define SS_B2 0.850233
#define SS_B3 2.382637

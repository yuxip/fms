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

//empirical cuts used by StFmsClusterFinder::categorise()
#define CAT_NTOWERS_PH1 5
#define CAT_EP1_PH2 2.1
#define CAT_EP0_PH2 7.0
#define CAT_SIGMAMAX_MIN_PH2 35.0
#define CAT_EP1_PH1 2.1
#define CAT_EP0_PH1 2.0
#define CAT_SIGMAMAX_MAX_PH1 10.0


//the following numbers are used by StFmsEventClusterer::OnePhotonFitParameters
#define PH1_START_NPH 1.0
#define PH1_DELTA_N 0.5
#define PH1_DELTA_X 0.5
#define PH1_DELTA_Y 0.5
#define PH1_DELTA_E 0.15

//these are used by StFmsEventClusterer::TwoPhotonFitParameters
#define PH2_START_NPH 2
#define PH2_START_FSIGMAMAX 2.2
#define PH2_RAN_LOW -0.1
#define PH2_RAN_HIGH 0.1
#define PH2_STEP_0 0
#define PH2_STEP_1 0.02
#define PH2_STEP_2 0.02
#define PH2_STEP_3 0.01
#define PH2_STEP_4 0.01
#define PH2_STEP_5 0.01
#define PH2_STEP_6 0.1
#define PH2_MAXTHETA_F 2.8
#define PH2_LOWER_NPH 1.5
#define PH2_LOWER_XF 0.2
#define PH2_LOWER_YF 0.2
#define PH2_LOWER_XMAX_F 18.0
#define PH2_LOWER_XMAX_POW 0.8
#define PH2_LOWER_XMAX_LIMIT 0.5
#define PH2_LOWER_5_F -1.0
#define PH2_LOWER_6_F 0.95
#define PH2_UPPER_NPH 2.5
#define PH2_UPPER_XF 0.2
#define PH2_UPPER_YF 0.2
#define PH2_UPPER_XMIN_F 0.085
#define PH2_UPPER_XMIN_P0 60.0
#define PH2_UPPER_XMIN_LIMIT 3.5
#define PH2_UPPER_5_F 1.0
#define PH2_UPPER_6_F 1.05
#define PH2_3_LIMIT_LOWER 0.9
#define PH2_3_LIMIT_UPPER 1.1

//StFmsEventClusterer::GlobalPhotonFitParameters
#define GL_LOWER_1 0.5
#define GL_UPPER_DELTA_MAXN 0.5
#define GL_0_DLOWER 1.25
#define GL_0_DUPPER 1.25
#define GL_1_DLOWER 1.25
#define GL_1_DUPPER 1.25
#define GL_2_DLOWER 0.3
#define GL_2_DUPPER 0.3


#ifndef Geom_
#define Geom_

#include <map>
#include <vector>

#include <TObject.h>

class fmsDetectorPosition_st;
class StFmsDbMaker;

namespace PSUGlobals {//$NMSPC
/**
 Enumeration for FPD and FMS subdetectors.
 
 The numerical values are defined for compatibility with the detector ID
 values expected by StFmsHit i.e. these are the detector ID numbers defined by
 StFmsChannelGeometry in the database.
 */
enum StFmsDetectorId {
  // FPD: forward pion detector, with preshower and shower-maximum detector
  kFpdNorth = 0,
  kFpdSouth = 1,
  kFpdNorthPreshower = 2,
  kFpdSouthPreshower = 3,
  kFpdNorthShowerMaxVertical = 4,
  kFpdSouthShowerMaxVertical = 5,
  kFpdNorthShowerMaxHorizontal = 6,
  kFpdSouthShowerMaxHorizontal = 7,
  // FMS: forward meson spectrometer, with four subregions
  kFmsNorthLarge = 8,
  kFmsSouthLarge = 9,
  kFmsNorthSmall = 10,
  kFmsSouthSmall = 11,
  // FHC: forward hadronic calorimeter (future detector)
  kFhcNorth = 12,
  kFhcSouth = 13,
  kFmsInvalidDetectorId = -1
};
/**
 Wrapper around FMS detector geometry database information
 
 Provides convenience functions to access fmsDetectorPosition_st and its
 members via detector ID
 */
class Geom : public TObject {
 public:
  /** Constructor */
  Geom();
  /** Destructor */
  ~Geom();
  /**
   Initialise geometry from the FMS database
   
   If the argument is NULL, attempt to locate an StFmsDbMaker in the current
   chain and use that.
   Return true if the geometry is initialised, false if it is not.
   */
  Bool_t initialize(StFmsDbMaker* fmsDbMaker);
  /**
   Return the z position of a detector in cm
  
   Gives the database value fmsDetectorPosition_st::zoffset
   Returns 0 if the argument is an invalid detector ID or the geometry is
   uninitialized
   */
  Float_t z(Int_t detectorId) const;
  /**
   Return the x coordinate offset of a detector in cm
  
   Gives the database value fmsDetectorPosition_st::xoffset
   Returns 0 if the argument is an invalid detector ID or the geometry is
   uninitialized
  */
  Float_t xOffset(Int_t detectorId) const;
  /**
   Return the y coordinate offset of a detector in cm
  
   Gives the database value fmsDetectorPosition_st::yoffset
   Returns 0 if the argument is an invalid detector ID or the geometry is
   uninitialized
  */
  Float_t yOffset(Int_t detectorId) const;
  /**
   Return [x, y] tower widths in cm
   
   Database values fmsDetectorPosition_st::xwidth and ::ywidth
   Returns [0, 0] if the argument is an invalid detector ID or the geometry is
   uninitialized
   */
  std::vector<Float_t> towerWidths(Int_t detectorId) const;
  /**
   Return the position information of a detector
   
   Returns NULL if the argument is an invalid detector ID or the geometry is
   uninitialized
   */
  const fmsDetectorPosition_st* find(Int_t detectorId) const;

 private:
  typedef std::map<int, fmsDetectorPosition_st*> Table;
  Table mPositions;  ///< Detector ID: position information pairs
  ClassDef(Geom, 3)
};  // class Geom
}  // namespace PSUGlobals
#endif

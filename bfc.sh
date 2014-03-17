#!/bin/sh

# Run BFC macro
# Chain options are SL14a with ry2013_2 geometry
# http://www.star.bnl.gov/devcgi/dbProdOptionRetrv.pl
root4star -q -b bfc.C'(1,"DbV20140222 pp2013a mtd btof fmsDb fmsDat fgt fgtPoint VFPPVnoCTB beamline BEmcChkStat Corr4 OSpaceZ2 OGridLeak3D -hitfilt","st_physics_14158051_raw_1030004.daq")'
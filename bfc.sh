#!/bin/sh

# Run BFC macro
# Chain options are SL14a with ry2013_2 geometry
# http://www.star.bnl.gov/devcgi/dbProdOptionRetrv.pl
#root4star -q -b bfc.C'(5,"DbV20130608 pp2013a fmsDb fmsDat fmsPoint","st_physics_14158051_raw_1030004.daq",0,0,"bfcQa.root")'
root4star -q -b bfc.C'(5,"sdt20110923 pp2011a fmsDb fmsDat fmsPoint","st_fms_12098008_raw_2150001.daq",0,0,"bfcQa.root")'

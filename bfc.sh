#!/bin/sh

# Run BFC macro
# Chain options are SL14a with ry2013_2 geometry
# http://www.star.bnl.gov/devcgi/dbProdOptionRetrv.pl
root4star -q -b bfc.C'(100,"DbV20140222 pp2013a fmsDb fmsDat fmsPoint","st_physics_14158051_raw_1030004.daq",0,0,"bfcQa.root")'

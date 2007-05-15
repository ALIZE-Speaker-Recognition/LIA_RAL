#!/bin/bash

# This will test your TrainWorld.exe app in (LIA_RAL_TOPDIR)/LIA_SpkDet/TrainWorld
# compare the resulting model "wld" with wmd.validate (models are seaved under XML format)
../TrainWorld.exe --config TrainWorld1.cfg
../TrainWorld.exe --config TrainWorld2.cfg

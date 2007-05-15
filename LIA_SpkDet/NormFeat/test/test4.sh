#!/bin/bash

# This will test your NormFeat.exe app in (LIA_RAL_TOPDIR)/LIA_SpkDet/NormFeat
# compare the resulting test1.norm.prm with test1.validate.prm or simply check the output so that each mean and covariance in each dimension are equal to 0 and 1 respectively
# this test checks that stats can be loaded from external file
../NormFeat.exe --config NormFeat.cfg --mode info
../NormFeat.exe --config NormFeat.cfg --externalStatsFilename test1.stat
../NormFeat.exe --config NormFeat.cfg --loadFeatureFileExtension .norm.prm --mode info

#!/bin/bash

# This will test your NormFeat.exe app in (LIA_RAL_TOPDIR)/LIA_SpkDet/NormFeat
# compare the resulting test1.norm.prm with test1.validate.prm or simply check the output so that each mean and covariance in each dimension are equal to 0 and 1 respectively
../NormFeat.exe --config NormFeat.cfg 
../NormFeat.exe --config NormFeat.cfg --gaussianize 
../NormFeat.exe --config NormFeat.cfg --inputFeatureFileName test1.norm --mode info
../NormFeat.exe --config FeatMap.cfg

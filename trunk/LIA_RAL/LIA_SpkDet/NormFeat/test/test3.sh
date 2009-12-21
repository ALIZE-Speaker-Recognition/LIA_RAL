#!/bin/bash

# This will test your NormFeat.exe app in (LIA_RAL_TOPDIR)/LIA_SpkDet/NormFeat
# In this example 6 output files are generated
# .0. is the real data, normalized
# .1. to .5. are generated using randomly moved mean
../NormFeat.exe --config NormFeat3.cfg 
../NormFeat.exe --config NormFeat3.cfg --loadFeatureFileExtension .norm.prm --mode info --inputFeatureFilename test1
../NormFeat.exe --config NormFeat3.cfg --loadFeatureFileExtension .norm.prm --mode info --inputFeatureFilename test1.1
../NormFeat.exe --config NormFeat3.cfg --loadFeatureFileExtension .norm.prm --mode info --inputFeatureFilename test1.2
../NormFeat.exe --config NormFeat3.cfg --loadFeatureFileExtension .norm.prm --mode info --inputFeatureFilename test1.3
../NormFeat.exe --config NormFeat3.cfg --loadFeatureFileExtension .norm.prm --mode info --inputFeatureFilename test1.4
../NormFeat.exe --config NormFeat3.cfg --loadFeatureFileExtension .norm.prm --mode info --inputFeatureFilename test1.5

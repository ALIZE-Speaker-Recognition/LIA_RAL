#!/bin/bash

# This will test your ComputeTest.exe app in (LIA_RAL_TOPDIR)/LIA_SpkDet/ComputeTest
# compare the resulting test1.res with test1.validate.res
../ComputeTest.exe --config ./ComputeTest.cfg
../ComputeTest.exe --config ./ComputeTestByLabel.cfg


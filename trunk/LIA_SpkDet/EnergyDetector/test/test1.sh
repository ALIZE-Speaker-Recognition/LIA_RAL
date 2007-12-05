#!/bin/bash

# This will test your EnergyDetector.exe app in (LIA_RAL_TOPDIR)/LIA_SpkDet/EnergyDetector
# compare the resulting test1.enr.lbl with test1.validate.enr.lbl
../EnergyDetector.exe --config ./EnergyDetector.cfg --debug --verbose

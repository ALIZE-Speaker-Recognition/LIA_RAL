#!/bin/bash

BIN=../ShiftedDeltaFeat
CFG=../cfg/ShiftedDeltaFeat.cfg

$BIN --config $CFG --inputFeatureFilename lid07700  | tee test.log
$BIN --config $CFG --inputFeatureFilename lid08800 --debug true >> test.log


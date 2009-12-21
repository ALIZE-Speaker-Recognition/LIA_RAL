#!/bin/bash

# compare the resulting test1.sym with test1.sym.ref
../GmmTokenizer --config GmmTokenizer.cfg --verbose true
../GmmTokenizer --config GmmTokenizer.cfg --confusionMatrix --verbose true --matrixOutputName mce_matrix.mat
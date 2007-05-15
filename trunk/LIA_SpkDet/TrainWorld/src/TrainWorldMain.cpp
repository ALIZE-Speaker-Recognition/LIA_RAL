#include <iostream>
#include "liatools.h"
#include "TrainWorld.h"

int main(int argc, char* argv[])
{
  ConfigChecker cc;
  cc.addIntegerParam("verboseLevel",false,true,"level of the berose information 0=no verbose, 1=normal, 2=more");
  cc.addStringParam("inputFeatureFilename",false, true,"feature filename or filename of a text file with the list of feature filenames");
  cc.addStringParam("inputStreamList",false, true,"filename of a text file with the filename of input streams");
  cc.addStringParam("weightStreamList",false,true,"filename of a text file with the weight of each input stream - default=equal weights");
  cc.addStringParam("outputWorldFilename",true,true,"output worldmodel filename");                            
  cc.addStringParam("inputWorldFilename",false,true,"if set, the init is based on a model get from this file, else frrom scratch");
  cc.addStringParam("saveInitModel",false,true,"if set (default), save the initial model");
  cc.addStringParam("labelSelectedFrames",true,true,"the segments with this label are used for training the worldmodel");
  cc.addFloatParam("baggedFrameProbability",true,true,"defines the % of frames taken for each iterations");
  cc.addFloatParam("baggedFrameProbabilityInit",false,true,"defines the % of frames taken BY COMPONENT for the initializing of the mixture- mandatory if init from scratch");
  cc.addIntegerParam("baggedMinimalLength",false,true,"minimum length for selected segments in bagged (default=3)");
  cc.addIntegerParam("baggedMaximalLength",false,true,"maximal length for selected segments in bagged (default=7)");
  cc.addFloatParam("initVarianceFlooring",true,true,"variance control parameters - relative to global data variance - initial value (moved during the it)");
  cc.addFloatParam("initVarianceCeiling",true,true,"variance control parameters - relative to global data variance - initial value (moved during the it)"); 
  cc.addFloatParam("finalVarianceFlooring",true,true,"variance control parameters - relative to global data variance - final value");
  cc.addFloatParam("finalVarianceCeiling",true,true,"variance control parameters - relative to global data variance - final value"); 
  cc.addIntegerParam("nbTrainIt",true,true,"number of it, the ceiling and flooring are moved and the baggedFrameProbability is used"); 
  cc.addBooleanParam("normalizeModel",false,true,"if set to true,  normalize the world (at each iteration)");
  cc.addBooleanParam("normalizeModelMeanOnly",false,true,"used only if normalizeModel is On, says if only mean parameters should be normalized"); 
  cc.addIntegerParam("normalizeModelNbIt",false,true,"used only if noramlizeModelMeanOnly is set, nb of normalization it");
  cc.addBooleanParam("use01",false,true,"if set at true, don't compute the global mean and cov but uses 0 mean and 1 cov");
  cc.addBooleanParam("componentReduction",false,true,"if set reduce the number of components at each it, selecting the best weights until targetDistribCount (default false)");
  cc.addIntegerParam("targetMixtureDistribCount",false,true,"final number of components if componentReduction is selected"); 


  try {
      CmdLine cmdLine(argc, argv);
      if (cmdLine.displayHelpRequired()){
        cout <<"TrainWorld.exe"<<endl<<"This program is used for training a world model from scratch or from a model"
             <<endl<<cc.getParamList()<<endl;
        return 0;  
      }
      if (cmdLine.displayVersionRequired()){
        cout <<"Version 2-beta"<<endl;
      } 
      Config tmp;
      cmdLine.copyIntoConfig(tmp);
      Config config;
      if (tmp.existsParam("config")) config.load(tmp.getParam("config"));
      cmdLine.copyIntoConfig(config);
      cc.check(config);
      debug=config.getParam_debug();
      if (config.existsParam("verbose"))verbose=config.getParam("verbose").toBool();else verbose=false;
      if (verbose) verboseLevel=1;else verboseLevel=0;
      if (config.existsParam("verboseLevel"))verboseLevel=config.getParam("verboseLevel").toLong();
      if (verboseLevel>0) verbose=true;
      trainWorld(config);
    }
catch (alize::Exception& e) {cout << e.toString() << endl << cc.getParamList()<< endl;}
return 0;     
}

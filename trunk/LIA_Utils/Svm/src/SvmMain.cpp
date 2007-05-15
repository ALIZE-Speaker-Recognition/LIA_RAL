#include <iostream>

#include "Svm.h"
#include "alize.h"
#include "liatools.h"
int main(int argc, char* argv[])
{
  using namespace std;
  using namespace alize;
  ConfigChecker cc;
  cc.addStringParam("mode",true,true,"train | predict | tnormpredict");
  cc.addStringParam("inputFilename",true, true,"train: positive instance | predict: test instance");
  cc.addStringParam("inputBCKList",false, true,"negative exemples to use every time world.vect");  
  cc.addStringParam("inputSVMModel",false, true,"in test mode");  
  cc.addStringParam("outputFilename",false, true,"train: model file name | predict: results");  
  cc.addStringParam("vectorFilesPath",true, true,"Path to instances");  
  cc.addStringParam("modelFilesPath",true, true,"Path to models");    
  cc.addStringParam("vectorFilesExtension",true, false,"ext instances");  
  cc.addStringParam("modelFilesExtension",true, false,"ext models");      
  cc.addIntegerParam("vsize",false, true,"vector sizes");
  cc.addFloatParam("targetPenalty",false, true,"penalty for error on target speaker class");  
  cc.addFloatParam("C",false, true,"C parameter");      
 
  try {
        CmdLine cmdLine(argc, argv);
        if (cmdLine.displayHelpRequired()){
          cout <<"Svm.exe"<<endl<<"Binding to libSvm for Speaker Verification in large Scale Evaluation as NIST SRE" <<endl<<cc.getParamList()<<endl;
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
        debug=config.getParam_debug();
        cc.check(config);
        if (config.existsParam("verbose"))verbose=config.getParam("verbose").toBool();else verbose=false;
        if (verbose) verboseLevel=1;else verboseLevel=0;
        if (config.existsParam("verboseLevel"))verboseLevel=config.getParam("verboseLevel").toLong();
        if (verboseLevel>0) verbose=true;       
          
        String mode=config.getParam("mode");
        if (verbose) cout << "(Svm) Mode " << mode << endl;
        if (mode=="train") svmTrain(config);
        else if (mode=="predict") svmPredict(config);
        else if(mode=="tnormpredict") svmPredictTnorm(config);
        else throw Exception("No mode",__FILE__,__LINE__);
        }
  catch (alize::Exception & e) {cout << e.toString () << endl << cc.getParamList() << endl;}
  return 0;
}

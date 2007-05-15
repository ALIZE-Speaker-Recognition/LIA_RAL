#include <iostream>

#include "ReadModel.h"

int main(int argc, char* argv[])
{
  using namespace std;
  using namespace alize;
  ConfigChecker cc;
  cc.addStringParam("inputModelFilename",true, true,"model filename in relative path with extension eg: ./gmm/world.gmm");
  cc.addBooleanParam("outputWeightVector",false,true,"if set to true, output weight vectors to STDOUT");
  cc.addStringParam("outputWeightVectorFormat",false,true,"LIA or SVMLight");
  try {
        CmdLine cmdLine(argc, argv);
        if (cmdLine.displayHelpRequired()){
          cout <<"ReadModel.exe"<<endl<<"This program is used to read ALIZE model file to STDOUT" <<endl<<cc.getParamList()<<endl;
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
        config.setParam("loadMixtureFileFormat","RAW");
        if (config.existsParam("outputWeightVector")) {outputWeightVector(config);}
        else {readModel(config);}
        }
  catch (alize::Exception & e) {cout << e.toString () << endl << cc.getParamList() << endl;}
  return 0;
}

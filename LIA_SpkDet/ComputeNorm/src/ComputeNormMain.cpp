#include <iostream>

#include "ComputeNorm.h"
#include <liatools.h>


int main(int argc, char* argv[])
{
    using namespace std;
    using namespace alize;
  ConfigChecker cc;
  cc.addStringParam("testNistFile",true, true,"target_seg score file");
  cc.addStringParam("normType",true, true,"tnorm|znorm|ztnorm");
  cc.addStringParam("tnormNistFile",false,true,"imp_seg score file");
  cc.addStringParam("znormNistFile",false,true,"target_imp score file");
  cc.addStringParam("ztnormNistFile",false,true,"imp_imp score file");
  try {
      CmdLine cmdLine(argc, argv);
      if (cmdLine.displayHelpRequired()){
	cout << "************************************" << endl;
	cout << "********** ComputeNorm.exe *************" << endl;
	cout << "************************************" << endl;
	cout << endl;
	cout << "Apply Z-T-ZT Norm to a Score File" << endl;
        cout <<endl<<cc.getParamList()<<endl;
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
      ComputeNorm(config);
    }
catch (alize::Exception& e) {cout << e.toString() << endl << cc.getParamList()<< endl;}
return 0;     
}



#include <iostream>

#include "BNgram.h"


int main(int argc, char* argv[]){
  try{

    ConfigChecker cc;
    cc.addIntegerParam("maxOrder",true,true,"maximum order of the input ngrams (from order 1 to maxOrder)");

    
    CmdLine cmdLine(argc, argv);
    if (cmdLine.displayHelpRequired()){
      cout <<"BNgram.exe"<<endl<<"This program is a toolset for ngram computation"<<endl<<cc.getParamList()<<endl;
	return 0;
      }
      if (cmdLine.displayVersionRequired()){
	cout <<"Version alpha"<<endl;
      }
      Config tmp;
      cmdLine.copyIntoConfig(tmp);
      Config config;
      if (tmp.existsParam("config")) config.load(tmp.getParam("config"));
      cmdLine.copyIntoConfig(config);
      cc.check(config);
      debug=config.getParam_debug();
      if (config.existsParam("verbose")) verbose=config.getParam("verbose").toBool();
      else verbose=false;
     testBNgram(config);
    }
    catch (alize::Exception& e){
      cout << e.toString() << endl;
    }
  return 0;
}

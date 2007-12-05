#include <iostream>
#include "ReadFeatFile.h"
#include "alize.h"


int main (int argc, char *argv[])
{
  using namespace std;
  using namespace alize;

  ConfigChecker cc;
  cc.addStringParam("inputFile",true, true,"input feature file name in relative path: ./foo.prm (./ mandatory to avoid loadFeatureFile[Path|Extension] options)");
  cc.addStringParam("loadFeatureFileFormat",false,true,"ALIZE option, see ALIZE doc");
  try {
        CmdLine cmdLine(argc, argv);
        if (cmdLine.displayHelpRequired()){
          cout <<"ReadFeatFile.exe"<<endl<<"This program is used to read feature files and print their contents to STDOUT" <<endl<<cc.getParamList()<<endl;
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
        readFeatFile(config);
        }
  catch (alize::Exception & e) {cout << e.toString () << endl << cc.getParamList() << endl;}
  return 0;
}

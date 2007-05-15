#include <iostream>
#include "Hist.h"

int main (int argc, char *argv[])
{
  using namespace std;
  using namespace alize;
  ConfigChecker cc;
  cc.addStringParam("dataFile",true, true,"input data file: text file with a single column");
  cc.addStringParam("outFile",true, true,"output histo file");
  cc.addStringParam("format",false,true,"could be gnuplot(default) or txt");
  cc.addIntegerParam("nbBins",true,true,"could be gnuplot(default) or txt");  
  try {
        CmdLine cmdLine(argc, argv);
        if (cmdLine.displayHelpRequired()){
          cout << "Hist: Build an histogram on a txt data file (1 column)" << endl;
          cout << "To plot foo.hist : $ gnuplot> plot \"foo.hist\" with lines" << endl;
          cout <<cc.getParamList()<<endl;
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
        Hist(config);
        }
  catch (alize::Exception & e) {cout << e.toString () << endl << cc.getParamList() << endl;}
  return 0;
}

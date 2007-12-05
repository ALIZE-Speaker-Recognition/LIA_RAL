#include <iostream>
#include "alize.h"
#include "liatools.h"
#include "PolyExpand.h"
using namespace std;
using namespace alize;


int main(int argc, char* argv[]) {
	ConfigChecker cc;
try {
	cc.addStringParam("inputFeatureFilename",true,true,"List of features");
	cc.addStringParam("computeR",false,true,"accumulateStatMode");
	cc.addStringParam("normalize",false,true,"the R matrix");
	cc.addStringParam("vectorFilesPath",true,true,"path to instance");
	cc.addStringParam("vectorFilesExtension",true,true,"ext of vectors");
	cc.addStringParam("format",true,true,"outputFileFormat SVMLight");
	cc.addStringParam("exType",false,true,"1/0/-1 defines if positive or negative exemple");
	CmdLine cmdLine(argc, argv);
	if (cmdLine.displayHelpRequired()){
	cout <<"PolyExpand.exe"<<endl<<"This program is used for expanding features into polynomial space of order 3"
	   <<endl<<cc.getParamList()<<endl;
	return 0;
	}
	Config tmp;
	cmdLine.copyIntoConfig(tmp);
	Config config;
	if (tmp.existsParam("config")) config.load(tmp.getParam("config"));
	cmdLine.copyIntoConfig(config);
	debug=config.getParam_debug();
	if (cmdLine.displayHelpRequired()) {cout << cc.getParamList() << endl;}
	cc.check(config);
	PolyExpand(config);
	}
catch (Exception& e) {cout << e.toString() << cc.getParamList() << endl;}
return 0;
}

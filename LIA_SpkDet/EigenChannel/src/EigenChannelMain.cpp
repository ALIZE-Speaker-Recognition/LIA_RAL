#include <iostream>
#include <liatools.h>
#include <EigenChannel.h>

using namespace alize;
using namespace std;

int main(int argc, char* argv[]) {
	ConfigChecker cc;
	try {
		cc.addStringParam("ndxFilename",true,true,"NDX of multiple GMM speaker recordings");
		cc.addStringParam("inputWorldFilename",true,true,"the world model file");
		cc.addIntegerParam("nbIt",true,true,"number of ml it");	
		cc.addStringParam("channelMatrix",true,true,"filename to save Channel Matrix ");		
		cc.addStringParam("initChannelMatrix",false,true,"init Channel Matrix");	
		cc.addBooleanParam("loadAccs",false,true,"if true do not compute UBM stats, load matrices");		
		cc.addIntegerParam("computeLLK",false,true,"optional: nb of files where LLK is computed");				
		cc.addIntegerParam("channelMatrixRank",true,true,"final rank of channel matrix");	
		cc.addFloatParam("regulationFactor",true,true,"map tau");
		cc.addStringParam("saveMatrixFormat",true,true,"matrix format: DB (binary) or DT (ascii)");		  
		cc.addStringParam("loadMatrixFormat",true,true,"matrix format: DB (binary) or DT (ascii)");				
		
		CmdLine cmdLine(argc, argv);
		Config tmp;
		cmdLine.copyIntoConfig (tmp);
		Config config;
		if (tmp.existsParam("config")) config.load(tmp.getParam("config"));
		cmdLine.copyIntoConfig(config);
		cc.check(config);
		debug=config.getParam_debug();	
		if (config.existsParam("verbose"))verbose=config.getParam("verbose").toBool();else verbose=false;
		if (verbose) verboseLevel=1;else verboseLevel=0;
		if (config.existsParam("verboseLevel"))verboseLevel=config.getParam("verboseLevel").toLong();
		if (verboseLevel>0) verbose=true;		
		if (cmdLine.displayHelpRequired()) {cout << cc.getParamList() << endl;}	
		EigenChannel(config);	
		}
	catch (Exception& e) {cout << e.toString() << cc.getParamList() << endl;}
if (debug) {
}
return 0;
}

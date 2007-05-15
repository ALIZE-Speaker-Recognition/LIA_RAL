#include <iostream>
#include "Scoring.h"

int main(int argc, char* argv[])
{
using namespace std;
using namespace alize;

	ConfigChecker cc;
	cc.addStringParam("mode",true,true,"NIST,ETF,MDTM,leaveMaxOutTnorm");
	cc.addStringParam("inputFile",true,true,"fooIn.nist");
	cc.addStringParam("outputFile",true,true,"fooOut.nist");
	cc.addStringParam("adaptationMode",false,true,"u|n");
	cc.addStringParam("segTypeTest",false,true,"<10sec|...|1side>");
	cc.addStringParam("trainTypeTest",false,true,"<10sec|...|1side>");
	cc.addIntegerParam("threshold",false,true,"threshold on scores to accept or reject test");
	cc.addIntegerParam("hard",false,false,"in ETF accept one speaker per test (identification)");
	try {
		CmdLine cmdLine(argc, argv);
		if (cmdLine.displayVersionRequired()){cout <<"Version 2-beta"<<endl;} 
		if (cmdLine.displayHelpRequired()){
			cout << "************************************" << endl;
			cout << "********** Scoring.exe *************" << endl;
			cout << "************************************" << endl;
			cout << endl;
			cout << "Apply a decision on scores according to a thershold and Format test Files from ComputeTest into either NIST04 format, ETF format or MDTM format" << endl;
			cout <<endl<<cc.getParamList()<<endl;
		return 0;  
		}
		else {
			Config tmp;
			cmdLine.copyIntoConfig(tmp);
			Config config;
			if (tmp.existsParam("config")) config.load(tmp.getParam("config"));
			cmdLine.copyIntoConfig(config);
			cc.check(config);
			config.setParam("minLLK","-200");
			if (config.existsParam("warpScores")) {WarpScores(config);}
			else {Scoring(config);}
		}  
	}
	catch (alize::Exception& e) {cout << e.toString() << endl << cc.getParamList() << endl;}
return 0;
}

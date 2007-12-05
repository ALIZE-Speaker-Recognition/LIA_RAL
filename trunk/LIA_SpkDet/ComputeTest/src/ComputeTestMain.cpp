#include <iostream>
#include "liatools.h"
#include "ComputeTest.h"

int main(int argc, char* argv[]) {
	using namespace std;
	using namespace alize;

		ConfigChecker cc;
		cc.addStringParam("config", false, true, "default config filename");
		cc.addStringParam("ndxFilename",true,true, "NDX file listing all verification tests to achieve: first column: test File, all others: models");	
		cc.addStringParam("outputFilename",true,true, "output scores in this file: 'gender - test model - scores'");	
		cc.addStringParam("inputWorldFilename",true,true,"model repsresenting the H0 hypothesis to get the LLR");
		cc.addStringParam("labelSelectedFrames",true,true,"the segments with this label are used for training the worldmodel");		
		cc.addStringParam("computeLLKWithTopDistribs",true,true, "PARTIAL/COMPLETE: will compute LLK with topdistribs. COMPLETE: add world LLK to client LLK, PARTIAL: just keeps the topDistrib LLK");
		cc.addIntegerParam("topDistribsCount ",false,true,"Number of distrib to approximate complete LLK");
		cc.addStringParam("gender",true,true, "M/F: will output the gender in score file");	
		cc.addStringParam("fileLLR",false,true, "will output a score for the entire file (default true)");
		cc.addStringParam("segmentLLR",false,true, "will output a score for each speech segment (default false)");
		cc.addStringParam("windowLLR",false,true, "windowLLR: will output a score for each set of windowLLRSize frames. windowLLRDec gives the shift of the window (default false)");
		cc.addIntegerParam("windowLLRSize",false,true, "if windowLLR is set, gives the size of the window (default 30)");
		cc.addIntegerParam("windowLLRDec",false,true, "if windowLLR is set, gives the shift of the window (default windowLLRSize)");
		cc.addBooleanParam("byLabelModel",false,true, "if the parameter is present, we work with a model by client and by  cluster (default false)");
		cc.addBooleanParam("histoMode",false,true, "if the parameter is present, entropy of LR distrib is used as score (default false)");
		cc.addStringParam("channelCompensation",false,true,"if set to NAP or FA, launch the corresponding channel compensation (default not set)");
	try {
		CmdLine cmdLine(argc, argv);
		if (cmdLine.displayHelpRequired()){
		cout << "************************************" << endl;
		cout << "********** ComputeTest.exe **********" << endl;
		cout << "************************************" << endl;
		cout << endl;
		cout << "LLR computation for an NDX (NIST format) File" << endl;
		cout << "" << endl;
		cout << cc.getParamList()<<endl;
		return 0;  
		}
		if (cmdLine.displayVersionRequired()){cout <<"ComputeTest - for computing the LLR using nixt style ndx files"<<endl;} 
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
		if (config.existsParam("byLabelModel"))        // if the parameter is present, we work with a model by client and by  cluster 
		ComputeTestByLabel(config);
		if (config.existsParam("histoMode"))        // if the parameter is present, entropy of LR distrib is used as score
		ComputeTestHisto(config);
		if (config.existsParam("channelCompensation")) {
			if (config.getParam("channelCompensation")=="FA") 
				ComputeTestFA(config);
			else
				ComputeTestNAP(config);			
		}
		else 
		ComputeTest(config);
	}
catch (alize::Exception& e) {cout << e.toString() << endl << cc.getParamList()<< endl;}
return 0;
}

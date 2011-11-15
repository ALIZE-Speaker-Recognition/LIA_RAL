/*
This file is part of LIA_RAL which is a set of software based on ALIZE
toolkit for speaker recognition. ALIZE toolkit is required to use LIA_RAL.

LIA_RAL project is a development project was initiated by the computer
science laboratory of Avignon / France (Laboratoire Informatique d'Avignon -
LIA) [http://lia.univ-avignon.fr <http://lia.univ-avignon.fr/>]. Then it
was supported by two national projects of the French Research Ministry:
	- TECHNOLANGUE program [http://www.technolangue.net]
	- MISTRAL program [http://mistral.univ-avignon.fr]

LIA_RAL is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation, either version 3 of
the License, or any later version.

LIA_RAL is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with LIA_RAL.
If not, see [http://www.gnu.org/licenses/].

The LIA team as well as the LIA_RAL project team wants to highlight the
limits of voice authentication in a forensic context.
The "Person Authentification by Voice: A Need of Caution" paper
proposes a good overview of this point (cf. "Person
Authentification by Voice: A Need of Caution", Bonastre J.F.,
Bimbot F., Boe L.J., Campbell J.P., Douglas D.A., Magrin-
chagnolleau I., Eurospeech 2003, Genova].
The conclusion of the paper of the paper is proposed bellow:
[Currently, it is not possible to completely determine whether the
similarity between two recordings is due to the speaker or to other
factors, especially when: (a) the speaker does not cooperate, (b) there
is no control over recording equipment, (c) recording conditions are not
known, (d) one does not know whether the voice was disguised and, to a
lesser extent, (e) the linguistic content of the message is not
controlled. Caution and judgment must be exercised when applying speaker
recognition techniques, whether human or automatic, to account for these
uncontrolled factors. Under more constrained or calibrated situations,
or as an aid for investigative purposes, judicious application of these
techniques may be suitable, provided they are not considered as infallible.
At the present time, there is no scientific process that enables one to
uniquely characterize a persones voice or to identify with absolute
certainty an individual from his or her voice.]

Copyright (C) 2004-2010
Laboratoire d'informatique d'Avignon [http://lia.univ-avignon.fr]
LIA_RAL admin [alize@univ-avignon.fr]
Jean-Francois Bonastre [jean-francois.bonastre@univ-avignon.fr]
*/

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
		cc.addStringParam("channelCompensation",true,true,"if set to NAP or JFA, launch the corresponding channel compensation (default not set)");
	

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
		if (config.getParam("channelCompensation")=="JFA"){ 
			if(config.getParam("Scoring") == "FrameByFrame")		ComputeTestJFA(config);
			else if(config.getParam("Scoring") == "DotProduct")		ComputeTestDotProduct(config);
		}
		else if (config.getParam("channelCompensation")=="LFA")		ComputeTestLFA(config);
		else if(config.getParam("channelCompensation")=="NAP")
			ComputeTestNAP(config);
		else
			ComputeTest(config);
	}
catch (alize::Exception& e) {cout << e.toString() << endl << cc.getParamList()<< endl;}
return 0;
}

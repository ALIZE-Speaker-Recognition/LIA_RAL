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
#include "EigenChannel.h"

using namespace std;
using namespace alize;


int main(int argc, char* argv[]) {
	ConfigChecker initCc,cc;
	try{
		// Insertion of config compatibility rules
		CmdLine cmdLine(argc, argv);

		//Add list of parameters always mandatory 
		initCc.addStringParam("eigenChannelMode",true,true,"compute EigenChannel matrix for JFA or LFA paradigm");
		initCc.addStringParam("ndxFilename",true,true,"NDX of multiple GMM speaker recordings");
		initCc.addStringParam("inputWorldFilename",true,true,"the world model file");
		initCc.addIntegerParam("nbIt",true,true,"number of ml it");
		initCc.addBooleanParam("loadInitChannelMatrix",true,true,"save initialisation EigenChannel Matrix");	
		initCc.addStringParam("eigenChannelMatrix",true,true,"filename to save EigenChannel Matrix ");					
		initCc.addStringParam("saveMatrixFormat",true,true,"matrix format: DB (binary) or DT (ascii)");		  
		initCc.addStringParam("saveMatrixFilesExtension",true,true,"extension of the matrices");	
		initCc.addStringParam("loadMatrixFilesExtension",true,true,"extension of the matrices for loading");
		initCc.addStringParam("matrixFilesPath",true,true,"directory to store the matrices");
		initCc.addStringParam("randomInitLaw",true,true,"random law to initialize the matrix, could be uniform or normal");	
		initCc.addStringParam("loadMatrixFormat",true,true,"matrix format: DB (binary) or DT (ascii)");	

		// Check existing parameters to create the appropriate ConfigChecker
		Config tmpConfig;
		cmdLine.copyIntoConfig(tmpConfig);
		Config initConfig;
		if (tmpConfig.existsParam("config")) initConfig.load(tmpConfig.getParam("config"));

		cmdLine.copyIntoConfig(initConfig);
		initCc.check(initConfig);
		
		//If Cosine scoring is selected require WCCN parameter
		if(initConfig.getParam("eigenChannelMode") == "JFA"){
			cc.addIntegerParam("eigenChannelNumber",true,true, "final rank of EigenChannel matrix");
		}
		else{
			cc.addIntegerParam("eigenChannelRank",false,true,"final rank of EigenChannel matrix");
		}


		// Optionnal
		cc.addStringParam("initEigenChannelMatrix",false,true,"init EigenChannel Matrix");
		cc.addStringParam("eigenVoiceMatrix",false,true,"name of the EigenVoice Matrix");
		cc.addBooleanParam("loadAccs",false,true,"if true do not compute UBM stats, load matrices");
		cc.addBooleanParam("checkLLK",false,true,"if true do compute the likelihood of training data after each iteration");
		cc.addBooleanParam("saveInitChannelMatrix",false,true,"if true save the matrix used for initialisation");
		cc.addBooleanParam("saveAllECMatrices",false,true,"if true save the matrices after each iteration");
		cc.addIntegerParam("computeLLK",false,true,"optional: nb of files where LLK is computed");	
		

		// Insertion of config compatibility rules
	     if (cmdLine.displayHelpRequired()){
			cout << "****************************************" << endl;
			cout << "********** EigenChannel.exe ************" << endl;
			cout << "****************************************" << endl;
			cout << endl;
			cout << "Evaluate EigenChannel Matrix from sessions data" << endl;
			cout <<endl<<cc.getParamList()<<endl;
			return 0;  
		}
		if (cmdLine.displayVersionRequired()){
		cout <<"Version 3.0 ALIZE Package"<<endl;
		} 

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

		if(config.getParam("eigenChannelMode") == "JFA")		EigenChannelJFA(config);
		else if(config.getParam("eigenChannelMode") == "LFA")	EigenChannelLFA(config);
		else{
			cout<<"Error : wrong eigenChannelMode parameter, please chose JFA or LFA"<<endl;
		}
	}
	catch (Exception& e) {cout << e.toString() << cc.getParamList() << endl;}
if (debug) {
}
return 0;
}

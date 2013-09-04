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
#include <alize.h>
#include <liatools.h>
#include <PLDA.h>

using namespace alize;
using namespace std;

int main(int argc, char* argv[]) {
	ConfigChecker initCc,cc;

	try{
		// Insertion of config compatibility rules
		CmdLine cmdLine(argc, argv);

		// Needed params
		initCc.addStringParam("backgroundNdxFilename",true,true,"NDX of multiple GMM speaker recordings");
		initCc.addStringParam("loadVectorFilesPath",true,true,"directory of vectors to load");
		initCc.addStringParam("loadMatrixFormat",true,true,"format of matrices to load, DB for binary, DT for ASCII");
		initCc.addStringParam("saveMatrixFormat",true,true,"format of matrices to save, DB for binary, DT for ASCII");
		initCc.addStringParam("matrixFilesPath",true,true,"directory to store matrices");
		initCc.addStringParam("saveMatrixFilesExtension",true,true,"extension of matrix files to save");
		initCc.addStringParam("vectorFilesExtension",true,true,"extension of vector files");
		initCc.addIntegerParam("pldaEigenVoiceNumber",true,true,"Rank of the EigenVoice matrix");
		initCc.addIntegerParam("pldaEigenChannelNumber",true,true,"Rank of the EigenChannel matrix");
		initCc.addIntegerParam("pldaNbIt",true,true,"number of EM it");	
		initCc.addBooleanParam("pldaLoadInitMatrices",true,true,"If true, load existing initialization matrices");
		
		// Check existing parameters to create the appropriate ConfigChecker
		Config tmpConfig;
		cmdLine.copyIntoConfig(tmpConfig);
		Config initConfig;
		if (tmpConfig.existsParam("config")) initConfig.load(tmpConfig.getParam("config"));

		if (cmdLine.displayHelpRequired()){
			cout << "****************************************" << endl;
			cout << "********** PLDA.exe ************" << endl;
			cout << "****************************************" << endl;
			cout << endl;
			cout << "Estimate PLDA model from development data" << endl;
			cout <<endl<<initCc.getParamList()<<endl;
			return 0;  
		}

		cmdLine.copyIntoConfig(initConfig);
		initCc.check(initConfig);

		// If initialize matrices from scratch
		if(initConfig.getParam("pldaLoadInitMatrices").toBool()){
			cc.addStringParam("pldaEigenVoiceMatrixInit",true,true,"file name of the EigenVoice matrix to load for initialization");
			cc.addStringParam("pldaEigenChannelMatrixInit",true,true,"file name of the EigenChannel matrix to load for initialization");
			cc.addStringParam("pldaSigmaMatrixInit",true,true,"file name of the noise matrix to load for initialization");
			cc.addStringParam("pldaMeanVecInit",true,true,"file name of the mean vector to load for initialization");	
			cc.addStringParam("loadMatrixFilesExtension",true,true,"extension of matrix files to load");
		}
		else{	// If load matrices from existing files
			cc.addStringParam("randomInitLaw",true,true,"random law to use for matrices initialization");
		}

		cc.addStringParam("pldaOriginalMean",false,true,"file name of the mean vector to save");
		cc.addStringParam("pldaEigenVoiceMatrix",false,true,"file name of the EigenVoice matrix to save");
		cc.addStringParam("pldaEigenChannelMatrix",false,true,"file name of the EigenChannel matrix to save");
		cc.addStringParam("pldaSigmaMatrix",false,true,"file name of the noise matrix of the noise matrix to save");
		cc.addStringParam("pldaMinDivMean",false,true,"file name of the mean estimate after minimum divergence");
		
		if (cmdLine.displayHelpRequired()){
			cout << "****************************************" << endl;
			cout << "********** PLDA.exe ************" << endl;
			cout << "****************************************" << endl;
			cout << endl;
			cout << "Estimate PLDA model from development data" << endl;
			cout <<endl<<cc.getParamList()<<endl;
			return 0;  
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

		PLDA(config);
	}
	catch (Exception& e) {cout << e.toString() << cc.getParamList() << endl;}
if (debug) {
}
return 0;
}

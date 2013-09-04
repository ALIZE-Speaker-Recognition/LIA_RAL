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

#include <IvNorm.h>

using namespace alize;
using namespace std;

int main(int argc, char* argv[]) {
	ConfigChecker initCc, cc;
	try{
		// Insertion of config compatibility rules
		CmdLine cmdLine(argc, argv);

		//Add list of parameters mandatory 
		initCc.addStringParam("config", false, true, "default config filename");
		initCc.addBooleanParam("ivNormLoadParam",true,true,"load existing normalization parameters");
		initCc.addStringParam("saveMatrixFilesExtension",true,true,"extension of matrix files to save");
	
		initCc.addBooleanParam("LDA",false,true,"apply Linear Discriminant Analysis");
		initCc.addIntegerParam("ivNormIterationNb",true,true,"apply Eigen Factor Radial normalization");

		initCc.addStringParam("loadMatrixFilesExtension",true,true,"extension of matrix files to load");
		initCc.addStringParam("saveMatrixFormat",true,true,"format to save matrices, DB for binary | DT for ASCII");
		initCc.addStringParam("loadMatrixFormat",true,true,"format to load matrices, DB for binary | DT for ASCII");
		initCc.addStringParam("matrixFilesPath",true,true,"path to matrices");

		initCc.addStringParam("loadVectorFilesExtension",true,true,"extension of vector files to load");
		initCc.addStringParam("saveVectorFilesExtension",true,true,"extension of vector files to save");
		initCc.addStringParam("ivNormEfrMatrixBaseName",true,true,"Eigen Factor Radial matrix name");
		initCc.addStringParam("ivNormEfrMeanBaseName",true,true,"Eigen Factor Radial mean vector name");
		
		


		// Check existing parameters to create the appropriate ConfigChecker
		Config tmpConfig;
		cmdLine.copyIntoConfig(tmpConfig);
		Config initConfig;
		if (tmpConfig.existsParam("config")) initConfig.load(tmpConfig.getParam("config"));

		cmdLine.copyIntoConfig(initConfig);
		initCc.check(initConfig);

		// List of parameters required if ivNormLoadParam
		if(initConfig.getParam("ivNormLoadParam").toBool())
			if(!initConfig.existsParam("inputVectorFilename") && !initConfig.existsParam("ndxFilename"))
				cc.addStringParam("inputVectorFilename",true,true,"list of vectors to normalize");
		else{		// if data, requires parameter to know if we save the matrices or not
			initCc.addStringParam("loadVectorFilesPath",true,true,"path to vectors");
			cc.addStringParam("backgroundNdxFilename",true,true,"list of files for normalization parameters estimation");
			cc.addStringParam("inputVectorFilename",false,true,"input list of vectors to normalize");
			cc.addStringParam("ndxFilename",false,true,"input list of vectors to normalize as test Ndx file");
		}

		if(initConfig.existsParam("inputVectorFilename") || initConfig.existsParam("ndxFilename"))
			cc.addStringParam("saveVectorFilesPath",true,true,"path to save normalized vectors");

		if(initConfig.existsParam("ndxFilename")){
			cc.addStringParam("testVectorFilesPath",true,true,"path to test vectors");
			cc.addStringParam("targetIdList",true,true,"Index for model training");
		}

		if(initConfig.getParam("ivNormIterationNb").toLong() != 0 ){
			cc.addStringParam("ivNormEfrMode",true,true,"normalization to apply, EFR | sphNorm");
			cc.addStringParam("ivNormEfrMatrixBaseName",false,true,"root of EFR matrices name to save");
			cc.addStringParam("ivNormEfrMeanBaseName",false,true,"root of EFR mean vector name to save");
		}
		if(initConfig.existsParam("LDA")){
			cc.addIntegerParam("ldaRank",true,true,"rank of the Linear Discriminant Analysis matrix");
			cc.addStringParam("ldaMatrix",true,true,"filename of the Linear Discriminant Analysis matrix");
		}
			


		if (cmdLine.displayHelpRequired()){
			cout << "****************************************" << endl;
			cout << "********** IvNorm.exe ************" << endl;
			cout << "****************************************" << endl;
			cout << endl;
			cout << "Normalize i-vectors" << endl;
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
		
		IvNorm(config);

	}
	catch (Exception& e) {cout << e.toString() << cc.getParamList() << endl;}
if (debug) {
}
return 0;
}

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
#include "liatools.h"
#include "IvTest.h"

using namespace std;
using namespace alize;

int main(int argc, char* argv[]) {
		ConfigChecker initCc,cc;
try {
		// Insertion of config compatibility rules
		CmdLine cmdLine(argc, argv);

		//Add list of parameters always mandatory 
		initCc.addStringParam("config", false, true, "default config filename");
		initCc.addStringParam("outputFilename",true,true, "output scores in this file: 'gender - test model - scores'");	
		initCc.addStringParam("gender",true,true, "M/F: will output the gender in score file");
		initCc.addStringParam("targetIdList",true,false, "Index file describing model enrollment segments");
		initCc.addStringParam("outputScoreFormat",true,true, "Format of the score file to write, ascii || binary");

		// Check existing parameters to create the appropriate ConfigChecker
		Config tmpConfig;
		cmdLine.copyIntoConfig(tmpConfig);
		Config initConfig;
		if (tmpConfig.existsParam("config")) initConfig.load(tmpConfig.getParam("config"));

		cmdLine.copyIntoConfig(initConfig);
		initCc.check(initConfig);

		cc.addStringParam("loadVectorFilesExtension",true,true,"extension of vector files");
		cc.addStringParam("loadVectorFilesPath",true,true,"path to vectors");
		cc.addStringParam("testVectorFilesPath",true,true,"path to test vectors");
		cc.addStringParam("scoring",true,true,"cosine / mahalanobis / 2cov / plda");
		cc.addBooleanParam("ivNorm",true,true,"true/false, Normalize i-vector before scoring");
		cc.addStringParam("matrixFilesPath",true,true,"directory for matrix storage");
		cc.addStringParam("saveMatrixFilesExtension",true,true,"extension for matrix files to save");
		cc.addStringParam("loadMatrixFilesExtension",true,true,"extension for matrix files to load");
		cc.addStringParam("loadMatrixFormat",true,true,"format to load matrices, DB for binary | DT for ASCII");
		cc.addStringParam("saveMatrixFormat",true,true,"format to save matrices, DB for binary | DT for ASCII");


		//If Cosine scoring is selected require WCCN parameter
		if(initConfig.getParam("scoring") == "cosine"){
			cc.addStringParam("wccn",true,true, "true/false, apply Within Class Covariance Scoring");
			if(initConfig.getParam("wccn").toBool())
				cc.addStringParam("loadWccnMatrix",true,true, "true/false, load the WCCN matrix as it exists");
		}

		//If Mahalanobis scoring is selected require WCCN parameter
		if(initConfig.getParam("scoring") == "mahalanobis"){
			cc.addStringParam("mahalanobisMatrix",true,true, "Mahalanobis matrix to load or to save if it does not exist");
		}
		
		// Choose the mode of file list to load: inputVectorFilename or ndxFilename and optional targetIdList
		if(!initConfig.existsParam("ndxFilename") && !initConfig.existsParam("inputVectorFilename")){
			cc.addStringParam("ndxFilename",true,true, "NDX file listing all verification tests to achieve: first column: test File, all others: models");
		}
		if(initConfig.existsParam("ndxFilename"))
			cc.addStringParam("ndxFilename",false,true, "NDX file listing all verification tests to achieve: first column: test File, all others: models");
		if(initConfig.existsParam("inputVectorFilename"))
			cc.addStringParam("inputVectorFilename",false,true, "List of files use as test-segments. If no targetIdList is given, this list of file is also use as model list");

		// if normalize ivectors then check the required options
		if(initConfig.getParam("ivNorm").toBool()){
			cc.addBooleanParam("ivNormLoadParam",true,true,"load existing normalization parameters");
			cc.addStringParam("backgroundNdxFilename",true,true,"list of files for normalization parameters estimation");

			cc.addBooleanParam("LDA",false,true,"apply Linear Discriminant Analysis");
			cc.addIntegerParam("ivNormIterationNb",true,true,"apply Eigen Factor Radial normalization");

			if(initConfig.existsParam("saveNormVect")){
				cc.addBooleanParam("saveNormVect",true,true,"true/false, save enrollment and test i-vectors after normalization");
				cc.addStringParam("saveVectorFilesPath",true,true,"path to save normalized vectors");
			}

			if(initConfig.getParam("ivNormIterationNb").toLong() != 0 ){
				cc.addStringParam("ivNormEfrMode",true,true,"normalization to apply, EFR | sphNorm");
				cc.addStringParam("ivNormEfrMatrixBaseName",false,true,"root of EFR matrices name to save");
				cc.addStringParam("ivNormEfrMeanBaseName",false,true,"root of EFR mean vector name to save");
			}
			if(initConfig.existsParam("LDA")&& initConfig.getParam("LDA").toBool()){
				cc.addIntegerParam("ldaRank",true,true,"rank of the Linear Discriminant Analysis matrix");
				cc.addStringParam("ldaMatrix",true,true,"filename of the Linear Discriminant Analysis matrix");
			}
		}

		if(initConfig.getParam("scoring") == "plda"){
			cc.addIntegerParam("iVectSize",true,true,"vector size for PLDA training");
			cc.addIntegerParam("pldaEigenVoiceNumber",true,true,"rank of the PLDA eigenvoice matrix");
			cc.addIntegerParam("pldaEigenChannelNumber",true,true,"rank of the PLDA eigenchannel matrix");
			cc.addBooleanParam("pldaLoadModel",true,true,"if true, use existing PLDA matrices");
				
			cc.addStringParam("pldaEigenVoiceMatrix",true,true,"file name of the EigenVoice matrix to load");
			cc.addStringParam("pldaEigenChannelMatrix",true,true,"file name of the EigenChannel matrix to load");
			cc.addStringParam("pldaSigmaMatrix",true,true,"file name of the noise matrix to load");
			cc.addStringParam("pldaMeanVec",true,true,"file name of the mean vector to load");	
			cc.addStringParam("pldaMinDivMean",true,true,"file name of the mean vector to load for minimum divergence");

			if(initConfig.getParam("pldaLoadModel").toBool()){	
				cc.addStringParam("loadMatrixFilesExtension",true,true,"extension of matrix files to load");
			}
			else{
				cc.addBooleanParam("pldaLoadInitMatrices",true,true,"if true, initializae PLDA model from existing matrices");
				cc.addStringParam("backgroundNdxFilename",true,true,"list of files for normalization parameters estimation");
				cc.addStringParam("saveMatrixFormat",true,true,"format to save matrices, DB for binary | DT for ASCII");
				cc.addStringParam("saveMatrixFilesExtension",true,true,"extension of matrix files to save");

				if(initConfig.getParam("pldaLoadInitMatrices").toBool()){
					cc.addStringParam("pldaEigenVoiceMatrixInit",true,true,"file name of the EigenVoice matrix to load for initialization");
					cc.addStringParam("pldaEigenChannelMatrixInit",true,true,"file name of the EigenChannel matrix to load for initialization");
					cc.addStringParam("pldaSigmaMatrixInit",true,true,"file name of the noise matrix to load for initialization");
					cc.addStringParam("pldaMeanVecInit",true,true,"file name of the mean vector to load for initialization");	
					cc.addStringParam("loadMatrixFilesExtension",true,true,"extension of matrix files to load");
				}
				else{
					cc.addIntegerParam("pldaNbIt",true,true,"number of iterations for PLDA training");
				}
			}
		}
		

		
		if (cmdLine.displayHelpRequired()){
		cout << "************************************" << endl;
		cout << "********** IvTest.exe **********" << endl;
		cout << "************************************" << endl;
		cout << endl;
		cout << "Score computation for an NDX (NIST format) File" << endl;
		cout << "" << endl;
		cout << cc.getParamList()<<endl;
		return 0;  
		}

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

		IvTest(config);
}
catch (alize::Exception& e) {cout << e.toString() << endl << cc.getParamList()<< endl;}
return 0;
}

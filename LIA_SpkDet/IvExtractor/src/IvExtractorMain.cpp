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
#include "IvExtractor.h"

int main(int argc, char* argv[]){

    ConfigChecker cc;
    cc.addStringParam("config", false, true, "default config filename");
    cc.addIntegerParam("verboseLevel",false,true,"level of the berose information 0=no verbose, 1=normal, 2=more");
    cc.addStringParam("inputFeatureFilename",false, true,"feature filename or filename of a text file with the list of feature filenames");
    cc.addStringParam("targetIdList",true,true,"The file with the list of models to train. A line is composed by client_id file1 file2 ...");
    cc.addStringParam("inputWorldFilename",false,true,"if set, the init is based on a model get from this file, else from scratch");
    cc.addStringParam("labelSelectedFrames",false,true,"The segments with this label are used for training the worldmodel (if UseIdForSelectedFrame is not used)"); 
    cc.addBooleanParam("info",false,false,"If info is requested, just info on the train set is outputed");
    cc.addStringParam("saveVectorFilesPath",true,true,"Directory to store i-vectors after extraction");
	cc.addStringParam("mode",false,true,"Mode of i-vector extraction, classic (default) | ubmWeight | eigenDecomposition");
	cc.addBooleanParam("loadUbmWeightParam",false,true,"If true, load pre-computed parameters for ubmWeight approximation");
	

try{
    CmdLine cmdLine(argc, argv);
    if (cmdLine.displayHelpRequired()){
      cout <<"IvExtractor.exe"<<endl<<"This program is used for extracting i-vectors"
	   <<endl<<cc.getParamList()<<endl;
      return 0;  
    }
    if (cmdLine.displayVersionRequired()){
    cout <<"Version 1.0"<<endl;
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
    bool train=true;                                 				// Training the target models
    if (config.existsParam("info"))                  				// Information on the data available for the targets and output a file with good amount of data
      train=false;

	String mode = "classic";
	if(config.existsParam("mode"))
		mode = config.getParam("mode");
    if (train){
//		if (config.existsParam("byLabelModel"))    					// if the paramter is present, for each client, we train one model for each cluster 
//			TrainTargetByLabel(config);                  			// (a cluster is a set of segments with the same label)
		if(mode == "classic")
			IvExtractor(config); 
		else if(mode == "ubmWeight")
			IvExtractorUbmWeigth(config);
		else if(mode == "eigenDecomposition")
			IvExtractorEigenDecomposition(config);
		else{
			ExtractTargetDataInfo(config);
		}
    }
    else   ExtractTargetDataInfo(config);   
    
    return 0;
  }
  catch(alize::Exception & e){ cout <<"IvExtractor "<< e.toString() << endl << cc.getParamList()<<endl;}
}

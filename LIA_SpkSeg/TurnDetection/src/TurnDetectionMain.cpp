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
#include "TurnDetection.h"

int main(int argc, char* argv[]){
  using namespace std;
  using namespace alize;

  ConfigChecker cc;

  cc.addStringParam("config",true,true,"defines the default config filename");

  cc.addFloatParam("frameLength",true,true,"defines length of a frame (default=10ms)");
  cc.addIntegerParam("nbTrainIt",true,true,"defines the number of it, the ceiling and flooring are moved and the baggedFrameProbability is used"); 
  cc.addFloatParam("baggedFrameProbabilityInit",false,true,"defines the % of frames taken by component for the initializing of the mixture - mandatory if init from scratch");
  cc.addFloatParam("baggedFrameProbability",true,true,"defines the % of frames taken for each iterations");
  cc.addFloatParam("initVarianceFlooring",true,true,"defines the variance control parameters - relative to global data variance - initial value (moved during the it)");
  cc.addFloatParam("initVarianceCeiling",true,true,"defines the variance control parameters - relative to global data variance - initial value (moved during the it)"); 
  cc.addFloatParam("finalVarianceFlooring",true,true,"defines the variance control parameters - relative to global data variance - final value");
  cc.addFloatParam("finalVarianceCeiling",true,true,"defines the variance control parameters - relative to global data variance - final value"); 

  // For the MAP Algorithm
//  cc.addStringParam("MAPAlgo",true,true,"defines the adaptation method (MAPConst,MAPConst2,MAPOccDep,MLLR)");
//  cc.addBooleanParam("meanAdapt",false,true,"if set to true, the gaussians means are adapted");
//  cc.addBooleanParam("varAdapt",false,true,"if set to true, the gaussians variances are adapted");
//  cc.addBooleanParam("weightAdapt",false,true,"if set to true, the gaussians weights are adapted");
//  cc.addFloatParam("MAPAlphaMean",false,true,"defines the value alpha for the adaptation of the means for MAP (MAPConst, MAPConst2)");
//  cc.addFloatParam("MAPAlphaVar",false,true,"defines the value alpha for the adaptation of the variances for MAP (MAPConst, MAPConst2)");
//  cc.addFloatParam("MAPAlphaWeight",false,true,"defines the value alpha for the adaptation of the weights for MAP (MAPConst, MAPConst2)");
//  cc.addFloatParam("MAPRegFactorMean",false,true,"defines the value alpha for the adaptation of the means for MAP (MAPOccDep, MAPModelBased, MLLR)");
//  cc.addFloatParam("MAPRegFactorVar",false,true,"defines the value alpha for the adaptation of the variances for MAP (MAPOccDep, MAPModelBased, MLLR)");
//  cc.addFloatParam("MAPRegFactorWeight",false,true,"defines the value alpha for the adaptation of the weights for MAP (MAPOccDep, MAPModelBased, MLLR)");

//  cc.addBooleanParam("normalizeModel",false,true,"if set to true, normalize the model (at each iteration)");
//  cc.addBooleanParam("normalizeModelMeanOnly",false,true,"if set to true, normalize only the means of the model");
//  cc.addIntegerParam("normalizeModelNbIt",false,true,"defines the number of it (only if normalizeModelMeanOnly is set to true)"); 
//  cc.addBooleanParam("normalizeClient",false,true,"if set to true, normalize the client (at each iteration)");

  cc.addStringParam("labelFilesPath",false,true,"defines the path where to load original labelFiles");
  cc.addStringParam("labelFilesExtension",false,true,"defines the extension to load original labelFiles");
//  cc.addStringParam("labelSelectedFrames",false,true,"only the frames from segments with this label will be used");

  // For the TurnDetection Algorithm
  cc.addStringParam("listFileToSegment",true,true,"defines the list of the files to segment");
  cc.addStringParam("outputFilesPath",true,true,"defines the path where will be to store the produced files");
  cc.addStringParam("fileRefPath",false,true,"defines the path of the reference files");
  cc.addStringParam("saveSegmentationExtension",true,true,"defines the extension to save the files of the label");
  cc.addStringParam("clusteringCrit",false,true,"defines the criteria (BIC, GLR, DGLR or DELTABIC) for the clustering (default=DGLR)");
  cc.addFloatParam("clusteringCritThresh",false,true,"defines the threshold value for the clustering (default=0.0)");
  cc.addIntegerParam("winSize",false,true,"defines the size of the window for the clustering (default=50)");
  cc.addIntegerParam("winStep",false,true,"defines the step of the window for the clustering (default=5)");
  cc.addFloatParam("alpha",false,true,"factor to find the maxima of the criterion value curve : alpha * standard deviation (default=0.7)");

  try {
      CmdLine cmdLine(argc, argv);
      if (cmdLine.displayHelpRequired()){
        cout << "*******************************" << endl;
        cout << "********* TurnDetection *******" << endl;
        cout << "*******************************" << endl;
        cout << endl;
        cout << "This program is a speaker turn detection program." << endl;
        cout << "This program is provided but not used in the LIA ......" << endl;
        cout << "" << endl;
        cout << cc.getParamList()<<endl;
        return 0;  
      }
      if (cmdLine.displayVersionRequired()){
        cout << "TurnDetection - Version 2" << endl;
        return 0;  
      } 
      Config tmp;
      cmdLine.copyIntoConfig(tmp);
      Config config;
      if (tmp.existsParam("config")) config.load(tmp.getParam("config"));
      cmdLine.copyIntoConfig(config);
      cc.check(config); 
      debug=config.getParam_debug();
      if (config.existsParam("verbose")) verbose=config.getParam("verbose").toBool(); else verbose=false;
      if (verbose) verboseLevel=1; else verboseLevel=0;
      if (config.existsParam("verboseLevel")) verboseLevel=config.getParam("verboseLevel").toLong();
      if (verboseLevel>0) verbose=true;
      launchTurnDetectionProcess(config);
    }
  catch (alize::Exception& e) { cout << e.toString() << endl << cc.getParamList() << endl; }
  return 0;     
}

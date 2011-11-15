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

#include "Svm.h"
#include "alize.h"
#include "liatools.h"
int main(int argc, char* argv[])
{
  using namespace std;
  using namespace alize;
  ConfigChecker cc;

  // params 
  cc.addStringParam("mode",true,true,"train | predict | tnormpredict");
  cc.addStringParam("inputFilename",true, true,"train: positive instance | predict: test instance");
  cc.addStringParam("vectorFilesPath",true, true,"Path to instances");  
  cc.addStringParam("modelFilesPath",true, true,"Path to models");    
  cc.addStringParam("vectorFilesExtension",true, false,"ext instances");  
  cc.addStringParam("modelFilesExtension",true, false,"ext models");      

  // optional params
  cc.addStringParam("inputBCKList",false, true,"negative exemples to use every time world.vect");  
  cc.addStringParam("inputSVMModel",false, true,"in test mode");  
  cc.addIntegerParam("vsize",false, true,"vector sizes");
  cc.addFloatParam("targetPenalty",false, true,"penalty for error on target speaker class");  
  cc.addFloatParam("C",false, true,"C parameter");      

  // To modify --> The param have utility only in the predic and tnormpredic step to generate results in a 
  // file. This option must be reprogrammed to be true only in the case of mode predict or tnormpredict
  cc.addStringParam("outputFilename",false, true,"train: model file name | predict: results");  
 
  try {
        CmdLine cmdLine(argc, argv);
        if (cmdLine.displayHelpRequired()){
          cout <<"Svm.exe"<<endl<<"Binding to libSvm for Speaker Verification in large Scale Evaluation as NIST SRE" <<endl<<cc.getParamList()<<endl;
          return 0;  
        }
        if (cmdLine.displayVersionRequired()){
          cout <<"Version 2.0"<<endl;
        } 
        Config tmp;
        cmdLine.copyIntoConfig(tmp);
        Config config;
        if (tmp.existsParam("config")) config.load(tmp.getParam("config"));
        cmdLine.copyIntoConfig(config);
        debug=config.getParam_debug();
        cc.check(config);
        if (config.existsParam("verbose"))verbose=config.getParam("verbose").toBool();else verbose=false;
        if (verbose) verboseLevel=1;else verboseLevel=0;
        if (config.existsParam("verboseLevel"))verboseLevel=config.getParam("verboseLevel").toLong();
        if (verboseLevel>0) verbose=true;       
          
        String mode=config.getParam("mode");
        if (verbose) cout << "(Svm) Mode " << mode << endl;
        if (mode=="train") svmTrain(config);
        else if (mode=="predict") svmPredict(config);
        else if(mode=="tnormpredict") svmPredictTnorm(config);
        else throw Exception("No mode",__FILE__,__LINE__);
        }
  catch (alize::Exception & e) {cout << e.toString () << endl << cc.getParamList() << endl;}
  return 0;
}

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
#include "TrainWorld.h"

int main(int argc, char* argv[])
{
  ConfigChecker cc;
  cc.addIntegerParam("verboseLevel",false,true,"level of the verbose information 0=no verbose, 1=normal, 2=more");
  cc.addStringParam("inputFeatureFilename",false, true,"feature filename or filename of a text file with the list of feature filenames");
  cc.addStringParam("inputStreamList",false, true,"filename of a text file with the filename of input streams");
  cc.addStringParam("weightStreamList",false,true,"filename of a text file with the weight of each input stream - default=equal weights");
  cc.addStringParam("outputWorldFilename",true,true,"output worldmodel filename");                            
  cc.addStringParam("inputWorldFilename",false,true,"if set, the init is based on a model get from this file, else frrom scratch");
  cc.addStringParam("saveInitModel",false,true,"if set (default), save the initial model");
  cc.addStringParam("labelSelectedFrames",true,true,"the segments with this label are used for training the worldmodel");
  cc.addFloatParam("baggedFrameProbability",true,true,"defines the % of frames taken for each iterations");
  cc.addFloatParam("baggedFrameProbabilityInit",false,true,"NOT LONGER USED IN TRAINWORLD !! deprecated and remplaced by nbFrameToSelect (defines the % of frames taken BY COMPONENT for the initializing of the mixture- mandatory if init from scratch)");
  cc.addIntegerParam("nbFrameToSelect",false,true,"Defines the number of frames selected to initialise one component, default=50 ");
  cc.addIntegerParam("baggedMinimalLength",false,true,"minimum length for selected segments in bagged (default=3)");
  cc.addIntegerParam("baggedMaximalLength",false,true,"maximal length for selected segments in bagged (default=7)");
  cc.addFloatParam("initVarianceFlooring",true,true,"variance control parameters - relative to global data variance - initial value (moved during the it)");
  cc.addFloatParam("initVarianceCeiling",true,true,"variance control parameters - relative to global data variance - initial value (moved during the it)"); 
  cc.addFloatParam("finalVarianceFlooring",true,true,"variance control parameters - relative to global data variance - final value");
  cc.addFloatParam("finalVarianceCeiling",true,true,"variance control parameters - relative to global data variance - final value"); 
  cc.addIntegerParam("nbTrainIt",true,true,"number of it, the ceiling and flooring are moved and the baggedFrameProbability is used"); 
  cc.addBooleanParam("normalizeModel",false,true,"if set to true,  normalize the world (at each iteration)");
  cc.addBooleanParam("normalizeModelMeanOnly",false,true,"used only if normalizeModel is On, says if only mean parameters should be normalized"); 
  cc.addIntegerParam("normalizeModelNbIt",false,true,"used only if noramlizeModelMeanOnly is set, nb of normalization it");
  cc.addBooleanParam("use01",false,true,"if set at true, don't compute the global mean and cov but uses 0 mean and 1 cov");
  cc.addBooleanParam("componentReduction",false,true,"if set reduce the number of components at each it, selecting the best weights until targetDistribCount (default false)");
  cc.addIntegerParam("targetMixtureDistribCount",false,true,"final number of components if componentReduction is selected"); 
  cc.addIntegerParam("initRand",false,true,"initialisation of the random generator for bagged set of data (default=0)"); 


  try {
      CmdLine cmdLine(argc, argv);
      if (cmdLine.displayHelpRequired()){
        cout <<"TrainWorld.exe"<<endl<<"This program is used for training a world model from scratch or from a model"
             <<endl<<cc.getParamList()<<endl;
        return 0;  
      }
      if (cmdLine.displayVersionRequired()){
        cout <<"Version 2-beta"<<endl;
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
      trainWorld(config);
    }
catch (alize::Exception& e) {cout << e.toString() << endl << cc.getParamList()<< endl;}
return 0;     
}

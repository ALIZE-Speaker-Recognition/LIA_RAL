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

#include <liatools.h>
#include "alize.h"
#include "EnergyDetector.h"

//---------------------------------------------------------------------------------------------------------------
// Energy detector
// Input : Xlist	Format: File1 File2 on separate or same lines

int main (int argc, char *argv[])
{
  using namespace std;
  using namespace alize;

    ConfigChecker cc;
    cc.addStringParam("config", false, true, "default config filename");
    cc.addStringParam("inputFeatureFilename",true,true,"input feature - could be a simple feature file or a list of filename");
    cc.addStringParam("labelFilesPath",false,true,"path of input and output labels");
    cc.addStringParam("saveLabelFileExtension",true,true,"output labelFile extension");
    cc.addStringParam("labelOutputFrames",true,true,"output label tag");
    cc.addFloatParam("alpha",true,true,"threshold giving the percentage of data of the middle Gaussian to take into account");
    cc.addStringParam("segmentalMode",true,true,"file to have a file by file behaviour");

    // gestion du thresholdMode obligatoire
    cc.addStringParam("thresholdMode",true,true,"this parameter must be set to select 3 top gaussian. It's a default parameter");

  try
  {
    CmdLine cmdLine (argc, argv);
    if (cmdLine.displayHelpRequired ()){	// --help
	cout << "************************************" << endl;
	cout << "********** EnergyDetector.exe *********" << endl;
	cout << "************************************" << endl;
	cout << endl;
	cout << "Speech/Non Speech Detector." << endl;
	cout<<cc.getParamList()<<endl;
      }
    else if (cmdLine.displayVersionRequired ())
      cout << "Version 2beta" << endl;
    else{
      // copy parameters from command line into config
      Config tmp;
      cmdLine.copyIntoConfig (tmp);
      Config config;
      if (tmp.existsParam("config")) config.load(tmp.getParam("config"));
      cmdLine.copyIntoConfig(config);
      cc.check(config);
      debug=config.getParam_debug();
      if (config.existsParam("verbose")) verbose=config.getParam("verbose").toBool();
      else verbose=false;
      if (verbose) verboseLevel=1;else verboseLevel=0;
      if (config.existsParam("verboseLevel"))verboseLevel=config.getParam("verboseLevel").toLong();
      if (verboseLevel>0) verbose=true;
	      
	bool energy=true;
	if (config.existsParam("turnDetection")) energy=false;                 // choice between energy detection or GLR turn detection
	String inputFeatureFileName =config.getParam("inputFeatureFilename");  // input feature - could be a simple feature file or a list of filename
	String extOutput=".lbl";                                               // the extension of the output files    
	if (config.existsParam("saveLabelFileExtension")) extOutput=config.getParam("saveLabelFileExtension");   
	String pathOutput="./";                                                // the path of the output files    
	if (config.existsParam("labelFilesPath")) pathOutput=config.getParam("labelFilesPath");
	XLine inputFeatureFileNameList;                                        // The (feature) input filename list - list of files to process
	    try{                                            
	      XList inputFileNameXList(inputFeatureFileName,config);               // Read the filename list file if it is really a list
	      inputFeatureFileNameList=inputFileNameXList.getAllElements();        // And put the filename in a list if the file is a list of feature filenames
	    }
	    catch(FileNotFoundException& e){                                       // It was a simple feature file and not a filename list
	      inputFeatureFileNameList.addElement(inputFeatureFileName);           // add the filename in the list
	    }
      String *file;
      while ((file=inputFeatureFileNameList.getElement())!= NULL){         // Loop on each feature file
	String & fileName=*file;                                           // current input file basename 
	SegServer segServer;                                               // Create the segment server for dealing with selected/unselected segments
	if (energy){                                                       // The energy detector
	  SegCluster &outputSeg=energyDetector(config,segServer,fileName); 
	  outputLabelFile(outputSeg,pathOutput+fileName+extOutput,config); 
	}
      }                                                                    // end of file loop
    }
  }
  catch (alize::Exception & e) {cout << e.toString () << endl << cc.getParamList() << endl;}
return 0;
}


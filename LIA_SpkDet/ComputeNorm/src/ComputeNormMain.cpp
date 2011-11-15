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

#include "ComputeNorm.h"
#include "liatools.h"

using namespace std;
using namespace alize;

int main(int argc, char* argv[])
{
  ConfigChecker cc;
  cc.addStringParam("testNistFile",true, true,"target_seg score file");
  cc.addStringParam("normType",true, true,"tnorm|znorm|ztnorm|tznorm, the normalization method (ztnorm is outputing both tnorm and ztnorm scores, tznorm both tznorm and znorm scores)");
  cc.addStringParam("tnormNistFile",false,true,"imp_seg score file, impostor scores used for tnorm and ztnorm");
  cc.addStringParam("znormNistFile",false,true,"target_imp score file, impostor scores used for znorm and ztnorm");
  cc.addStringParam("ztnormNistFile",false,true,"imp_imp score file, impostor scores used for ztnorm and tznorm");
  cc.addStringParam("outputFileBaseName",true, true,"the output file(s) basename");
  cc.addStringParam("znormFilesExtension",false, true,"znorm output file extension (default='.znorm'");
  cc.addStringParam("ztnormFilesExtension",false, true,"ztnorm output file extension (default='.znorm'");
  cc.addStringParam("tnormFilesExtension",false, true,"tnorm output file extension (default='.tnorm'");
  cc.addStringParam("tznormFilesExtension",false, true,"tznorm output file extension (default='.tznorm'");
  cc.addStringParam("cohortFilePath",false, true,"cohort files path, for selectTargetDependentCohortInFile");
  cc.addStringParam("cohortFileExt",false, true,"cohort files extension, for selectTargetDependentCohortInFile");
  cc.addIntegerParam("maxIdNb",false, true,"Max target speakers - use to fix the max number of znorm score distributions (default=1000)");
  cc.addIntegerParam("maxSegNb",false, true,"Max test segments - use to fix the max number of tnorm score distributions (default=1000)");
  cc.addIntegerParam("maxScoreDistribNb",false, true,"Max scores per distribution - use to fix the max number of score in a distribution (default=1000)");
  cc.addIntegerParam("fieldGender",false, true,"The field for gender in the nist file format (default=0)");
  cc.addIntegerParam("fieldName",false, true,"The field for gender in the nist file format (default=1)");
  cc.addIntegerParam("fieldDecision",false, true,"The field for gender in the nist file format (default=2)");
  cc.addIntegerParam("fieldSeg",false, true,"The field for gender in the nist file format (default=3)");
  cc.addIntegerParam("fieldLLR",false, true,"The field for gender in the nist file format (default=4)");
  cc.addStringParam("impostorIDList",false, true,"If the option is set, it limits the used impostor scores to the one of the given list");
  cc.addIntegerParam("meanMode",false,true,"Score distrib mean computation mode. 0 classical, 1 median (default 0)"); 
  cc.addFloatParam("percentH",false,true,"% of higest scores discarded (default 0)");  
  cc.addFloatParam("percentL",false,true,"% of lowest scores discarded (default 0)");  
  
  try {
      CmdLine cmdLine(argc, argv);
      if (cmdLine.displayHelpRequired()){
	cout << "************************************" << endl;
	cout << "********** ComputeNorm.exe *************" << endl;
	cout << "************************************" << endl;
	cout << endl;
	cout << "Apply Z-T-ZT or ZT Norm to a Score File" << endl<< "ztnorm includes tnorm"<<endl;
        cout <<endl<<cc.getParamList()<<endl;
        return 0;  
      }
      if (cmdLine.displayVersionRequired()){
        cout <<"Version 3"<<endl;
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
      // Initialize the default values for the parameters TODO : add this option in the addParam functions...
 	  if(!config.existsParam("znormFilesExtension")) config.setParam("znormFilesExtension",".znorm");
 	  if(!config.existsParam("tnormFilesExtension")) config.setParam("tnormFilesExtension",".tnorm");
	  if(!config.existsParam("ztnormFilesExtension")) config.setParam("ztnormFilesExtension",".ztnorm");
	  if(!config.existsParam("tznormFilesExtension")) config.setParam("tznormFilesExtension",".tznorm");
	  if(!config.existsParam("maxIdNb")) config.setParam("maxIdNb","1000");
	  if(!config.existsParam("maxSegNb")) config.setParam("maxSegNb","1000");
	  if(!config.existsParam("maxScoreDistribNb")) config.setParam("maxScoreDistribNb","1000");
	  if(!config.existsParam("fieldGender")) config.setParam("fieldGender","0");
	  if(!config.existsParam("fieldName")) config.setParam("fieldName","1");
	  if(!config.existsParam("fieldDecision")) config.setParam("fieldDecision","2");	  	  	  	  	  	  	   	  
	  if(!config.existsParam("fieldSeg")) config.setParam("fieldSeg","3");	  	  	  	  	  	  	   	  
	  if(!config.existsParam("fieldLLR")) config.setParam("fieldLLR","4");	  
	  if(!config.existsParam("meanMode")) config.setParam("meanMode","0");
	  if(!config.existsParam("percentH")) config.setParam("percentH","0.0");
	  if(!config.existsParam("percentL")) config.setParam("percentL","0.0");
	  
	  
	  	  	  	  	  	  	   	  

      // start the prog!
      ComputeNorm(config);
    }
catch (alize::Exception& e) {cout << e.toString() << endl << cc.getParamList()<< endl;}
return 0;     
}



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
 
#include "SequenceExtractor.h"


int main(int argc, char* argv[])
{    
    try    
    {    
      // String inputTreeFilename=config.getParam("decoderFilename");  // The complete filename for the decoder tree
      // String outputFilename=config.getParam("outputFilename");      // The complete output filename
      // String inputFilename=config.getParam("inputFilename");        // The complete input filename 
    
      ConfigChecker cc;
      cc.addStringParam("config", false, true, "default config filename");
      cc.addStringParam("decoderFilename",true,true,"filename of the decoder tree");
      cc.addStringParam("inputFilename",true,true,"filename of the input symbol stream");
      cc.addStringParam("outputFilename",true,true,"filename of the output symbol file");
      cc.addStringParam("labelFilename",false,true,"if segment processing is needed, name of the label file");
      cc.addStringParam("labelSelectedFrames",false,true,"if segment processsing is needed, label of the selected segments");        
      cc.addFloatParam("frameLength",false,true,"If segment processing is required, it gives the length in s of an input symbol");
      cc.addBooleanParam("overlapMode",false,true,"if set to true (default=false), find overlap sequences. i.e decode one sequance by input symbol");
      cc.addStringParam("ngramOutputFilename",false,true,"If set, computes the ngram on the output seq and save it in this file");
      cc.addIntegerParam("ngramOutputOrder",false,true,"If output ngram is needed, gives the max order of the ngram");
      CmdLine cmdLine(argc, argv);
      if (cmdLine.displayHelpRequired())
	{ 
	  cout <<"SequenceDecoder.exe"<<endl<<"This program uses in input a decoder tree and an input stream"
	       <<endl<<"It decodes the input stream and outputs the finded sequences"<<endl<<cc.getParamList()<<endl;
	  cout <<"HelpRequired"<<endl; 
	  return 0;  
	}
      if (cmdLine.displayVersionRequired())
	{
	  cout <<"Version Beta "<<endl;
	}
      Config tmp;
      cmdLine.copyIntoConfig(tmp);
      Config config;
      if (tmp.existsParam("config")) config.load(tmp.getParam("config"));
      cmdLine.copyIntoConfig(config);
      cc.check(config);
      
      
      debug=config.getParam_debug();
      if (config.existsParam("verbose")) verbose=config.getParam("verbose").toBool();
      else verbose=false;
      sequenceDecoder(config);
	
     
    }
    catch (alize::Exception& e)
    {
        cout << e.toString() << endl;
    }

    return 0;
}

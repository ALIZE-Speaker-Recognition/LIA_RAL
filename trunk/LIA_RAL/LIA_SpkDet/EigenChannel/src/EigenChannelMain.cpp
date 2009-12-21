// EigenChannelMain.cpp
// 
// This file is a part of Mistral Package and LIA Software 
// LIA_SpkDet, based on Mistral_Ral toolkit 
// LIA_SpkDet is a free, open tool for speaker recognition
// LIA_SpkDet is a development project initiated and funded by the LIA lab.
//
// See mistral.univ-avignon.fr 
// 
// ALIZE is needed for LIA_SpkDet
// for more information about ALIZE, see http://alize.univ-avignon.fr
//
// Copyright (C) 2004 - 2005 - 2006 - 2007 -2008
//  Laboratoire d'informatique d'Avignon [www.lia.univ-avignon.fr]
//  Jean-Francois Bonastre [jean-francois.bonastre@univ-avignon.fr]
//      
// Mistral is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
// You should have received a copy of the GNU General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// The LIA team as well as the ALIZE project want to highlight the limits of voice authentication
// in a forensic context. 
// The following paper proposes a good overview of this point:
// [Bonastre J.F., Bimbot F., Boe L.J., Campbell J.P., Douglas D.A., Magrin-chagnolleau I.,
//  Person  Authentification by Voice: A Need of Caution,
//  Eurospeech 2003, Genova]
// The conclusion of the paper of the paper is proposed bellow:
// [Currently, it is not possible to completely determine whether the
//  similarity between two recordings is due to the speaker or to other
//  factors, especially when: (a) the speaker does not cooperate, (b) there
//  is no control over recording equipment, (c) recording conditions are not 
//  known, (d) one does not know whether the voice was disguised and, to a
//  lesser extent, (e) the linguistic content of the message is not
//  controlled. Caution and judgment must be exercised when applying speaker
//  recognition techniques, whether human or automatic, to account for these
//  uncontrolled factors. Under more constrained or calibrated situations,
//  or as an aid for investigative purposes, judicious application of these
//  techniques may be suitable, provided they are not considered as infallible.
//  At the present time, there is no scientific process that enables one to
//  uniquely characterize a person=92s voice or to identify with absolute
//  certainty an individual from his or her voice.]
//
// Contact Jean-Francois Bonastre (jean-francois.bonastre@lia.univ-avignon.fr) for
// more information about the licence or the use of LIA_SpkDet
// First version 15/07/2004
// New version 23/02/2005
// 
// Last review 4 nov 2008
//

#include <iostream>
#include <liatools.h>
#include <EigenChannel.h>

using namespace alize;
using namespace std;

int main(int argc, char* argv[]) {
	ConfigChecker cc;
	try {
		// Needed params
		cc.addStringParam("ndxFilename",true,true,"NDX of multiple GMM speaker recordings");
		cc.addStringParam("inputWorldFilename",true,true,"the world model file");
		cc.addIntegerParam("nbIt",true,true,"number of ml it");	
		cc.addStringParam("channelMatrix",true,true,"filename to save Channel Matrix ");					
		cc.addIntegerParam("channelMatrixRank",true,true,"final rank of channel matrix");	
		cc.addFloatParam("regulationFactor",true,true,"map tau");
		cc.addStringParam("saveMatrixFormat",true,true,"matrix format: DB (binary) or DT (ascii)");		  
		cc.addStringParam("loadMatrixFormat",true,true,"matrix format: DB (binary) or DT (ascii)");				

		// Optionnal
		cc.addStringParam("initChannelMatrix",false,true,"init Channel Matrix");	
		cc.addBooleanParam("loadAccs",false,true,"if true do not compute UBM stats, load matrices");		
		cc.addIntegerParam("computeLLK",false,true,"optional: nb of files where LLK is computed");	
		

		// Insertion of config compatibility rules
		CmdLine cmdLine(argc, argv);
	     if (cmdLine.displayHelpRequired()){
			cout << "****************************************" << endl;
			cout << "********** EigenChannel.exe ************" << endl;
			cout << "****************************************" << endl;
			cout << endl;
			cout << "Evaluate Channel Matrix from speakers data" << endl;
	  	   cout <<endl<<cc.getParamList()<<endl;
       	 return 0;  
	      }
      	if (cmdLine.displayVersionRequired()){
      	  cout <<"Version 2.0 Mistral Package"<<endl;
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
		EigenChannel(config);	
		}
	catch (Exception& e) {cout << e.toString() << cc.getParamList() << endl;}
if (debug) {
}
return 0;
}

// ReadFeatFileMain.cpp
//
// This file is a part of Mistral Software MISTRAL_SpkDet, based on ALIZE library 
// MISTRAL_SpkDet  is a free, open tool for speaker recognition
// MISTRAL_SpkDet is a development project initiated and funded by the ANR (agence nationale 
// de la recherche).
// See http://mistral.univ-avignon.fr
// 
// ALIZE is needed for MISTRAL_SpkDet
//
// Copyright (C) 2004, 2005, 2006, 2007, 2008
//  Laboratoire d'informatique d'Avignon [www.lia.univ-avignon.fr]
//  Jean-Francois Bonastre [jean-francois.bonastre@lia.univ-avignon.fr]
//      
// MISTRAL_SpkDet is free software; you can redistribute it and/or
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
// Contact Jean-Francois Bonastre (jean-francois.bonastre@lia.univ-avignon.fr) for
// more information about the licence or the use of LIA_SpkDet
// First version 15/07/2004
// This version 28/01/2008 - (edited by e.c / eric.charton@univ-avignon.fr)

#include <iostream>
#include "ReadFeatFile.h"
#include "alize.h"

/*!
	This program is used to read feature files and print their contents to STDOUT

	@param argc	
	@param argv[]  params for config checker
	@return 		send back a 0 value in all cases
*/
int main (int argc, char *argv[])
{
  using namespace std;
  using namespace alize;

  ConfigChecker cc;
  /*!
	For this tool, you must give following parameters: inputFeatureFileName, loadFeatureFileFormat, loadFeatureFileExtension
  */
  /*! 
	Sample command : ./ReadFeatFile --inputFeatureFileName ocpy_A --loadFeatureFileFormat SPRO4 --loadFeatureFileExtension .sph --featureFilesPath ./

  */
  cc.addStringParam("inputFeatureFileName",true, true,"input feature file name");
  cc.addStringParam("loadFeatureFileFormat",false,true,"ALIZE and MISTRAL_RAL option, see doc. Give a feature format ie spro, raw, etc");
  cc.addStringParam("loadFeatureFileExtension",false,true,"ALIZE MISTRAL_RAL option, see doc. Give a feature ext ie .prm, .raw, etc");
  cc.addStringParam("loadFeatureFilesPath",false,true,"ALIZE MISTRAL_RAL option, see doc. Give the feature path");
  try {
        CmdLine cmdLine(argc, argv);
	   if (cmdLine.displayVersionRequired()){
          cout <<"Version 1.0M"<<endl;
        } 
        if (cmdLine.displayHelpRequired()){
          cout <<"ReadFeatFile.exe"<<endl<<"This program is used to read feature files and print their contents to STDOUT" <<endl<<cc.getParamList()<<endl;
          return 0;  
        }
        
        Config tmp;
        cmdLine.copyIntoConfig(tmp);
        Config config;
        if (tmp.existsParam("config")) config.load(tmp.getParam("config"));
        cmdLine.copyIntoConfig(config);
        cc.check(config);
        readFeatFile(config);
        }
  catch (alize::Exception & e) {cout << e.toString () << endl << cc.getParamList() << endl;}
  return 0;
}

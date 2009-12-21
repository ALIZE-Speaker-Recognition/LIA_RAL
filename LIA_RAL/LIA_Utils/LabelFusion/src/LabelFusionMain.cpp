// LabelFusionMain.cpp
// This file is a part of LIA Softwares LIA_SpkDet and LIA_SpkSeg, based on ALIZE toolkit 
// LIA_SpkDet is a free, open tool for speaker recognition
// LIA_SpkSeg is a free, open tool for speaker segmentation
// This project is a development project initiated and funded by the LIA lab.
// See www.lia.univ-avignon.fr
// 
// ALIZE is needed for LIA_SpkDet and LIA_SpkSeg
// for more information about ALIZE, see http://www.lia.univ-avignon.fr/heberges/ALIZE/
// First version March 26 2005
//
// Copyright (C) 2005
//  Laboratoire d'informatique d'Avignon [www.lia.univ-avignon.fr]
// Main authors
//  Jean-Francois Bonastre [jean-francois.bonastre@lia.univ-avignon.fr]
//      
// LIA_SpkDet and LIA_SpkSeg are free software; you can redistribute it and/or
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
// more information about the licence or the use of LIA_SpkDet or LIA_SpkSeg
#include <iostream>
#include "LabelFusion.h"


int main(int argc, char* argv[])
{
using namespace std;
using namespace alize;
ConfigChecker cc;
cc.addStringParam("config",false, true,"default config file");
cc.addStringParam("saveLabelFileExtension",true, true,"produced labelFiles extension");
cc.addStringParam("labelOutputPath",false,true,"path where to store produced labelFiles");
cc.addStringParam("labelInputPath",false,true,"path where to load original labelFiles");
cc.addIntegerParam("winLength",false,true,"window of the morphological filter in frames");
cc.addFloatParam("selectThreshold",false,true,"applying filter when selected/unselected matches this threshold");
cc.addStringParam("mode",false,true,"decide the mode: fusion(2 files) or morphing (one file) default=fusion");
cc.addStringParam("labelOneFilename",false,true,"if fusion is selected, name of the first feature file");
cc.addStringParam("labelTwoFilename",false,true,"if fusion is selected, name of the second feature file");
cc.addStringParam("labelFilename",false,true,"if morphing is selected, name of the feature file");

  try {
    CmdLine cmdLine(argc, argv);
    if (cmdLine.displayHelpRequired()){
      cout <<"LabelFusion.exe"<<endl<<"This program is used to merge labelFiles and to apply a morphological window" <<endl<<cc.getParamList()<<endl;
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
    bool fusion=true;
    if (config.existsParam("mode")) fusion=(config.getParam("mode")=="fusion");
    if (fusion) labelFusion(config);
    	else labelMorphing(config);	
  }
catch (alize::Exception & e) {cout << e.toString () << endl << cc.getParamList() << endl;}
return 0;
}

// NormFeatMain.cpp
// This file is a part of LIA Software LIA_SpkDet, based on ALIZE toolkit 
// LIA_SpkDet  is a free, open tool for speaker recognition
// LIA_SpkDet is a development project initiated and funded by the LIA lab.
// See www.lia.univ-avignon.fr
// 
// ALIZE is needed for LIA_SpkDet
// for more information about ALIZE, see http://www.lia.univ-avignon.fr/heberges/ALIZE/
//
// Copyright (C) 2004
//  Laboratoire d'informatique d'Avignon [www.lia.univ-avignon.fr]
//  Jean-Francois Bonastre [jean-francois.bonastre@lia.univ-avignon.fr]
//      
// LIA_SpkDet is free software; you can redistribute it and/or
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
// Augmented version 20/03/2006 (JFB) - adding window based norm/cleaning
#include <iostream>
#include "NormFeat.h"
#include <alize.h>
#include <liatools.h>

int main (int argc, char *argv[])
{
  using namespace std; 
  using namespace alize;

    ConfigChecker cc;
    cc.addStringParam("config", false, true, "default config filename");
    cc.addStringParam("mode", true, true, "<norm: mean, var, warp normalization | info: save stats of file | featmap: feature mapping | FA EigenChannel Compensation >");
    cc.addStringParam("inputFeatureFilename",true,true,"input feature - could be a simple feature file or a list of filename");
    cc.addStringParam("labelSelectedFrames",true,true,"Only the frames from segments with this label  will be used");
    cc.addBooleanParam("fileMode",false,true,"One normalisation for all the selected data, could be applied after segmental/window modes if both are selected");  
    cc.addBooleanParam("segmentalMode",false,true,"if set to true, one normalisation by segment is computed");  
    cc.addBooleanParam("windowMode",false,true,"if set to true, the normalisation is computed by windows");  
    cc.addFloatParam("windowDuration",false,true,"Length of the window (for window mode) in seconds (default=5)");
    cc.addFloatParam("windowComp",false,true,"compensation factor for window mode (default =1.0 - no compensation");
    cc.addBooleanParam("writeAllFeatures",false,true,"if set to true,values for all the input frames are outputed (default true)");  
    cc.addStringParam("featureServerMode",true,true,"FEATURE_WRITABLE to write normalized features");  
    cc.addStringParam("outputInfoFilename",false,true,"The complete output filename for the info file");
    cc.addBooleanParam("cmsOnly",false,true,"If cmsOnly is set, only remove means from features (default false)");
    cc.addBooleanParam("varOnly",false,true,"If cmsOnly is set, reduce only variance from features (default false)");
    cc.addBooleanParam("warp",false,true,"If warp is set, uses featureWarping (default false)");
    cc.addIntegerParam("warpBinCount",false,true,"If warp is set, this parameter fixes the number of bins in the data histo (default=40)");
    cc.addStringParam("warpGaussHistoFilename",false,true,"for warping, if this parameter is set, it loads the destination histo, default it creates it (N(0,1))");
    cc.addBooleanParam("warpAdd01Norm",false,true,"for warping, add a global N(0,1) normalisation after warping (default false)");  
    cc.addStringParam("externalStatsFilename",false,true,"filename containing external stats to apply for normalization, usually obtained with info mode (opt.)");
    cc.addFloatParam("frameLength",false,true,"length of a frame, by default 10ms");
  try
  {
    CmdLine cmdLine (argc, argv);
    if (cmdLine.displayHelpRequired ()){	// --help
      cout << "NormFeat" << endl;
      cout << "NormFeat.exe --config <foo.cfg> --inputFeatureFileName <foo.prm> --mode <norm|info|featMap>"<< endl;
      cout<<cc.getParamList()<<endl;
      }
    else if (cmdLine.displayVersionRequired ())	// --version
      cout << "Version 2.0" << endl;
    else{
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
      String mode=config.getParam("mode");
      if ( mode == "norm"){normFeat(config);}
      else if ( mode == "info"){infoFeat(config);}
      else if (mode == "featMap"){featMap(config);}
      else if (mode == "NAP"){normFeatNAP(config);}      
      else if (mode == "FA"){normFeatFA(config);}      
      else throw Exception("Error: Mode unknown!",__FILE__,__LINE__);
    }
  }
  catch (alize::Exception & e) {cout << e.toString () << endl << cc.getParamList() << endl;}
  #ifdef NDEBUG 
  cout<<"*** Objects created and destroyed **"<<Object::getCreationCounter()<<"-"<<Object::getDestructionCounter()<<endl;    
  #endif
  return 0;
}

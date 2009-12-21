// LIA_UTILS
//
// This file is a part of Mistral Software LIA_Utils, based on Mistral_Ral toolkit 
// Mistral_Ral  is a free, open tool for speaker recognition
// Mistral_Ral is a development project initiated and funded by the LIA lab.
//
// See mistral.univ-avignon.fr 
// 
// ALIZE is needed for Mistral_SpkDet
// for more information about ALIZE, see http://alize.univ-avignon.fr
//
// Copyright (C) 2004 - 2005 - 2006 - 2007 -2008
//  Laboratoire d'informatique d'Avignon [www.lia.univ-avignon.fr]
//  Jean-Francois Bonastre [jean-francois.bonastre@lia.univ-avignon.fr]
//      
// Mistral and Mistral_SpkDet is free software; you can redistribute it and/or
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
//
// Contact Jean-Francois Bonastre (jean-francois.bonastre@lia.univ-avignon.fr) for
// more information about the licence or the use of LIA_SpkDet
//
// Reviewed on 5 nov 2008 by eric
//
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

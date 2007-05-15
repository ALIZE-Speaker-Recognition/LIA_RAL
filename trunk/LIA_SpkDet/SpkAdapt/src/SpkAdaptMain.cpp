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
//Author : Alexandre PRETI.
#include <iostream>
#include "SpkAdapt.h"
#define MANDATORY true
#define OPTIONAL false
#define ARG_REQUIRED true
#include "liatools.h"

int main(int argc, char *argv[])
{
  try
  {
    ConfigChecker cc;
    cc.addIntegerParam("verboseLevel", false, true,
      "level of the verbose information 0=no verbose, 1=normal, 2=more");
    cc.addStringParam("config", OPTIONAL, ARG_REQUIRED,
      "Train config file name ");
    cc.addStringParam("configTest", OPTIONAL, ARG_REQUIRED,
      "Unsupervised adaptation config file name ");
    cc.addStringParam("targetIdList", MANDATORY, ARG_REQUIRED,
      "targets list");
    cc.addStringParam("inputWorldFilename", MANDATORY, ARG_REQUIRED,
      "world model");
    cc.addStringParam("testsIdList", MANDATORY, ARG_REQUIRED, "tests list");
    cc.addBooleanParam("NISTprotocol", MANDATORY, ARG_REQUIRED,
      "true to follow unsupervised NIST protocol, false to follow BATCH protocol");

    cc.addBooleanParam("WMAP", false, false,
      "NO MORE USED, Choice of WMAP (one Gaussian) to compute trial feature server weights, need other param : MUtarget,MUimp,SIGMAtarget and SIGMAIMP;TAR weight and IMPweight, thrMin and thrMax");
    cc.addBooleanParam("WMAPGMM", false, false,
      "Choice of WMAP GMM to compute trial feature server weights, param needed : two gmm of scores, TARScoresModel and NONScoresModel. Need also a value for Priors update, initPriorImp (10) and initPriorTar(1) ");
    cc.addBooleanParam("TNORM", false, false,
      " TNORM scores before computing WMAP, need gmm of Tnormed scores and a param : impScoreFile, score file for TNORM and testsNames, the list of trial segments ");
    cc.addBooleanParam("FAST", false, false,
      "For computing LLR quicker, use top ten info files");
    cc.addBooleanParam("FromResFile", false, false,
      " Avoid to compute LLR for WMAP computing, search LLR in a score file, param needed: InputResFilename");
    cc.addBooleanParam("TrainConfident", false, false,
      " A try?! : compute LLR of a percentage of train data on a model learnt on a percentage of train data AverageIt times and average the results, need AverageIt, SelectedTest and SelectedTrain. result used for WMAP of train data on client model ");
    cc.addBooleanParam("Oracle", false, false,
      " For Oracle (supervised) experiments, need targetTests to perform adaptation on target tests only");
    cc.addBooleanParam("REGRESS", false, false,
      "Use logistic regression for adaptation weights computing, need BETA and THETA");
    cc.addIntegerParam("MaxMixturesCount", false, true,
      "max mixtures stored in the mixture server (in memory)");
    cc.addIntegerParam("LLKthreshold", false, true,
      "thresholds the results of LLR in order to avoid problem in WMAP computing");
    cc.addStringParam("computeLLKWithTopDistribs", true, true,
      "PARTIAL/COMPLETE: will compute LLK with topdistribs. COMPLETE: add world LLK to client LLK, PARTIAL: just keeps the topDistrib LLK");
    cc.addIntegerParam("topDistribsCount ", false, true,
      "Number of distrib to approximate complete LLK");
    cc.addStringParam("InfoExtension", false, true,
      "Extension for top ten info files");
    cc.addStringParam("InfoPath", false, true, "Path for top ten info files");
    cc.addFloatParam("OptimalScore", false, true,
      "Used for updating  priors, if LLR is superior add one for target count, else add one for impostors count ");
    Config tmp;
    CmdLine cmdLine(argc, argv);
    if (cmdLine.displayHelpRequired())
      {				// --help
	cout << "************************************" << endl;
	cout << "********** SpkAdapt.exe *********" << endl;
	cout << "************************************" << endl;
	cout << endl;
	cout <<
	  "Unsupervised adaptation Train target, update model using test trial information following NIST or BATCH protocol"
	  << endl;
	cout << "" << endl;
	cout << endl;
	cout << "Command Line example: " << endl;
	cout <<
	  "SpkAdapt.exe --config cfg/TrainTarget.cfg --configTest cfg/TrainTest.cfg --inputWorldFilename world --targetIdList ./lst/id.lst --testsIdList ./lst/tests.list --outputLLRFilename file --NISTprotocol true"
	  << endl;
	cout << cc.getParamList() << endl;
      }
    else
      {
	cmdLine.copyIntoConfig(tmp);
	Config config(tmp.getParam("config"));
	cmdLine.copyIntoConfig(config);
	Config configTest(tmp.getParam("configTest"));
	if (config.existsParam("verbose"))
	  verbose = config.getParam("verbose").toBool();
	else
	  verbose = false;
	if (verbose)
	  verboseLevel = 1;
	else
	  verboseLevel = 0;
	if (config.existsParam("verboseLevel"))
	  verboseLevel = config.getParam("verboseLevel").toLong();
	if (verboseLevel > 0)
	  verbose = true;
	TrainTargetAdapt(config, configTest);	// Training the target models

      }
  }
  catch(alize::Exception & e)
  {
    cout << e.toString() << endl;
  }

  return 0;
}

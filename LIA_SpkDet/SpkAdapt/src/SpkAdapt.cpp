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
#if !defined(ALIZE_SpkAdapt_cpp)
#define ALIZE_SpkAdapt_cpp

#include <iostream>
#include <fstream>		// pour outFile
#include <cstdio>		// pour printf()
#include <cassert>		// pour le debug pratique
#include <cmath>
#include <liatools.h>
#include "SpkAdapt.h"
#include "FileInfo.h"
#include "UnsupervisedTools.h"

using namespace alize;
using namespace std;



// Training of client Speakers
// Input: Xlist Format: ID_Client Seg1 Seg2 ..
// Output: ALIZE_MixtureServer (binaire) + GMM / Client (binary)
int TrainTargetAdapt(Config & config, Config & configTest)
{
  String outputLLRFileName = config.getParam("outputLLRFilename");
  ofstream outLLR(outputLLRFileName.c_str(), ios::out | ios::trunc);
  String fullFileName;
  bool NIST = false;
  if (config.existsParam("NISTprotocol"))
    NIST = config.getParam("NISTprotocol").toBool();
  bool becomeCurrent = true;
  String inputClientListFileName = config.getParam("targetIdList");
  String inputTestsListFileName = config.getParam("testsIdList");
  String inputWorldFilename = config.getParam("inputWorldFilename");
  String outputSERVERFilename = config.getParam("mixtureServer");
  double LLRRatio = 0.0, LLRRatio2 = 0.0;
  bool regress = false;
  if (configTest.existsParam("REGRESS"))
    regress = configTest.getParam("REGRESS").toBool();
  bool wmap = false;
  if (configTest.existsParam("WMAP"))
    wmap = configTest.getParam("WMAP").toBool();
  bool wmapGmm = false;
  if (configTest.existsParam("WMAPGMM"))
    wmapGmm = configTest.getParam("WMAPGMM").toBool();
  bool Fast = false;
  if (configTest.existsParam("FAST"))
    Fast = configTest.getParam("FAST").toBool();
  int countTests;
  bool saveEmptyModel = false;
  if (config.existsParam("saveEmptyModel"))
    saveEmptyModel = config.getParam("saveEmptyModel").toBool();

  bool tnorm = false;
  if (configTest.existsParam("TNORM"))
    tnorm = configTest.getParam("TNORM").toBool();

  bool oracle = false;
  if (configTest.existsParam("Oracle"))
    oracle = configTest.getParam("Oracle").toBool();
  bool confident = false;
  if (configTest.existsParam("TrainConfident"))
    confident = configTest.getParam("TrainConfident").toBool();
  bool fromResFile = false;
  if (configTest.existsParam("FromResFile"))
    fromResFile = configTest.getParam("FromResFile").toBool();

  String inputResFilename = configTest.getParam("InputResFilename");


  ObjectRefVector stockTnorm;	//To stock Tnrom parameters, not need to be declared if Tnorm=false!
  DoubleVector decision;	//To stock the LK for the client models with the impostor data !!!!!!!CAPACITE A CHANGER!!!!!!!!
  DoubleVector tmp;		//For NIST protocol 
  DoubleVector NbFramesSelected;
  ObjectRefVector FeatServ;	//To stock the feature servers
  ObjectRefVector SegServ;	//To stock the seg clusters 

  if (tnorm)
    {				//load Tnorm parameters
      String testFile = (config.getParam("impScoreFile"));
      String TestNames = (config.getParam("testNames"));
      loadTnormParam(TestNames, testFile, stockTnorm, config);
    }

  // label for selected frames - Only the frames associated with this label, in the label files, will be used
  bool fixedLabelSelectedFrame = true;
  String labelSelectedFrames;
  if (config.existsParam("useIdForSelectedFrame"))	// the ID of each speaker is used as labelSelectedFrame ?
    fixedLabelSelectedFrame =
      (config.getParam("useIdForSelectedFrame") == "false");
  if (fixedLabelSelectedFrame)	// the label is decided by the command line and is unique for the run
    labelSelectedFrames = config.getParam("labelSelectedFrames");
  bool saveCompleteServer = false;
  debug = config.getParam_debug();
  if (config.existsParam("verbose"))
    verbose = true;
  else
    verbose = false;

  XList inputClientList(inputClientListFileName, config);	// read the Id + filenames for each client
  XLine *linep;
  inputClientList.getLine(0);
  XList inputTestsList(inputTestsListFileName, configTest);
  XLine *lineTests;
  inputTestsList.getLine(0);
  MixtureServer ms(config);
  Config configScores(config);	//Later for GMM on scores
  configScores.setParam("VectSize", "1");
  MixtureServer msScores(configScores);
  StatServer ss(config, ms);
  StatServer ssScores(configScores, msScores);
  XLine testsToCompute;

  if (verbose)
    cout << "TrainTarget - Load world model [" << inputWorldFilename << "]"
      << endl;
  MixtureGD & world = ms.loadMixtureGD(inputWorldFilename);
  //      load TAR and NON GMM
  MixtureGD & non =
    msScores.loadMixtureGD(configTest.getParam("NONScoresModel"));
  MixtureGD & tar =
    msScores.loadMixtureGD(configTest.getParam("TARScoresModel"));
  //temporate mixtures used in computeMAPmodelFromEMones, create here to accelerate
  /*MixtureGD & tmpMixture = ms.createMixtureGD(world.getDistribCount());       
     MixtureGD & auxMixture = ms.createMixtureGD(world.getDistribCount()); */
  String idAux = (ms.createMixtureGD(world.getDistribCount())).getId();
  String idTmp = (ms.createMixtureGD(world.getDistribCount())).getId();

  // *********** Target loop ***************** 
  while ((linep = inputClientList.getLine()) != NULL)
    {				// linep gives the XLine with the Id of a given client and the list of files
      MixtureGD & world = ms.getMixtureGD(ms.getMixtureIndex(inputWorldFilename));	//   load needed for each client because of the mixture server reset 
      MixtureGD & tmpMixture = ms.getMixtureGD(ms.getMixtureIndex(idAux));
      MixtureGD & auxMixture = ms.getMixtureGD(ms.getMixtureIndex(idTmp));
      String *id = linep->getElement();	// Get the Client ID (id)
      XLine featureFileListp = linep->getElements();	// Get the list of feature file for the client (end of the line)
      if (verbose)
	cout << "Train model [" << *id << "]" << endl;
      if (!fixedLabelSelectedFrame)
	{			// the ID is used as label for selecting the frame
	  labelSelectedFrames = *id;
	  if (debug)
	    cout << *id << " is used for label selected frames" << endl;
	}

      FeatureServer & fs = *new FeatureServer(config, featureFileListp);	// Reading the features (from several files)
      SegServer & segmentsServer = *new SegServer();	// Create the segment server for managing the segments/clusters
      LabelServer labelServer;	// Create the lable server, for indexing the segments/clusters
      initializeClusters(featureFileListp, segmentsServer, labelServer, config);	// Reading the segmentation files for each feature input file
      verifyClusterFile(segmentsServer, fs, config);	// Verify if the segments ending before the end of the feature files...
      MixtureGD & clientMixture = ms.duplicateMixture(world, DUPL_DISTRIB);
      MixtureGD & clientMixtureEM = ms.duplicateMixture(world, DUPL_DISTRIB);
      long codeSelectedFrame = labelServer.getLabelIndexByString(labelSelectedFrames);	// Get the index of the cluster with in interest audio segments
      if (codeSelectedFrame == -1)
	{			// No data for this model !!!!!!!!!!!!!!
	  cout << " WARNING - NO DATA FOR TRAINING [" << *id << "]";
	}
      else
	{
	  SegCluster & selectedSegments = segmentsServer.getCluster(codeSelectedFrame);	// Gives the cluster of the selected/used segments                                   
	  countTests = 0;	//reset the counter of tests segments (  train data included)
	  if (!Fast)
	    {
	      //add train data to objectrefvector
	      FeatServ.addObject(fs);
	      SegServ.addObject(selectedSegments);

	    }
	  adaptModel(config, ss, ms, fs, selectedSegments, world,
	    clientMixture);
	  adaptModelEM(config, ss, ms, fs, selectedSegments, world,
	    clientMixtureEM);
	  ms.setMixtureId(clientMixtureEM, *id);
	  testsToCompute.addElement(*id);
	  //compute the LLR with the client train data and store it to compute probability
	  if (confident)
	    {			//Method to set a WMAP weight for the train data
	      LLRRatio =
		computeLLRForTrain(configTest, ss, ms, fs, selectedSegments,
		world);
	      if (tnorm)
		TnormalizeScore(*id, LLRRatio, stockTnorm);	//should add the train data in Tnorm computing
	    }
	  else
	    {
	      //Fast LLR is not used for L(data train/ client model)
	      LLRRatio = computeLLR(ss, fs, world, clientMixture, selectedSegments);	//NO more needed if proba is fixed to 1

	    }
	  decision.addValue(LLRRatio);	//L(data train/ client model) always the first element of the doublevector

	  countTests++;
	  NbFramesSelected.addValue(SegClusterFrame(selectedSegments));	//stock the number of frames selected (used for fuseModels)
	  MixtureGD & SAVclientMixture =
	    ms.duplicateMixture(clientMixture, DUPL_DISTRIB);

	/***********************************************************************************/
	  /*  Adaptation   */
	/***********************************************************************************/

	  // ***********  tests loop ***************** 
	  while ((lineTests = inputTestsList.getLine()) != NULL)
	    {
	      String *idclient = lineTests->getElement();
	      if (*idclient == *id)
		{		//file format : column client model name, column gender, column tests in one word (i.e : toto_A)
		  MixtureGD & testMixtureEM =
		    ms.duplicateMixture(world, DUPL_DISTRIB);
		  String & idTest = lineTests->getElement(2, becomeCurrent);	// Get the Tests ID (id)
		  XLine featureFileListTests = lineTests->getElements();	// Get the list of feature file for the tests (end of the line)
		  if (verbose)
		    cout << "Train model for test [" << idTest << "]" << endl;
		  if (!fixedLabelSelectedFrame)
		    {		// the ID is used as label for selecting the frame
		      labelSelectedFrames = idTest;
		      if (debug)
			cout << idTest << " is used for label selected frames"
			  << endl;
		    }
		  FeatureServer & fsTests =
		    *new FeatureServer(configTest, featureFileListTests);
		  SegServer & segmentsServerTests = *new SegServer();	// Create the segment server for managing the segments/clusters
		  LabelServer labelServerTests;	// Create the lable server, for indexing the segments/clusters
		  initializeClusters(featureFileListTests, segmentsServerTests, labelServerTests, configTest);	// Reading the segmentation files for each feature input file
		  verifyClusterFile(segmentsServerTests, fsTests, configTest);	// Verify if the segments ending before the end of the feature files...
		  long codeSelectedFrame = labelServerTests.getLabelIndexByString(labelSelectedFrames);	// Get the index of the cluster with in interest audio segments
		  if (codeSelectedFrame == -1)
		    {		// No data for this model !!!!!!!!!!!!!!
		      cout << " WARNING - NO DATA FOR TRAINING [" << idTest
			<< "]";
		    }
		  else
		    {
		      SegCluster & selectedSegmentsTests = segmentsServerTests.getCluster(codeSelectedFrame);	// Gives the cluster of the selected/used segments                                   
		      //add tests data to objectrefvector
		      if (!Fast)
			{
			  FeatServ.addObject(fsTests);
			  SegServ.addObject(selectedSegmentsTests);

			}

		      NbFramesSelected.addValue(SegClusterFrame(selectedSegmentsTests));	//stock the number of frames selected (used for fuseModels)
		      if (Fast)
			{
			  testsToCompute.addElement(idTest);

			  int index = ms.getMixtureIndex(idTest);
			  if (index == -1)
			    {
			      adaptModelEM(configTest, ss, ms, fsTests, selectedSegmentsTests, world, testMixtureEM);	//create the EM model for the test 
			      ms.setMixtureId(testMixtureEM, idTest);
			    }	// Set the client model Id
			  else
			    testMixtureEM = ms.getMixtureGD(index);
			}
		      //Compute LLR and write it on a res file
		      LLRRatio = 0.0;
		      fullFileName = getFullFileName(idTest, configTest);
		      //TAKE CARE : by using it you should have 2 GMM of scores and a TNORM RES file according to the system used for computing the RES file 
		      if (fromResFile)
			{
			  LLRRatio =
			    searchLLRFromResFile(*id, idTest,
			    inputResFilename, configTest);
			}
		      else
			{
			  //Choice wether to use FastLLR or LLR (FastLLR is used when the top component are already computed for this Feature Server)
			  if (FileExists(fullFileName) == false
			    || Fast == false)
			    {	//first time seeing this test
			      LLRRatio = computeLLR(configTest, ss, fsTests, world, SAVclientMixture, selectedSegmentsTests, fullFileName);	//For proba, take the decision on the original target model
			    }
			  else
			    LLRRatio = computeFastLLR(ss, fsTests, world, SAVclientMixture, selectedSegmentsTests, fullFileName, configTest);	//For proba, take the decision on the original target model
			}
		      if (tnorm)
			TnormalizeScore(idTest, LLRRatio, stockTnorm);

		      if (NIST)
			{
			  LLRRatio2 = 0.0;
			  if (FileExists(fullFileName) == false
			    || Fast == false)
			    {	//first time seeing this test
			      LLRRatio2 =
				computeLLR(configTest, ss, fsTests, world,
				clientMixture, selectedSegmentsTests,
				fullFileName);
			    }
			  else
			    LLRRatio2 = computeFastLLR(ss, fsTests, world, clientMixture, selectedSegmentsTests, fullFileName, configTest);	//for LLR output 
			  if (tnorm)
			    TnormalizeScore(idTest, LLRRatio2, stockTnorm);
			  outLLR << "M " << *id << " 0 " << idTest << " " << LLRRatio2 << endl;	// Output the result
			}
		      else
			outLLR << "M " << *id << " 0 " << idTest << " " << LLRRatio << endl;	// Output the result
		      if (oracle)
			Oracle(*id, idTest, LLRRatio, configTest, tar, non, ssScores);	//if Oracle == true , DO NOT SET regress, wmap or wmapgmm to true!!!
		      decision.addValue(LLRRatio);
		      countTests++;
				  /**********************************NIST protocol*******************************/
		      if (NIST)	//compute the client model now
			{
			  tmp = decision;
			  //compute FS weights
			  if (regress)
			    expandLLR(decision, configTest);
			  else if (wmap)
			    WMAP(decision, configTest);
			  else if (wmapGmm)
			    WMAPGMM(decision, configTest, tar, non, ssScores);
			  if (!confident)
			    decision[0] = 1;	//SET a weight of 1 for the train data
			  //create the new client model 
			  if (Fast)
			    computeMAPmodelFromEMones(config, ss, ms,
			      NbFramesSelected, world, clientMixture,
			      auxMixture, tmpMixture, decision,
			      testsToCompute);
			  else
			    {
			      clientMixture = world;	//reset clientMixture to MAP from world
			      adaptModel(config, ss, ms, FeatServ, SegServ,
				world, clientMixture, decision);
			    }
			  decision = tmp;
			  tmp.clear();
			}
		    }
		  if (Fast)
		    {		//delete object created with new, not nedded anymore
		      delete & segmentsServerTests;
		      delete & fsTests;
		    }
		}
	    }			//end tests loop
	  if (!NIST)
	    {
	      //compute the probability associated with each feature server 
	      if (regress)
		expandLLR(decision, configTest);
	      else if (wmap)
		WMAP(decision, configTest);
	      else if (wmapGmm)
		WMAPGMM(decision, configTest, tar, non, ssScores);
	      if (!confident)
		decision[0] = 1;	//SET a weight of 1 for the train data
	      //create the new client model 
	      if (Fast)
		computeMAPmodelFromEMones(config, ss, ms, NbFramesSelected,
		  world, clientMixture, auxMixture, tmpMixture, decision,
		  testsToCompute);
	      else
		{
		  clientMixture = world;	//reset clientMixture to MAP from world
		  adaptModel(config, ss, ms, FeatServ, SegServ, world,
		    clientMixture, decision);
		}
	    }
	  if (verbose)
	    cout << "Save client model [" << *id << "]" << endl;
	  clientMixture.save(*id, config);	// Save the client model
	  inputTestsList.rewind();	//rewind the xlist, needed when changing target
	  if (!saveCompleteServer)
	    {

	      //delete all client Mixtures 
	      ms.deleteMixture(clientMixture);
	      ms.deleteMixture(clientMixtureEM);
	      ms.deleteMixture(SAVclientMixture);
	      NbFramesSelected.clear();
	      decision.clear();
	      if (!Fast)
		{
		  FeatServ.deleteAllObjects();
		  SegServ.deleteAllObjects();
		}
	      testsToCompute.reset();
	      if (Fast)
		{		//delete object created with new, not nedded anymore
		  delete & fs;
		  delete & segmentsServer;

		}
	    }
	  //reset mixtureServer to free memory (reload needed for several mixtures)
	  if (ms.getMixtureCount() >
	    (unsigned long) configTest.getParam("MaxMixturesCount").toLong())
	    {
	      if (verbose)
		cout << "Reset MixtureServer and reload  of the used mixtures"
		  << endl;
	      ms.reset();
	      MixtureGD & world = ms.loadMixtureGD(inputWorldFilename);	//load needed for each client because of the mixture server reset 
	      //temporate mixtures used in computeMAPmodelFromEMones, create here to accelerate
	      idAux = (ms.createMixtureGD(world.getDistribCount())).getId();
	      idTmp = (ms.createMixtureGD(world.getDistribCount())).getId();


	    }
	}
    }				// end of the the target loop 

//When end delete all mixture in the server
  ms.deleteMixtures(0, ms.getMixtureCount());
  ms.deleteUnusedDistribs();
  if (tnorm)
    stockTnorm.deleteAllObjects();;
  outLLR.close();



  return 0;
}




#endif //!defined(ALIZE_TrainTarget_cpp)

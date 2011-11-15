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

#if !defined(ALIZE_SpkAdapt_cpp)
#define ALIZE_SpkAdapt_cpp

#include <iostream>
#include <fstream>        // pour outFile
#include <cstdio>        // pour printf()
#include <cassert>        // pour le debug pratique
#include <cmath>
#include <liatools.h>
#include "SpkAdapt.h"
#include "FileInfo.h"
#include "UnsupervisedTools.h"

using namespace alize;
using namespace std;

/************************************************************************/
/*
UNSUPERVISED ADAPTATION PROCESS

NEED THE TRAIN TARGET CONFIG FILE + NDX
NEED A CONFIG FOR TRIALS ADAPTATION (SEE CONFIG CHECKER) + NDX, SYNTAX (LINE  =     MODEL GENDER, TRIAL)  
-- IF WMAPGMM MODE : NEED TWO GMM OF SCORES
-- IF SCORES BY CLIENT MODE : NEED A SET OF TAR GMM LEARNT ON ANOTHER DATABASE AND 
A SET OF NON GMM LEARN BY COMPUTING LLR OF SOME IMPOSTORS ON CLIENTS MODELS
-- SEGMODE NOT DEBUGGED

FAST MODE : FAST LLR GIVES ABOUT THE SAME RESULTS AS CLASSICAL LLR (10-3)

*/
/*************************************************************************/




int TrainTargetAdapt(Config & config, Config & configTest)
{
   //LOAD PARAM FROM CONFIG
  String outputLLRFileName  =     config.getParam("outputLLRFilename");
  ofstream outLLR(outputLLRFileName.c_str(), ios::out | ios::trunc);
  String fullFileName,fullFileNameTrain;
	
  // CHOICE OF THE PROTOCOL BATCH OR NIST
   bool NIST  =     true; 
  if (config.existsParam("NISTprotocol"))        
    NIST  =     config.getParam("NISTprotocol").toBool();
  
  bool becomeCurrent  =     true,
       selected	      =     true;
  
  //LOAD LISTS
  String inputClientListFileName    =     config.getParam("targetIdList");
  String inputTestsListFileName     =     config.getParam("testsIdList");
  String inputWorldFilename         =     config.getParam("inputWorldFilename");
  String outputSERVERFilename       =     config.getParam("mixtureServer");
  
  double LLRRatio   =     0.0, 
         LLRRatio2  =     0.0,
         LLRTrain   =     0.0;
         //previousLLR =    -10000;        //NEW
  
  //FOR LOGISTIC REGRESSION ON SCORES
  bool regress  =     false;
  if (configTest.existsParam("REGRESS"))
    regress  =     configTest.getParam("REGRESS").toBool();
  
  //SIMPLE WMAP CURVE (1 GAUSSIAN)
  bool wmap  =     false;
  if (configTest.existsParam("WMAP"))
    wmap  =     configTest.getParam("WMAP").toBool();
  
  //WMAP USING GMM
  bool wmapGmm  =     false;
  if (configTest.existsParam("WMAPGMM"))
    wmapGmm  =     configTest.getParam("WMAPGMM").toBool();
  
  //FAST MODE USE TOP DISTRIB FOR LLR COMPUTATION AND FUSEMODEL FOR UPDATE OF MODEL
  bool Fast  =     false;
  if (configTest.existsParam("FAST"))
    Fast  =     configTest.getParam("FAST").toBool();
  int countTests;
  bool saveEmptyModel  =     false;
  if (config.existsParam("saveEmptyModel"))
    saveEmptyModel  =     config.getParam("saveEmptyModel").toBool();
  
  //TNORM SCORES, USES TNORM SCORES FOR WMAP
  bool tnorm  =     false;
  if (configTest.existsParam("TNORM"))
    tnorm  =     configTest.getParam("TNORM").toBool();
  
  //ZNORM SCORES, USES TNORM SCORES FOR WMAP
  bool znorm  =  false, ztnorm = false;
  if (configTest.existsParam("ZNORM"))
    znorm  =     configTest.getParam("ZNORM").toBool();
  
  //SUPERVISED MODE, NEED A LIST OF TRUE TARGET CLIENT (targetTests)
  bool oracle  =     false;
  if (configTest.existsParam("Oracle"))
    oracle  =     configTest.getParam("Oracle").toBool();
  
  //CROSS VALIDATION OF TEST DATA
  bool cross  =     false;
  if (configTest.existsParam("CrossValid"))
    cross  =     configTest.getParam("CrossValid").toBool();
  
  //USES A SCORES FILE FOR LLR OF TESTS ON CLIENT MODEL, A WAY OF USING ANOTHER DECISION SYSTEM
  String inputResFilename;
  bool fromResFile  =     false;
  if (configTest.existsParam("FromResFile")){
    fromResFile  =     configTest.getParam("FromResFile").toBool();
    inputResFilename =     configTest.getParam("InputResFilename");
  }
  //NEEDEED FOR STORING A GMM AFTER EACH ADAPTATION ITERATION
  bool saveAfterIt  =     false;
  if (configTest.existsParam("SaveAfterIt"))
    saveAfterIt  =     configTest.getParam("SaveAfterIt").toBool();
  
  //GMM OF TAR AND NON FOR EACH CLIENT INSTEAD OF GLOBAL ONES  
  bool ScoresByTarget  =     false;
  if (configTest.existsParam("ScoresByTarget"))
    ScoresByTarget  =     configTest.getParam("ScoresByTarget").toBool();
  
   //TRY USING TRAIN DATA TO CONTROLL ADAPTATION
  bool assess  =     false;
  if (configTest.existsParam("Assess"))
    assess  =     configTest.getParam("Assess").toBool();
  
  //FIXED PRIORS FROM THE CONFIG OR UPDATE DURING PROCESSING
     bool fixedPriors  =     true;
   if (configTest.existsParam("fixedPriors"))
    fixedPriors  =     configTest.getParam("fixedPriors").toBool(); 
   
  //USE MAP ADAPTATION FOR CLIENT AND TESTS MODELS INSTEAD OF ML.
   bool map  =     false;
   if (configTest.existsParam("MAP"))
    map  =     configTest.getParam("MAP").toBool(); 
   
   unsigned long resetNb =    0;
   if (configTest.existsParam("ResetNbAdapt"))
    resetNb  =     configTest.getParam("ResetNbAdapt").toLong(); 
   
  //To stock Tnorm parameters, not need to be declared if Tnorm =    false!
  ObjectRefVector stockTnorm; 
   
   //To stock Znorm parameters, not need to be declared if Znorm =    false!
  ObjectRefVector stockZnorm;
   
  //NEW //To stock the LK for the client models with the impostor data !!!!!!!CAPACITE A CHANGER!!!!!!!!
  DoubleVector decision;//,decisionOnlyTests; 
   
  DoubleVector tmp;        //For NIST protocol 
  DoubleVector NbFramesSelected, NbFramesSelectedForTests;
  ObjectRefVector FeatServ;    //To stock the feature servers
  ObjectRefVector SegServ;    //To stock the seg clusters 
   
  if (tnorm)
  {                //load Tnorm parameters
      String testFile  =     (config.getParam("impScoreFile"));
      String TestNames  =     (config.getParam("testNames"));
      loadTnormParam(TestNames, testFile, stockTnorm, config);
  }
  
  String impCohortFile ;
  if(znorm) {
      impCohortFile  =     (configTest.getParam("impCohortFile"));
  }
  
  // label for selected frames - Only the frames associated with this label, in the label files, will be used
  bool fixedLabelSelectedFrame  =     true;
  String labelSelectedFrames;
  // the ID of each speaker is used as labelSelectedFrame ?
  if (config.existsParam("useIdForSelectedFrame"))    
    fixedLabelSelectedFrame  =     (config.getParam("useIdForSelectedFrame")  ==     "false");
  
  if (fixedLabelSelectedFrame)    // the label is decided by the command line and is unique for the run
    labelSelectedFrames  =     config.getParam("labelSelectedFrames");
  
  debug  =     config.getParam_debug();
  if (config.existsParam("verbose"))
    verbose  =     true;
  else
    verbose  =     false;
  
  // read the Id + filenames for each client
  XList inputClientList(inputClientListFileName, config);    
  XLine *linep;
  inputClientList.getLine(0);
  XList inputTestsList(inputTestsListFileName, configTest);
  XLine *lineTests;
  inputTestsList.getLine(0);
  MixtureServer ms(config);
  
  //CREATE CONFIG + MIXTURE SERVER FOR GMM OF SCORES
  Config configScores(config);    
  configScores.setParam("vectSize", "1");
  configScores.setParam("mixtureDistribCount", configTest.getParam("mixtureDistribCountForGMMScores"));
  MixtureServer msScores(configScores);
  StatServer ss(config, ms);
  StatServer ssScores(configScores, msScores);
  
  XLine testsToCompute;//,onlyTests;//NEW
  
  if (verbose)
    cout << "TrainTarget - Load world model [" << inputWorldFilename << "]"  << endl;
  
  MixtureGD & world  =     ms.loadMixtureGD(inputWorldFilename);
  
  // load TAR and NON GMM
  MixtureGD & non  =        msScores.createMixtureGD();
  MixtureGD & tar  =        msScores.createMixtureGD();
   
  if(!ScoresByTarget){
    non   =      msScores.loadMixtureGD(configTest.getParam("NONScoresModel"));
    tar   =      msScores.loadMixtureGD(configTest.getParam("TARScoresModel"));
  }

  
  /****************ADAPTED TNORM STUFF*********************/
 double shift =    0.0;//used for TNorm shift computation
 bool tnormAdapt  =     false;
 if (configTest.existsParam("TnormAdapt"))
    tnormAdapt  =     configTest.getParam("TnormAdapt").toBool(); 
 Matrix <double> tnormA;
 if(tnormAdapt && tnorm) {    //compute Adaptative TNORM
    //load Matrix of LLRs
    tnormA.load(configTest.getParam("TnormAMatrix"),configTest);
 }
 unsigned long session_nb =    1;
 /****************END ADAPTED TNORM STUFF*********************/  
 
 String idAux  =     (ms.createMixtureGD(world.getDistribCount())).getId();
 String idTmp  =     (ms.createMixtureGD(world.getDistribCount())).getId();

  // *********** Target loop ***************** 
       // linep gives the XLine with the Id of a given client and the list of files
  while ((linep  =     inputClientList.getLine()) !=     NULL)
  {   
      //   load needed for each client because of the mixture server reset 	    
      MixtureGD & world  =     ms.getMixtureGD(ms.getMixtureIndex(inputWorldFilename));    
      MixtureGD & tmpMixture  =     ms.getMixtureGD(ms.getMixtureIndex(idAux));
      MixtureGD & auxMixture  =     ms.getMixtureGD(ms.getMixtureIndex(idTmp));
      // Get the Client ID (id)	    
      String *id  =     linep->getElement();   
      // Get the list of feature file for the client (end of the line)	    
      XLine featureFileListp  =     linep->getElements();  
	    
      if (verbose)
          cout << "Train model [" << *id << "]" << endl;
      
      if (!fixedLabelSelectedFrame)
      {            // the ID is used as label for selecting the frame
          labelSelectedFrames  =     *id;
          if (debug)
              cout << *id << " is used for label selected frames" << endl;
      }
      //STUFF FOR TRAIN DATA
      // Reading the features (from several files)
      FeatureServer & fs          =     *new FeatureServer(config, featureFileListp);
      
      // Create the segment server for managing the segments/clusters
      SegServer & segmentsServer  =     *new SegServer();   
      
      // Create the lable server, for indexing the segments/clusters
      LabelServer labelServer;
      
      // Reading the segmentation files for each feature input file
      initializeClusters(featureFileListp, segmentsServer, labelServer, config);
      
      // Verify if the segments ending before the end of the feature files...
      verifyClusterFile(segmentsServer, fs, config); 
      
      MixtureGD & clientMixture  =     ms.duplicateMixture(world, DUPL_DISTRIB);
      //MixtureGD & mixtureOnlyTests  =     ms.duplicateMixture(world, DUPL_DISTRIB);//NEW
      MixtureGD & clientMixtureEM  =     ms.duplicateMixture(world, DUPL_DISTRIB);
      
      // Get the index of the cluster with in interest audio segments
      long codeSelectedFrame  =     labelServer.getLabelIndexByString(labelSelectedFrames);    
      if (codeSelectedFrame  ==     -1)
      {            // No data for this model !!!!!!!!!!!!!!
          cout << " WARNING - NO DATA FOR TRAINING [" << *id << "]";
      }
      else
      {
	  // Gives the cluster of the selected/used segments                                   
          SegCluster & selectedSegments  =     segmentsServer.getCluster(codeSelectedFrame);    
          countTests  =     0;    //reset the counter of tests segments (  train data included)
          if (!Fast)
          {
              //add train data to objectrefvector
              FeatServ.addObject(fs);
              SegServ.addObject(selectedSegments);
          }
	  
          adaptModel(config, ss, ms, fs, selectedSegments, world, clientMixture);
          if(!map){
              adaptModelEM(config, ss, ms, fs, selectedSegments, world, clientMixtureEM);
              ms.setMixtureId(clientMixtureEM, *id);
          }
          else ms.setMixtureId(clientMixture, *id);
		  
          testsToCompute.addElement(*id);
	  
          //compute the LLR with the client train data and store it to compute probability
          fullFileNameTrain  =     getFullFileName(*id, configTest);
	  
          //Fast LLR is not used for L(data train/ client model)
	  //NO more needed if proba is fixed to 1, tweak to init LKVector
          LLRTrain  =     computeLLR(configTest,ss, fs, world, clientMixture, selectedSegments,fullFileNameTrain); 
	  
	  //L(data train/ client model) always the first element of the doublevector
          decision.addValue(LLRTrain);    
          countTests++;
	  
	  //stock the number of frames selected (used for fuseModels)
          if(!map) NbFramesSelected.addValue(SegClusterFrame(selectedSegments)); 
		  
          MixtureGD & SAVclientMixture  =      ms.duplicateMixture(clientMixture, DUPL_DISTRIB);
	  
	  /*
           MEMORY OF 2 MODELS, INIT  = train model
          */
	  MixtureGD & clientMixture_2  =      ms.duplicateMixture(clientMixture, DUPL_DISTRIB);
	  MixtureGD & clientMixture_3  =      ms.duplicateMixture(clientMixture, DUPL_DISTRIB);
	  
	  
          //If NonByTarget is set the global Non model is replaced by the non model idclient.non.gmm (should be created before) 
          if(ScoresByTarget){
              String TARListFilename  =     configTest.getParam("TARListFilename");
              non  =     msScores.loadMixtureGD(((*id)+".non"));
              if (verbose) cout << "load "<<(*id)+".non"<<" NON file"<<endl;
              String idTar  =     selectNearestTarModel(TARListFilename, fullFileNameTrain,configTest, ss, fs, world,  selectedSegments,ms);
              if (verbose) cout <<"TAR SELECTED : "<< idTar<<endl;
              tar  =     msScores.loadMixtureGD(((idTar)+".tar"));
          }
	  
	  //COMPUTE AND STORE THE TNORM PARAMETERS FOR THIS CLIENT MODEL
	  if(znorm)
	  {   if(tnorm && znorm) ztnorm = true;
	      computeAndStoreZnormParam(ss, impCohortFile, (*id), clientMixture,stockZnorm, world, configTest, ztnorm, stockTnorm); 
	  }
	  
          /***********************************************************************************/
          /*******************************  ADAPTATION  PROCESS ******************************/
          /************!!!!!!INCLUDING THE CURRENT SEGMENT IS NOT PERFORMED !!!!!!!!!!!*******/
          /***********************************************************************************/
          shift =    0.0;
          // ***********  tests loop ***************** 
          while ((lineTests  =     inputTestsList.getLine()) !=     NULL)
          {
	      //file format : column client model name, column gender, column tests in one word (i.e : toto_A)
              String *idclient  =     lineTests->getElement();
              if (*idclient  ==     *id)
              {        
                  MixtureGD & testMixtureEM   =     ms.duplicateMixture(world, DUPL_DISTRIB);
		      
		  // Get the Tests ID (id)
                  String & idTest             =     lineTests->getElement(2, becomeCurrent);
		      
		  // Get the list of feature file for the tests (end of the line)
                  XLine featureFileListTests  =     lineTests->getElements();    
                  
		  if (verbose)
                      cout << "Train model for test [" << idTest << "]" << endl;
		  
		  // the ID is used as label for selecting the frame
                  if (!fixedLabelSelectedFrame)
                  {        
                      labelSelectedFrames  =     idTest;
                      if (debug)
                          cout << idTest << " is used for label selected frames" << endl;
                  }
		  
                  //TEST DATA STUFF
                  FeatureServer & fsTests  =       *new FeatureServer(configTest, featureFileListTests);
            
                  // Create the segment server for managing the segments/clusters
                  SegServer & segmentsServerTests  =     *new SegServer();    
            
                  // Create the lable server, for indexing the segments/clusters
                  LabelServer labelServerTests; 
		  
                  // Reading the segmentation files for each feature input file
                  initializeClusters(featureFileListTests, segmentsServerTests, labelServerTests, configTest);
		  
                  // Verify if the segments ending before the end of the feature files...
		  verifyClusterFile(segmentsServerTests, fsTests, configTest);  
		  
                  // Get the index of the cluster with in interest audio segments
		  long codeSelectedFrame  =     labelServerTests.getLabelIndexByString(labelSelectedFrames);    
                  if (codeSelectedFrame  ==     -1)
                  {        // No data for this model !!!!!!!!!!!!!!
                      cout << " WARNING - NO DATA FOR TRAINING [" << idTest << "]";
                  }
                  else
                  { 
                      SegCluster & allFramesTest          =     segmentsServerTests.getCluster(codeSelectedFrame);
                      SegCluster & selectedSegmentsTests  =     segmentsServerTests.createCluster(1,"","");
                
                      if(cross)
                      {
                          /* int indx =    -1;
                          if((indx  =     findClusterInRefvector(SegServ,idTest))! =    -1) copyCluster((static_cast <SegCluster & >(SegServ.getObject(indx))),selectedSegmentsTests);
			  else{*/
                          crossValid(configTest, ss,  ms, fsTests, allFramesTest,world,testMixtureEM,selectedSegmentsTests,idTest);
                          ms.setMixtureId(testMixtureEM, idTest);
              
			  /* selectedSegmentsTests.setSourceName(idTest); //set name of the segcluster for findClusterInRefVector function
			  //cout <<"SOURCE NAME : " <<selectedSegmentsTests.sourceName()<<endl;
			  SegServ.addObject(selectedSegmentsTests); //stock selected cluster in RefVector to avoid recalculation (take care it could be memory consumming)
			  }*/
                      }
                      else {//copy allFrameTests in selectedSegmensTests
                          copyCluster(allFramesTest,selectedSegmentsTests);
			  segmentsServerTests.remove(allFramesTest);
		      }
                     
		      
		      if (Fast)
                      {
                           
                           //onlyTests.addElement(idTest);//NEW
                          int index  =     ms.getMixtureIndex(idTest); //IF STOCKED (ALREADY COMPUTED) LOAD MODEL
                          if (index  ==     -1)
                          { 
                              String fullMixtureName =    getFullMixtureName(idTest,config);
                              
			      //IF TEST MIXTURE NOT STORED ON THE HARD DRIVE
			      if (FileExists(fullMixtureName)  ==     false)
			      {
                                  //create the EM model for the test 
				  if(!map) 
				      adaptModelEM(config, ss, ms, fsTests, selectedSegmentsTests, world, testMixtureEM);    
				  
                                  else adaptModel(config, ss, ms, fsTests, selectedSegmentsTests, world, testMixtureEM);
					  
                                  ms.setMixtureId(testMixtureEM, idTest);
				  
				  //SAVE test mixture 
                                  testMixtureEM.save(idTest,config);        
                              }
                              else testMixtureEM =     ms.loadMixtureGD(idTest);
                         }    // Set the client model Id
                         else
                             testMixtureEM  =     ms.getMixtureGD(index);
                      }
                      //Compute LLR and write it on a res file
                      LLRRatio  =     0.0;
                      fullFileName  =     getFullFileName(idTest, configTest);
		      
                      //TAKE CARE : by using it you should have 2 GMM of scores and a TNORM RES file according to the system used for computing the RES file 
                      if(!tnormAdapt){
			  
                          //TEST SHIFT TNORM			      
                          if (fromResFile)        
                          {
                              LLRRatio  =      searchLLRFromResFile(*id, idTest,  inputResFilename, configTest);
                          }
                          else
                          {
                              //Choice wether to use FastLLR or LLR (FastLLR is used when the top component are already computed for this Feature Server)
                          
			      //first time seeing this test
			      if (FileExists(fullFileName)  ==     false  || Fast  ==     false || cross)
                              {   
                                  //For proba, take the decision on the original target model
			          LLRRatio  =     computeLLR(configTest, ss, fsTests, world, SAVclientMixture, selectedSegmentsTests, fullFileName);    
                              }
			      //For proba, take the decision on the original target model
                              else
                                  LLRRatio  =     computeFastLLR(ss, fsTests, world, SAVclientMixture, selectedSegmentsTests, fullFileName, configTest);    
                          }
                
                          if (tnorm) normalizeScore(idTest, LLRRatio, stockTnorm,shift);
			  if (znorm) normalizeScore(*id, LLRRatio, stockZnorm,shift);	  
                     }
            
                     
                     
		     LLRRatio2  =     0.0;
		     if (FileExists(fullFileName)  ==     false  || Fast  ==     false || cross)
		     {    
			 //first time seeing this test
			 LLRRatio2  =    computeLLR(configTest, ss, fsTests, world,clientMixture, selectedSegmentsTests,fullFileName);
		     }
	             else   LLRRatio2  =     computeFastLLR(ss, fsTests, world, clientMixture, selectedSegmentsTests, fullFileName, configTest);    //for LLR output

		     if(tnormAdapt && tnorm)
		     {
                
		         shift =    findNearestLLRInMatrix(tnormA, session_nb, LLRRatio2);
			 cout <<"SHIFT  =     " <<shift<<endl;
			 LLRRatio2 -=    shift;
		     }
                
		     //if (tnorm)   TnormalizeScore(idTest, LLRRatio2, stockTnorm,shift);
			 
		     if (tnorm) normalizeScore(idTest, LLRRatio2, stockTnorm,shift);
		     if (znorm) normalizeScore(*id, LLRRatio2, stockZnorm,shift);
				
		    /*
		    TEST TRIAL ON MEMORY 
		    */
		    double L2 = computeFastLLR(ss, fsTests, world, clientMixture_2, selectedSegmentsTests, fullFileName, configTest);    
                    double L3 = computeFastLLR(ss, fsTests, world, clientMixture_3, selectedSegmentsTests, fullFileName, configTest);    
		     
		    /*
		     TRAIN DATA ON MEMORY, take care TNORM not available for train data
                    */	
		    double Ltrain   = computeFastLLR(ss, fs, world, clientMixture, selectedSegments, fullFileNameTrain, configTest);    
                    double Ltrain_2 = computeFastLLR(ss, fs, world, clientMixture_2, selectedSegments, fullFileNameTrain, configTest);    
                    double Ltrain_3 = computeFastLLR(ss, fs, world, clientMixture_3, selectedSegments, fullFileNameTrain, configTest);    
		     
		     
		    //DISPLAY STATS
		    cout <<" Model : "<< *id << " , train data : "
		         <<" LLR test adapt -1 : " << Ltrain
		         <<" LLR test adapt -2 : " << Ltrain_2
		         <<" LLR test adapt -3 : " << Ltrain_3 <<endl;
		 
		     
		     
		    if (tnorm) normalizeScore(idTest, L2, stockTnorm,shift);
		    if (tnorm) normalizeScore(idTest, L3, stockTnorm,shift);
			    
		    //DISPLAY STATS
		    cout <<" Model : "<< *id << " , current test : " << idTest
		         <<" LLR test adapt -1 : " << LLRRatio2
		         <<" LLR test adapt -2 : " << L2
		         <<" LLR test adapt -3 : " << L3
		         <<" LLR test train model : " << LLRRatio <<endl;
		
		    
		    
		    //TEST FOR ADAPTATION ERRORS, N-1 equal current client model
		   /* if(LLRRatio2 > (L2 + 0.24)) //ERROR OCCURS AT N-1
		    {   
			//IF ERRORS OUTPUT LLR ON MODEL N-1 AND SET TRIAL N-1 WEIGHT TO 0
			cout << "Found adaptation error, Model : "<<*id<<" Test : " << testsToCompute.getElement(countTests-1, false)<<endl; 
			LLRRatio2 = L2;
			clientMixture =  clientMixture_2;
			decision[countTests-1] = configTest.getParam("meanImp").toDouble();
		    }
		    
		    if(LLRRatio2 > (L3 + 0.24) && (L2-L3 > 0.24)) //ERROR OCCURS AT N-2
		    {   
			//IF ERRORS OUTPUT LLR ON MODEL N-2 AND SET TRIAL N-2 WEIGHT TO 0
			cout << "Found adaptation error, Model : "<<*id<<" Test : " << testsToCompute.getElement(countTests-2, false)<<endl; 
			LLRRatio2 = L3;
			clientMixture =  clientMixture_3;
			clientMixture_2 =  clientMixture_3;
			decision[countTests-2] = configTest.getParam("meanImp").toDouble();
		    }
		    */
		    
		    // Output the result
		    outLLR << "M " << *id << " 0 " << idTest << " " << LLRRatio2 << endl;    
               
		    
                    /**********************DEBUG******************************************/
                    //if(LLRRatio < 0) LLRRatio  =     configTest.getParam("meanImp").toDouble();
                    /**********************END DEBUG******************************************/
            
                    if(tnormAdapt)
                        LLRRatio =    LLRRatio2;
		    
		    //if Oracle == true , DO NOT SET regress, wmap or wmapgmm to true!!!
                    if (oracle)
                        Oracle(*id, idTest, LLRRatio, configTest, tar, non, ssScores);
		    
		    //IF TRIAL SELECTED FOR ADAPTATION (DEFAULT : ALWAYS ON)
		    if(selected) //VALID FOR THE CURRENT, SHOULD BE UPDATED FOR PAST TRIALS
		    {
                        //Take the WMAP decision on adapted client model (LLRRatio2) or one session record (LLRRatio)		    
                        decision.addValue(LLRRatio);
			
			 //add tests data to objectrefvector
                        if (!Fast)
                        {
                            FeatServ.addObject(fsTests);
                            SegServ.addObject(selectedSegmentsTests);
		        }    
			countTests++; 
			
			//stock the number of frames selected (used for fuseModels)
		        if(!map) NbFramesSelected.addValue(SegClusterFrame(selectedSegmentsTests));    
                        testsToCompute.addElement(idTest);
			
			selected = true;
		    }
              
                   /**********************************NIST protocol*******************************/
	           //SAVE LLRs because WMAP will modify them 	   
	           tmp  =     decision;    
	           //compute FS weights
	           if (regress)
		       expandLLR(decision, configTest);
	           else if (wmap)
		       WMAP(decision, configTest);
	           else if (wmapGmm){
		       if(fixedPriors)
		           WMAPGMMFixedPriors(decision, configTest, tar, non, ssScores);
		       else                        
		           WMAPGMM(decision, configTest, tar, non, ssScores);
		    }
		
		    /**************************for debug***************/
		    /*   cout <<"Reset weights > 0.5 (for test)"<<endl;
		    for(unsigned long e =    1; e<decision.size();e++){
		    //if(decision[e]>0.5 || decision[e]<0.1)
		    if (decision[e]<0.8)
		    decision[e] =    0;
	     
		    }*/
		    /***********************end***********************/
    
		    cout <<"Set weight of train data to 1"<<endl;
		    decision[0]  =     1;    //SET a weight of 1 for the train data
      
		    //create the new client model 
		    if (Fast)
		    {
		        if(assess){
			    assessAdaptation(ss,clientMixture ,selectedSegments, configTest, fs, world, tmp, fullFileNameTrain, countTests, LLRTrain);
		        }
	
		        //STUFF FOR DELETING PREVIOUS ADAPTATION, BACK TO TRAIN MODEL
		        if(resetNb && resetNb < (unsigned long)countTests)
		        {   
		            //new model  =     init client
		            //RESET
		            testsToCompute.reset();
		            NbFramesSelected.clear();
		            decision.clear();
		            //RELOAD WITH CLIENT ONLY

		            if(!map) NbFramesSelected.addValue(SegClusterFrame(selectedSegments));
		            testsToCompute.addElement(*id);
		            decision.addValue(1);    
		        }
		   
		        //UPDATE MEMORY MODELS
		        clientMixture_3 = clientMixture_2; 
		        clientMixture_2 = clientMixture;
		        //EXIT the oldest clientMixture_3
		   
		   
		        //Create adapted clientmodel
		        if(!map)
		            computeMAPmodelFromEMones(configTest, ss, ms, NbFramesSelected, world, clientMixture, auxMixture, tmpMixture, decision, testsToCompute);
		        else fuseMAP(config, ss,  ms, world, clientMixture, auxMixture,tmpMixture, decision, testsToCompute);
	       
		        if(saveAfterIt){
		            /* SAVE THE TARGET MODEL AFTER EACH ITERATION        */                
		            cout << "Save client model [" << *id << "]" <<" after IT : "<< countTests-1 << endl;
		            char  toto[256];
		            sprintf(toto,"%d",(countTests-1));         
		            clientMixture.save((*id)+"."+toto, config);    // Save the client model
		        } 

			    
                       }
		    else 
		    {
			clientMixture  =     world;    //reset clientMixture to MAP from world
			adaptModel(config, ss, ms, FeatServ, SegServ, world, clientMixture, decision);
		    }
              
		    //TEST SHIFT TNORM
		    if(tnormAdapt)
		    {    
			double toto =    0.0;
			for (unsigned long e =    0; e < NbFramesSelected.size(); e++)
                               toto  +=    NbFramesSelected[e]*decision[e];
			cout <<"Total Frames used for learning the adapted model  =     " << toto<<endl;
                
			/*if(0.0 <toto && toto <9500.0)          shift  =     0.0;
			else if (9500.0 <toto && toto <17000.0)  shift  =     0.10869;    
			else if (17000.0 <toto && toto <25500.0) shift  =     0.15848;    
			else if (25500.0 <toto && toto <34000.0) shift  =     0.18503;    
			else if (34000.0 <toto && toto <42400.0) shift  =     0.20503;    
			else if (42400.0 <toto && toto <50000.0) shift  =     0.22061;    
			else if ( 50000.0<toto && toto <59400.0) shift  =     0.23417;    
			else if ( 59400.0<toto && toto <67900.0) shift  =     0.24945;    
			else if ( 67900.0<toto && toto <76400.0) shift  =     0.26698;    
			else if (76400.0 <toto && toto <85000.0) shift  =     0.28624;    
			else if (toto>85000.0)                   shift  =     0.28624;    */
                
			if (9500.0 < toto && toto < 17000.0)        session_nb =    2;    
			else if ( 17000.0 < toto && toto < 25500.0) session_nb =    3;    
			else if ( 25500.0 < toto && toto < 34000.0) session_nb =    4;    
			else if ( 34000.0 < toto && toto < 42400.0) session_nb =    5;    
			else if ( 42400.0 < toto && toto < 50000.0) session_nb =    6;    
			else if ( 50000.0 < toto && toto < 59400.0) session_nb =    7;    
			else if ( 59400.0 < toto && toto < 67900.0) session_nb =    8;    
			else if ( 67900.0 < toto && toto < 76400.0) session_nb =    9;    
			else if ( 76400.0 < toto && toto < 85000.0) session_nb =    10;    
			else if (toto > 85000.0)                    session_nb =    10;
                
                
		    }
                
		    decision  =     tmp; //RESTORE LLRs
		    tmp.clear();
		      
	      }
              if (Fast)
              {        //delete object created with new, not nedded anymore
                  delete & segmentsServerTests;
                  delete & fsTests;
              }
	      //UNTIL NOW LOAD EACH TEST MIXTURE
              if (ms.getMixtureCount() > (unsigned long)configTest.getParam("MaxMixturesCount").toLong())
	      { 
	      
	          if(verbose && verboseLevel>1)
	              cout <<"Exceeding capacity of the mixtureServer ("<<(unsigned long)configTest.getParam("MaxMixturesCount").toLong()<<")"<<endl;
                
	          ms.deleteMixture(testMixtureEM); //delete current test mixture
	          ms.deleteUnusedDistribs();
                
              }
            
            }
        
        }            //end tests loop

        if (verbose)
            cout << "Save client model [" << *id << "]" << endl;
        clientMixture.save(*id, config);    // Save the client model
        inputTestsList.rewind();    //rewind the xlist, needed when changing target
    
	//delete all client Mixtures 
	ms.deleteMixture(clientMixture);
	ms.deleteMixture(clientMixtureEM);
	ms.deleteMixture(SAVclientMixture);
	ms.deleteUnusedDistribs();

	NbFramesSelected.clear();
	decision.clear();
	
	/***************************NEW STUFF*********************  */    
        /* NbFramesSelectedForTests.clear();
	decisionOnlyTests.clear(); 
	ms.deleteMixture(mixtureOnlyTests);
        onlyTests.reset();
	previousLLR =    -10000;*/
        /*************************** END NEW STUFF*********************  */
	
	if (!Fast)
        {
           FeatServ.deleteAllObjects();
           SegServ.deleteAllObjects();
        }
	
	testsToCompute.reset();
	
	if (Fast)
        {        //delete object created with new, not nedded anymore
           delete & fs;
           delete & segmentsServer;
        }
        
        //reset mixtureServer to free memory (reload needed for several mixtures)
        if (ms.getMixtureCount() > (unsigned long)configTest.getParam("MaxMixturesCount").toLong())
        {
            if (verbose)
                cout << "Reset MixtureServer and reload  of the used mixtures" << endl;
	    
	    ms.reset();
            MixtureGD & world  =     ms.loadMixtureGD(inputWorldFilename);    //load needed for each client because of the mixture server reset 
            
	    //temporate mixtures used in computeMAPmodelFromEMones, create here to accelerate
            idAux  =     (ms.createMixtureGD(world.getDistribCount())).getId();
            idTmp  =     (ms.createMixtureGD(world.getDistribCount())).getId();
        }
    }
}                // end of the the target loop 

//When end delete all mixture in the server
  ms.deleteMixtures(0, ms.getMixtureCount());
  ms.deleteUnusedDistribs();

  if (tnorm)
      stockTnorm.deleteAllObjects();
  if(znorm)
      stockZnorm.deleteAllObjects();
  outLLR.close();

  return 0;
}





#endif //!defined(ALIZE_TrainTarget_cpp)

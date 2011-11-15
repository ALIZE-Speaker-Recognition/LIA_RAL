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

#if !defined(ALIZE_TrainTarget_cpp)
#define ALIZE_TrainTarget_cpp

#include <iostream>
#include <fstream>  // pour outFile
#include <cstdio>   // pour printf()
#include <cassert> // pour le debug pratique
#include <cmath>
#include <liatools.h>
#include "TrainTarget.h"

using namespace alize;
using namespace std;

// Information on the quantity of data available by client
// could output a list with the selected files for a defined quantity of data
// Uses the same input file than traintarget and output a new list
int InfoTarget(Config& config)
{
  String inputClientListFileName = config.getParam("targetIdList");
  bool fixedLabelSelectedFrame;
  String labelSelectedFrames;
  if (config.existsParam("useIdForSelectedFrame"))      // the ID of each speaker is used as labelSelectedFrame
    fixedLabelSelectedFrame=false;
  else{                                                // the label is decided by the command line and is unique for the run
    labelSelectedFrames=config.getParam("labelSelectedFrames");
    if (verbose) cout << "Computing on" << labelSelectedFrames << " label" << endl;
    fixedLabelSelectedFrame=true;
  }
  unsigned long maxFrame=config.getParam("maxFrame").toLong();
  String outputFilename=config.getParam("outputFilename");
  
  
  ofstream outputFile(outputFilename.c_str(),ios::out| ios::trunc);
  try{
    XList inputClientList(inputClientListFileName,config);          // read the Id + filenames for each client
    XLine * linep;
    if (verbose) cout << "InfoTarget" << endl;
    // *********** Target loop ***************** 
    while ((linep=inputClientList.getLine()) != NULL){             // linep gives the XLine with the Id of a given client and the list of files
      String *id=linep->getElement();                              // Get the Client ID (id)
      outputFile<<*id;
      String currentFile="";
      XLine featureFileListp=linep->getElements();	           // Get the list of feature file for the client (end of the line)
      if (verbose) cout << "Info model ["<<*id<<"]"<<endl;   
      if (!fixedLabelSelectedFrame){                                // the ID is used as label for selecting the frame
	labelSelectedFrames=*id;
	if (debug) cout <<*id<<" is used for label selected frames"<<endl;
      }
      // label files reading - It creates, for each file and each label, a cluster of segments - will be integrated witth the featre s - asap
      SegServer segmentsServer;                                    // Reading the segmentation files for each feature input file
      LabelServer labelServer;                                     
      initializeClusters(featureFileListp,segmentsServer,labelServer,config);           // Reading the segmentation files for each feature input file
      unsigned long codeSelectedFrame=labelServer.getLabelIndexByString(labelSelectedFrames);            // Get the index of the cluster with in interest audio segments
      SegCluster& selectedSegments=segmentsServer.getCluster(codeSelectedFrame); // Gives the cluster of the selected/used segments   
      Seg *seg;                                                                  // Will give the current segment
      unsigned long frameCount=0;
      selectedSegments.rewind();                                                 // at the begin of the selected segments list
      while(((seg=selectedSegments.getSeg())!=NULL) && (frameCount<maxFrame)){   // For each of the selected segments until the amount of data is get
	frameCount+=seg->length();
	cout << seg->sourceName()<<" "<<seg->begin()<<" "<<seg->length()<<" Total time="<<frameCount<<endl;
	if (seg->sourceName()!=currentFile){                                     
	  outputFile<<" "<<seg->sourceName();
	  currentFile=seg->sourceName();
	}
      }                                                                          // end of the initial Train Iteration loop
      outputFile<<endl;
      if (verbose) cout << "Save info client ["<<*id<<"]" << endl;
    }                                                                            // end of the the target loop 
  } // fin try

  catch (Exception& e)
    { 
      cout << e.toString().c_str() << endl;
    }
  return 0;
}
 

// Training of client Speakers
// Input: Xlist	Format: ID_Client Seg1 Seg2 ..
// Output: ALIZE_MixtureServer (binaire) + GMM / Client (binary)
 int TrainTarget(Config& config)
{
  String inputClientListFileName = config.getParam("targetIdList");
  String inputWorldFilename = config.getParam("inputWorldFilename");
  String outputSERVERFilename = "";
  if (config.existsParam("mixtureServer")) outputSERVERFilename =config.getParam("mixtureServer");
  bool initByClient=false;                                              // In this case, the init model is read from the file
  if (config.existsParam("initByClient")) initByClient=config.getParam("initByClient").toBool();
  bool saveEmptyModel=false;
  if (config.existsParam("saveEmptyModel")) saveEmptyModel=config.getParam("saveEmptyModel").toBool();
  // label for selected frames - Only the frames associated with this label, in the label files, will be used
  bool fixedLabelSelectedFrame=true;
  String labelSelectedFrames;
  if (config.existsParam("useIdForSelectedFrame"))    // the ID of each speaker is used as labelSelectedFrame ?
    fixedLabelSelectedFrame=(config.getParam("useIdForSelectedFrame").toBool()==false);  
  if (fixedLabelSelectedFrame)                        // the label is decided by the command line and is unique for the run
    labelSelectedFrames=config.getParam("labelSelectedFrames");
  bool modelData=false;
  if (config.existsParam("useModelData")) modelData=config.getParam("useModelData").toBool();
  String initModelS=inputWorldFilename;
  if (modelData) if (config.existsParam("initModel")) initModelS=config.getParam("initModel"); // Use a specific model for Em init
  bool outputAdaptParam=false;
   if (config.existsParam("superVector")) outputAdaptParam=true;
  bool NAP=false;
  Matrix <double> ChannelMatrix;
  if (config.existsParam("NAP")) {
     if (verbose) cout<< "Removing channel effect with NAP from " << config.getParam("NAP") << " of size: ["; 
    NAP=true; // enable NAP
    ChannelMatrix.load(config.getParam("NAP"),config); //get Channel Matrix from args and load in a Matrix object
     if (verbose) cout << ChannelMatrix.rows() << "," <<ChannelMatrix.cols() << "]" << endl;
    }
    

  bool saveCompleteServer=false;
 
  try{
    XList inputClientList(inputClientListFileName,config);          // read the Id + filenames for each client
    XLine * linep;
    inputClientList.getLine(0);
    MixtureServer ms(config);
    StatServer ss(config, ms);
    if (verbose) cout << "TrainTarget - Load world model [" << inputWorldFilename<<"]"<<endl;
    MixtureGD& world = ms.loadMixtureGD(inputWorldFilename);
    MixtureGD& initModel =ms.loadMixtureGD(initModelS);
   
    if (verbose) cout <<"Use["<<initModelS<<"] for initializing EM"<<endl;
    
    // *********** Target loop ***************** 
    while ((linep=inputClientList.getLine()) != NULL){             // linep gives the XLine with the Id of a given client and the list of files
      String *id=linep->getElement();                              // Get the Client ID (id)
      XLine featureFileListp=linep->getElements();	           // Get the list of feature file for the client (end of the line)
      if (verbose) cout << "Train model ["<<*id<<"]"<<endl;   
      if (!fixedLabelSelectedFrame){                                // the ID is used as label for selecting the frame
	labelSelectedFrames=*id;
	if (verbose) cout <<*id<<" is used for label selected frames"<<endl;
      }
      FeatureServer fs(config,featureFileListp);                                            // Reading the features (from several files)
      SegServer segmentsServer;                                                             // Create the segment server for managing the segments/clusters
      LabelServer labelServer;                                                              // Create the lable server, for indexing the segments/clusters
      initializeClusters(featureFileListp,segmentsServer,labelServer,config);               // Reading the segmentation files for each feature input file
      verifyClusterFile(segmentsServer,fs,config);                                          // Verify if the segments ending before the end of the feature files...
      MixtureGD & adaptedMixture = ms.duplicateMixture(world,DUPL_DISTRIB);                 // Creating final as a copy of the world model
      MixtureGD & clientMixture= ms.duplicateMixture(world,DUPL_DISTRIB);
      if (initByClient){                                                                   // During trainig data statistic estimation by EM,
	clientMixture= ms.loadMixtureGD(*id);                                               // the client model is used for initalization
	adaptedMixture=clientMixture;
      }
      long codeSelectedFrame=labelServer.getLabelIndexByString(labelSelectedFrames);        // Get the index of the cluster with in interest audio segments
      if (codeSelectedFrame==-1){                                                           // No data for this model !!!!!!!!!!!!!!
	cout << " WARNING - NO DATA FOR TRAINING ["<<*id<<"]";
	if (saveEmptyModel){
	  cout <<" World model is returned"<<endl;                                    // In this case, the client model is the world model
	  if (verbose) cout << "Save client model ["<<*id<<"]" << endl;
	  adaptedMixture.save(*id, config);                                           // Save the client model
	}
      }
      else{
	SegCluster& selectedSegments=segmentsServer.getCluster(codeSelectedFrame); // Gives the cluster of the selected/used segments                                   
	if (!initByClient) ms.setMixtureId(clientMixture,*id);                                        // Set the client model Id
	if (modelData) modelBasedadaptModel(config,ss,ms,fs,selectedSegments,world,clientMixture,initModel);          // EM algo with MAP criterion
	else adaptModel(config,ss,ms,fs,selectedSegments,world,clientMixture);          // EM algo with MAP criterion
        if (NAP) {
          if (verbose) cout << "NAP on SVs" << endl;
          computeNap(clientMixture,ChannelMatrix);
        }
        if (outputAdaptParam) {
            RealVector<double> v;
            getSuperVector(v,world,clientMixture,config);   
           String out=config.getParam("vectorFilesPath")+*id+config.getParam("vectorFilesExtension");        
            Matrix <double> vv=(Matrix<double>)v;
            vv.save(out,config);     
          }
        if (!outputAdaptParam) {
            if (verbose) cout << "Save client model ["<<*id<<"]" << endl;
            clientMixture.save(*id, config);                                           // Save the client model
        }
	if (!saveCompleteServer){
	  long tid=ms.getMixtureIndex(*id);      // TO BE SUPPRESSED BY
	  ms.deleteMixtures(tid,tid);            // ADDING a delete on a mixture pointor
	  ms.deleteUnusedDistribs();
	  }
      }
    }                                                                              // end of the the target loop 
    
    // Save the complete mixture server
    // TODO
  } // fin try
   


  catch (Exception& e)
    { 
      cout << e.toString().c_str() << endl;
    }
  return 0;
}


 
// Training of client Speakers
// The same than TrainTarget but train simultaneoulsy 1 model for each cluster (set of segments with the same label)
// found in the input files labels.
// One option in order to save the n models as a modification of the world model - save disk space
 int TrainTargetByLabel(Config& config)
{
  String inputClientListFileName = config.getParam("targetIdList");
  String inputWorldFilename = config.getParam("inputWorldFilename");
  String outputSERVERFilename = config.getParam("mixtureServer");
  // label for selected frames - Only the frames associated with this label, in the label files, will be used
  //bool fixedLabelSelectedFrame;
  bool initByClient=false;
  bool aprioriWorld=true;
  if (config.existsParam("initByClient")) initByClient=true;
  if (config.existsParam("aprioriClient")){
    aprioriWorld=false;
    initByClient=true;
  }
  bool saveCompleteServer=false;
  bool outputAdaptParam=false;
  if (config.existsParam("outputAdaptParam")) outputAdaptParam=config.getParam("outputAdaptParam").toBool();
  
  try{
    XList inputClientList(inputClientListFileName,config);          // read the Id + filenames for each client
    XLine *linep;
    inputClientList.getLine(0);
    MixtureServer ms(config);
    StatServer ss(config, ms);
    if (verbose) cout << "TrainTarget - by label opption - Load world model [" << inputWorldFilename<<"]"<<endl;
    MixtureGD& world = ms.loadMixtureGD(inputWorldFilename);
    // *********** Target loop ***************** 
    while ((linep=inputClientList.getLine()) != NULL){             // linep gives the XLine with the Id of a given client and the list of files
      String clientId=(*linep->getElement());                      // Get the Client ID (clientId)
      XLine featureFileListp=linep->getElements();	           // Get the list of feature file for the client (end of the line)
      FeatureServer fs(config,featureFileListp);                   // Reading the features (from several files)
      if (verbose) cout << "Train label models for client ["<<clientId<<"]"<<endl;   
      MixtureGD &clientGModel=ms.createMixtureGD();
      if (initByClient) {
          if (verbose) cout << "Load client model [" << clientId <<"]"<<endl;
          clientGModel = ms.loadMixtureGD(clientId); //not necessary to load client model
      }
      SegServer segmentsServer;                                                             // Create the segment server for managing the segments/clusters
      LabelServer labelServer;                                                              // Create the lable server, for indexing the segments/clusters
      initializeClusters(featureFileListp,segmentsServer,labelServer,config);               // Reading the segmentation files for each feature input file
      verifyClusterFile(segmentsServer,fs,config);                                          // Verify if the segments ending before the end of the feature files...
      for (unsigned long codeSelectedFrame=0;codeSelectedFrame<segmentsServer.getClusterCount();codeSelectedFrame++){ // For each cluster
	String clientIdByLabel=clientId+"_"+labelServer.getLabel(codeSelectedFrame).getString(); // Build the model name for the client and the label
	if (verbose) cout << "Train labeldependent model ["<<clientIdByLabel<<"]"<<endl;   
	SegCluster& selectedSegments=segmentsServer.getCluster(codeSelectedFrame);          // Gives the cluster of the selected/used segments
	MixtureGD & clientMixture = ms.duplicateMixture(world,DUPL_DISTRIB);       // Creating clientMixture as a copy of the world model
	ms.setMixtureId(clientMixture,clientIdByLabel);                                     // Set the client model Id
	if (initByClient)                                                                   // During trainig data statistic estimation by EM,
	  clientMixture=clientGModel;                                                       // the global client model is used for initalization
	if (aprioriWorld)                                                                   // EM algo with MAP criterion
	  adaptModel(config,ss,ms,fs,selectedSegments,world,clientMixture);          // A priori info is the world model  
	else adaptModel(config,ss,ms,fs,selectedSegments,clientGModel,clientMixture);// A priori info is the client model-by default initByClient is also set
        if (!outputAdaptParam) {
            if (verbose) cout << "Save client model ["<<clientIdByLabel<<"]" << endl;
            clientMixture.save(clientIdByLabel, config);                                           // Save the client model
        }
	if (!saveCompleteServer){
	  long tid=ms.getMixtureIndex(clientIdByLabel);      // TO BE SUPPRESSED BY
	  ms.deleteMixtures(tid,tid);            // ADDING a delete on a mixture pointor
	  ms.deleteUnusedDistribs();
	  }
      }      
      if (!saveCompleteServer){
	long tid=ms.getMixtureIndex(clientId);      // TO BE SUPPRESSED BY
	ms.deleteMixtures(tid,tid);                 // ADDING a delete on a mixture pointor
	ms.deleteUnusedDistribs();
      }                                                                   // end of the the label loop fr a speaker
    } // end of the the target loop 
    
    // Save the complete mixture server
    // TODO
  } // fin try
  

  
  catch (Exception& e)
    { 
      cout << e.toString().c_str() << endl;
    }
  return 0;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
int TrainTargetFA(Config& config)
{
  String inputClientListFileName = config.getParam("targetIdList");
  String inputWorldFilename = config.getParam("inputWorldFilename");
  String outputSERVERFilename = "";
  if (config.existsParam("mixtureServer")) outputSERVERFilename =config.getParam("mixtureServer");
  bool initByClient=false;                                              // In this case, the init model is read from the file
  if (config.existsParam("initByClient")) initByClient=config.getParam("initByClient").toBool();
  bool saveEmptyModel=false;
  if (config.existsParam("saveEmptyModel")) saveEmptyModel=config.getParam("saveEmptyModel").toBool();
  // label for selected frames - Only the frames associated with this label, in the label files, will be used
  bool fixedLabelSelectedFrame=true;
  String labelSelectedFrames;
  if (config.existsParam("useIdForSelectedFrame"))    // the ID of each speaker is used as labelSelectedFrame ?
    fixedLabelSelectedFrame=(config.getParam("useIdForSelectedFrame").toBool()==false);  
  if (fixedLabelSelectedFrame)                        // the label is decided by the command line and is unique for the run
    labelSelectedFrames=config.getParam("labelSelectedFrames");
  bool modelData=false;
  if (config.existsParam("useModelData")) modelData=config.getParam("useModelData").toBool();
  String initModelS=inputWorldFilename;
  if (modelData) if (config.existsParam("initModel")) initModelS=config.getParam("initModel"); // Use a specific model for Em init
  bool outputAdaptParam=false;
  if (config.existsParam("superVectors")) outputAdaptParam=true;
  Matrix <double> ChannelMatrix;
  if (verbose) cout<< "EigenMAP and Eigenchannel with [" << config.getParam("initChannelMatrix") << "] of size: ["; 
  ChannelMatrix.load(config.getParam("initChannelMatrix"),config); //get Channel Matrix from args and load in a Matrix object
  if (verbose) cout << ChannelMatrix.rows() << "," <<ChannelMatrix.cols() << "]" << endl;
  bool varAdapt=false;
  if (config.existsParam("FAVarAdapt")) varAdapt=true;
  bool saveCompleteServer=false;
 
  try{
    XList inputClientList(inputClientListFileName,config);          // read the Id + filenames for each client
    XLine * linep;
    inputClientList.getLine(0);
    MixtureServer ms(config);
    StatServer ss(config, ms);
    if (verbose) cout << "(TrainTarget) Factor Analysis - Load world model [" << inputWorldFilename<<"]"<<endl;
    MixtureGD& world = ms.loadMixtureGD(inputWorldFilename);      
    if (verbose) cout <<"(TrainTarget) Use["<<initModelS<<"] for initializing EM"<<endl;
    
    // *********** Target loop ***************** 
    while ((linep=inputClientList.getLine()) != NULL){             // linep gives the XLine with the Id of a given client and the list of files

      String *id=linep->getElement();                              // Get the Client ID (id)
      XLine featureFileListp=linep->getElements();	           // Get the list of feature file for the client (end of the line)
      if (verbose) cout << "(TrainTarget) Train model ["<<*id<<"]"<<endl;   
      FeatureServer fs(config,featureFileListp);                                            // Reading the features (from several files)
      SegServer segmentsServer;                                                             // Create the segment server for managing the segments/clusters
      LabelServer labelServer;                                                              // Create the lable server, for indexing the segments/clusters
      initializeClusters(featureFileListp,segmentsServer,labelServer,config);               // Reading the segmentation files for each feature input file
      verifyClusterFile(segmentsServer,fs,config);                                          // Verify if the segments ending before the end of the feature files...
      MixtureGD & adaptedMixture = ms.duplicateMixture(world,DUPL_DISTRIB);                 // Creating final as a copy of the world model
      MixtureGD & clientMixture= ms.duplicateMixture(world,DUPL_DISTRIB);
      long codeSelectedFrame=labelServer.getLabelIndexByString(labelSelectedFrames);        // Get the index of the cluster with in interest audio segments
      if (codeSelectedFrame==-1){                                                           // No data for this model !!!!!!!!!!!!!!
	cout << " WARNING - NO DATA FOR TRAINING ["<<*id<<"]";
	if (saveEmptyModel){
	  cout <<" World model is returned"<<endl;                                    // In this case, the client model is the world model
	  if (verbose) cout << "Save client model ["<<*id<<"]" << endl;
	  adaptedMixture.save(*id, config);                                           // Save the client model
	}
      }
      else{
	SegCluster& selectedSegments=segmentsServer.getCluster(codeSelectedFrame); // Gives the cluster of the selected/used segments                                   
        /// **** Factor Analysis Stuff
        XList faNdx;
        faNdx.addLine()=featureFileListp; 
        FactorAnalysisStat FA(faNdx,fs,config); // give all features to FA stats
        
        //FA.computeAndAccumulateGeneralFAStats(selectedSegments,fs,config);    
        for(int i=0;i<config.getParam("nbTrainIt").toLong();i++){
          if (verbose) cout << "------ Iteration ["<<i<<"] ------"<<endl;
          FA.computeAndAccumulateGeneralFAStats(selectedSegments,fs,config);                
          /*if (!varAdapt) FA.getTrueSpeakerModel(clientMixture,linep->getElement(1));
          else FA.getFactorAnalysisModel(clientMixture,linep->getElement(1));
          if (verbose) cout << "LLK for model["<<*id<<"] at it["<<i-1<<"]="<<FA.getLLK(selectedSegments,clientMixture,fs,config) << endl; */
          FA.estimateAndInverseL(config);
          FA.substractSpeakerStats();
          FA.getXEstimate();
          FA.substractChannelStats(); 
          FA.getYEstimate();    
      }      
      MixtureGD & sessionMixture= ms.duplicateMixture(world,DUPL_DISTRIB);
      bool saveSessionModel=false;
      if (config.existsParam("saveSessionModel")) saveSessionModel=true;
      if (saveSessionModel) FA.getSessionModel(sessionMixture,linep->getElement(1));
      if (!varAdapt) FA.getTrueSpeakerModel(clientMixture,linep->getElement(1)); // basically compute M_s_h=M+Dy_s and get a model
      else FA.getFactorAnalysisModel(clientMixture,linep->getElement(1)); // get FA variance adapted model
      if (verbose) cout << "Final LLK for model["<<*id<<"]="<<FA.getLLK(selectedSegments,clientMixture,fs,config) << endl;    

      /// **** End of FA
        if (!outputAdaptParam) {
            if (verbose) cout << "Save client model ["<<*id<<"]" << endl;
            clientMixture.save(*id, config);                                           // Save the client model
            if (saveSessionModel) {
              String sessionfile=*id+".session";
              if (verbose) cout << "Save session model ["<<sessionfile<<"]" << endl;              
              sessionMixture.save(sessionfile,config);   
            }              
        }
	if (!saveCompleteServer){
	  long tid=ms.getMixtureIndex(*id);      // TO BE SUPPRESSED BY
	  ms.deleteMixtures(tid,tid);            // ADDING a delete on a mixture pointor
	  ms.deleteUnusedDistribs();
	  }
      }
    }    
  } // fin try
catch (Exception& e) {cout << e.toString().c_str() << endl;}
  return 0;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
int TrainTargetJFA(Config& config)
{
	String inputClientListFileName = config.getParam("targetIdList");
	String inputWorldFilename = config.getParam("inputWorldFilename");
	String outputSERVERFilename = "";
	if (config.existsParam("mixtureServer")) outputSERVERFilename =config.getParam("mixtureServer");
	bool initByClient=false;                                              // In this case, the init model is read from the file
	if (config.existsParam("initByClient")) initByClient=config.getParam("initByClient").toBool();
	bool saveEmptyModel=false;
	if (config.existsParam("saveEmptyModel")) saveEmptyModel=config.getParam("saveEmptyModel").toBool();
	// label for selected frames - Only the frames associated with this label, in the label files, will be used
	bool fixedLabelSelectedFrame=true;
	String labelSelectedFrames;
	if (config.existsParam("useIdForSelectedFrame"))    // the ID of each speaker is used as labelSelectedFrame ?
		fixedLabelSelectedFrame=(config.getParam("useIdForSelectedFrame").toBool()==false);  
	if (fixedLabelSelectedFrame)                        // the label is decided by the command line and is unique for the run
		labelSelectedFrames=config.getParam("labelSelectedFrames");
	bool modelData=false;
	if (config.existsParam("useModelData")) modelData=config.getParam("useModelData").toBool();
	String initModelS=inputWorldFilename;
	if (modelData) if (config.existsParam("initModel")) initModelS=config.getParam("initModel"); // Use a specific model for Em init
	bool outputAdaptParam=false;
	if (config.existsParam("superVectors")) outputAdaptParam=true;
 
try{
	XList inputClientList(inputClientListFileName,config);          // read the Id + filenames for each client
	XLine * linep;
	inputClientList.getLine(0);
	MixtureServer ms(config);
	StatServer ss(config, ms);
	if (verbose) cout << "(TrainTarget) Joint Factor Analysis - Load world model [" << inputWorldFilename<<"]"<<endl;
	MixtureGD& world = ms.loadMixtureGD(inputWorldFilename);      
	if (verbose) cout <<"(TrainTarget) Use["<<initModelS<<"] for initializing EM"<<endl;
	
	//LOAD JFA MAtrices
	Matrix<double> U, V; 
	DoubleVector D;
	
	//Initialise EC matrix
	if(config.existsParam("eigenChannelMatrix")){
		String uName = config.getParam("matrixFilesPath") + config.getParam("eigenChannelMatrix") + config.getParam("loadMatrixFilesExtension");
 		U.load (uName, config);
		if (verboseLevel >=1) cout << "(TrainTargetJFA) Init EC matrix from "<< config.getParam("eigenChannelMatrix") <<"  from EigenChannel Matrix: "<<", rank: ["<<U.rows() << "] sv size: [" << U.cols() <<"]"<<endl;
	}
	else{
		unsigned long sS = world.getVectSize() * world.getDistribCount();
		U.setDimensions(1,sS);
		U.setAllValues(0.0);
		if (verboseLevel >=1) cout << "(TrainTargetJFA) Init EC matrix to 0"<<endl;
	}
	
	//Initialise EV matrix
	if(config.existsParam("eigenVoiceMatrix")){
		String vName = config.getParam("matrixFilesPath") + config.getParam("eigenVoiceMatrix") + config.getParam("loadMatrixFilesExtension");
		V.load (vName, config);
		if (verboseLevel >=1) cout << "(TrainTargetJFA) Init EV matrix from "<< config.getParam("eigenVoiceMatrix") <<"  from EigenVoice Matrix: "<<", rank: ["<<V.rows() << "] sv size: [" << V.cols() <<"]"<<endl;
	}
	else{
		unsigned long sS = world.getVectSize() * world.getDistribCount();
		V.setDimensions(1,sS);
		V.setAllValues(0.0);
		if (verboseLevel >=1) cout << "(TrainTargetJFA) Init EV matrix to 0"<<endl;
	}
	
	//Initialise D matrix
	if(config.existsParam("DMatrix")){
		String dName = config.getParam("matrixFilesPath") + config.getParam("DMatrix") + config.getParam("loadMatrixFilesExtension");
		Matrix<double> tmpD(dName, config);
		
		if( (tmpD.rows() != 1) || ( tmpD.cols() != world.getVectSize()*world.getDistribCount() ) ){
			throw Exception("Incorrect dimension of D Matrix",__FILE__,__LINE__);
		}
		else{
			D.setSize(world.getVectSize()*world.getDistribCount());
			D.setAllValues(0.0);
			for(unsigned long i=0; i<world.getVectSize()*world.getDistribCount(); i++){
				D[i] = tmpD(0,i);
			}
			if (verboseLevel >=1) cout << "(TrainTargetJFA) Init D matrix from "<<config.getParam("DMatrix")<<endl;
		}
	}
	else{
		unsigned long sS = world.getVectSize() * world.getDistribCount();
		D.setSize(sS);
		D.setAllValues(0.0);
		if (verboseLevel >1) cout << "(TrainTargetJFA) Init D matrix to 0"<<endl;
	}
	
	// *********** Target loop ***************** 
	while ((linep=inputClientList.getLine()) != NULL){             	// linep gives the XLine with the Id of a given client and the list of files

		String *id=linep->getElement();                              		// Get the Client ID (id)
		XLine featureFileListp=linep->getElements();	           	// Get the list of feature file for the client (end of the line)
		if (verbose) cout << "(TrainTarget) Train model ["<<*id<<"]"<<endl;
	
		XList ndx; ndx.addLine() = featureFileListp;
		JFAAcc jfaAcc(ndx,config,"TrainTarget");

		//Charger les matrices V, U et D a partir des objets matrice existant.
		jfaAcc.loadEV(V, config); jfaAcc.loadEC(U, config); jfaAcc.loadD(D);  

		//Initialise VU matrix
		jfaAcc.initVU();

		FeatureServer fs(config,featureFileListp);                                            			// Reading the features (from several files)
		SegServer segmentsServer;                                                             				// Create the segment server for managing the segments/clusters
		LabelServer labelServer;                                                              				// Create the lable server, for indexing the segments/clusters
		initializeClusters(featureFileListp,segmentsServer,labelServer,config);         		// Reading the segmentation files for each feature input file
		verifyClusterFile(segmentsServer,fs,config);                                     				// Verify if the segments ending before the end of the feature files...

		MixtureGD & adaptedMixture = ms.duplicateMixture(world,DUPL_DISTRIB);                 	// Creating final as a copy of the world model
		MixtureGD & clientMixture= ms.duplicateMixture(world,DUPL_DISTRIB);
		long codeSelectedFrame=labelServer.getLabelIndexByString(labelSelectedFrames);   	// Get the index of the cluster with in interest audio segments
		if (codeSelectedFrame==-1){                                                           					// No data for this model !!!!!!!!!!!!!!
			cout << " WARNING - NO DATA FOR TRAINING ["<<*id<<"]";
			if (saveEmptyModel){
				cout <<" World model is returned"<<endl;                                    				// In this case, the client model is the world model
				if (verbose) cout << "Save client model ["<<*id<<"]" << endl;
				adaptedMixture.save(*id, config);                                           					// Save the client model
			}
		}			

		else{
			SegCluster& selectedSegments=segmentsServer.getCluster(codeSelectedFrame); // Gives the cluster of the selected/used segments                                   

			//Compute the JFA statistics
			jfaAcc.computeAndAccumulateJFAStat(selectedSegments,fs,config);

			//Estimate X and Y in one time for each speaker
			jfaAcc.storeAccs();
			jfaAcc.estimateVUEVUT(config);

			jfaAcc.estimateAndInverseL_VU(config);

			jfaAcc.substractMplusDZ(config);

			jfaAcc.estimateYX();
			//Reinitialise the accumulators
			jfaAcc.resetTmpAcc();
			jfaAcc.restoreAccs();
			
			//Split X and Y estimates
			jfaAcc.splitYX();

			//Substract speaker and channel statistics M + VUYX
			jfaAcc.substractMplusVUYX();		
			//Estimate Z for each speaker
			jfaAcc.estimateZ();
			//Reinitialise the accumulators
			jfaAcc.resetTmpAcc();
			jfaAcc.restoreAccs();

			bool varAdapt = false;
			if((config.existsParam("varAdapt")) && ( config.getParam("varAdapt").toBool() )){
				varAdapt = true;
			}

			DoubleVector clientSV(jfaAcc.getSvSize(), jfaAcc.getSvSize());
			clientSV.setSize(jfaAcc.getSvSize());
			DoubleVector clientModel(jfaAcc.getSvSize(), jfaAcc.getSvSize());
			clientModel.setSize(jfaAcc.getSvSize());

			bool saveMixture = true;
			if((config.existsParam("saveMixture")) && !( config.getParam("saveMixture").toBool() ))	saveMixture = false;
			bool saveSuperVector = true;
			if((config.existsParam("saveSuperVector")) && !( config.getParam("saveSuperVector").toBool() ))	saveSuperVector = false;
			bool saveX = false;
			bool saveY = false;
			bool saveZ = false;


			if(config.existsParam("saveX"))			saveX = config.getParam("saveX").toBool();
			if(config.existsParam("saveY"))			saveY = config.getParam("saveY").toBool();
			if(config.existsParam("saveZ"))			saveZ = config.getParam("saveZ").toBool();
			String xExtension = ".x"; String yExtension = ".y"; String zExtension = ".z";
			if(config.existsParam("xExtension"))	xExtension = config.getParam("xExtension");
			if(config.existsParam("yExtension"))	xExtension = config.getParam("yExtension");
			if(config.existsParam("zExtension"))	xExtension = config.getParam("zExtension");

			jfaAcc.getVYplusDZ(clientSV, 0);
			jfaAcc.getMplusVYplusDZ(clientModel, 0);
			
			//WARNING !!!!! only the SuperVector model is divided by the UBM Co-Variance.
			for(unsigned long i=0; i<jfaAcc.getSvSize(); i++){
				clientSV[i] *= jfaAcc.getUbmInvVar()[i];
			}
			
			//Create the ClientMixture to save if required
			if(saveMixture){
				svToModel(clientModel, clientMixture);
				clientMixture.save(*id, config);
			}

			if(saveSuperVector){
				String svPath=config.getParam("vectorFilesPath");
				String svExt=config.getParam("vectorFilesExtension"); 
				String svFile=svPath+*id+svExt; 
				((Matrix<double>)clientSV).save(svFile,config);   
			}

			String svPath=config.getParam("vectorFilesPath");

			if(saveX){
				String xFile=svPath+*id+xExtension;
				jfaAcc.saveX(xFile,config);
			}
			if(saveY){
				String yFile=svPath+*id+yExtension;
				jfaAcc.saveY(yFile,config);
			}
			if(saveZ){
				String zFile=svPath+*id+zExtension;
				jfaAcc.saveZ(zFile,config);
			}

			long tid=ms.getMixtureIndex(*id);
			ms.deleteMixtures(tid,tid);
			ms.deleteUnusedDistribs();
		}
	}
} // fin try
catch (Exception& e) {cout << e.toString().c_str() << endl;}
return 0;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//int TrainTargetIVector(Config& config){
//
//	String inputClientListFileName = config.getParam("targetIdList");
//	String inputWorldFilename = config.getParam("inputWorldFilename");
//	String outputSERVERFilename = "";
//	if (config.existsParam("mixtureServer")) outputSERVERFilename =config.getParam("mixtureServer");
//	bool initByClient=false;                                              // In this case, the init model is read from the file
//	if (config.existsParam("initByClient")) initByClient=config.getParam("initByClient").toBool();
//	bool saveEmptyModel=false;
//	if (config.existsParam("saveEmptyModel")) saveEmptyModel=config.getParam("saveEmptyModel").toBool();
//	// label for selected frames - Only the frames associated with this label, in the label files, will be used
//	bool fixedLabelSelectedFrame=true;
//	String labelSelectedFrames;
//	if (config.existsParam("useIdForSelectedFrame"))    // the ID of each speaker is used as labelSelectedFrame ?
//		fixedLabelSelectedFrame=(config.getParam("useIdForSelectedFrame").toBool()==false);  
//	if (fixedLabelSelectedFrame)                        // the label is decided by the command line and is unique for the run
//		labelSelectedFrames=config.getParam("labelSelectedFrames");
//	bool modelData=false;
//	if (config.existsParam("useModelData")) modelData=config.getParam("useModelData").toBool();
//	String initModelS=inputWorldFilename;
//	if (modelData) if (config.existsParam("initModel")) initModelS=config.getParam("initModel"); // Use a specific model for Em init
//	bool outputAdaptParam=false;
//	if (config.existsParam("superVectors")) outputAdaptParam=true;
// 
//try{
//	XList inputClientList(inputClientListFileName,config);          // read the Id + filenames for each client
//	XLine * linep;
//	inputClientList.getLine(0);
//	MixtureServer ms(config);
//	StatServer ss(config, ms);
//	if (verbose) cout << "(TrainTarget) Joint Factor Analysis - Load world model [" << inputWorldFilename<<"]"<<endl;
//	MixtureGD& world = ms.loadMixtureGD(inputWorldFilename);      
//	if (verbose) cout <<"(TrainTarget) Use["<<initModelS<<"] for initializing EM"<<endl;
//	
//	//LOAD JFA MAtrices
//	Matrix<double> U, V; 
//	DoubleVector D;
//	
//	//Initialise EC matrix
//	unsigned long sS = world.getVectSize() * world.getDistribCount();
//	U.setDimensions(1,sS);
//	U.setAllValues(0.0);
//	if (verboseLevel >1) cout << "(TrainTargetJFA) Init EC matrix to 0"<<endl;
//	
//	//Initialise EV matrix
//	if(config.existsParam("eigenVoiceMatrix")){
//		String vName = config.getParam("matrixFilesPath") + config.getParam("eigenVoiceMatrix") + config.getParam("loadMatrixFilesExtension");
//		V.load (vName, config);
//		if (verboseLevel >=1) cout << "(TrainTargetJFA) Init EV matrix from "<< config.getParam("eigenVoiceMatrix") <<"  from EigenVoice Matrix: "<<", rank: ["<<V.rows() << "] sv size: [" << V.cols() <<"]"<<endl;
//	}
//	else{
//		unsigned long sS = world.getVectSize() * world.getDistribCount();
//		V.setDimensions(1,sS);
//		V.setAllValues(0.0);
//		if (verboseLevel >=1) cout << "(TrainTargetJFA) Init EV matrix to 0"<<endl;
//	}
//
//	D.setSize(sS);
//	D.setAllValues(0.0);
//	if (verboseLevel >1) cout << "(TrainTargetJFA) Init D matrix to 0"<<endl;
//	
//	// *********** Target loop ***************** 
//	while ((linep=inputClientList.getLine()) != NULL){             	// linep gives the XLine with the Id of a given client and the list of files
//
//		String *id=linep->getElement();                              		// Get the Client ID (id)
//		XLine featureFileListp=linep->getElements();	           	// Get the list of feature file for the client (end of the line)
//		if (verbose) cout << "(TrainTarget) Train model ["<<*id<<"]"<<endl;
//	
//		XList ndx; ndx.addLine() = featureFileListp;
//		JFAAcc jfaAcc(ndx,config);
//		
//		//Charger les matrices V, U et D a partir des objets matrice existant.
//		jfaAcc.loadEV(V, config); jfaAcc.loadEC(U, config); jfaAcc.loadD(D);  
//
//		//Initialise VU matrix
//		jfaAcc.initVU();
//
//		FeatureServer fs(config,featureFileListp);                                            			// Reading the features (from several files)
//		SegServer segmentsServer;                                                             				// Create the segment server for managing the segments/clusters
//		LabelServer labelServer;                                                              				// Create the lable server, for indexing the segments/clusters
//		initializeClusters(featureFileListp,segmentsServer,labelServer,config);         		// Reading the segmentation files for each feature input file
//		verifyClusterFile(segmentsServer,fs,config);                                     				// Verify if the segments ending before the end of the feature files...
//
//		MixtureGD & adaptedMixture = ms.duplicateMixture(world,DUPL_DISTRIB);                 	// Creating final as a copy of the world model
//		MixtureGD & clientMixture= ms.duplicateMixture(world,DUPL_DISTRIB);
//		long codeSelectedFrame=labelServer.getLabelIndexByString(labelSelectedFrames);   	// Get the index of the cluster with in interest audio segments
//		if (codeSelectedFrame==-1){                                                           					// No data for this model !!!!!!!!!!!!!!
//			cout << " WARNING - NO DATA FOR TRAINING ["<<*id<<"]";
//			if (saveEmptyModel){
//				cout <<" World model is returned"<<endl;                                    				// In this case, the client model is the world model
//				if (verbose) cout << "Save client model ["<<*id<<"]" << endl;
//				adaptedMixture.save(*id, config);                                           					// Save the client model
//			}
//		}			
//
//		else{
//			SegCluster& selectedSegments=segmentsServer.getCluster(codeSelectedFrame); // Gives the cluster of the selected/used segments                                   
//
//			//Compute the JFA statistics
//			jfaAcc.computeAndAccumulateJFAStat(selectedSegments,fs,config);
//
//			//Estimate X and Y in one time for each speaker
//			jfaAcc.storeAccs();
//			jfaAcc.estimateVUEVUT(config);
//			jfaAcc.estimateAndInverseL_VU(config);
//			jfaAcc.substractMplusDZ(config);
//			jfaAcc.estimateYX();
//			//Reinitialise the accumulators
//			jfaAcc.resetTmpAcc();
//			jfaAcc.restoreAccs();
//			
//			//Split X and Y estimates
//			jfaAcc.splitYX();
//
//			//Reinitialise the accumulators
//			jfaAcc.resetTmpAcc();
//			jfaAcc.restoreAccs();
//
//			//Save the I-Vector for the current speaker
//			DoubleVector clientIVect(jfaAcc.getSvSize(), jfaAcc.getSvSize());
//			String svPath=config.getParam("vectorFilesPath");
//			String svExt=config.getParam("vectorFilesExtension"); 
//			String svFile=svPath+*id+svExt; 
//			((Matrix<double>)clientIVect).save(svFile,config);
//
//		}
//	}
//} // fin try
//catch (Exception& e) {cout << e.toString().c_str() << endl;}
//return 0;
//}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
int TrainTargetLFA(Config& config)
{
	String inputClientListFileName = config.getParam("targetIdList");
	String inputWorldFilename = config.getParam("inputWorldFilename");
	String outputSERVERFilename = "";
	if (config.existsParam("mixtureServer")) outputSERVERFilename =config.getParam("mixtureServer");
	bool initByClient=false;                                              // In this case, the init model is read from the file
	if (config.existsParam("initByClient")) initByClient=config.getParam("initByClient").toBool();
	bool saveEmptyModel=false;
	if (config.existsParam("saveEmptyModel")) saveEmptyModel=config.getParam("saveEmptyModel").toBool();
	// label for selected frames - Only the frames associated with this label, in the label files, will be used
	bool fixedLabelSelectedFrame=true;
	String labelSelectedFrames;
	if (config.existsParam("useIdForSelectedFrame"))    // the ID of each speaker is used as labelSelectedFrame ?
		fixedLabelSelectedFrame=(config.getParam("useIdForSelectedFrame").toBool()==false);  
	if (fixedLabelSelectedFrame)                        // the label is decided by the command line and is unique for the run
		labelSelectedFrames=config.getParam("labelSelectedFrames");
	bool modelData=false;
	if (config.existsParam("useModelData")) modelData=config.getParam("useModelData").toBool();
	String initModelS=inputWorldFilename;
	if (modelData) if (config.existsParam("initModel")) initModelS=config.getParam("initModel"); // Use a specific model for Em init
	bool outputAdaptParam=false;
	if (config.existsParam("superVectors")) outputAdaptParam=true;
 
try{
	XList inputClientList(inputClientListFileName,config);          // read the Id + filenames for each client
	XLine * linep;
	inputClientList.getLine(0);
	MixtureServer ms(config);
	StatServer ss(config, ms);
	if (verbose) cout << "(TrainTarget) Latent Factor Analysis - Load world model [" << inputWorldFilename<<"]"<<endl;
	MixtureGD& world = ms.loadMixtureGD(inputWorldFilename);      
	if (verbose) cout <<"(TrainTarget) Use["<<initModelS<<"] for initializing EM"<<endl;
	
	//LOAD JFA MAtrices
	unsigned long svsize=world.getDistribCount()*world.getVectSize();
	Matrix<double> U, V; 
	DoubleVector D(svsize,svsize);
	
	//Initialise EC matrix
	if(config.existsParam("eigenChannelMatrix")){
		String uName = config.getParam("matrixFilesPath") + config.getParam("eigenChannelMatrix") + config.getParam("loadMatrixFilesExtension");
 		U.load (uName, config);
		if (verboseLevel >=1) cout << "(TrainTargetLFA) Init EC matrix from "<< config.getParam("eigenChannelMatrix") <<"  from EigenChannel Matrix: "<<", rank: ["<<U.rows() << "] sv size: [" << U.cols() <<"]"<<endl;
	}
	else{
		U.setDimensions(1,svsize);
		U.setAllValues(0.0);
		if (verboseLevel >1) cout << "(TrainTargetLFA) Init EC matrix to 0"<<endl;
	}
	
	V.setDimensions(1,svsize);
	V.setAllValues(0.0);
	if (verboseLevel >=1) cout << "(TrainTargetLFA) Init EV matrix to 0"<<endl;

	//Initialise the D matrix for MAP adaptation
	for(unsigned long i=0; i<world.getDistribCount(); i++){
		for(unsigned long j = 0; j<world.getVectSize(); j++){
			D[i*world.getVectSize()+j] = sqrt(1.0/(world.getDistrib(i).getCovInv(j)*config.getParam("regulationFactor").toDouble()));
		}
	}

	// *********** Target loop ***************** 
	while ((linep=inputClientList.getLine()) != NULL){             	// linep gives the XLine with the Id of a given client and the list of files

		String *id=linep->getElement();                              		// Get the Client ID (id)
		XLine featureFileListp=linep->getElements();	           	// Get the list of feature file for the client (end of the line)
		if (verbose) cout << "(TrainTargetLFA) Train model ["<<*id<<"]"<<endl;
	
		XList ndx; ndx.addLine() = featureFileListp;
		JFAAcc jfaAcc(ndx,config,"TrainTarget");
		
		//Charger les matrices V, U et D a partir des objets matrice existant.
		jfaAcc.loadEV(V, config); jfaAcc.loadEC(U, config); jfaAcc.loadD(D);  

		//Initialise VU matrix
		jfaAcc.initVU();

		FeatureServer fs(config,featureFileListp);											// Reading the features (from several files)
		SegServer segmentsServer;															// Create the segment server for managing the segments/clusters
		LabelServer labelServer;															// Create the lable server, for indexing the segments/clusters
		initializeClusters(featureFileListp,segmentsServer,labelServer,config);				// Reading the segmentation files for each feature input file
		verifyClusterFile(segmentsServer,fs,config);                                     	// Verify if the segments ending before the end of the feature files...

		MixtureGD & adaptedMixture = ms.duplicateMixture(world,DUPL_DISTRIB);               // Creating final as a copy of the world model
		MixtureGD & clientMixture= ms.duplicateMixture(world,DUPL_DISTRIB);
		long codeSelectedFrame=labelServer.getLabelIndexByString(labelSelectedFrames);   	// Get the index of the cluster with in interest audio segments
		if (codeSelectedFrame==-1){                                                         	// No data for this model !!!!!!!!!!!!!!
			cout << " WARNING - NO DATA FOR TRAINING ["<<*id<<"]";
			if (saveEmptyModel){
				cout <<" World model is returned"<<endl;                                    // In this case, the client model is the world model
				if (verbose) cout << "Save client model ["<<*id<<"]" << endl;
				adaptedMixture.save(*id, config);                                           // Save the client model
			}
		}			

		else{
			SegCluster& selectedSegments=segmentsServer.getCluster(codeSelectedFrame); // Gives the cluster of the selected/used segments                                   

			//Compute the JFA statistics
			jfaAcc.computeAndAccumulateJFAStat(selectedSegments,fs,config);

			//Estimate X and Y in one time for each speaker
			jfaAcc.storeAccs();
			jfaAcc.estimateVUEVUT(config);
			jfaAcc.estimateAndInverseL_VU(config);
			jfaAcc.substractMplusDZ(config);
			jfaAcc.estimateYX();
			//Reinitialise the accumulators
			jfaAcc.resetTmpAcc();
			jfaAcc.restoreAccs();
			
			//Split X and Y estimates
			jfaAcc.splitYX();

			//Substract speaker and channel statistics M + VUYX
			jfaAcc.substractMplusVUYX();		
			//Estimate Z for each speaker
			double tau = config.getParam("regulationFactor").toLong();
			jfaAcc.estimateZMAP(tau);
			//Reinitialise the accumulators
			jfaAcc.resetTmpAcc();
			jfaAcc.restoreAccs();

			bool varAdapt = false;
			if((config.existsParam("varAdapt")) && ( config.getParam("varAdapt").toBool() )){
				varAdapt = true;
			}

			DoubleVector clientModel(jfaAcc.getSvSize(), jfaAcc.getSvSize());
			clientModel.setSize(jfaAcc.getSvSize());

			jfaAcc.getMplusVYplusDZ(clientModel, 0);
			
			//Create the ClientMixture
			svToModel(clientModel, clientMixture);
			clientMixture.save(*id, config);

			long tid=ms.getMixtureIndex(*id);
			ms.deleteMixtures(tid,tid);
			ms.deleteUnusedDistribs();
		}
	}
} // fin try
catch (Exception& e) {cout << e.toString().c_str() << endl;}
return 0;
}

#endif //!defined(ALIZE_TrainTarget_cpp)

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

#if !defined(ALIZE_TrainWorld_cpp)
#define ALIZE_TrainWorld_cpp

#include <iostream>
#include <cmath>
#include "liatools.h"
#include "TrainWorld.h"

using namespace alize;
using namespace std;

void featureStream(Config &config,String filename,FeatureServer *&fs,SegServer *&segServ,SegCluster *&segCluster,String labelSelectedFrames){
  fs=new FeatureServer(config,filename);
  try{   
  	// TODO Test if the stream is not empty
  	
  	segServ=new SegServer;                                                                      // Create the segment server for managing the segments/clusters
  	LabelServer labelServer;                                                                    // Create the label server, for indexing the segments/clusters
 	initializeClusters(filename,*segServ,labelServer,config);                                   // Reading the segmentation files for each feature input file
 	verifyClusterFile(*segServ,*fs,config);                                                     // Verify if the segments ending before the end of the feature files...
	unsigned long codeSelectedFrame=labelServer.getLabelIndexByString(labelSelectedFrames);     // Get the index of the cluster with in interest audio segments
 	segCluster=&(segServ->getCluster(codeSelectedFrame));
  }
  catch (Exception& e){
    cout << e.toString() << endl;
  }
}
void reserveMem(FeatureServer** &fsTab,SegServer** &segServTab,SegCluster** &segTab,double *&weightTab,unsigned long nbStream){
  fsTab=new FeatureServer*[nbStream];
  segServTab=new SegServer*[nbStream];
  segTab=new SegCluster*[nbStream];
  weightTab=new double [nbStream];
  for (unsigned long i=0;i<nbStream;i++)weightTab[i]=1/(double) nbStream;
}
void freeMem(FeatureServer** &fsTab,SegServer** &segServTab,SegCluster** &segTab,double *&weightTab,unsigned long nbStream){
  for (unsigned long i=0;i<nbStream;i++){
    delete fsTab[i]; 
    delete segServTab[i];
  }
  delete [] fsTab;
  delete [] segServTab;
  delete [] segTab;
  delete [] weightTab;
}

//-------------------------------------------------------------------------
int trainWorld(Config& config){
  if (verbose) cout << "Begin world model training"<<endl;   
  try{    
    // Reading the data, one or multiple separate streams
    unsigned long nbStream=0;                                                             // Number of Streams
    FeatureServer **fsTab=NULL;                                                           // Array of FeatureServer (address) - one by input stream
    SegServer     **segServTab=NULL;                                                      // Array of segment server (address)- one by input stream
    SegCluster    **segTab=NULL;                                                          // Array of selected segments cluster(address) - one by stream
    double         *weightTab=NULL;                                                       // Array of weight of each stream. i.e influence of a stream on the final model
    String outputWorldFilename = config.getParam("outputWorldFilename");                  // output worldmodel file filename                            
    bool fileInit=config.existsParam("inputWorldFilename");                               // if a inputWorlFilename is given, init by file, else from scratch
    bool saveInitModel=true;
    if (config.existsParam("saveInitModel")) saveInitModel=config.getParam("saveInitModel").toBool();
    String inputWorldFilename="";
    if (fileInit) inputWorldFilename=config.getParam("inputWorldFilename");                // if file init, the initial model filename
    String labelSelectedFrames =config.getParam("labelSelectedFrames");                    // label for selected frames
    TrainCfg trainCfg(config);                                                             // Get the training algo params

    // Reading the data
    if(config.existsParam("inputStreamList")){// We want to work on separated list 
      XList tmp(config.getParam("inputStreamList"),config);                                        // Each data set influence will be balanced during training
      XLine & listInputFilename=tmp.getAllElements();                                              // Read the list of (list) filenames in tmp -> listInputFilename
      nbStream=listInputFilename.getElementCount();
      if (nbStream==0) throw Exception("TrainWorld error:no input stream" , __FILE__, __LINE__);
      reserveMem(fsTab,segServTab,segTab,weightTab,nbStream);
      for (unsigned i=0;i<nbStream;i++)
		featureStream(config,listInputFilename.getElement(i),fsTab[i],segServTab[i],segTab[i],labelSelectedFrames);
      if (config.existsParam("weightStreamList")){ // Read the weight of each stream, text file
		XList tmpW(config.getParam("weightStreamList"),config);
		XLine & listW=tmpW.getAllElements();                                              // Read the list of (list) filenames in tmp -> listInputFilename
		if (listW.getElementCount()!=nbStream) throw Exception("TrainWorld error: number of weigths differs than number of input streams" , __FILE__, __LINE__);
		for (unsigned i=0;i<nbStream;i++) weightTab[i]=listW.getElement(i).toDouble();
      }
    }
    else{ // Only one input stream, no stream list
      nbStream=1;
      reserveMem(fsTab,segServTab,segTab,weightTab,nbStream);
      featureStream(config,config.getParam("inputFeatureFilename"),fsTab[0],segServTab[0],segTab[0],labelSelectedFrames);
    }
    unsigned long vectSize=fsTab[0]->getVectSize();                                          // size of the input vectors
    // Create stat server and mixture server
    MixtureServer ms(config);
    StatServer ss(config, ms);
    if (debug || verbose) cout << "Stream mode, nb Stream="<<nbStream<<endl;
    if (debug|| (verboseLevel>2)){
      for (unsigned long i=0;i<nbStream;i++){
		cout <<"Stream["<<i<<"]"<<endl;
		segTab[i]->rewind(); 
		Seg *seg;                                                                            // Reset to the first segment
		while((seg=segTab[i]->getSeg())!=NULL)                                         // For each of the selected segments
	  		cout << "File["<<seg->sourceName()<<"] Segment begin["<<
	    		seg->begin()<<"] length["<<seg->length()<<"] index in the feature server["<<fsTab[i]->getFirstFeatureIndexOfASource(seg->sourceName())<<"]"<<endl;
      }
    }  
    // Global mean and variance matrices initialisation (computed from dataa or set to 0,1)
    bool use01=false;
    if (config.existsParam("use01")) use01=config.getParam("use01").toBool();
    if (verbose){ if (use01) cout<<"Use 0 mean, 1 cov "<<endl; else cout << "Compute global mean and cov"<<endl;}
    DoubleVector globalMean;
    DoubleVector globalCov;
    if (!use01){
      FrameAccGD globalFrameAcc;
      unsigned long nbFrame=computeMeanCov(config,fsTab,segTab,nbStream,globalMean,globalCov);                             // Compute the global mean and covariance
      if (verboseLevel>1){
	cout <<"global mean and cov of training data, number of frame= ["<<nbFrame<<"]"<<endl;
	for (unsigned i=0; i < vectSize; i++)cout << "mean[" << i << "=" << globalMean[i] << "]\tcov[" << globalCov[i] << "]" << endl;
      }
    }
    else initialize01(vectSize,globalMean,globalCov);
    MixtureGD &world=ms.createMixtureGD();
    if (fileInit){                                                                                        // Load or initialize the initial model
      if (verbose) cout << "Load initial world model ["<<inputWorldFilename<<"]" << endl;                 
      world=ms.loadMixtureGD(inputWorldFilename);                                                         // Load
    } 
    else{ 
      if (verbose) cout <<"World model init from scratch"<<endl;
      mixtureInit(ms,fsTab,segTab,weightTab,nbStream,world,globalCov,config,trainCfg);                             // Initialize    
      if (saveInitModel) world.save(outputWorldFilename+"init", config);
    }
    MixtureGD *newWorld=&world; // TODO Verify and suppress...
    trainModelStream(config,ms,ss,fsTab,segTab,weightTab,nbStream,globalMean,globalCov,newWorld,trainCfg);
    if (verbose) cout << "Save world model ["<<outputWorldFilename<<"]" << endl;
    newWorld->save(outputWorldFilename, config);                                          
    // Cleaning the memory
    freeMem(fsTab,segServTab,segTab,weightTab,nbStream);
  }
  catch (Exception& e){
    cout << e.toString() << endl;
  }
  return 0;
}


#endif // !defined(ALIZE_TrainWorld_cpp)

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

#if !defined(ALIZE_UnsupervisedTools_cpp)
#define ALIZE_UnsupervisedTools_cpp

#include <iostream>
#include <fstream>		// pour outFile
#include <cstdio>		// pour printf()
#include <cassert>		// pour le debug pratique
#include <cmath>
#include <liatools.h>
#include <DoubleSquareMatrix.h>
#include <RealVector.h>
#include "FileInfo.h"
#include "UnsupervisedTools.h"

// Methods of class WindowLLR
void WindowLLR::_initMem(){
    _freeMem();
    _idxA= new ULongVector(_size,_size);
    _accLlrA = new DoubleVector(_nClient,_nClient);
    _llrM= new Matrix <double>(_size,_nClient);
    for (unsigned long idxC=0;idxC<_nClient;idxC++){
	for (unsigned long idxF=0;idxF<_size;idxF++)
	    (*_llrM)(idxF,idxC)=0;
	(*_accLlrA)[idxC]=0;
    }
   
    _bIdx=0;
    _count=0;
}
void WindowLLR::_freeMem(){
    if (_llrM) {
	delete _llrM;
	delete _accLlrA;
	delete _idxA;
    }
    _llrM=NULL;
}
WindowLLR::WindowLLR(Config &config){
    _set=false;
    _size=0;
    _dec=0;
    _bIdx=0;
    _count=0;
    _nClient=0;
    _llrM=NULL;
    if (config.existsParam("windowLLR")) _set=config.getParam("windowLLR").toBool();
    if (_set){
	if (config.existsParam("windowLLRSize")) _size=config.getParam("windowLLRSize").toLong();
	else _size=30;
	if (config.existsParam("windowLLRDec")) _dec=config.getParam("windowLLRDec").toLong();
	else _dec=_size;	
	_nClient=1;
	_initMem();
    }
}
WindowLLR::~WindowLLR(){
    _freeMem();
}
void WindowLLR::showConfig(){
    if (_set) cout<<"windowLLR mode size["<<_size<<"] dec["<<_dec<<"]"<<endl; 
}
unsigned long WindowLLR::wCount(){
    return (_count);
}
void WindowLLR::dec(unsigned long idxFrame){
    if (_count<_size){       //window is not full
	_count++;
	unsigned long eIdx=(_bIdx+_count-1)%_size;
	(*_idxA)[eIdx]=idxFrame;
    }
    else{// window is full, real dec (shift the window, step _dec frame)
	for (unsigned long wIdx=0;wIdx<_dec;wIdx++){
	    //suppress the begin value
	    for (unsigned long cIdx=0;cIdx<_nClient;cIdx++)
		(*_accLlrA)[cIdx]-=(*_llrM)(_bIdx,cIdx);	    	    
	    _bIdx=(_bIdx+1)%_size;
	}
	_count-=(_dec-1);
	(*_idxA)[(_bIdx+_count-1)%_size]=idxFrame;	
    }   
}
void WindowLLR::accLLR(unsigned long clientIdx,double llr){
    (*_llrM)((_bIdx+_count-1)%_size,clientIdx)=llr;
    (*_accLlrA)[clientIdx]+=llr;
}
double WindowLLR::getLLR(unsigned long clientIdx){
    return (*_accLlrA)[clientIdx]/(double)_size;
}   
bool WindowLLR::isEnd(){
    return (wCount()==_size);
}


//-------------------------------------------------------------------------------------------------------


//-- Accumulate the occupation for the selected frames and a given model

void accumulateStatLK(StatServer & ss, FeatureServer & fs, MixtureStat & acc,
  unsigned long idxBeginFrame, unsigned long nbFrames, Config & config)
{
    fs.seekFeature(idxBeginFrame);	// go to the frame in the buffer (and load it if needed)
    for (unsigned long n = 0; n < nbFrames; n++)
    {
        Feature f;
        if (fs.readFeature(f) == false)
	    cout << "No more features" << endl;
        acc.computeAndAccumulateLLK(f);
    }
}

// one a Segment
void
accumulateStatLK(StatServer & ss, FeatureServer & fs, MixtureStat & acc,
  Seg * seg, Config & config)
{
    unsigned long begin = seg->begin() + fs.getFirstFeatureIndexOfASource(seg->sourceName());	// Find the index of the first frame of the file in the buffer
    accumulateStatLK(ss, fs, acc, begin, seg->length(), config);
}

// One on Cluster
void accumulateStatLK(StatServer & ss, FeatureServer & fs, MixtureStat & acc,
  SegCluster & selectedSegments, Config & config)
{
    Seg *seg;			// reset the reader at the begin of the input stream
    selectedSegments.rewind();
    while ((seg = selectedSegments.getSeg()) != NULL)	// For each of the selected segments
        accumulateStatLK(ss, fs, acc, seg, config);
}






// Estimation client model from different Feature server with a weight associated
// Using EM and MAP 
void adaptModel(Config & config, StatServer & ss, MixtureServer & ms,
  ObjectRefVector & FeatServ, ObjectRefVector & ClusterSeg,
  MixtureGD & aprioriModel, MixtureGD & clientMixture,
  DoubleVector & decision)
{
    MAPCfg mapCfg(config);
    if (verbose)
        mapCfg.showConfig(cout);
    if (verboseLevel > 1)
        cout << "Mean LLK Init = " << meanLikelihood(ss, FeatServ, ClusterSeg, clientMixture, decision, config) << endl;
    
    for (unsigned long trainIt = 0; trainIt < mapCfg.getNbTrainIt(); trainIt++)
    {				// Begin the initial adaptation loop (with bagged frames)
	// Create a statistic accumulator using the curent model
        MixtureStat & emAcc = ss.createAndStoreMixtureStat(clientMixture);	
        emAcc.resetEM();
        double llkPreviousIt = 0;
        for (unsigned long nbFs = 0; nbFs < FeatServ.size(); nbFs++)
	{

	    SegServer segServer;	// Create a local segment server 
	    // Create the cluster for describing the selected frames
	    SegCluster & baggedFramesCluster = segServer.createCluster(1, "", "");	
	    baggedSegments((static_cast < SegCluster & >(ClusterSeg.getObject(nbFs))),baggedFramesCluster, mapCfg.getBaggedFrameProbability());

	    if (verboseLevel > 2)
	        cout <<"Accumulate statistics on the feature server weighted by : " << decision[nbFs] << endl;
	    
	    llkPreviousIt += accumulateStatEM(ss, (static_cast < FeatureServer & >(FeatServ.getObject(nbFs))), emAcc, baggedFramesCluster, decision[nbFs], config);	// Accumulate the EM statistics

	}
	
        clientMixture = emAcc.getEM();	// Get the EM estimate   
        unsigned long frameCount = (unsigned long) emAcc.getEMFeatureCount();
        cout << "Total Frames :" << frameCount << endl;
        llkPreviousIt = llkPreviousIt / (double) frameCount;
        if (verbose)
	    cout << "ML (partial) estimate it[" << trainIt <<"] (take care, it corresponds to the previous it,0 means init likelihood) = "
	  << llkPreviousIt << endl;
	
        computeMAP(ms, aprioriModel, clientMixture, frameCount, config);	// Bayesian Adaptation client=MAP(aprioriModel,client)
        if (mapCfg.getNormalizeModel())
	    normalizeMixture(clientMixture, mapCfg, config);	// Normalize/fit the model if needed
        ss.deleteMixtureStat(emAcc);

        if (verboseLevel > 2)
	    cout << "Likelihood on all frames =" << meanLikelihood(ss, FeatServ,
	  ClusterSeg, clientMixture, decision, config) << endl;

    }

  if (verboseLevel > 1)
      cout << "Final likelihood on all frames =" << meanLikelihood(ss, FeatServ,ClusterSeg, clientMixture, decision, config) << endl;

  if (debug)
      cout << "adaptModel nb distrib:" << ms.getDistribCount() << "nb mixt:" << ms.getMixtureCount() << endl;
}










//------------------------------------------------------------------------
// ** New training algo based on a true EM/ML estimate of the training data before to apply MAP
void modelBasedadaptModelEM(Config & config, StatServer & ss,
  MixtureServer & ms, FeatureServer & fs, SegCluster & selectedSegments,
  FeatureServer & fsTests, SegCluster & selectedSegmentsTests,
  MixtureGD & aprioriModel, MixtureGD & clientMixture, MixtureGD & initModel)
{
    MAPCfg mapCfg(config);
    if (verbose)
    {
        cout << "Model adaptation based on true EM/ML estimate of training data"<< endl;

    }
    MixtureServer msTmp(config);
    MixtureGD & data = msTmp.duplicateMixture(initModel, DUPL_DISTRIB);

    /*unsigned long totalFrameCount1 = totalFrame(selectedSegments);
     unsigned long totalFrameCount2 = totalFrame(selectedSegmentsTests);
     unsigned long totalFrameCount = totalFrameCount1 + totalFrameCount2; */
    
    if (verboseLevel > 1)
    cout << "Mean LLK Init = " << meanLikelihood(ss, fs, data,selectedSegments, config) << endl;
    
    for (unsigned long emIt = 0; emIt < mapCfg.getNbEmIt(); emIt++)
    {	
	// begin the true EM/ML estimate of the adpatation data 
	    
	// Create a statistic accumulator using the curent model
        MixtureStat & emAcc = ss.createAndStoreMixtureStat(data);	
        SegServer segServer;	// Create a local segment server
	    
	// Create the cluster for describing the selected frames    
        SegCluster & baggedFramesCluster = segServer.createCluster(1, "", "");	
        baggedSegments(selectedSegments, baggedFramesCluster,
	mapCfg.getBaggedFrameProbability());
        SegServer segServerTests;	// Create a local segment server
	    
	// Create the cluster for describing the selected frames
        SegCluster & baggedFramesClusterTests = segServerTests.createCluster(1, "", "");	
        baggedSegments(selectedSegmentsTests, baggedFramesClusterTests,
	mapCfg.getBaggedFrameProbability());
        emAcc.resetEM();
	    
	// Accumulate the EM statistics on the first feature server
        double llkPreviousIt = accumulateStatEM(ss, fs, emAcc, baggedFramesCluster, config);
	    
	// Accumulate the EM statistics on the second feature server
        llkPreviousIt += accumulateStatEM(ss, fsTests, emAcc, baggedFramesClusterTests, config);	
        data = emAcc.getEM();	// Get the EM estimate         
        unsigned long frameCount = (unsigned long) emAcc.getEMFeatureCount();
        llkPreviousIt = llkPreviousIt / (double) frameCount;
	
        if (verbose)
	    cout << "ML (partial) estimate it[" << emIt <<"] (take care, it corresponds to the previous it,0 means init likelihood) = "
	  << llkPreviousIt << endl;
	
        ss.deleteMixtureStat(emAcc);
    }
    // Begin the estimation of the statistic using the EM/ML model of the adaptation data
    unsigned long modelNbComp = aprioriModel.getDistribCount();
    
    // Begin the estimation of the statistic using the EM/ML model of the adaptation data
    
    // Complete log likelihood of the adaptation data given the apriori model
    unsigned long vectSize = fs.getVectSize();	
    for (unsigned long idxModel = 0; idxModel < modelNbComp; idxModel++)
    {	
	// Initialize the client mixture
        DistribGD & c = clientMixture.getDistrib(idxModel);
	    
        for (unsigned long idxC = 0; idxC < vectSize; idxC++)
	{
	    c.setMean(0.00, idxC);
	    c.setCov(0.00, idxC);
	}
    }
    DoubleVector apProbaTot(modelNbComp, modelNbComp);
    apProbaTot.setAllValues(0.0);
    for (unsigned long idxData = 0; idxData < data.getDistribCount(); idxData++)
    {
        if (debug)
	cout << "Distrib Data[" << idxData << "]" << endl;
        DistribGD & d = data.getDistrib(idxData);
        double totLk = 0.0;	// Likelihood of the current data component given the apriori model
        DoubleVector apProba(modelNbComp, modelNbComp);
        apProba.setAllValues(0.0);
	
        for (unsigned long idxModel = 0; idxModel < modelNbComp; idxModel++)
	{
	    if (debug)
	        cout << "Distrib A Priori model[" << idxModel << "]" << endl;
	    
	    DistribGD & m = aprioriModel.getDistrib(idxModel);
	    apProba[idxModel] = aprioriModel.weight(idxModel) * likelihoodGD(d, m);
	    totLk += apProba[idxModel];
	}
        for (unsigned long idxModel = 0; idxModel < modelNbComp; idxModel++)
	{
	    DistribGD & c = clientMixture.getDistrib(idxModel);
	    apProba[idxModel] /= totLk;
		
	    for (unsigned long idxC = 0; idxC < vectSize; idxC++)
	    {
	        c.setMean(c.getMean(idxC) +(d.getMean(idxC) * apProba[idxModel] * data.weight(idxData)),
		idxC);
	        c.setCov(c.getCov(idxC) +((d.getMean(idxC) * d.getMean(idxC)) * apProba[idxModel] *
		  data.weight(idxData)), idxC);
	    }
	    
	    apProbaTot[idxModel] += apProba[idxModel] * data.weight(idxData);
	}
    }
  for (unsigned long idxModel = 0; idxModel < modelNbComp; idxModel++)
    {
      DistribGD & c = clientMixture.getDistrib(idxModel);
      for (unsigned long idxC = 0; idxC < vectSize; idxC++)
	{
	  c.setMean(c.getMean(idxC) / apProbaTot[idxModel], idxC);
	  //        c.setCov(c.getCov(idxC),idxC);
	}
    }
  for (unsigned long idxModel = 0; idxModel < modelNbComp; idxModel++)
    {
      DistribGD & c = clientMixture.getDistrib(idxModel);
      c.computeAll();

    }
  //


}

//-------------------------------------------------------------------------




void modelBasedadaptModelEM(Config & config, StatServer & ss,
  MixtureServer & ms, FeatureServer & fs, SegCluster & selectedSegments,
  MixtureGD & aprioriModel, MixtureGD & clientMixture, MixtureGD & initModel)
{
  MAPCfg mapCfg(config);
  if (verbose)
  {
      cout << "Model adaptation based on true EM/ML estimate of training data"
	<< endl;

  }
  MixtureServer msTmp(config);
  MixtureGD & data = msTmp.duplicateMixture(initModel, DUPL_DISTRIB);

  //  unsigned long totalFrameCount = totalFrame(selectedSegments);
  //if (verboseLevel>1) cout << "Mean LLK Init = " << meanLikelihood(ss,fs,data,selectedSegments,config)<< endl;    
  for (unsigned long emIt = 0; emIt < mapCfg.getNbEmIt(); emIt++)
  {
      // begin the true EM/ML estimate of the adpatation data 
	  
      // Create a statistic accumulator using the curent model  
      MixtureStat & emAcc = ss.createAndStoreMixtureStat(data);	
      SegServer segServer;	// Create a local segment server
	  
      // Create the cluster for describing the selected frames
      SegCluster & baggedFramesCluster = segServer.createCluster(1, "", "");
	  
      baggedSegments(selectedSegments, baggedFramesCluster,mapCfg.getBaggedFrameProbability());
      emAcc.resetEM();
	  
      // Accumulate the EM statistics
      double llkPreviousIt = accumulateStatEM(ss, fs, emAcc, baggedFramesCluster, config);
	  
      data = emAcc.getEM();	// Get the EM estimate         
      unsigned long frameCount = (unsigned long) emAcc.getEMFeatureCount();
      llkPreviousIt = llkPreviousIt / (double) frameCount;
	  
      if (verbose)
	  cout << "ML (partial) estimate it[" << emIt << "] (take care, it corresponds to the previous it,0 means init likelihood) = "
	  << llkPreviousIt << endl;
      
      ss.deleteMixtureStat(emAcc);
    }
    unsigned long modelNbComp = aprioriModel.getDistribCount();
    // Begin the estimation of the statistic using the EM/ML model of the adaptation data
    
    // Complete log likelihood of the adaptation data given the apriori model
    unsigned long vectSize = fs.getVectSize();	
    
    for (unsigned long idxModel = 0; idxModel < modelNbComp; idxModel++)
    {				// Initialize the client mixture
        DistribGD & c = clientMixture.getDistrib(idxModel);
        for (unsigned long idxC = 0; idxC < vectSize; idxC++)
	{
	    c.setMean(0.00, idxC);
	    c.setCov(0.00, idxC);
	}
    }
    
    DoubleVector apProbaTot(modelNbComp, modelNbComp);
    apProbaTot.setAllValues(0.0);
    for (unsigned long idxData = 0; idxData < data.getDistribCount(); idxData++)
    {
        if (debug)
	    cout << "Distrib Data[" << idxData << "]" << endl;
	
        DistribGD & d = data.getDistrib(idxData);
        double totLk = 0.0;	// Likelihood of the current data component given the apriori model
	
        DoubleVector apProba(modelNbComp, modelNbComp);
        apProba.setAllValues(0.0);
        for (unsigned long idxModel = 0; idxModel < modelNbComp; idxModel++)
	{
	    if (debug)
	        cout << "Distrib A Priori model[" << idxModel << "]" << endl;
	    DistribGD & m = aprioriModel.getDistrib(idxModel);
	    apProba[idxModel] = aprioriModel.weight(idxModel) * likelihoodGD(d, m);
	    totLk += apProba[idxModel];
	}
      for (unsigned long idxModel = 0; idxModel < modelNbComp; idxModel++)
      {
	  DistribGD & c = clientMixture.getDistrib(idxModel);
	  apProba[idxModel] /= totLk;
		
	  for (unsigned long idxC = 0; idxC < vectSize; idxC++)
	    {
	        c.setMean(c.getMean(idxC) +(d.getMean(idxC) * apProba[idxModel] * data.weight(idxData)),idxC);
		    
	        c.setCov(c.getCov(idxC) + ((d.getMean(idxC) * d.getMean(idxC)) * apProba[idxModel] * data.weight(idxData)), idxC);
	    }
	    
	  apProbaTot[idxModel] += apProba[idxModel] * data.weight(idxData);
      }
  }
  
  for (unsigned long idxModel = 0; idxModel < modelNbComp; idxModel++)
  {
      DistribGD & c = clientMixture.getDistrib(idxModel);
	  
      for (unsigned long idxC = 0; idxC < vectSize; idxC++)
	{
	    c.setMean(c.getMean(idxC) / apProbaTot[idxModel], idxC);
	    //        c.setCov(c.getCov(idxC),idxC);
	}
  }
  
  for (unsigned long idxModel = 0; idxModel < modelNbComp; idxModel++)
  {
      DistribGD & c = clientMixture.getDistrib(idxModel);
      c.computeAll();
  }


}


// adapt model function to compute MAP from a given EM estimate
void adaptModelMAP(Config & config, StatServer & ss, MixtureServer & ms,
  MixtureGD & aprioriModel, MixtureGD & clientMixture, MAPCfg & mapCfg,
  unsigned long &frameCount)
{
  
    for (unsigned long trainIt = 0; trainIt < mapCfg.getNbTrainIt(); trainIt++)
    {	
	// Begin the initial adaptation loop (with bagged frames)
	    
	// Bayesian Adaptation client=MAP(aprioriModel,client)
        computeMAP(ms, aprioriModel, clientMixture, frameCount, config);	
        
	if (mapCfg.getNormalizeModel())
	    normalizeMixture(clientMixture, mapCfg, config);	// Normalize/fit the model if needed
     
    }
 
    if (debug)
        cout << "adaptModel nb distrib:" << ms.getDistribCount() << "nb mixt:" << ms.getMixtureCount() << endl;
}

void adaptModelMAP(Config & config, StatServer & ss, MixtureServer & ms,
  MixtureGD & aprioriModel, MixtureGD & clientMixture,
  unsigned long &frameCount)
{
  MAPCfg mapCfg(config);
  adaptModelMAP(config, ss, ms, aprioriModel, clientMixture, mapCfg,frameCount);
}




// Adapt model function to compute an EM estimate 

void adaptModelEM(Config & config, StatServer & ss, MixtureServer & ms,
  FeatureServer & fs, SegCluster & selectedSegments, MixtureGD & aprioriModel,
  MixtureGD & clientMixture, MAPCfg & mapCfg)
{
    unsigned long frameCount;	//A.P.
    if (verbose)
        cout << "Model adaptation based on true EM/ML estimate of training data" << endl;

    if (verboseLevel > 1)
    cout << "Mean LLK Init = " << meanLikelihood(ss, fs, clientMixture,
      selectedSegments, config) << endl;

    for (unsigned long trainIt = 0; trainIt < mapCfg.getNbTrainIt(); trainIt++)
    {	
	// Begin the initial adaptation loop (with bagged frames)
	    
	// Create a statistic accumulator using the curent model    
        MixtureStat & emAcc = ss.createAndStoreMixtureStat(clientMixture);	
        SegServer segServer;	// Create a local segment server 
	    
	// Create the cluster for describing the selected frames    
        SegCluster & baggedFramesCluster = segServer.createCluster(1, "", "");	
        baggedSegments(selectedSegments, baggedFramesCluster,mapCfg.getBaggedFrameProbability());
        emAcc.resetEM();
	    
	// Accumulate the EM statistics
        double llkPreviousIt = accumulateStatEM(ss, fs, emAcc, baggedFramesCluster, config);	
        clientMixture = emAcc.getEM();
	    
        // Get the EM estimate       
        frameCount = (unsigned long) emAcc.getEMFeatureCount();
        llkPreviousIt = llkPreviousIt / (double) frameCount;
	    
        if (verbose)
	    cout << "ML (partial) estimate it[" << trainIt <<"] (take care, it corresponds to the previous it,0 means init likelihood) = "
	  << llkPreviousIt << endl;
	
        ss.deleteMixtureStat(emAcc);
	
        if (verboseLevel > 2)
	    cout << "Likelihood on all frames= " << meanLikelihood(ss, fs, clientMixture, selectedSegments, config) << endl;

    }
    if (verboseLevel == 2)
        cout << "Final likelihood on all frames= " << meanLikelihood(ss, fs,clientMixture, selectedSegments, config) << endl;

}

void adaptModelEM(Config & config, StatServer & ss, MixtureServer & ms,
  FeatureServer & fs, SegCluster & selectedSegments, MixtureGD & aprioriModel,
  MixtureGD & clientMixture)
{
    MAPCfg mapCfg(config);
    adaptModelEM(config, ss, ms, fs, selectedSegments, aprioriModel,clientMixture, mapCfg);
}

// Adapt model function to compute an EM estimate from 2 feature server

void adaptModelEM(Config & config, StatServer & ss, MixtureServer & ms,
  FeatureServer & fs, SegCluster & selectedSegments, FeatureServer & fs2,
  SegCluster & selectedSegments2, MixtureGD & aprioriModel,
  MixtureGD & clientMixture, MAPCfg & mapCfg)
{
    unsigned long frameCount;	//A.P.
    if (verbose)
        cout << "Model adaptation based on true EM/ML estimate of training data" << endl;
    if (verboseLevel > 1)
        cout << "Mean LLK Init = " << meanLikelihood(ss, fs, clientMixture,selectedSegments, config) << endl;

    for (unsigned long trainIt = 0; trainIt < mapCfg.getNbTrainIt(); trainIt++)
    {	
	// Begin the initial adaptation loop (with bagged frames)
        double llkPreviousIt = 0;
	    
        // Create a statistic accumulator using the curent model
	MixtureStat & emAcc = ss.createAndStoreMixtureStat(clientMixture);
        SegServer segServer;	// Create a local segment server 
        SegServer segServer2;	// Create a local segment server 
	    
	// Create the cluster for describing the selected frames
        SegCluster & baggedFramesCluster = segServer.createCluster(1, "", "");	
	    
	// Create the cluster for describing the selected frames
        SegCluster & baggedFramesCluster2 = segServer2.createCluster(1, "", "");
	    
        baggedSegments(selectedSegments, baggedFramesCluster,mapCfg.getBaggedFrameProbability());
        baggedSegments(selectedSegments2, baggedFramesCluster2,	mapCfg.getBaggedFrameProbability());
	    
        emAcc.resetEM();
	// Accumulate the EM statistics on the 1st FS
        llkPreviousIt += accumulateStatEM(ss, fs, emAcc, baggedFramesCluster, config);	
	    
	// Accumulate the EM statistics on the 2nd FS
        llkPreviousIt += accumulateStatEM(ss, fs2, emAcc, baggedFramesCluster2, config);
	
	// Get the EM estimate 
        clientMixture = emAcc.getEM();
             
        frameCount = (unsigned long) emAcc.getEMFeatureCount();

        llkPreviousIt = llkPreviousIt / (double) frameCount;
        if (verbose)
	    cout << "ML (partial) estimate it[" << trainIt <<"] (take care, it corresponds to the previous it,0 means init likelihood) = "
	  << llkPreviousIt << endl;
	
        ss.deleteMixtureStat(emAcc);
	
        if (verboseLevel > 2)
	    cout << "Likelihood on all frames= " << meanLikelihood(ss, fs, clientMixture, selectedSegments, config) << endl;

    }
    if (verboseLevel == 2)
        cout << "Final likelihood on all frames= " << meanLikelihood(ss, fs,
      clientMixture, selectedSegments, config) << endl;

}


void adaptModelEM(Config & config, StatServer & ss, MixtureServer & ms,
  FeatureServer & fs, SegCluster & selectedSegments, FeatureServer & fs2,
  SegCluster & selectedSegments2, MixtureGD & aprioriModel,
  MixtureGD & clientMixture)
{
    MAPCfg mapCfg(config);
    adaptModelEM(config, ss, ms, fs, selectedSegments, fs2, selectedSegments2,
    aprioriModel, clientMixture, mapCfg);
}







double computeLLR(StatServer & ss, FeatureServer & fsTests, MixtureGD & world,
  MixtureGD & clientMixture, SegCluster & selectedSegmentsTests)
{
    ss.resetLLK(world);		// Reset the world LLK accumulator
    ss.resetLLK(clientMixture);	// ss.resetLLK(tabClientLine.getClientModel(i));                                   // Reset client LLK accumulator
    Seg *seg;			// reset the reader at the begin of the input stream
    selectedSegmentsTests.rewind();
	
    while ((seg = selectedSegmentsTests.getSeg()) != NULL)
    {	
	// For each of the selected segments
        unsigned long idxBeginFrame = seg->begin() + fsTests.getFirstFeatureIndexOfASource(seg->sourceName());
        fsTests.seekFeature(idxBeginFrame);
        Feature f;
	    
        for (unsigned long idxFrame = 0; idxFrame < seg->length(); idxFrame++)
	{			// For each frame of the segment
	    fsTests.readFeature(f);
	    // Determine the top components and compute wrld LLK
	    ss.computeAndAccumulateLLK(world, f, DETERMINE_TOP_DISTRIBS);	
	    ss.computeAndAccumulateLLK(clientMixture, f, USE_TOP_DISTRIBS);
	}
    }

    double LLKWorld  =  ss.getMeanLLK(world);	// Take the world LLK
    double LLKClient =  ss.getMeanLLK(clientMixture);	// Get the mean LLK 
    double LLRClient =  LLKClient - LLKWorld;	// Compute the LLR
    
    if ((verbose) && (verboseLevel == 2))
    {
        cout << "LLKWorld => " << LLKWorld << endl;
        cout << "LLKClient => " << LLKClient << endl;
        cout << "LLR => " << LLRClient << endl;
    }

    return LLRClient;
}



double computeLLR(Config & config, StatServer & ss, FeatureServer & fsTests,
  MixtureGD & world, MixtureGD & clientMixture,
  SegCluster & selectedSegmentsTests, String & idTest)
{
  cout << "compute LLR" << endl;
  FileInfo FI(idTest);		//create file to write top components info
  ss.resetLLK(world);		// Reset the world LLK accumulator
  ss.resetLLK(clientMixture);	// ss.resetLLK(tabClientLine.getClientModel(i));                                   
  Seg *seg;			// reset the reader at the begin of the input stream
  selectedSegmentsTests.rewind();
  RefVector <DoubleVector> stockLKVect(fsTests.getFeatureCount());
  while ((seg = selectedSegmentsTests.getSeg()) != NULL)
  {   
      // For each of the selected segments
      unsigned long idxBeginFrame =seg->begin() +fsTests.getFirstFeatureIndexOfASource(seg->sourceName());
      fsTests.seekFeature(idxBeginFrame);
      Feature f;
      for (unsigned long idxFrame = 0; idxFrame < seg->length(); idxFrame++)
      {			// For each frame of the segment
	  fsTests.readFeature(f);
	  // Determine the top components and compute wrld LLK
	  ss.computeAndAccumulateLLK(world, f, DETERMINE_TOP_DISTRIBS);	
	  
	  //STOCK the vector with the top selected component for reuse 
	  const LKVector & lkv = ss.getTopDistribIndexVector();
	      
	  //Stock all LKVectors for writing once in the file	
	  RealVector <real_t> & tmp= *new DoubleVector(config.getParam("topDistribsCount").toLong()+2,config.getParam("topDistribsCount").toLong()+2);
	  
	  for (unsigned long i=0;i<(unsigned long)config.getParam("topDistribsCount").toLong();i++){
	      tmp[i]=lkv.getArray()[i].idx;
	  }
	  
	  tmp[config.getParam("topDistribsCount").toLong()]=(lkv.sumNonTopDistribLK);
	  tmp[config.getParam("topDistribsCount").toLong()+1]=(lkv.sumNonTopDistribWeights);
	  stockLKVect.addObject(tmp);
	  ss.computeAndAccumulateLLK(clientMixture, f, USE_TOP_DISTRIBS);
     }
  }
  for (unsigned long i=0;i<stockLKVect.size();i++){			//write top info 
      FI.writeTopInfo(stockLKVect.getObject(i), config);
  }
  
  FI.close();
  stockLKVect.deleteAllObjects();
  
  double LLKWorld  =  ss.getMeanLLK(world);	// Take the world LLK
  double LLKClient =  ss.getMeanLLK(clientMixture);	// Get the mean LLK 
  double LLRClient =  LLKClient - LLKWorld;	// Compute the LLR
  
  if ((verbose) && (verboseLevel == 2))
  {
      cout << "LLKWorld => " << LLKWorld << endl;
      cout << "LLKClient => " << LLKClient << endl;
      cout << "LLR => " << LLRClient << endl;
  }

  return LLRClient;
}





double computeFastLLR(StatServer & ss, FeatureServer & fsTests,
  MixtureGD & world, MixtureGD & clientMixture,
  SegCluster & selectedSegmentsTests, String & idTest, Config & config)
{
  if (debug) cout << "FAST LLR COMPUTATION" << endl;
  FileInfo FI(idTest);		//load file to read top components info
  ss.resetLLK(world);		// Reset the world LLK accumulator
  ss.resetLLK(clientMixture);	// ss.resetLLK(tabClientLine.getClientModel(i));                                   // Reset client LLK accumulator
  Seg *seg;			// reset the reader at the begin of the input stream
  selectedSegmentsTests.rewind();
  unsigned long id = 0;		//id of the frame
  while ((seg = selectedSegmentsTests.getSeg()) != NULL)	// For each of the selected segments
    {
      unsigned long idxBeginFrame =seg->begin() +fsTests.getFirstFeatureIndexOfASource(seg->sourceName());
      fsTests.seekFeature(idxBeginFrame);
      Feature f;
      for (unsigned long idxFrame = 0; idxFrame < seg->length(); idxFrame++)	// For each frame of the segment
	{
	  fsTests.readFeature(f);
	  FI.loadTopInfo(ss, id, config);
	  ss.computeAndAccumulateLLK(world, f, USE_TOP_DISTRIBS);	// uses the top components and compute wrld LLK
	  id++;			
	  ss.computeAndAccumulateLLK(clientMixture, f, USE_TOP_DISTRIBS);
	}
    }
  FI.close();
  double LLKWorld = ss.getMeanLLK(world);	// Take the world LLK
  double LLKClient = ss.getMeanLLK(clientMixture);	// Get the mean LLK 
  double LLRClient = LLKClient - LLKWorld;	// Compute the LLR
  if ((verbose) && (verboseLevel == 2))
    {
      cout << "LLKWorld => " << LLKWorld << endl;
      cout << "LLKClient => " << LLKClient << endl;
      cout << "LLR => " << LLRClient << endl;
    }

  return LLRClient;
}


double computeLLRGD(Config & config, MixtureGD & clientMixture,
  MixtureGD & world, MixtureGD & dataTest)
{
  unsigned long topDistribs = 128;	//config.getParam ("topDistribsCount");
  TabWeight TabWeightData(dataTest, (unsigned long) 512);
  TabWeight TabWeightClient(clientMixture, topDistribs);
	
  double LLKClient = likelihoodGD(dataTest, clientMixture, TabWeightData, TabWeightClient);
  double LLKWorld  = likelihoodGD(dataTest, world, TabWeightClient, TabWeightClient);	//LLR between models 
  double LLRClient = LLKClient - LLKWorld;	// Compute the LLR
  if ((verbose) && (verboseLevel == 2))
  {
      cout << "LLKWorld => " << LLKWorld << endl;
      cout << "LLKClient => " << LLKClient << endl;
      cout << "LLR => " << LLRClient << endl;
  }

  return LLRClient;
}



void expandLLR(DoubleVector & decision, Config & configTest)
{
  cout <<
    " Choice of logistic regression to compute Feature Server Probabilities"
    << endl;
  double THETA = configTest.getParam("THETA").toDouble();
  double BETA = configTest.getParam("BETA").toDouble();
  for (unsigned long e = 0; e < decision.size(); e++)
  {
      cout << " LLR num[" << e << "] = " << decision[e] << endl;
      decision[e] =
	exp(THETA + BETA * decision[e]) / (1 + exp(THETA +
	  BETA * decision[e]));
      cout << " After logistic regression [" << e << "] = " << decision[e] <<
	endl;
    }

}









void WMAP(DoubleVector & decision, Config & configTest)
{
  cout << " Choice of WMAP to compute Feature Server Probabilities" << endl;
  double den = 0.0, num = 0.0;
  double pi = 3.14159;
  //From TARGET and IMPOSTOR score distributions

  double SIGMAclient = configTest.getParam("SIGMAtarget").toDouble();
  double SIGMAimp    = configTest.getParam("SIGMAimp").toDouble();
  double MUclient    = configTest.getParam("MUtarget").toDouble();
  double MUimp       = configTest.getParam("MUimp").toDouble();
  double poidsImp    = configTest.getParam("IMPweight").toDouble();
  double poidsTar    = configTest.getParam("TARweight").toDouble();
  double seuilMin    = configTest.getParam("thrMin").toDouble();
  double seuilMax    = configTest.getParam("thrMax").toDouble();
	
  for (unsigned long e = 0; e < decision.size(); e++)
  {
      if (e == 0)
	decision[e] = 1;	//Fixed proba for the train data = 1
      else
      {
	  cout << " LLR num[" << e << "] = " << decision[e] << endl;
	  //calcul de la proba de llr de x sachant que l'accès est client. P(client)=0.1 et calcul de la proba de llr de x sachant que l'accès est imposteur. P(imp)=0.9
	  if (decision[e] < seuilMin)
	      decision[e] = poidsTar;	//set a priori proba for LLR < -0.5 PB to Fix.
	  
	  else if (decision[e] > seuilMax)
	      decision[e] = 1;	//For train data pb to Fix 
	  
	  else
	  {
	      num =
		(1 / (SIGMAclient * sqrt(2 * pi))) * exp(-0.5 *
		((decision[e] - MUclient) / SIGMAclient) * ((decision[e] -
		    MUclient) / SIGMAclient)) * poidsTar;
	      num *= 100000;	//for precision
	      den =
		(1 / (SIGMAimp * sqrt(2 * pi))) * exp(-0.5 * ((decision[e] -
		    MUimp) / SIGMAimp) * ((decision[e] -
		    MUimp) / SIGMAimp)) * poidsImp;
	      den *= 100000;	//for precision
	      num = 1 + floor(num);
	      den = 1 + floor(den);
	      decision[e] = num / (num + den);
	  }
	  cout << " After WMAP [" << e << "] = " << decision[e] << endl;
      }
  }
}

String getFullFileName(String & id, Config & c)
{
  String ext = c.getParam("InfoExtension");
  String path = c.getParam("InfoPath");
  String fullFileName = path + id + ext;
  return fullFileName;
}
String getFullMixtureName(String & id, Config & c)
{
  String ext = c.getParam("loadMixtureFileExtension");
  String path = c.getParam("mixtureFilesPath");
  String fullFileName = path + id + ext;
  return fullFileName;
}


bool FileExists(String & fullFileName)
{
  bool result;
  ifstream outFile(fullFileName.c_str(), ios::binary);
  if (outFile.is_open())
    result= true;
  else
    result= false;
  outFile.close();
  return result;

}



// WMAPGMM with Fixed priors by config : TARweight and IMPweight

void WMAPGMMFixedPriors(DoubleVector & decision, Config & configTest,
  MixtureGD & tar, MixtureGD & non, StatServer & ss)
{
  cout << " Choice of WMAP GMM to compute Feature Server Probabilities" <<
    endl;
  double poidsImp = configTest.getParam("IMPweight").toDouble();
  double poidsTar = configTest.getParam("TARweight").toDouble();
  Feature f(1);
  double llkTar = 0.0;
  double llkNon = 0.0;
  for (unsigned long e = 0; e < decision.size(); e++)
  {
     /* if (e == 0)
	decision[0] = 1;	//the first LLR is the train data on the target model, we give WMAP=1. 
      else
	{*/
      f[0] = decision[e];	// to compute the LK between the score and the GMM learnt on scores we set a feature = value of the score. 
      if (e == 0)
          decision[0] = 1;	//train data : weight=1  
      else
      {
          llkTar = (ss.computeLLK(tar, f));
          llkNon = (ss.computeLLK(non, f));
          if (llkNon < configTest.getParam("LLKthreshold").toLong())
          {
              if (debug)
                  cout << "Flooring LLK at " << configTest.
              getParam("LLKthreshold").toLong() << endl;
              llkNon = configTest.getParam("LLKthreshold").toLong();
          }
          if (llkTar < configTest.getParam("LLKthreshold").toLong())
          {
              if (debug)
                  cout << "Flooring LLK at " << configTest.getParam("LLKthreshold").toLong() << endl;
                  llkTar = configTest.getParam("LLKthreshold").toLong();
	  }
	  
          llkTar = exp(llkTar);
          llkNon = exp(llkNon);
	  decision[e] =llkTar * poidsTar / (llkTar * poidsTar + llkNon * poidsImp);

       }
       if (verbose && verboseLevel >0)
           cout << "LLR : " << f[0] << " , After WMAP GMM [" << e << "] = "<< decision[e] << endl;



  }
}


//WMAP GMM with adaptative priors

void WMAPGMM(DoubleVector & decision, Config & configTest, MixtureGD & tar,
  MixtureGD & non, StatServer & ss)
{
  cout << " Choice of WMAP GMM to compute Feature Server Probabilities" <<endl;

  DoubleVector priorImp(decision.size(), decision.size());
  DoubleVector priorTar(decision.size(), decision.size());
  computePriors(decision, priorImp, priorTar, configTest);

  Feature f(1);
  double llkTar = 0.0;
  double llkNon = 0.0;
  for (unsigned long e = 0; e < decision.size(); e++)
  {
      /*if(e==0) decision[0]=1;               //the first LLR is the train data on the target model, we give WMAP=1. 
         else{        */
      f[0] = decision[e];	// to compute the LK between the score and the GMM learnt on scores we set a feature = value of the score. 
	  
      llkTar = (ss.computeLLK(tar, f));
      llkNon = (ss.computeLLK(non, f));
	  
      if (llkNon < configTest.getParam("LLKthreshold").toLong())
      {
	  if (debug)
	      cout << "Flooring LLK at " << configTest.getParam("LLKthreshold").
	      toLong() << endl;
	  llkNon = configTest.getParam("LLKthreshold").toLong();
      }
      if (llkTar < configTest.getParam("LLKthreshold").toLong())
      {
	  if (debug)
	      cout << "Flooring LLK at " << configTest.getParam("LLKthreshold").
	      toLong() << endl;
	  llkTar = configTest.getParam("LLKthreshold").toLong();
      }
      
      llkTar = exp(llkTar);
      llkNon = exp(llkNon);
      
      decision[e] = llkTar * priorTar[e] / (llkTar * priorTar[e] + llkNon * priorImp[e]);

      if (verbose && verboseLevel >0)
	cout << "LLR : " << f[0] << " , After WMAP GMM [" << e << "] = " <<
	  decision[e] << endl;
          //}

  }
}



void computePriors(DoubleVector & decision, DoubleVector & priorImp,
  DoubleVector & priorTar, Config & configTest)
{
  //the prior computation include the current trial
  double initPriorTar = configTest.getParam("initPriorTar").toDouble();
  double initPriorImp = configTest.getParam("initPriorImp").toDouble();
  priorTar.setAllValues(initPriorTar);
  priorImp.setAllValues(initPriorImp);
  double optiScore = configTest.getParam("OptimalScore").toDouble();
  for (unsigned long e = 1; e < decision.size(); e++)
    {
      if (decision[e] > optiScore)
	{			//target trial
	  initPriorTar = initPriorTar + 1;
	  //cout <<"prioTAR "<< initPriorTar<<endl;
	}
      else
	{			//Imp trial
	  initPriorImp = initPriorImp + 1;
	  // cout <<"prioIMP "<< initPriorImp<<endl;
	}

      priorTar[e] = initPriorTar / (initPriorTar + initPriorImp);
      priorImp[e] = 1 - priorTar[e];
      if (verbose && verboseLevel > 2)
	cout << "Priors after test : " << e << ", target = " << priorTar[e] <<
	  " impostors = " << priorImp[e] << endl;


    }



}





//fusing EM models and do a MAP  : for faster execution

void computeMAPmodelFromEMones(Config & config, StatServer & ss,
  MixtureServer & ms, DoubleVector & nbFramesSelected,
  MixtureGD & aprioriModel, MixtureGD & clientMixture, MixtureGD & aux,
  MixtureGD & tmp, DoubleVector & decision, XLine & testsToCompute)
{
  if (verbose && verboseLevel > 1)
    cout << "Fusing..." << endl;
  int i = 0;
  int nbModels = decision.size();	//number of EM models
  int index = ms.getMixtureIndex(testsToCompute.getElement(0)); //get the first model , its name is stocked in the XLine
  if (index == -1) aux =ms.loadMixtureGD(testsToCompute.getElement(0));
  else aux = ms.getMixtureGD(ms.getMixtureIndex(testsToCompute.getElement(0)));	  //IF STOCKED (ALREADY COMPUTED) LOAD MODEL

  unsigned long auxNbFrame = (unsigned long) (nbFramesSelected[0] * decision[0]);	//the elements in the ObjectRefvector follow the XLine order 
  //i.e. the objectRefvector and the XLine are reseted after each client                  

  for (i = 0; i < nbModels - 1; i++)
  {
      if (debug)
      {
	  cout << "Fusing models : " << testsToCompute. getElement(i)<<"  NbFrames = " <<auxNbFrame << " with " << testsToCompute.getElement(i +
	  1) <<" NbFrames = "<< (unsigned long) (nbFramesSelected[i + 1] * decision[i + 1])<< endl;
	     
      }
     
      index=ms.getMixtureIndex(testsToCompute.getElement(i + 1));
      if(index==-1) fuseModels(aux, auxNbFrame,ms.loadMixtureGD(testsToCompute.getElement(i + 1)),(unsigned long) (nbFramesSelected[i + 1] * decision[i + 1]), tmp);
      else fuseModels(aux, auxNbFrame,ms.getMixtureGD(ms.getMixtureIndex(testsToCompute.getElement(i + 1))),(unsigned long) (nbFramesSelected[i + 1] * decision[i + 1]), tmp);
      aux = tmp;
      auxNbFrame +=(unsigned long) (nbFramesSelected[i + 1] * decision[i + 1]);
      
      /*CA PLANTE, PKOI???
      ms.deleteMixtures(ms.getMixtureIndex(testsToCompute.getElement(i + 1)),ms.getMixtureIndex(testsToCompute.getElement(i + 1)));
      ms.deleteUnusedDistribs();*/
      
    }
    
   /* double regFactorAdapted=(auxNbFrame*config.getParam("MAPRegFactorMean").toDouble())/(nbFramesSelected[0] * decision[0]);
    if (verbose && verboseLevel >1) cout <<"New REG Factor = "<<regFactorAdapted<<endl;
    Config configMAP(config);
    char  regValue[256];
    sprintf(regValue,"%f",regFactorAdapted);
    configMAP.setParam("MAPRegFactorMean",regValue);
    adaptModelMAP(configMAP, ss, ms, aprioriModel, tmp, auxNbFrame);*/
    
    
    adaptModelMAP(config, ss, ms, aprioriModel, aux, auxNbFrame);
    clientMixture = aux;

}


//retrun the nb of frames in a cluster
double SegClusterFrame(SegCluster & SegC)
{
  double total = 0;
  SegC.rewind();
  Seg *p;
  while ((p = SegC.getSeg()) != NULL)
    total += p->length();
  return total;
}




class Norm:public Object
{
  public:String idTest;
  double mu;
  double sigma;
  String getClassName() const
  {
    return "Norm";
  };

};

//load impostors scores for TNORM computation from imp_seg.res file
//Must concatenate imp_seg.res and imp_imp.res if you want to do ZTNORM.

void loadTnormParam(String & inputTestListFileName, String & testFileTnorm,
  ObjectRefVector & stockTnorm, Config & config)
{

  int fieldId = 3;
  int fieldScore = 4;

  double nbTests = 0;
  //load imp_seg_male.res


  String idTest;

  XList testsTnorm(testFileTnorm, config);	// read the Id for each test
  XList clientsTnorm(inputTestListFileName, config);	//NDX target_seg
  XLine *line, *linep;

  while ((line = clientsTnorm.getLine()))
    {
      idTest = line->getElement(0);	//id client
      //norm tnorm=new norm();

      double accumScore = 0.0;
      double accumScore_2 = 0.0;
      double add = 0.0;
      while ((linep = testsTnorm.getLine()))
	{
	  if (linep->getElement(fieldId) == idTest)
	    {
	      add = (linep->getElement(fieldScore)).toDouble();
	      accumScore = accumScore + add;
	      accumScore_2 = accumScore_2 + (add * add);
	      nbTests++;
	    }

	}
      Norm & tnorm = *new Norm();	//don't forget to delete! use refVector.deleteAllObjects()
      testsTnorm.rewind();
      tnorm.idTest = (idTest);
      tnorm.mu = ((accumScore / nbTests));
      tnorm.sigma = sqrt((accumScore_2 / (nbTests) - (tnorm.mu * tnorm.mu)));
      stockTnorm.addObject(tnorm);
      if (debug)
	cout << "Test : " << tnorm.
	  idTest << " , TNORM Param : mu = " << tnorm.
	  mu << " sigma = " << tnorm.sigma << endl;
      nbTests = 0;

    }
}


//Do normalization on a score, look for test name in refVector to load mu and std variables.
void normalizeScore(String & test, double &decision, ObjectRefVector & stockNorm)
{

  if (verbose && verboseLevel > 2)
    cout << " Applying NORM, score =  " << decision << endl;

  for (unsigned long i = 0; i < stockNorm.size(); i++)
    {
      if (static_cast < Norm & >(stockNorm.getObject(i)).idTest == test)
	decision =
	  (decision - static_cast <
	  Norm & >(stockNorm.getObject(i)).mu) / static_cast <
	  Norm & >(stockNorm.getObject(i)).sigma;
    }
  if (verbose && verboseLevel > 1)
    cout << " Normed score =  " << decision << endl;

}

//Do normalization on a score, look for test name in refVector to load mu and std variables.
void normalizeScore(String & test, double &decision, ObjectRefVector & stockNorm,double &shift)
{

  if (verbose && verboseLevel > 2)
    cout << " Applying NORM, score =  " << decision << endl;

  for (unsigned long i = 0; i < stockNorm.size(); i++)
    {
      if (static_cast < Norm & >(stockNorm.getObject(i)).idTest == test)
	decision =
	  (decision - (static_cast < Norm & >(stockNorm.getObject(i)).mu + shift)) / static_cast <
	  Norm & >(stockNorm.getObject(i)).sigma;
    }
  if (verbose && verboseLevel > 1)
    cout << " Normed score =  " << decision << endl;

}


void computeAndStoreZnormParam(StatServer &ss, String & inputImpListFileName, String &idclient, MixtureGD &clientMixture,
  ObjectRefVector & stockZnorm, MixtureGD &world, Config & config, bool &ztnorm, ObjectRefVector & stockTnorm ){
	  
    // read the Id for each impostor test , format : one column with impostor tests names 
    XList testsZnorm(inputImpListFileName, config);
    XLine *line;
    String fullFileName;
    String *idImp = NULL;	  
    double LLRRatio = 0.0, nbTests = 0.0, accumScore = 0.0 , accumScore_2 = 0.0;
	  
	  
    //Loop on each impostor tests
    while ((line = testsZnorm.getLine()))
    {
	line->getElement();   //impostor name LOST because not used
	    
	idImp = line->getElement();
        fullFileName = getFullFileName(*idImp, config);
	
        if (verbose && verboseLevel > 2)
            cout << " IMPOSTOR SEGMENT TESTED   [" << *idImp << "]" << endl;
        String labelSelectedFrames =  config.getParam("labelSelectedFrames");
	
        //IMPOSTOR DATA STUFF
        FeatureServer fs(config, *idImp);
	
	// Create the segment server for managing the segments/clusters
        SegServer segmentsServer;
	
        // Create the lable server, for indexing the segments/clusters	
        LabelServer labelServer;
	
	// Reading the segmentation files for each feature input file
        initializeClusters(*line, segmentsServer, labelServer, config);
	
	// Verify if the segments ending before the end of the feature files...
	verifyClusterFile(segmentsServer, fs, config);
	
	// Get the index of the cluster with in interest audio segments
        long codeSelectedFrame = labelServer.getLabelIndexByString(labelSelectedFrames);
	
        if (codeSelectedFrame == -1)
        {		// No data for this model !!!!!!!!!!!!!!
            cout << " WARNING - NO DATA FOR [" << *idImp<< "]";
        }
        else
        { 
	    SegCluster& selectedSegments = segmentsServer.getCluster(codeSelectedFrame);
            
	    //TEST IF TOP TEN INFO FILE ARE ALREADY COMPUTED AND STORED ON THE HARD DRIVE	
	    if (FileExists(fullFileName)  ==  false  )
	    {   
	        LLRRatio  =     computeLLR(config, ss, fs, world, clientMixture, selectedSegments, fullFileName);    
	    }
		
	    else
	    {
	        LLRRatio  =     computeFastLLR(ss, fs, world, clientMixture, selectedSegments, fullFileName, config);    
            } 
	    
	    //IF ZTNORM is used scores must be tnormed before ZNORM computation
	    if(ztnorm)
	    {	
		if (verbose && verboseLevel > 1) 
		    cout <<" ZTNorm activated : T-norm score before Z-norm"<<endl;
		normalizeScore(*idImp, LLRRatio, stockTnorm);
		
	    }
		    
	    accumScore = accumScore + LLRRatio;
            accumScore_2 = accumScore_2 + (LLRRatio * LLRRatio);
	    nbTests++;
        }
	
    }
    Norm & znorm = *new Norm();	//don't forget to delete! use refVector.deleteAllObjects()
    testsZnorm.rewind();
    znorm.idTest = idclient;
    znorm.mu = ((accumScore / nbTests));
    znorm.sigma = sqrt((accumScore_2 / (nbTests) - (znorm.mu * znorm.mu)));
    stockZnorm.addObject(znorm);
    if (verbose && verboseLevel > 1)
        cout << "Test : " << znorm.idTest << " , TNORM Param : mu = " << znorm.mu << " sigma = " << znorm.sigma << endl;
	
}

//Reset the vector which contains LLR or WMAP weights ( do not reset the first value : LLR(train model/train data).
void resetWeights(DoubleVector & decision)
{

  //set all weights to 0 : TEST should be equal to baseline
  cout << "Reset all WMAP weights" << endl;
  for (unsigned long e = 1; e < decision.size(); e++)
    {
      decision[e] = 0;		//FOR DEBUGGING
    }
}

 //look for a true target trials in the targetTests file, if ok set score to 1 or WMAP else set to 0.
 // DO NOT SET wmap, regress or wmapgmm to true when using Oracle.

void Oracle(String & idTar, String & idTest, double &score, Config & config,
  MixtureGD & tar, MixtureGD & non, StatServer & ss)
{

  bool wmap = config.getParam("wmapOracleType").toBool();
  bool one = config.getParam("classicalOracleType").toBool();;
  String model, test;
  XLine *line;
  bool find = false;
  String targetTestsList = config.getParam("targetTests");
  XList targetTests(targetTestsList, config);

  while ((line = targetTests.getLine()))
    {
      model = line->getElement(0);	//id client in file
      test = line->getElement(2);	//id test in file

      if (idTar == model && test == idTest)	//test = target test
	{
	  if (wmap)
	    {
	      cout << "Find true target trial, set weight to WMAP computation"
		<< endl;
	      //HERE : leave the WMAP weight
	      DoubleVector tmp;
	      tmp.addValue(0);	//The first element is not evaluated by WMAP.
	      tmp.addValue(score);
	      WMAPGMMFixedPriors(tmp, config, tar, non, ss);
	      score = tmp[1];
	    }
	  else if (one)
	    {
	      cout << "Find true target trial, set weight to 1" << endl;
	      score = 1;	//true Oracle with weight = 1  
	    }


	  find = true;
	}



    }

  if (!find)
    {
      cout << "Impostor trial, set weight to 0" << endl;
      score = 0;		//test != test target, set proba to 0;  
    }



}

//Jack knife on data, WARNING ONLY ONE EM AND MAP ITERATION IS DONE
void crossValid(Config & configTest, StatServer & ss,  MixtureServer & ms, FeatureServer & fs, SegCluster & selectedSegments,
  MixtureGD & aprioriModel,MixtureGD &bestModel,  SegCluster& selectedSegmentsBagged, String & idTest)
{
  MixtureGD & clientMixture = ms.duplicateMixture(aprioriModel, DUPL_DISTRIB);
  MixtureGD & clientMixtureEM = ms.duplicateMixture(aprioriModel, DUPL_DISTRIB);
  double LLR = 0.0,previousLLR=100000;
  MAPCfg mapCfg(configTest);
  double trainSelected = configTest.getParam("SelectedTrain").toDouble();
  cout << "Compute Model on " << (trainSelected *  100) << "% of data" << endl;
 unsigned long baggedMinimalLength =3;
  if (configTest.existsParam("baggedMinimalLength"))
    baggedMinimalLength = configTest.getParam("baggedMinimalLength").toLong();
    unsigned long baggedMaximalLength =7;
  if (configTest.existsParam("baggedMaximalLength"))
    baggedMaximalLength = configTest.getParam("baggedMaximalLength").toLong();
	
  for (int It = 0; It < configTest.getParam("AverageIt").toLong(); It++)			//compute AverageIt times
    {	
	  if (debug) cout << "Iteration "<<It<<endl;
	  if (verboseLevel > 1)
	    cout << "Mean LLK Init = " << meanLikelihood(ss, fs, clientMixture,   selectedSegments, configTest) << endl;
	     MixtureStat & emAcc = ss.createAndStoreMixtureStat(clientMixture);	// Create a statistic accumulator using the curent model
	      SegCluster & baggedSelected =  (selectedSegments.getServer()).createCluster(1, "", "");	// Create the cluster for describing the selected frames
	      SegCluster & baggedUnselected = (selectedSegments.getServer()).createCluster(1, "", "");	// Create the cluster for describing the selected frames
	      baggedSegments(selectedSegments, baggedSelected, baggedUnselected,trainSelected,baggedMinimalLength,baggedMaximalLength);	//Train the client model on trainSelected Frame
	      emAcc.resetEM();
	      double llkPreviousIt = accumulateStatEM(ss, fs, emAcc, baggedSelected, configTest);	// Accumulate the EM statistics
	      clientMixture = emAcc.getEM();	// Get the EM estimate 
	      clientMixtureEM=clientMixture;  
	      unsigned long frameCount = (unsigned long) emAcc.getEMFeatureCount();
	      llkPreviousIt = llkPreviousIt / (double) frameCount;
	      if (verbose) cout << "ML (partial) estimate it[0] (take care, it corresponds to the previous it,0 means init likelihood) = " << llkPreviousIt << endl;
	      computeMAP(ms, aprioriModel, clientMixture, frameCount, configTest);	// Bayesian Adaptation client=MAP(aprioriModel,client)
	      if (mapCfg.getNormalizeModel())
		normalizeMixture(clientMixture, mapCfg, configTest);	// Normalize/fit the model if needed
	      ss.deleteMixtureStat(emAcc);
	      if (verboseLevel > 2)
		cout << "Likelihood on all frames= " << meanLikelihood(ss, fs, clientMixture, selectedSegments, configTest) << endl;
		 if (verboseLevel > 1)
		    cout << "Final likelihood on all frames= " << meanLikelihood(ss, fs, clientMixture, selectedSegments, configTest) << endl;
		  if(debug)  cout << "Compute LLR of  " << ((1-trainSelected) * 100) << "% of train data on the model learnt on " << (trainSelected* 100) << "%  data" << endl;
		  LLR = computeLLR(ss, fs, aprioriModel, clientMixture, baggedUnselected); //compute LLR with unselected data
		 if (debug) cout << "LLR for IT "<<It<<" = "<< LLR<<endl;
		 if(LLR < previousLLR) {
			  bestModel=clientMixtureEM;	// Get the EM estimate 
			//recopy cluster
			  copyCluster(baggedSelected,selectedSegmentsBagged);
		 }
		  (selectedSegments.getServer()).remove(baggedSelected);
		  (selectedSegments.getServer()).remove(baggedUnselected);
		  previousLLR=LLR;
		  clientMixture=aprioriModel;	
		
				
    }
  
    bestModel.save(idTest,configTest);
  
    //delete temporary mixture
    
    ms.deleteMixture(clientMixtureEM);
    ms.deleteMixture(clientMixture);
    ms.deleteUnusedDistribs();

}


//2 purposes : do not recalculate LLR of test data on a client model, and possibility to use RES file from a different recognition system to compute adaptation weights
double searchLLRFromResFile(String & idTar, String & test,
  String & inputResFilename, Config & config)
{
  double LLR = 0.0;
  cout << "Search from file : " << inputResFilename << " for the LLR : " <<
    idTar << " / " << test << endl;
  String idTest, idClient;
  XList listLLR(inputResFilename, config);
  XLine *line;
  while ((line = listLLR.getLine()))
    {

      idClient = line->getElement(1, 0);
      idTest = line->getElement(3, 0);

      if ((idClient == idTar) && (idTest == test))
	LLR = line->getElement(4, 0).toDouble();

    }
  if (LLR != 0.0)
    return LLR;
  else
    {
      cout << "LLR not found in :" << inputResFilename << " , return -1" <<
	endl;
      return -1;
    }

}






unsigned long adaptModelEMweightedFrames(String &labelSelectedFrames,XLine & featureFileName,StatServer &ss,MixtureGD &world, MixtureGD &tar, MixtureGD &non,MixtureGD &MixtureforLLR,MixtureGD &MixtureEMOutput, Config &config, FeatureServer &fs, String &fullFileName,String &idTest){
	unsigned long frameCount=0;
	WindowLLR windowLLR(config); // Initialize the windowLLR mode if requested
	if(verbose && verboseLevel >1) windowLLR.showConfig();
	if (windowLLR.isSet()) windowLLR.setNbClient(1); // set one to nb client for window LLR ???A ENLEVER????
	FileInfo FI(fullFileName);		//file to read or write top components info
	bool fileispresent=true;
	if (FileExists(fullFileName) == false ) fileispresent=false;
	DoubleVector FrameWeights(0,0); 
	MixtureGDStat &worldAcc=ss.createAndStoreMixtureGDStat(world);
	worldAcc.resetLLK();               
	MixtureGDStat &clientAccforLLR=ss.createAndStoreMixtureGDStat(MixtureforLLR);
	MixtureGDStat &clientAcc=ss.createAndStoreMixtureGDStat(MixtureEMOutput);
	clientAccforLLR.resetLLK(); 	
	SegServer segmentsServer;                                                             // Create the segment server for managing the segments/clusters
        LabelServer labelServer;                                                              // Create the lable server, for indexing the segments/clusters
        initializeClusters(featureFileName,segmentsServer,labelServer,config);                // Reading the segmentation files for each feature input file
        verifyClusterFile(segmentsServer,fs,config);                                          // Verify if the segments ending before the end of the feature files...
        long codeSelectedFrame=labelServer.getLabelIndexByString(labelSelectedFrames);        // Get the index of the cluster with in interest audio segments
        if (codeSelectedFrame==-1)                                                            // The file is empty
		cout << "ATTENTION, TEST FILE ["<<idTest<<"] is empty"<<endl;
        else{
		SegCluster& selectedSegments=segmentsServer.getCluster(codeSelectedFrame); // Gives the cluster of the selected/used segments   
		//compute LLR for each feature server test frames
		Seg *seg;			// reset the reader at the begin of the input stream
		selectedSegments.rewind();
		unsigned long idFrame=0;
		while ((seg = selectedSegments.getSeg()) != NULL)
			{				// For each of the selected segments
			unsigned long idxBeginFrame =seg->begin() +fs.getFirstFeatureIndexOfASource(seg->sourceName());
			fs.seekFeature(idxBeginFrame);
			Feature f;
			for (unsigned long idxFrame = 0; idxFrame < seg->length(); idxFrame++)
			{			// For each frame of the segment
				double llkw=0.0;
				double llkc=0.0;
				fs.readFeature(f);
				if(!fileispresent){
					llkw=worldAcc.computeAndAccumulateLLK(f,1.0,DETERMINE_TOP_DISTRIBS);    // Determine the top components and compute wrld LLK
					//STOCK the vector with the top selected component for reuse 
					const LKVector & lkv = ss.getTopDistribIndexVector();
					FI.writeTopInfo(lkv, config);
				}
				else {
					FI.loadTopInfo(ss, idFrame, config);
					llkw=worldAcc.computeAndAccumulateLLK(f,1.0,USE_TOP_DISTRIBS);       // Determine the top components and compute wrld LLK
					idFrame++;
				}
			   
				if (windowLLR.isSet()) windowLLR.dec(idxBeginFrame+idxFrame);                 
				llkc=clientAccforLLR.computeAndAccumulateLLK(f,1.0,USE_TOP_DISTRIBS);
				if (windowLLR.isSet()) windowLLR.accLLR(0,llkc-llkw);
				if (windowLLR.isSet() && windowLLR.isEnd()){				//fin window on adapte
					FrameWeights.addValue(windowLLR.getLLR(0));				//for each segments of window size the LLR is stocked
					
				}
				
			}
		}
		WMAPGMMFixedPriors(FrameWeights,config,tar,non,ss);								//compute WMAP for each LLR
		//Create EM model with a weight for each segment 
		clientAcc.resetEM();
		double llkAcc = 0.0;
		unsigned long cpt=0;
		selectedSegments.rewind();
		while ((seg = selectedSegments.getSeg()) != NULL)
			{				// For each of the selected segments
			unsigned long idxBeginFrame =seg->begin() +fs.getFirstFeatureIndexOfASource(seg->sourceName());
			fs.seekFeature(idxBeginFrame);
			Feature f;
			for (unsigned long idxFrame = 0; idxFrame < seg->length(); idxFrame++)
			{			// For each frame of the segment
				fs.readFeature(f);
				if (windowLLR.isSet()) windowLLR.dec(idxBeginFrame+idxFrame);  
				if (windowLLR.isSet() && windowLLR.isEnd())	cpt++;
				llkAcc += log(clientAcc.computeAndAccumulateEM(f, FrameWeights[cpt]))*FrameWeights[cpt];		//PB si un seg < 30 trames
						
			}
			
			    
	}
	MixtureEMOutput=clientAcc.getEM();
	frameCount = (unsigned long) clientAcc.getEMFeatureCount();
	if(verbose && verboseLevel >1) cout << "NBFRAMES  for Acc= "<< frameCount <<endl;
	if(verbose && verboseLevel >1) cout << "ML estimate Likelihood  = "<< (llkAcc/frameCount) <<endl;
			    
	}
	return frameCount;
}
	
	
	
String selectNearestTarModel(String &TARListFilename, String & fullFileName,Config &config, StatServer & ss, FeatureServer & fs,
  MixtureGD & world,  SegCluster & selectedSegments,MixtureServer &ms){
	DoubleVector LLR; 
	double LLRvalue= 0.0;
	XList TARXList(TARListFilename, config);	// read the Id + filenames for each client
	XLine *line;
	while ((line = TARXList.getLine()) != NULL){
		  String *idmodel = line->getElement();	// Get the TAR model ID 
		  MixtureGD & TARmodel = ms.loadMixtureGD(*idmodel);
		  if (FileExists(fullFileName) == false)
			    {	//first time seeing this test
			      LLRvalue = computeLLR(config, ss, fs, world, TARmodel, selectedSegments, fullFileName);	//For proba, take the decision on the original target model
			    }
		  else
			      LLRvalue = computeFastLLR(ss, fs, world, TARmodel, selectedSegments, fullFileName, config);	//For proba, take the decision on the original target model
		
	cout <<"Model : "<<*idmodel <<" LLR = " <<LLRvalue<<endl;
	LLR.addValue(LLRvalue);
	ms.deleteMixture(TARmodel);	
	ms.deleteUnusedDistribs();
	}
	
	unsigned long index=LLR.getIndexOfLargestValue();
	TARXList.rewind();
	return *(TARXList.getLine(index)).getElement();
	
  
		  
}
//Need all models in memory
//Not optimised with top ten info file
void computeLLRmatrix(DoubleMatrix &LLR, XLine &models, XList &features, Config &config, StatServer &ss, MixtureGD & world, MixtureServer &ms,String &labelSelectedFrames){
	String *idmodel;
	XLine *line;
	int j=0;
	int i=0;
	LLR.setDimensions(models.getElementCount(),models.getElementCount());
	double LLRvalue=0.0;
	while((idmodel=models.getElement())){			//loop on each target model
			
		while((line = features.getLine()) != NULL){		//loop on each target model
			
		      FeatureServer fs(config,*line);                                            // Reading the features (from several files)
		      SegServer segmentsServer;                                                             // Create the segment server for managing the segments/clusters
		      LabelServer labelServer;                                                              // Create the lable server, for indexing the segments/clusters
		      initializeClusters(*line,segmentsServer,labelServer,config);               // Reading the segmentation files for each feature input file
		      verifyClusterFile(segmentsServer,fs,config);                                          // Verify if the segments ending before the end of the feature files...
		      long codeSelectedFrame=labelServer.getLabelIndexByString(labelSelectedFrames);        // Get the index of the cluster with in interest audio segments
		      if (codeSelectedFrame==-1){                                                           // No data for this model !!!!!!!!!!!!!!
				cout << " WARNING - NO DATA FOR TRAINING ["<<*idmodel<<"]";
		      }
		      else{
				SegCluster& selectedSegments=segmentsServer.getCluster(codeSelectedFrame); // Gives the cluster of the selected/used segments                                   
				MixtureGD & TARmodel = ms.loadMixtureGD(*idmodel);
				if (debug) cout << "Model : "<<*idmodel<<" feature" << line->toString()<<endl;
				LLRvalue = computeLLR(ss, fs, world, TARmodel, selectedSegments);	
				LLR(i,j) = LLRvalue;
				j++;
				ms.deleteMixture(TARmodel);
				ms.deleteUnusedDistribs();
				 
			} 
                      		
			fs.reset();	
			labelServer.clear();
			segmentsServer.removeAllClusters();
			segmentsServer.removeAllSegs();
		}
		features.rewind();
		i++;
		j=0;
		
	}
	
	LLR.save("LLRmatrix.mat",config);
	
}


void addLineInXList(XLine &line,XList &list){
	unsigned long cpt=0;
	list.addLine();
	XLine *tmp=list.getLine();
	while(cpt!=line.getElementCount()){
		(*tmp).addElement(*line.getElement());
		cpt++;
	}
	list.rewind();
	line.rewind();
	
	
}



//weighted fusion of MAP models  

void fuseMAP(Config & config, StatServer & ss,  MixtureServer & ms, MixtureGD & aprioriModel, MixtureGD & clientMixture, MixtureGD & aux,
  MixtureGD & tmp, DoubleVector & decision, XLine & testsToCompute)
{
  if (verbose && verboseLevel > 1)
    cout << "MAP Fusion" << endl;
  int i = 0;
  int nbModels = decision.size();	//number of models
  int index = ms.getMixtureIndex(testsToCompute.getElement(0)); //get the first model , its name is stocked in the XLine
  if (index == -1) aux =ms.loadMixtureGD(testsToCompute.getElement(0));
  else aux = ms.getMixtureGD(ms.getMixtureIndex(testsToCompute.getElement(0)));	  //IF STOCKED (ALREADY COMPUTED) LOAD MODEL

  double auxWeight = decision[0];	//first weight, usually 1 for train data
  //i.e. the objectRefvector and the XLine are reseted after each client                  

  for (i = 0; i < nbModels - 1; i++)
    {
      if (debug){
	cout << "Fusing MAP models : " << testsToCompute.getElement(i)<<" with " << testsToCompute.getElement(i + 1) << endl;
	}

      index=ms.getMixtureIndex(testsToCompute.getElement(i + 1));
      if(index==-1) fuseMAPMeans(aux, auxWeight,ms.loadMixtureGD(testsToCompute.getElement(i + 1)),decision[i + 1], tmp);
      else fuseMAPMeans(aux, auxWeight,ms.getMixtureGD(ms.getMixtureIndex(testsToCompute.getElement(i + 1))),decision[i + 1], tmp);
      aux = tmp;
      auxWeight +=decision[i + 1];

    }
  clientMixture = tmp;

}



void fuseMAPMeans(const MixtureGD & model1,double &w1,const MixtureGD & model2,double &w2,MixtureGD &result){
  unsigned long vectSize=model1.getVectSize();  
  if (vectSize!=model2.getVectSize()) throw Exception("Feature vector size should be the same" , __FILE__, __LINE__);
  unsigned long nbDistrib=model1.getDistribCount();
  if (nbDistrib!=model2.getDistribCount()) throw Exception("Number of components should be the same" , __FILE__, __LINE__);
  for (unsigned long idx=0;idx<nbDistrib;idx++){
    DistribGD & d1=model1.getDistrib(idx);
    DistribGD & d2=model2.getDistrib(idx);
    DistribGD & dr=result.getDistrib(idx);

    for (unsigned long c=0;c<vectSize;c++){
	dr.setMean(((d1.getMean(c)*w1)+(d2.getMean(c)*w2))/(w1+w2),c);
	dr.setCov((d1.getCov(c)),c);// COPY of the model 1 covariance, should be equal to world (MAP mean only)
	    }
    dr.computeAll();
    result.weight(idx)=model1.weight(idx);// COPY of the model 1 weight, should be equal to world (MAP mean only)
   }
}


double computeWeightedLLR(Config & config, StatServer & ss, FeatureServer &fsTests, MixtureServer & ms, MixtureGD & aprioriModel, SegCluster & selectedSegmentsTests,MixtureGD & aux, DoubleVector & decision, XLine & testsToCompute, String &idTest, String &fullFileName){

if (verbose && verboseLevel > 1)
    cout << "Weighted SUM of LLR" << endl;
  int i = 0;
  int nbModels = decision.size();	//number of models
  double Wllr=0.0,weights=0.0,Wllrtot=0.0;
                  

  for (i = 0; i < nbModels - 1; i++)	//DON'T TAKE THE LAST MODEL because it is the current test model
    {
      if (debug){
	cout << "Weighted SUM of LLR model : " << testsToCompute.getElement(i)<<" with segment " << idTest << endl;
	}
      int index=ms.getMixtureIndex(testsToCompute.getElement(i));
      if(index==-1) aux=ms.loadMixtureGD(testsToCompute.getElement(i));
      else aux=ms.getMixtureGD(ms.getMixtureIndex(testsToCompute.getElement(i)));
      Wllr+=(computeFastLLR(ss, fsTests, aprioriModel, aux, selectedSegmentsTests, fullFileName, config)*decision[i]);
      weights+=decision[i];
	cout <<"LLR accumulated = "<<Wllr<< "sum weights = "<<weights;
    }
  Wllrtot=Wllr/weights;
    return Wllrtot;
	  
  }

  
void copyCluster(SegCluster &source, SegCluster &dest) {
	dest.removeAll(); // clean the destination
	source.rewind();
	Seg *p;
	while((p=source.getSeg())!=NULL)
		dest.add(*p);
}

int findClusterInRefvector(ObjectRefVector &SegServ,String & id){
	for (unsigned long i=0;i<SegServ.size();i++){
		if((static_cast <SegCluster & >(SegServ.getObject(i))).sourceName() == id)
			return i;
			
	}
	return -1;
}



void learnEMimpostorModels(Config & config, DoubleVector &NbFramesSelectedImp,MixtureServer &ms,StatServer &ss,MixtureGD &world){
XLine *line;
String ImpostorList = config.getParam("impCohortFile");
XList Imp(ImpostorList, config);
	
while ((line = Imp.getLine()))
    {
      	  MixtureGD & impMixtureEM = ms.duplicateMixture(world, DUPL_DISTRIB);
	  String  *idImp = line->getElement();	// Get the Tests ID (id)
	  String fullMixtureName=getFullMixtureName(*idImp,config);
	    
	  if (FileExists(fullMixtureName) == false) {
		  
	    	  XLine featureFileListImp = line->getElements();	// Get the list of feature file for the tests (end of the line)
		  if (verbose)
		    cout << "Train IMPOSTOR model  [" << *idImp << "]" << endl;
		  String labelSelectedFrames = config.getParam("labelSelectedFrames");
	
		   //IMPOSTOR DATA STUFF
		  FeatureServer fsImp(config, featureFileListImp);
		  SegServer segmentsServerImp;	// Create the segment server for managing the segments/clusters
		  LabelServer labelServerImp;	// Create the lable server, for indexing the segments/clusters
		  initializeClusters(featureFileListImp, segmentsServerImp, labelServerImp, config);	// Reading the segmentation files for each feature input file
		  verifyClusterFile(segmentsServerImp, fsImp, config);	// Verify if the segments ending before the end of the feature files...
		  long codeSelectedFrame = labelServerImp.getLabelIndexByString(labelSelectedFrames);	// Get the index of the cluster with in interest audio segments
		  if (codeSelectedFrame == -1)
		    {		// No data for this model !!!!!!!!!!!!!!
		      cout << " WARNING - NO DATA FOR TRAINING [" << *idImp
			<< "]";
		    }
		  else
		    { 
		        SegCluster& selectedSegmentsImp = segmentsServerImp.getCluster(codeSelectedFrame);
		        NbFramesSelectedImp.addValue(SegClusterFrame(selectedSegmentsImp));				//used for fuseModel
		        adaptModelEM(config, ss, ms, fsImp, selectedSegmentsImp, world, impMixtureEM);	//create the EM model for the impostor
		        ms.setMixtureId(impMixtureEM, *idImp);
		        impMixtureEM.save(*idImp,config);		//SAVE test mixture 
		        ms.deleteMixture(impMixtureEM);	
		        ms.deleteUnusedDistribs();
		    }
		    
	    }
		    
   }

}

//Learn all impostor models by classical MAP adaptation

void learnMAPimpostorModels(Config & config, DoubleVector &NbFramesSelectedImp,MixtureServer &ms,StatServer &ss,MixtureGD &world){
XLine *line;
String ImpostorList = config.getParam("impCohortFile");
XList Imp(ImpostorList, config);
	
while ((line = Imp.getLine()))
    {
      	  MixtureGD & impMixture = ms.duplicateMixture(world, DUPL_DISTRIB);
	  String  *idImp = line->getElement();	// Get the Tests ID (id)
	  String fullMixtureName=getFullMixtureName(*idImp,config);
	    
	   //IF IMP MODEL WAS NOT CREATED YET 
	  if (FileExists(fullMixtureName) == false) {
		  
	    	  XLine featureFileListImp = line->getElements();	// Get the list of feature file for the tests (end of the line)
		  if (verbose)
		    cout << "Train IMPOSTOR model  [" << *idImp << "]" << endl;
		  String labelSelectedFrames = config.getParam("labelSelectedFrames");
	
		   //IMPOSTOR DATA STUFF
		  
		  FeatureServer fsImp(config, featureFileListImp);
		  
		  // Create the segment server for managing the segments/clusters
		  SegServer segmentsServerImp;
		  
		  // Create the lable server, for indexing the segments/clusters
		  LabelServer labelServerImp;	
		  
		  // Reading the segmentation files for each feature input file
		  initializeClusters(featureFileListImp, segmentsServerImp, labelServerImp, config);
		  
		  // Verify if the segments ending before the end of the feature files...
		  verifyClusterFile(segmentsServerImp, fsImp, config);
		  
		  // Get the index of the cluster with in interest audio segments
		  long codeSelectedFrame = labelServerImp.getLabelIndexByString(labelSelectedFrames);
		  
		  if (codeSelectedFrame == -1)
		    {		// No data for this model !!!!!!!!!!!!!!
		      cout << " WARNING - NO DATA FOR TRAINING [" << *idImp
			<< "]";
		    }
		  else
		    { 
		        SegCluster& selectedSegmentsImp = segmentsServerImp.getCluster(codeSelectedFrame);
		        NbFramesSelectedImp.addValue(SegClusterFrame(selectedSegmentsImp));
			//create the EM model for the impostor
		        adaptModel(config, ss, ms, fsImp, selectedSegmentsImp, world, impMixture);	
		        ms.setMixtureId(impMixture, *idImp);
		        impMixture.save(*idImp,config);		//SAVE mixture 
		        ms.deleteMixture(impMixture);	
		        ms.deleteUnusedDistribs();
		    }
		    
	    }
		    
   }

}




double computeAdaptedTnorm(Config & configTest, MixtureServer & ms,StatServer &ss, MixtureGD &world, MixtureGD &auxMixture, MixtureGD &tmpMixture, double &nonorm_score, int &countTests, DoubleVector &decision, DoubleVector &NbFramesSelected, XLine &testsToCompute,
DoubleVector &NbFramesSelectedImp, FeatureServer &fsTests, SegCluster &selectedSegmentsTests, String &idTest){
	
//load modèle EM imposteur n-1 iteration
XLine *line;
String ImpostorList = configTest.getParam("ImpList");
XList Imp(ImpostorList, configTest);
//recopie nb trames
DoubleVector NbFramesSelectedForAdaptation=NbFramesSelected;
XLine testsToComputeForAdaptation;
double accumScore = 0.0;
double accumScore_2 = 0.0;
double add = 0.0;
int nbTests=0;
double mu=0.0, sigma=0.0, score=0.0;	
MixtureGD & impMixture = ms.duplicateMixture(world, DUPL_DISTRIB);

while ((line = Imp.getLine())){
	String *idImp = line->getElement();	// Get the Tests ID (id)    
	testsToComputeForAdaptation.addElement(*idImp);
	for (unsigned long e=1;e>testsToCompute.getElementCount();e++)
		testsToComputeForAdaptation.addElement(testsToCompute.getElement(e,false));
	
	NbFramesSelectedForAdaptation[0]=NbFramesSelectedImp[countTests-1];
	computeMAPmodelFromEMones(configTest, ss, ms, NbFramesSelectedForAdaptation, world, impMixture, auxMixture, tmpMixture, decision, testsToComputeForAdaptation);
	add= computeFastLLR(ss, fsTests, world, impMixture, selectedSegmentsTests,idTest, configTest);
	accumScore = accumScore + add;
	accumScore_2 = accumScore_2 + (add * add);
	nbTests++;
	testsToComputeForAdaptation.reset();
	impMixture=world;
}

mu = ((accumScore / nbTests));
sigma = sqrt((accumScore_2 / (nbTests) - (mu * mu)));    
score=(nonorm_score-mu)/sigma;

return score;	

}

//find the nearest LLR in a matrix of LLR learnt on different training duration 
double findNearestLLRInMatrix(Matrix <double> &Mat, unsigned long &index, double &LLR){
	
double delta=0.0;
double distance=0.0;
double previous=1000000;
unsigned long ind_sav=0;	
for (unsigned long i=0; i<Mat.rows(); i++){		//parcours la colonne index correspondant au nombre de session
	distance =LLR-Mat(i,index-1);
	if ((abs(distance))< abs(previous))
	{
		previous=distance;
		delta=Mat(i,index-1)-Mat(i,0);			//return the LLR shift between index session and 1 session
		if (delta<0)		//BAD SHIFT, SHOULD INCREASE
			delta=0;		//for the moment =0 instead of elimnating the locutor in matrix
		ind_sav=i;
	}
		
}
//if(debug)
cout << "LLR " << LLR <<" ; find nearest " <<  Mat(ind_sav,index-1)<< " for "<< index << " sessions"<<endl;
return delta;
}



/***************************************************************/
/* Use of the train data to assess the adaptation choice	*/
/********************************************************************/
void assessAdaptation(StatServer &ss,MixtureGD & adaptedMixture ,SegCluster &selectedSegments, 
Config &configTest, FeatureServer &fs, MixtureGD &world, DoubleVector &tmp, String &fullFileNameTrain, 
int countTests, double &tarScore){
		
double LLRtrain = 0.0,
       IntUp    = 0.0,
       IntDown  = 0.0;
	
//FROM CONFIG BUT COULD BE TAKEN FROM A GMM OF SCORES	
double stdTARscores = configTest.getParam("VarTarScores").toDouble();

LLRtrain = computeFastLLR(ss, fs, world, adaptedMixture, selectedSegments,fullFileNameTrain,configTest);

IntUp   = tarScore + stdTARscores;
IntDown	= tarScore - stdTARscores;
	
cout << "LLR for Train data on adapted model = "<<LLRtrain <<endl;
	
//OUT OF THE TAR DISTRIB
if (LLRtrain > IntUp || LLRtrain < IntDown )
{ 
    //should decrease the WMAP weight because of adaptation error	
    //score SHOULD BE a score for which the weight equal to 0 (mean of IMP distrib)!!! TO CHANGE!!!
    
    tmp[countTests-1]=configTest.getParam("meanImp").toDouble();

}



}





/*
Try to detect Baseline and unsupervised mode divergence on decision
case 1: baseline accept trials as a true target trial
	unsupervised refuses it

case 2 : baseline refuses trial
	 unsupervised accepts

The other cases (both systems agree) is not interesting

*/

void divDetect(double &LLRbaseline,double &LLRadapted, double &thresBaseline, double & thresAdapted){
	//case 1
	if(LLRbaseline > thresBaseline && LLRadapted < thresAdapted)
	{
		
	}
	
	//case 2
	if(LLRbaseline < thresBaseline && LLRadapted > thresAdapted)
	{
	    LLRadapted = LLRbaseline;
	}
}



#endif //!defined(ALIZE_UnsupervisedTools_cpp)

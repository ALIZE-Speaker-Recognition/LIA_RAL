// UnsupervisedTools.cpp
// This file is a part of LIA Softwares LIA_SpkDet and LIA_SpkSeg, based on ALIZE toolkit 
// LIA_SpkDet is a free, open tool for speaker recognition
// LIA_SpkSeg is a free, open tool for speaker segmentation
// This project is a development project initiated and funded by the LIA lab.
// See www.lia.univ-avignon.fr
// 
// ALIZE is needed for LIA_SpkDet and LIA_SpkSeg
// for more information about ALIZE, see http://www.lia.univ-avignon.fr/heberges/ALIZE/
// New version February 2005
//
// Copyright (C) 2004
//  Laboratoire d'informatique d'Avignon [www.lia.univ-avignon.fr]
// Main author:
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
//Author : Alexandre PRETI.
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
    cout << "Mean LLK Init = " << meanLikelihood(ss, FeatServ, ClusterSeg,
      clientMixture, decision, config) << endl;



  for (unsigned long trainIt = 0; trainIt < mapCfg.getNbTrainIt(); trainIt++)
    {				// Begin the initial adaptation loop (with bagged frames)
      MixtureStat & emAcc = ss.createAndStoreMixtureStat(clientMixture);	// Create a statistic accumulator using the curent model
      emAcc.resetEM();
      double llkPreviousIt = 0;
      for (unsigned long nbFs = 0; nbFs < FeatServ.size(); nbFs++)
	{

	  SegServer segServer;	// Create a local segment server 
	  SegCluster & baggedFramesCluster = segServer.createCluster(1, "", "");	// Create the cluster for describing the selected frames

	  baggedSegments((static_cast <
	      SegCluster & >(ClusterSeg.getObject(nbFs))),
	    baggedFramesCluster, mapCfg.getBaggedFrameProbability());

	  if (verboseLevel > 2)
	    cout <<
	      "Accumulate statistics on the feature server weighted by : " <<
	      decision[nbFs] << endl;
	  llkPreviousIt += accumulateStatEM(ss, (static_cast < FeatureServer & >(FeatServ.getObject(nbFs))), emAcc, baggedFramesCluster, decision[nbFs], config);	// Accumulate the EM statistics

	}
      clientMixture = emAcc.getEM();	// Get the EM estimate   
      unsigned long frameCount = (unsigned long) emAcc.getEMFeatureCount();
      cout << "Total Frames :" << frameCount << endl;
      llkPreviousIt = llkPreviousIt / (double) frameCount;
      if (verbose)
	cout << "ML (partial) estimate it[" << trainIt <<
	  "] (take care, it corresponds to the previous it,0 means init likelihood) = "
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
    cout << "Final likelihood on all frames =" << meanLikelihood(ss, FeatServ,
      ClusterSeg, clientMixture, decision, config) << endl;



  if (debug)
    cout << "adaptModel nb distrib:" << ms.
      getDistribCount() << "nb mixt:" << ms.getMixtureCount() << endl;
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
      cout << "Model adaptation based on true EM/ML estimate of training data"
	<< endl;

    }
  MixtureServer msTmp(config);
  MixtureGD & data = msTmp.duplicateMixture(initModel, DUPL_DISTRIB);

  /*unsigned long totalFrameCount1 = totalFrame(selectedSegments);
     unsigned long totalFrameCount2 = totalFrame(selectedSegmentsTests);
     unsigned long totalFrameCount = totalFrameCount1 + totalFrameCount2; */
  if (verboseLevel > 1)
    cout << "Mean LLK Init = " << meanLikelihood(ss, fs, data,
      selectedSegments, config) << endl;
  for (unsigned long emIt = 0; emIt < mapCfg.getNbEmIt(); emIt++)
    {				// begin the true EM/ML estimate of the adpatation data 
      MixtureStat & emAcc = ss.createAndStoreMixtureStat(data);	// Create a statistic accumulator using the curent model
      SegServer segServer;	// Create a local segment server
      SegCluster & baggedFramesCluster = segServer.createCluster(1, "", "");	// Create the cluster for describing the selected frames
      baggedSegments(selectedSegments, baggedFramesCluster,
	mapCfg.getBaggedFrameProbability());
      SegServer segServerTests;	// Create a local segment server
      SegCluster & baggedFramesClusterTests = segServerTests.createCluster(1, "", "");	// Create the cluster for describing the selected frames
      baggedSegments(selectedSegmentsTests, baggedFramesClusterTests,
	mapCfg.getBaggedFrameProbability());
      emAcc.resetEM();
      double llkPreviousIt = accumulateStatEM(ss, fs, emAcc, baggedFramesCluster, config);	// Accumulate the EM statistics on the first feature server
      llkPreviousIt += accumulateStatEM(ss, fsTests, emAcc, baggedFramesClusterTests, config);	// Accumulate the EM statistics on the second feature server
      data = emAcc.getEM();	// Get the EM estimate         
      unsigned long frameCount = (unsigned long) emAcc.getEMFeatureCount();
      llkPreviousIt = llkPreviousIt / (double) frameCount;
      if (verbose)
	cout << "ML (partial) estimate it[" << emIt <<
	  "] (take care, it corresponds to the previous it,0 means init likelihood) = "
	  << llkPreviousIt << endl;
      ss.deleteMixtureStat(emAcc);
    }
  // Begin the estimation of the statistic using the EM/ML model of the adaptation data
  unsigned long modelNbComp = aprioriModel.getDistribCount();
  // Begin the estimation of the statistic using the EM/ML model of the adaptation data
  unsigned long vectSize = fs.getVectSize();	// Complete log likelihood of the adaptation data given the apriori model
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
	  apProba[idxModel] =
	    aprioriModel.weight(idxModel) * likelihoodGD(d, m);
	  totLk += apProba[idxModel];
	}
      for (unsigned long idxModel = 0; idxModel < modelNbComp; idxModel++)
	{
	  DistribGD & c = clientMixture.getDistrib(idxModel);
	  apProba[idxModel] /= totLk;
	  for (unsigned long idxC = 0; idxC < vectSize; idxC++)
	    {
	      c.setMean(c.getMean(idxC) +
		(d.getMean(idxC) * apProba[idxModel] * data.weight(idxData)),
		idxC);
	      c.setCov(c.getCov(idxC) +
		((d.getMean(idxC) * d.getMean(idxC)) * apProba[idxModel] *
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
    {				// begin the true EM/ML estimate of the adpatation data 
      MixtureStat & emAcc = ss.createAndStoreMixtureStat(data);	// Create a statistic accumulator using the curent model
      SegServer segServer;	// Create a local segment server
      SegCluster & baggedFramesCluster = segServer.createCluster(1, "", "");	// Create the cluster for describing the selected frames
      baggedSegments(selectedSegments, baggedFramesCluster,
	mapCfg.getBaggedFrameProbability());
      emAcc.resetEM();
      double llkPreviousIt = accumulateStatEM(ss, fs, emAcc, baggedFramesCluster, config);	// Accumulate the EM statistics
      data = emAcc.getEM();	// Get the EM estimate         
      unsigned long frameCount = (unsigned long) emAcc.getEMFeatureCount();
      llkPreviousIt = llkPreviousIt / (double) frameCount;
      if (verbose)
	cout << "ML (partial) estimate it[" << emIt <<
	  "] (take care, it corresponds to the previous it,0 means init likelihood) = "
	  << llkPreviousIt << endl;
      ss.deleteMixtureStat(emAcc);
    }
  unsigned long modelNbComp = aprioriModel.getDistribCount();
  // Begin the estimation of the statistic using the EM/ML model of the adaptation data
  unsigned long vectSize = fs.getVectSize();	// Complete log likelihood of the adaptation data given the apriori model
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
	  apProba[idxModel] =
	    aprioriModel.weight(idxModel) * likelihoodGD(d, m);
	  totLk += apProba[idxModel];
	}
      for (unsigned long idxModel = 0; idxModel < modelNbComp; idxModel++)
	{
	  DistribGD & c = clientMixture.getDistrib(idxModel);
	  apProba[idxModel] /= totLk;
	  for (unsigned long idxC = 0; idxC < vectSize; idxC++)
	    {
	      c.setMean(c.getMean(idxC) +
		(d.getMean(idxC) * apProba[idxModel] * data.weight(idxData)),
		idxC);
	      c.setCov(c.getCov(idxC) +
		((d.getMean(idxC) * d.getMean(idxC)) * apProba[idxModel] *
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


// adapt model function to compute MAP from a given EM estimate
void adaptModelMAP(Config & config, StatServer & ss, MixtureServer & ms,
  MixtureGD & aprioriModel, MixtureGD & clientMixture, MAPCfg & mapCfg,
  unsigned long &frameCount)
{
  if (verbose)
    mapCfg.showConfig(cout);
  // if (verboseLevel>1) cout << "Mean LLK Init = " << meanLikelihood(ss,fs,clientMixture,selectedSegments,config)<< endl;    
  for (unsigned long trainIt = 0; trainIt < mapCfg.getNbTrainIt(); trainIt++)
    {				// Begin the initial adaptation loop (with bagged frames)
      computeMAP(ms, aprioriModel, clientMixture, frameCount, config);	// Bayesian Adaptation client=MAP(aprioriModel,client)
      if (mapCfg.getNormalizeModel())
	normalizeMixture(clientMixture, mapCfg, config);	// Normalize/fit the model if needed

      // if (verboseLevel>2) cout << "Likelihood on all frames= "<< meanLikelihood(ss,fs,clientMixture, selectedSegments,config) << endl;        
    }
  // if (verboseLevel>1) cout << "Final likelihood on all frames= "<< meanLikelihood(ss,fs,clientMixture, selectedSegments,config) << endl;                               
 
  if (debug)
    cout << "adaptModel nb distrib:" << ms.
      getDistribCount() << "nb mixt:" << ms.getMixtureCount() << endl;
}
void adaptModelMAP(Config & config, StatServer & ss, MixtureServer & ms,
  MixtureGD & aprioriModel, MixtureGD & clientMixture,
  unsigned long &frameCount)
{
  MAPCfg mapCfg(config);
  adaptModelMAP(config, ss, ms, aprioriModel, clientMixture, mapCfg,
    frameCount);
}




// Adapt model function to compute an EM estimate 

void adaptModelEM(Config & config, StatServer & ss, MixtureServer & ms,
  FeatureServer & fs, SegCluster & selectedSegments, MixtureGD & aprioriModel,
  MixtureGD & clientMixture, MAPCfg & mapCfg)
{
  unsigned long frameCount;	//A.P.
  if (verbose)
    cout << "Model adaptation based on true EM/ML estimate of training data"
      << endl;

  if (verboseLevel > 1)
    cout << "Mean LLK Init = " << meanLikelihood(ss, fs, clientMixture,
      selectedSegments, config) << endl;

  for (unsigned long trainIt = 0; trainIt < mapCfg.getNbTrainIt(); trainIt++)
    {				// Begin the initial adaptation loop (with bagged frames)
      MixtureStat & emAcc = ss.createAndStoreMixtureStat(clientMixture);	// Create a statistic accumulator using the curent model
      SegServer segServer;	// Create a local segment server 
      SegCluster & baggedFramesCluster = segServer.createCluster(1, "", "");	// Create the cluster for describing the selected frames
      baggedSegments(selectedSegments, baggedFramesCluster,
	mapCfg.getBaggedFrameProbability());
      emAcc.resetEM();
      double llkPreviousIt = accumulateStatEM(ss, fs, emAcc, baggedFramesCluster, config);	// Accumulate the EM statistics
      clientMixture = emAcc.getEM();
      // Get the EM estimate       
      frameCount = (unsigned long) emAcc.getEMFeatureCount();

      llkPreviousIt = llkPreviousIt / (double) frameCount;
      if (verbose)
	cout << "ML (partial) estimate it[" << trainIt <<
	  "] (take care, it corresponds to the previous it,0 means init likelihood) = "
	  << llkPreviousIt << endl;
      ss.deleteMixtureStat(emAcc);


      if (verboseLevel > 2)
	cout << "Likelihood on all frames= " << meanLikelihood(ss, fs,
	  clientMixture, selectedSegments, config) << endl;

    }
  if (verboseLevel == 2)
    cout << "Final likelihood on all frames= " << meanLikelihood(ss, fs,
      clientMixture, selectedSegments, config) << endl;

}

void adaptModelEM(Config & config, StatServer & ss, MixtureServer & ms,
  FeatureServer & fs, SegCluster & selectedSegments, MixtureGD & aprioriModel,
  MixtureGD & clientMixture)
{
  MAPCfg mapCfg(config);
  adaptModelEM(config, ss, ms, fs, selectedSegments, aprioriModel,
    clientMixture, mapCfg);
}

// Adapt model function to compute an EM estimate from 2 feature server

void adaptModelEM(Config & config, StatServer & ss, MixtureServer & ms,
  FeatureServer & fs, SegCluster & selectedSegments, FeatureServer & fs2,
  SegCluster & selectedSegments2, MixtureGD & aprioriModel,
  MixtureGD & clientMixture, MAPCfg & mapCfg)
{
  unsigned long frameCount;	//A.P.
  if (verbose)
    cout << "Model adaptation based on true EM/ML estimate of training data"
      << endl;
  if (verboseLevel > 1)
    cout << "Mean LLK Init = " << meanLikelihood(ss, fs, clientMixture,
      selectedSegments, config) << endl;

  for (unsigned long trainIt = 0; trainIt < mapCfg.getNbTrainIt(); trainIt++)
    {				// Begin the initial adaptation loop (with bagged frames)
      double llkPreviousIt = 0;
      MixtureStat & emAcc = ss.createAndStoreMixtureStat(clientMixture);	// Create a statistic accumulator using the curent model
      SegServer segServer;	// Create a local segment server 
      SegServer segServer2;	// Create a local segment server 
      SegCluster & baggedFramesCluster = segServer.createCluster(1, "", "");	// Create the cluster for describing the selected frames
      SegCluster & baggedFramesCluster2 = segServer2.createCluster(1, "", "");	// Create the cluster for describing the selected frames
      baggedSegments(selectedSegments, baggedFramesCluster,
	mapCfg.getBaggedFrameProbability());
      baggedSegments(selectedSegments2, baggedFramesCluster2,
	mapCfg.getBaggedFrameProbability());
      emAcc.resetEM();
      llkPreviousIt += accumulateStatEM(ss, fs, emAcc, baggedFramesCluster, config);	// Accumulate the EM statistics on the 1st FS
      llkPreviousIt += accumulateStatEM(ss, fs2, emAcc, baggedFramesCluster2, config);	// Accumulate the EM statistics on the 2nd FS
      clientMixture = emAcc.getEM();
      // Get the EM estimate       
      frameCount = (unsigned long) emAcc.getEMFeatureCount();

      llkPreviousIt = llkPreviousIt / (double) frameCount;
      if (verbose)
	cout << "ML (partial) estimate it[" << trainIt <<
	  "] (take care, it corresponds to the previous it,0 means init likelihood) = "
	  << llkPreviousIt << endl;
      ss.deleteMixtureStat(emAcc);


      if (verboseLevel > 2)
	cout << "Likelihood on all frames= " << meanLikelihood(ss, fs,
	  clientMixture, selectedSegments, config) << endl;

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
    {				// For each of the selected segments
      unsigned long idxBeginFrame =
	seg->begin() +
	fsTests.getFirstFeatureIndexOfASource(seg->sourceName());
      fsTests.seekFeature(idxBeginFrame);
      Feature f;
      for (unsigned long idxFrame = 0; idxFrame < seg->length(); idxFrame++)
	{			// For each frame of the segment
	  fsTests.readFeature(f);
	  ss.computeAndAccumulateLLK(world, f, DETERMINE_TOP_DISTRIBS);	// Determine the top components and compute wrld LLK

	  ss.computeAndAccumulateLLK(clientMixture, f, USE_TOP_DISTRIBS);
	}
    }

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



double computeLLR(Config & config, StatServer & ss, FeatureServer & fsTests,
  MixtureGD & world, MixtureGD & clientMixture,
  SegCluster & selectedSegmentsTests, String & idTest)
{
  cout << "compute LLR" << endl;
  FileInfo FI(idTest);		//create file to write top components info
  ss.resetLLK(world);		// Reset the world LLK accumulator
  ss.resetLLK(clientMixture);	// ss.resetLLK(tabClientLine.getClientModel(i));                                   // Reset client LLK accumulator
  Seg *seg;			// reset the reader at the begin of the input stream
  selectedSegmentsTests.rewind();
  while ((seg = selectedSegmentsTests.getSeg()) != NULL)
    {				// For each of the selected segments
      unsigned long idxBeginFrame =
	seg->begin() +
	fsTests.getFirstFeatureIndexOfASource(seg->sourceName());
      fsTests.seekFeature(idxBeginFrame);
      Feature f;
      for (unsigned long idxFrame = 0; idxFrame < seg->length(); idxFrame++)
	{			// For each frame of the segment
	  fsTests.readFeature(f);
	  ss.computeAndAccumulateLLK(world, f, DETERMINE_TOP_DISTRIBS);	// Determine the top components and compute wrld LLK
	  //STOCK the vector with the top selected component for reuse 
	  const LKVector & lkv = ss.getTopDistribIndexVector();
	  FI.writeTopInfo(lkv, config);
	  ss.computeAndAccumulateLLK(clientMixture, f, USE_TOP_DISTRIBS);
	}
    }
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





double computeFastLLR(StatServer & ss, FeatureServer & fsTests,
  MixtureGD & world, MixtureGD & clientMixture,
  SegCluster & selectedSegmentsTests, String & idTest, Config & config)
{
  cout << "FAST LLR COMPUTATION" << endl;
  FileInfo FI(idTest);		//load file to read top components info
  ss.resetLLK(world);		// Reset the world LLK accumulator
  ss.resetLLK(clientMixture);	// ss.resetLLK(tabClientLine.getClientModel(i));                                   // Reset client LLK accumulator
  Seg *seg;			// reset the reader at the begin of the input stream
  selectedSegmentsTests.rewind();
  unsigned long id = 0;		//id of the frame
  while ((seg = selectedSegmentsTests.getSeg()) != NULL)	// For each of the selected segments
    {
      unsigned long idxBeginFrame =
	seg->begin() +
	fsTests.getFirstFeatureIndexOfASource(seg->sourceName());
      fsTests.seekFeature(idxBeginFrame);
      Feature f;
      for (unsigned long idxFrame = 0; idxFrame < seg->length(); idxFrame++)	// For each frame of the segment
	{
	  fsTests.readFeature(f);
	  FI.loadTopInfo(ss, id, config);
	  ss.computeAndAccumulateLLK(world, f, USE_TOP_DISTRIBS);	// uses the top components and compute wrld LLK
	  id++;			//id = index dans la XLIST
	  ss.computeAndAccumulateLLK(clientMixture, f, USE_TOP_DISTRIBS);
	}
    }
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
  double LLKClient =
    likelihoodGD(dataTest, clientMixture, TabWeightData, TabWeightClient);
  double LLKWorld = likelihoodGD(dataTest, world, TabWeightClient, TabWeightClient);	//LLR between models 
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
  double SIGMAimp = configTest.getParam("SIGMAimp").toDouble();
  double MUclient = configTest.getParam("MUtarget").toDouble();
  double MUimp = configTest.getParam("MUimp").toDouble();
  double poidsImp = configTest.getParam("IMPweight").toDouble();
  double poidsTar = configTest.getParam("TARweight").toDouble();
  double seuilMin = configTest.getParam("thrMin").toDouble();
  double seuilMax = configTest.getParam("thrMax").toDouble();
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

bool FileExists(String & fullFileName)
{
  ifstream outFile(fullFileName.c_str(), ios::binary);
  if (outFile.is_open())
    return true;
  else
    return false;


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
      if (e == 0)
	decision[0] = 1;	//the first LLR is the train data on the target model, we give WMAP=1. 
      else
	{
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
		    cout << "Flooring LLK at " << configTest.
		      getParam("LLKthreshold").toLong() << endl;
		  llkTar = configTest.getParam("LLKthreshold").toLong();
		}
	      llkTar = exp(llkTar);
	      llkNon = exp(llkNon);
	      decision[e] =
		llkTar * poidsTar / (llkTar * poidsTar + llkNon * poidsImp);
	    }
	  if (verbose)
	    cout << "LLR : " << f[0] << " , After WMAP GMM [" << e << "] = "
	      << decision[e] << endl;
	}

    }
}


//WMAP GMM with adaptative priors

void WMAPGMM(DoubleVector & decision, Config & configTest, MixtureGD & tar,
  MixtureGD & non, StatServer & ss)
{
  cout << " Choice of WMAP GMM to compute Feature Server Probabilities" <<
    endl;



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
      decision[e] =
	llkTar * priorTar[e] / (llkTar * priorTar[e] + llkNon * priorImp[e]);


      if (verbose)
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
  aux = ms.getMixtureGD(ms.getMixtureIndex(testsToCompute.getElement(0)));	//get the first model , its name is stocked in the XLine

  unsigned long auxNbFrame = (unsigned long) (nbFramesSelected[0] * decision[0]);	//the elements in the ObjectRefvector follow the XLine order 
  //i.e. the objectRefvector and the XLine are reseted after each client                  

  for (i = 0; i < nbModels - 1; i++)
    {
      if (debug)
	cout << "Fusing models : " << testsToCompute.
	  getElement(i) << " with " << testsToCompute.getElement(i +
	  1) << endl;
      fuseModels(aux, auxNbFrame,
	ms.getMixtureGD(ms.getMixtureIndex(testsToCompute.getElement(i + 1))),
	(unsigned long) (nbFramesSelected[i + 1] * decision[i + 1]), tmp);
      aux = tmp;
      auxNbFrame +=
	(unsigned long) (nbFramesSelected[i + 1] * decision[i + 1]);
    }
  adaptModelMAP(config, ss, ms, aprioriModel, tmp, auxNbFrame);
  clientMixture = tmp;
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

//load impostors scores for TNORM computation
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


//Do Tnormalization on a score
void TnormalizeScore(String & test, double &decision,
  ObjectRefVector & stockTnorm)
{

  if (verbose && verboseLevel > 2)
    cout << " Applying TNORM, score =  " << decision << endl;

  for (unsigned long i = 0; i < stockTnorm.size(); i++)
    {
      if (static_cast < Norm & >(stockTnorm.getObject(i)).idTest == test)
	decision =
	  (decision - static_cast <
	  Norm & >(stockTnorm.getObject(i)).mu) / static_cast <
	  Norm & >(stockTnorm.getObject(i)).sigma;
    }
  if (verbose && verboseLevel > 2)
    cout << " T-Normed score =  " << decision << endl;

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

  bool wmap = false;
  bool one = true;
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
	      //score=1;       true Oracle with weight = 1    
	      //HERE : leave the WMAP weight

	      DoubleVector tmp;
	      tmp.addValue(0);	//The first element is not evaluated by WMAP.
	      tmp.addValue(score);
	      WMAPGMM(tmp, config, tar, non, ss);
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

//Jack knife on train data
double computeLLRForTrain(Config & configTest, StatServer & ss,
  MixtureServer & ms, FeatureServer & fs, SegCluster & selectedSegments,
  MixtureGD & aprioriModel)
{
  MixtureGD & clientMixture = ms.duplicateMixture(aprioriModel, DUPL_DISTRIB);
  double LLR = 0.0;
  MAPCfg mapCfg(configTest);
  double trainSelected = configTest.getParam("SelectedTrain").toDouble();
  double testSelected = configTest.getParam("SelectedTest").toDouble();
  cout << "Compute Client Model on " << (trainSelected *
    100) << "% of train data" << endl;


  if (verboseLevel > 1)
    cout << "Mean LLK Init = " << meanLikelihood(ss, fs, clientMixture,
      selectedSegments, configTest) << endl;
  for (unsigned long trainIt = 0; trainIt < mapCfg.getNbTrainIt(); trainIt++)
    {				// Begin the initial adaptation loop (with bagged frames)
      MixtureStat & emAcc = ss.createAndStoreMixtureStat(clientMixture);	// Create a statistic accumulator using the curent model
      SegServer segServer;	// Create a local segment server 
      SegCluster & baggedFramesCluster = segServer.createCluster(1, "", "");	// Create the cluster for describing the selected frames
      baggedSegments(selectedSegments, baggedFramesCluster, trainSelected);	//Train the client model on trainSelected Frame
      emAcc.resetEM();
      double llkPreviousIt = accumulateStatEM(ss, fs, emAcc, baggedFramesCluster, configTest);	// Accumulate the EM statistics
      clientMixture = emAcc.getEM();	// Get the EM estimate   
      unsigned long frameCount = (unsigned long) emAcc.getEMFeatureCount();
      llkPreviousIt = llkPreviousIt / (double) frameCount;
      if (verbose)
	cout << "ML (partial) estimate it[" << trainIt <<
	  "] (take care, it corresponds to the previous it,0 means init likelihood) = "
	  << llkPreviousIt << endl;
      computeMAP(ms, aprioriModel, clientMixture, frameCount, configTest);	// Bayesian Adaptation client=MAP(aprioriModel,client)
      if (mapCfg.getNormalizeModel())
	normalizeMixture(clientMixture, mapCfg, configTest);	// Normalize/fit the model if needed
      ss.deleteMixtureStat(emAcc);
      if (verboseLevel > 2)
	cout << "Likelihood on all frames= " << meanLikelihood(ss, fs,
	  clientMixture, selectedSegments, configTest) << endl;
    }
  if (verboseLevel > 1)
    cout << "Final likelihood on all frames= " << meanLikelihood(ss, fs,
      clientMixture, selectedSegments, configTest) << endl;

  cout << "Compute LLR of  " << (testSelected *
    100) << "% of train data on the client model learnt on " << (trainSelected
    * 100) << "% train data" << endl;



  for (int It = 0; It < configTest.getParam("AverageIt").toLong(); It++)
    {				//Compute LLR "It" times to have a better estimation (average)

      SegServer segServer;	// Create a local segment server 
      SegCluster & baggedFramesCluster = segServer.createCluster(1, "", "");	// Create the cluster for describing the selected frames
      baggedSegments(selectedSegments, baggedFramesCluster, testSelected);	//Train the client model on trainSelected Frame

      double accum =
	computeLLR(ss, fs, aprioriModel, clientMixture, baggedFramesCluster);

      if (verbose && verboseLevel > 2)
	cout << "LLR It[" << It << "] = " << accum << endl;
      LLR += accum;
    }

  LLR = LLR / (configTest.getParam("AverageIt").toLong());
  if (verbose && verboseLevel > 1)
    cout << "LLR of train data = " << LLR << endl;

  return LLR;
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
#endif //!defined(ALIZE_UnsupervisedTools_cpp)

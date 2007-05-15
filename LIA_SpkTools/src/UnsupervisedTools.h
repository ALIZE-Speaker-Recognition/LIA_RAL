// UnsupervisedTools.h
// This file is a part of LIA Softwares LIA_SpkDet and LIA_SpkSeg, based on ALIZE toolkit 
// LIA_SpkDet is a free, open tool for speaker recognition
// LIA_SpkSeg is a free, open tool for speaker segmentation
// This project is a development project initiated and funded by the LIA lab.
// See www.lia.univ-avignon.fr
// 
// ALIZE is needed for LIA_SpkDet and LIA_SpkSeg
// for more information about ALIZE, see http://www.lia.univ-avignon.fr/heberges/ALIZE/
// First Version 07/15/2004
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
//Author: Alexandre PRETI.

#if !defined(ALIZE_UnsupervisedTools_h)
#define ALIZE_UnsupervisedTools_h

#include <alize.h>
#include "liatools.h"

using namespace alize;
using namespace std;



//Accumulate statistics on a frame, a segment,a cluster.
void accumulateStatLK(StatServer & ss, FeatureServer & fs, MixtureStat & acc, unsigned long idxBeginFrame, unsigned long nbFrames, Config & config);
void accumulateStatLK(StatServer & ss, FeatureServer & fs, MixtureStat & acc, Seg * seg, Config & config);
void accumulateStatLK(StatServer & ss, FeatureServer & fs, MixtureStat & acc, SegCluster & selectedSegments, Config & config);




// The main function for estimate a client model by bayesian adaptattion of a world model
// Using EM and MAP 
//void adaptModel(Config& config,StatServer &ss,MixtureServer &ms,FeatureServer &fs,SegCluster& selectedSegments,MixtureGD &world,MixtureGD &clientMixture);//A.P.
void adaptModelEM(Config& config,StatServer &ss,MixtureServer &ms,FeatureServer &fs,SegCluster& selectedSegments,MixtureGD &world,MixtureGD &clientMixture);//A.P.
void adaptModelEM(Config& config,StatServer &ss,MixtureServer &ms,FeatureServer &fs,SegCluster& selectedSegments,FeatureServer &fs2,SegCluster& selectedSegments2,MixtureGD &aprioriModel,MixtureGD &clientMixture,MAPCfg &mapCfg);
void adaptModelEM(Config& config,StatServer &ss,MixtureServer &ms,FeatureServer &fs,SegCluster& selectedSegments,FeatureServer &fs2,SegCluster& selectedSegments2,
		MixtureGD &aprioriModel,MixtureGD &clientMixture);
void adaptModelMAP(Config& config,StatServer &ss,MixtureServer &ms,MixtureGD &world,MixtureGD &clientMixture, unsigned long &frameCount);
void adaptModelMAP(Config & config, StatServer & ss, MixtureServer & ms, MixtureGD & aprioriModel, MixtureGD & clientMixture, MAPCfg & mapCfg,  unsigned long &frameCount);
//void adaptModel(Config& config,StatServer &ss,MixtureServer &ms,FeatureServer &fs,SegCluster& selectedSegments,MixtureGD &world,MixtureGD &clientMixture, MAPCfg &cfg);
// ** New training algo based on a true EM/ML estimate of the training data before to apply MAP
void modelBasedadaptModel(Config& config,StatServer &ss,MixtureServer &ms,FeatureServer &fs,SegCluster& selectedSegments,
			  MixtureGD &aprioriModel,MixtureGD &clientMixture, MixtureGD &initModel);
void modelBasedadaptModelEM(Config& config,StatServer &ss,MixtureServer &ms,FeatureServer &fs,SegCluster& selectedSegments,
			  MixtureGD &aprioriModel,MixtureGD &clientMixture, MixtureGD &initModel);
void modelBasedadaptModelEM(Config& config,StatServer &ss,MixtureServer &ms,FeatureServer &fs,SegCluster& selectedSegments,FeatureServer &fsTests,SegCluster& selectedSegmentsTests,MixtureGD &aprioriModel,MixtureGD &clientMixture, MixtureGD &initModel);
void adaptModel(Config& config,StatServer &ss,MixtureServer &ms,ObjectRefVector &FeatServ,ObjectRefVector & ClusterSeg,
		MixtureGD &aprioriModel,MixtureGD &clientMixture,DoubleVector &decision);








//Compute LLR 
//basic
double computeLLR(StatServer & ss, FeatureServer & fsTests, MixtureGD & world,
MixtureGD & clientMixture, SegCluster & selectedSegmentsTests);
//basic and write top information 
double computeLLR(Config & config, StatServer & ss, FeatureServer & fsTests,  MixtureGD & world, MixtureGD & clientMixture,
SegCluster & selectedSegmentsTests, String & idTest);
//Fast by using top information
double computeFastLLR(StatServer & ss, FeatureServer & fsTests, MixtureGD & world, MixtureGD & clientMixture, SegCluster & selectedSegmentsTests,
String & idTest, Config & config);
//modelbased
double computeLLRGD(Config & config, MixtureGD & clientMixture, MixtureGD & world,  MixtureGD & dataTest);

//compute weights used for adaptation for each tests
//logistic regression, take parameter from the config
void expandLLR(DoubleVector & decision, Config & configTest);
//WMAP using gaussian formulae
void WMAP(DoubleVector & decision, Config & configTest);
//WMAP on GMM on scores
void WMAPGMM(DoubleVector & decision, Config & configTest, MixtureGD & tar,  MixtureGD & non, StatServer & ss);
void WMAPGMMFixedPriors(DoubleVector & decision, Config & configTest, MixtureGD & tar,
  MixtureGD & non, StatServer & ss);
String getFullFileName(String & id, Config & c);
bool FileExists(String & fullFileName);

//compute a MAP model from EM models (fast)
//takes a list of EM models, the number of frames each model was learnt with, and a weight for each model  
void computeMAPmodelFromEMones(Config & config, StatServer & ss,  MixtureServer & ms, DoubleVector & nbFramesSelected,
MixtureGD & aprioriModel, MixtureGD & clientMixture, MixtureGD & aux,  MixtureGD & tmp, DoubleVector & decision, XLine & testsToCompute);

double SegClusterFrame(SegCluster & SegC); //return the number of frames in a cluster

void computePriors(DoubleVector & decision, DoubleVector &priorImp, DoubleVector &priorTar,Config & configTest); //update Imp and Tar priors during the adaptation process

//load TNORM parameters, need trials scores on impostor models, a list of the trials name.
void loadTnormParam(String &inputClientListFileName,String &testFile,ObjectRefVector &stockTnorm,Config &config);

//T-Normalize a score

void TnormalizeScore(String &test, double &decision, ObjectRefVector &stockTnorm);

//Oracle TEST 

void resetWeights(DoubleVector & decision);
void Oracle(String &idTar, String & idTest, double &score, Config & config,MixtureGD & tar,
  MixtureGD & non, StatServer & ss);

//Compute a train model with a selected percentage of train data, and compute the LLR with the rest of the data (this is done several times for robustness)
double computeLLRForTrain(Config & configTest, StatServer & ss, MixtureServer & ms, FeatureServer & fs, SegCluster & selectedSegments, MixtureGD & aprioriModel);


//2 purposes : do not recalculate LLR of test data on a client model, and possibility to use RES file from a different recognition system to compute adaptation weights
double searchLLRFromResFile(String &idTar, String &test ,String &inputResFilename, Config & config);
#endif //!defined(ALIZE_UnsupervisedTools_h)

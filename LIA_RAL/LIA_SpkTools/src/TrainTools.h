// TrainTools.h
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

#if !defined(ALIZE_TrainTools_h)
#define ALIZE_TrainTools_h

#include <alize.h>
#include "liatools.h"

using namespace alize;
using namespace std;

class MAPCfg {
  String _method;
  double _r[3];
  bool _mean,_var,_weight;
  unsigned long _nbTrainIt; 
  unsigned long _nbEmIt; // used only for modelBasedadaptMode
  bool  _normalizeModel;
  bool _normalizeModelMeanOnly;
  unsigned long _normalizeModelNbIt;
  double _baggedFrameProbability;
public:
  MAPCfg(Config &);
  String &getMethod(){return _method;}
  void showConfig(ostream & st);
  bool getMeanAdapt() { return _mean;}
  bool getVarAdapt() { return _var;}
  bool getWeightAdapt(){ return _weight;}
  double getMeanReg() { return _r[0];}
  double getVarReg() { return _r[1];}
  double getWeightReg() { return _r[2];}
  double getMeanAlpha() { return _r[0];}
  double getVarAlpha() { return _r[1];}
  double getWeightAlpha() { return _r[2];}
  unsigned long getNbTrainIt(){return _nbTrainIt;}
  unsigned long getNbEmIt(){return _nbEmIt;} // used only for modelBasedadaptMode
  bool getNormalizeModel(){return  _normalizeModel;}
  bool getNormalizeModelMeanOnly(){return _normalizeModelMeanOnly;}
  unsigned long getNormalizeModelNbIt(){return _normalizeModelNbIt;}
  double getBaggedFrameProbability(){return _baggedFrameProbability;}
  void setMethod(String v) {_method=v;}
  void setMeanAdapt(bool v) {_mean=v;}
  void setVarAdapt(bool v) { _var=v;}
  void setWeightAdapt(bool v){ _weight=v;}
  void setMeanReg(double v) { _r[0]=v;}
  void setVarReg(double v) { _r[1]=v;}
  void setWeightReg(double v) { _r[2]=v;}
  void setMeanAlpha(double v) { _r[0]=v;}
  void setVarAlpha(double v) { _r[1]=v;}
  void setWeightAlpha(double v) { _r[2]=v;}  
  void setNbTrainIt(unsigned long v){_nbTrainIt=v;}
  void setNbEmIt(unsigned long v){_nbEmIt=v;} // used only for modelBasedadaptMode
  void setNormalizeModel(bool v){_normalizeModel=v;}
  void setNormalizeModelMeanOnly(bool v){_normalizeModelMeanOnly=v;}  
  void setNormalizeModelNbIt(unsigned long v){_normalizeModelNbIt=v;}
  void setBaggedFrameProbability(double v){_baggedFrameProbability=v;}
};


class TrainCfg{
  double _initVarianceFlooring; 
  double _initVarianceCeiling; 
  double _finalVarianceFlooring;
  double _finalVarianceCeiling; 
  unsigned long _nbTrainIt; 
  bool  _normalizeModel;
  bool _normalizeModelMeanOnly;
  unsigned long _normalizeModelNbIt;
  double _baggedFrameProbability;
  double _baggedFrameProbabilityInit;
  bool _componentReduction;
  unsigned long _targetDistribCount;
 public:
  TrainCfg(Config &);
  double getInitVarFloor(){return _initVarianceFlooring;}
  double getInitVarCeil(){return _initVarianceCeiling;}
  double getFinalVarFloor(){return _finalVarianceFlooring;}
  double getFinalVarCeil(){return _finalVarianceCeiling;}
  unsigned long getNbTrainIt(){return _nbTrainIt;}
  bool getNormalizeModel(){return  _normalizeModel;}
  bool getNormalizeModelMeanOnly(){return _normalizeModelMeanOnly;}
  unsigned long getNormalizeModelNbIt(){return _normalizeModelNbIt;}
  double getBaggedFrameProbability(){return _baggedFrameProbability;}
  double getBaggedFrameProbabilityInit(){return _baggedFrameProbabilityInit;}
  bool getComponentReduction(){return _componentReduction;}
  unsigned long getTargetDistribCount(){return _targetDistribCount;}
  void setInitVarFlooring (double v){_initVarianceFlooring=v;}
  void setInitVarCeiling (double v){_initVarianceCeiling=v;} 
  void setFinalVarFlooring (double v){_finalVarianceFlooring=v;}
  void setFinalVarCeiling (double v){_finalVarianceCeiling=v;} 
  void setNbTrainIt(unsigned long v){_nbTrainIt=v;}
  void setNormalizeModel(bool v){_normalizeModel=v;}
  void setNormalizeModelMeanOnly(bool v){_normalizeModelMeanOnly=v;}  
  void setNormalizeModelNbIt(unsigned long v){_normalizeModelNbIt=v;}
  void setBaggedFrameProbability(double v){_baggedFrameProbability=v;}
  void setBaggedFrameProbabilityInit(double v){_baggedFrameProbabilityInit=v;}
  void showConfig(ostream &);
};
// Mixture and ditrib tools
void fuseModels(const MixtureGD &,unsigned long,const MixtureGD &,unsigned long,MixtureGD &);
void copyMixture(DistribGD &,DistribGD &);
unsigned long selectComponent(bool selectCompA[],XList & distribL, MixtureGD &inputM);
unsigned long selectComponent(bool selectCompA[],double wFactor,MixtureGD &inputM);
unsigned long selectComponent(bool selectCompA[],unsigned long nbTop,MixtureGD &inputM);
double reduceModel(bool selectCompA[],MixtureGD &inputM,MixtureGD &outputM);
void normalizeWeights(MixtureGD &outputM);

// gaussian/mixture fusion using less likelihood loss criterion
void gaussianFusion(const DistribGD &g1,double w1,const DistribGD & g2,double w2, DistribGD &res,double &w);
void mixtureFusion(const MixtureGD &mixt,DistribGD &res,double &wres);

// normalizeMixture() normmalizes  the mixture in order to fit
// the data distribution
// Usually used with a mean=0, cov=1 distribution
void normalizeMixture(MixtureGD &mixt,const DoubleVector &meanSignal,const DoubleVector &covSignal,bool zeroOne,
		      unsigned long nbIt, bool meanOnly,Config &config);
void normalizeMixture(MixtureGD &mixt,const DoubleVector &meanSignal,const DoubleVector &covSignal,Config &config);
// for compatibility reasons only, same function but the target is always D(mean=0,std=1)
void normalizeMixture(MixtureGD &mixt,Config &config);
void normalizeMixture(MixtureGD & mixt, MAPCfg & mapCfg, Config & config);

double likelihoodLoss(const DistribGD &g1,double w1,const DistribGD &g2,double w2);
double *outputMAPTransformation(MixtureGD & world, MixtureGD & clientMixture,
  Config & config);
//-------------------------------------------------------------------------
// MAP computation functions 
// the main MAP functions
void computeMAP(MixtureServer &ms,const MixtureGD& world,MixtureGD &client,unsigned long frameCount,Config &config);
void computeMAP(MixtureServer &ms,const MixtureGD& world,MixtureGD &client,unsigned long frameCount,MAPCfg &cfg);

// Model manipulation tools
void copyMean(MixtureGD & mixtS, MixtureGD & mixtD);
void copyVar(MixtureGD & mixtS, MixtureGD & mixtD);
void copyWeight(MixtureGD & mixtS, MixtureGD & mixtD);

//-------------------------------------------------------------------------
double setItParameter(double begin, double end, int nbIt, int it);

//-------------------------------------------------------------------------
void varianceControl(MixtureGD& model,double flooring,double ceiling,const DoubleVector &covSignal);

// **********************************************************************
// Mixture Initialization stuff
//***********************************************************************
// cov and mean initialization
unsigned long computeMeanCov(Config &config,FeatureServer **fsTab,SegCluster ** segTab,unsigned long nbStream,DoubleVector &mean,DoubleVector &cov);
unsigned long computeMeanCov(Config &config,FeatureServer &fs,SegCluster &seg,DoubleVector &mean,DoubleVector &cov);
void initialize01(unsigned long vectSize,DoubleVector &mean,DoubleVector &cov);
// Mixture initialisation, based on random picking of frames
MixtureGD &mixtureInit(MixtureServer &ms,FeatureServer &fs,MixtureGD &world,
	               SegCluster &selectedSegments,const DoubleVector &globalCov, Config& config);
MixtureGD &mixtureInit(MixtureServer &ms,FeatureServer &fs,MixtureGD &world,
	               SegCluster &selectedSegments,const DoubleVector &globalCov, Config& config,TrainCfg & trainCfg);
MixtureGD &mixtureInit(MixtureServer &ms,FeatureServer **fsTab, SegCluster **segTab,double *weightTab,unsigned long nbStream,MixtureGD &world,
		       const DoubleVector &globalCov, Config& config,TrainCfg &trainCfg);
MixtureGD &mixtureInit(MixtureServer &ms,FeatureServer **fsTab, SegCluster **segTab,double*weightTab,unsigned long nbStream,MixtureGD &world,
		       const DoubleVector &globalCov, Config& config);
// The main function for estimate a client model by bayesian adaptattion of a world model
// Using EM and MAP 
void adaptModel(Config& config,StatServer &ss,MixtureServer &ms,FeatureServer &fs,SegCluster& selectedSegments,MixtureGD &world,MixtureGD &clientMixture);//A.P.
void adaptModel(Config& config,StatServer &ss,MixtureServer &ms,FeatureServer &fs,SegCluster& selectedSegments,MixtureGD &world,MixtureGD &clientMixture, MAPCfg &cfg);
// ** New training algo based on a true EM/ML estimate of the training data before to apply MAP
void modelBasedadaptModel(Config& config,StatServer &ss,MixtureServer &ms,FeatureServer &fs,SegCluster& selectedSegments,
			  MixtureGD &aprioriModel,MixtureGD &clientMixture, MixtureGD &initModel);

		

		

//-------------------------------------------------------------------------
// TrainModel is the main function for model training
// It works on an initialized model and on a clmuster of segments
void trainModel(Config& config,StatServer &ss,FeatureServer &fs,SegCluster& selectedSegments,
		DoubleVector & globalMean,DoubleVector &globalCov,MixtureGD &world);
void trainModel(Config& config,StatServer &ss,FeatureServer &fs,SegCluster& selectedSegments,
		DoubleVector & globalMean,DoubleVector &globalCov,MixtureGD &world,TrainCfg &trainCfg);
void trainModelStream(Config& config,MixtureServer &ms,StatServer &ss,FeatureServer **fsTab,SegCluster** segTab,double *weightTab,unsigned long nbStream,
		DoubleVector & globalMean,DoubleVector &globalCov,MixtureGD* &world,TrainCfg &trainCfg);
void trainModel(Config& config,MixtureServer &ms,StatServer &ss,FeatureServer **fsTab,SegCluster** segTab,double *weightTab,unsigned long nbStream,
		DoubleVector & globalMean,DoubleVector &globalCov,MixtureGD* &world);
		

//  Function for adapting a model by MLLR
DoubleSquareMatrix computeMLLR (MixtureGD &inM,MixtureGD& outM,unsigned long frameCount, Config &config);


#endif //!defined(ALIZE_TrainTools_h)

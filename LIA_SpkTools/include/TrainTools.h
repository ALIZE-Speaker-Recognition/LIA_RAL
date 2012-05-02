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

#if !defined(ALIZE_TrainTools_h)
#define ALIZE_TrainTools_h

#if defined(_WIN32)
#if defined(LIA_SPKTOOLS_EXPORTS)
#define LIA_SPKTOOLS_API __declspec(dllexport)
#else
#define LIA_SPKTOOLS_API __declspec(dllimport)
#endif
#else
#define LIA_SPKTOOLS_API
#endif

#include <alize.h>
#include "liatools.h"

using namespace alize;
using namespace std;

class LIA_SPKTOOLS_API MAPCfg {
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


class LIA_SPKTOOLS_API TrainCfg{
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
LIA_SPKTOOLS_API void fuseModels(const MixtureGD &,unsigned long,const MixtureGD &,unsigned long,MixtureGD &);
LIA_SPKTOOLS_API void copyMixture(DistribGD &,DistribGD &);
LIA_SPKTOOLS_API unsigned long selectComponent(bool selectCompA[],XList & distribL, MixtureGD &inputM);
LIA_SPKTOOLS_API unsigned long selectComponent(bool selectCompA[],double wFactor,MixtureGD &inputM);
LIA_SPKTOOLS_API unsigned long selectComponent(bool selectCompA[],unsigned long nbTop,MixtureGD &inputM);
LIA_SPKTOOLS_API double reduceModel(bool selectCompA[],MixtureGD &inputM,MixtureGD &outputM);
LIA_SPKTOOLS_API void normalizeWeights(MixtureGD &outputM);

// gaussian/mixture fusion using less likelihood loss criterion
LIA_SPKTOOLS_API void gaussianFusion(const DistribGD &g1,double w1,const DistribGD & g2,double w2, DistribGD &res,double &w);
LIA_SPKTOOLS_API void mixtureFusion(const MixtureGD &mixt,DistribGD &res,double &wres);

// normalizeMixture() normmalizes  the mixture in order to fit
// the data distribution
// Usually used with a mean=0, cov=1 distribution
LIA_SPKTOOLS_API void normalizeMixture(MixtureGD &mixt,const DoubleVector &meanSignal,const DoubleVector &covSignal,bool zeroOne,
		      unsigned long nbIt, bool meanOnly,Config &config);
LIA_SPKTOOLS_API void normalizeMixture(MixtureGD &mixt,const DoubleVector &meanSignal,const DoubleVector &covSignal,Config &config);
// for compatibility reasons only, same function but the target is always D(mean=0,std=1)
LIA_SPKTOOLS_API void normalizeMixture(MixtureGD &mixt,Config &config);
LIA_SPKTOOLS_API void normalizeMixture(MixtureGD & mixt, MAPCfg & mapCfg, Config & config);

LIA_SPKTOOLS_API double likelihoodLoss(const DistribGD &g1,double w1,const DistribGD &g2,double w2);
LIA_SPKTOOLS_API double *outputMAPTransformation(MixtureGD & world, MixtureGD & clientMixture,
  Config & config);
//-------------------------------------------------------------------------
// MAP computation functions 
// the main MAP functions
LIA_SPKTOOLS_API void computeMAP(MixtureServer &ms,const MixtureGD& world,MixtureGD &client,unsigned long frameCount,Config &config);
LIA_SPKTOOLS_API void computeMAP(MixtureServer &ms,const MixtureGD& world,MixtureGD &client,unsigned long frameCount,MAPCfg &cfg);

// Model manipulation tools
LIA_SPKTOOLS_API void copyMean(MixtureGD & mixtS, MixtureGD & mixtD);
LIA_SPKTOOLS_API void copyVar(MixtureGD & mixtS, MixtureGD & mixtD);
LIA_SPKTOOLS_API void copyWeight(MixtureGD & mixtS, MixtureGD & mixtD);

//-------------------------------------------------------------------------
LIA_SPKTOOLS_API double setItParameter(double begin, double end, int nbIt, int it);

//-------------------------------------------------------------------------
LIA_SPKTOOLS_API void varianceControl(MixtureGD& model,double flooring,double ceiling,const DoubleVector &covSignal);

// **********************************************************************
// Mixture Initialization stuff
//***********************************************************************
// cov and mean initialization
LIA_SPKTOOLS_API unsigned long computeMeanCov(Config &config,FeatureServer **fsTab,SegCluster ** segTab,unsigned long nbStream,DoubleVector &mean,DoubleVector &cov);
LIA_SPKTOOLS_API unsigned long computeMeanCov(Config &config,FeatureServer &fs,SegCluster &seg,DoubleVector &mean,DoubleVector &cov);
LIA_SPKTOOLS_API void initialize01(unsigned long vectSize,DoubleVector &mean,DoubleVector &cov);
// Mixture initialisation, based on random picking of frames
LIA_SPKTOOLS_API MixtureGD &mixtureInit(MixtureServer &ms,FeatureServer &fs,MixtureGD &world,
	               SegCluster &selectedSegments,const DoubleVector &globalCov, Config& config);
LIA_SPKTOOLS_API MixtureGD &mixtureInit(MixtureServer &ms,FeatureServer &fs,MixtureGD &world,
	               SegCluster &selectedSegments,const DoubleVector &globalCov, Config& config,TrainCfg & trainCfg);
LIA_SPKTOOLS_API MixtureGD &mixtureInit(MixtureServer &ms,FeatureServer **fsTab, SegCluster **segTab,double *weightTab,unsigned long nbStream,MixtureGD &world,
		       const DoubleVector &globalCov, Config& config,TrainCfg &trainCfg);
LIA_SPKTOOLS_API MixtureGD &mixtureInit(MixtureServer &ms,FeatureServer **fsTab, SegCluster **segTab,double*weightTab,unsigned long nbStream,MixtureGD &world,
		       const DoubleVector &globalCov, Config& config);
// The main function for estimate a client model by bayesian adaptattion of a world model
// Using EM and MAP 
LIA_SPKTOOLS_API void adaptModel(Config& config,StatServer &ss,MixtureServer &ms,FeatureServer &fs,SegCluster& selectedSegments,MixtureGD &world,MixtureGD &clientMixture);//A.P.
LIA_SPKTOOLS_API void adaptModel(Config& config,StatServer &ss,MixtureServer &ms,FeatureServer &fs,SegCluster& selectedSegments,MixtureGD &world,MixtureGD &clientMixture, MAPCfg &cfg);
// ** New training algo based on a true EM/ML estimate of the training data before to apply MAP
LIA_SPKTOOLS_API void modelBasedadaptModel(Config& config,StatServer &ss,MixtureServer &ms,FeatureServer &fs,SegCluster& selectedSegments,
			  MixtureGD &aprioriModel,MixtureGD &clientMixture, MixtureGD &initModel);

		

		

//-------------------------------------------------------------------------
// TrainModel is the main function for model training
// It works on an initialized model and on a clmuster of segments
LIA_SPKTOOLS_API void trainModel(Config& config,StatServer &ss,FeatureServer &fs,SegCluster& selectedSegments,
		DoubleVector & globalMean,DoubleVector &globalCov,MixtureGD &world);
LIA_SPKTOOLS_API void trainModel(Config& config,StatServer &ss,FeatureServer &fs,SegCluster& selectedSegments,
		DoubleVector & globalMean,DoubleVector &globalCov,MixtureGD &world,TrainCfg &trainCfg);
LIA_SPKTOOLS_API void trainModelStream(Config& config,MixtureServer &ms,StatServer &ss,FeatureServer **fsTab,SegCluster** segTab,double *weightTab,unsigned long nbStream,
		DoubleVector & globalMean,DoubleVector &globalCov,MixtureGD* &world,TrainCfg &trainCfg);
LIA_SPKTOOLS_API void trainModel(Config& config,MixtureServer &ms,StatServer &ss,FeatureServer **fsTab,SegCluster** segTab,double *weightTab,unsigned long nbStream,
		DoubleVector & globalMean,DoubleVector &globalCov,MixtureGD* &world);
		

//  Function for adapting a model by MLLR
LIA_SPKTOOLS_API Matrix<double> computeMLLR (MixtureGD &inM,MixtureGD& outM,unsigned long frameCount, Config &config);


#endif //!defined(ALIZE_TrainTools_h)

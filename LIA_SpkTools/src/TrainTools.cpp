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

#if !defined(ALIZE_TrainTools_cpp)
#define ALIZE_TrainTools_cpp

#include <iostream>
#include <fstream>  // pour outFile
#include <cstdio>   // pour printf()
#include <cassert> // pour le debug pratique
#include <cmath>
#include <liatools.h>
#include <DoubleSquareMatrix.h>

// Small stuff (command line, small functions
TrainCfg::TrainCfg(Config &config){
 _initVarianceFlooring  = config.getParam("initVarianceFlooring").toDouble();    // variance control parameters - relative to global data variance 
 _initVarianceCeiling   = config.getParam("initVarianceCeiling").toDouble();     // (1.0 is 100 % of the global vriance, 0.2 is 20% of the global variance
 _finalVarianceFlooring = config.getParam("finalVarianceFlooring").toDouble();   // Each parameter varies from initial value to final value 
 _finalVarianceCeiling  = config.getParam("finalVarianceCeiling").toDouble();    
 _nbTrainIt = config.getParam("nbTrainIt").toLong();                             // number of  it 
 _baggedFrameProbability= config.getParam("baggedFrameProbability").toDouble();  // Defines the percentage of frames used for one IT
 if (config.existsParam("baggedFrameProbabilityInit")) _baggedFrameProbabilityInit =config.getParam("baggedFrameProbabilityInit").toDouble();  // TODO : SUPRESS (deprecated)
 else _baggedFrameProbabilityInit=0;// Defines the percentage of frames BY COMPONENT used for the init from scratch
 if (config.existsParam("normalizeModel")) _normalizeModel=config.getParam("normalizeModel").toBool();// normalize the world (at each iteration) or not (default)
 else _normalizeModel=false;
 if (_normalizeModel){
   _normalizeModelNbIt=1;
   if (config.existsParam("normalizeModelMeanOnly")) 
     _normalizeModelMeanOnly=config.getParam("normalizeModelMeanOnly").toBool(); 
   else _normalizeModelMeanOnly=false;
   if (_normalizeModelMeanOnly) 
     _normalizeModelNbIt=config.getParam("normalizeModelNbIt").toLong();
 }

 _componentReduction=false;
 if (config.existsParam("componentReduction")) _componentReduction=config.getParam("componentReduction").toBool();
  _targetDistribCount=0;
 if (_componentReduction){
   _targetDistribCount=config.getParam("targetMixtureDistribCount").toLong();
 }
}
void TrainCfg::showConfig(ostream & st=cout){
  st << "Training parameters "<<endl;
  st <<"nb training iterations:"<<_nbTrainIt<<endl;
  st <<"init variance floor:"<<_initVarianceFlooring<<" final floor:"<<_finalVarianceFlooring<<endl;
  st <<"init variance ceil:"<<_initVarianceCeiling<<" final ceil:"<<_finalVarianceCeiling<<endl;
  if (_normalizeModel){
    st <<"normaliseModel on";
    if (_normalizeModelMeanOnly){
      st <<" Mean Only mode, nbIt:"<<_normalizeModelNbIt<<endl;
    }
    else st<<endl;
  }
  else st <<"normalise model off"<<endl;
  if (_componentReduction)
    st <<"Component number reduction, final["<<_targetDistribCount<<"]"<<endl;
}
MAPCfg::MAPCfg(Config &config){
  // MAP Algorithm 
  _mean=false;_var=false;_weight=false;
  if (config.existsParam("meanAdapt")) _mean=config.getParam("meanAdapt").toBool();
  if (config.existsParam("varAdapt")) _var=config.getParam("varAdapt").toBool();
  if (config.existsParam("weightAdapt")) _weight=config.getParam("weightAdapt").toBool();
  if (!_var && !_mean && !_weight) cerr <<"MAP, no adaptation selected !"<<endl;
  _method=config.getParam("MAPAlgo");
  if ((_method=="MAPConst")||(_method=="MAPConst2")){
    if (_mean)_r[0]=config.getParam("MAPAlphaMean").toDouble();//0.75 a priori probability of the initModel			
    if (_var)_r[1]=config.getParam("MAPAlphaVar").toDouble();
    if (_weight)_r[2]=config.getParam("MAPAlphaWeight").toDouble();
  }	  		
  else 
    if ((_method=="MAPOccDep")||(_method=="MAPModelBased")|| (_method=="MLLR")){
      if (_mean) _r[0]= config.getParam("MAPRegFactorMean").toDouble();// r factor, r=14			
      if (_var)  _r[1]= config.getParam("MAPRegFactorVar").toDouble();// r factor, r=14
      if (_weight)_r[2]= config.getParam("MAPRegFactorWeight").toDouble();// r factor, r=14 
    } else cerr <<"mapAlgo["<<_method<<"] unknown"<<endl; 
  _normalizeModel=false;
  if (config.existsParam("normalizeModel")) _normalizeModel=config.getParam("normalizeModel").toBool();// normalize the world (at each iteration) or not (default)
  if (_normalizeModel){
    _normalizeModelNbIt=1;
    if (config.existsParam("normalizeModelMeanOnly")) 
      _normalizeModelMeanOnly=config.getParam("normalizeModelMeanOnly").toBool(); 
    else _normalizeModelMeanOnly=false;
    if (_normalizeModelMeanOnly) 
      _normalizeModelNbIt=config.getParam("normalizeModelNbIt").toLong();
  }
  else _normalizeModel=false;	
  _nbTrainIt = 1;
  if (config.existsParam("nbTrainIt"))  _nbTrainIt=config.getParam("nbTrainIt").toLong();                      // number of  it 
  _nbEmIt = 1;
  if (config.existsParam("nbEmIt"))  _nbEmIt=config.getParam("nbEmIt").toLong();                      // number of EM/ML  it for modelbased map
  if (config.existsParam("baggedFrameProbability")) 
      _baggedFrameProbability= config.getParam("baggedFrameProbability").toDouble();
      else  _baggedFrameProbability=1.0;// Defines the percentage of frames used for one IT
}
void MAPCfg::showConfig(ostream & st=cout){
  st << "MAP Algo "<<_method<<endl;
  if (_mean) st<<"Mean adaption, param["<<_r[0]<<"]"<<endl;
  if (_var) st<<"Variance adaption, param["<<_r[1]<<"]"<<endl;
  if (_weight) st<<"Weight adaptation, param["<<_r[2]<<"]"<<endl;
  st <<"Nb training iterations["<<_nbTrainIt<<"]"<<endl;
  if (_normalizeModel){
    st <<"NormaliseModel on";
    if (_normalizeModelMeanOnly){
      st <<" Mean Only mode, nbIt["<<_normalizeModelNbIt<<"]"<<endl;
    }
    else st<<endl;
  }
  else st <<"Normalise model off"<<endl;
}

// Mixture and ditrib tools

void copyMixture(DistribGD &s,DistribGD &d){
  unsigned long vectSize=s.getVectSize();
  if (vectSize!=d.getVectSize()) throw Exception("vectSize != for copyMixture" , __FILE__, __LINE__);
  for (unsigned long c=0;c<vectSize;c++){
    d.setMean(s.getMean(c),c);
    d.setCov(s.getCov(c),c);
  }
  d.computeAll();
}
unsigned long selectComponent(bool selectCompA[],XList & distribL, MixtureGD &inputM){
  unsigned long nbInputDistrib=inputM.getDistribCount();  
  for (unsigned long idx=0;idx<nbInputDistrib;idx++) selectCompA[idx]=true;
  for (unsigned long idx=0;idx<nbInputDistrib;idx++) selectCompA[idx]=true;
  unsigned long nbOutputDistrib=nbInputDistrib;
  for (unsigned long idxS=0;idxS<distribL.getLineCount();idxS++){
    selectCompA[distribL.getLine(idxS).getElement(0).toLong()]=false;
    nbOutputDistrib--;
  } 
  return nbOutputDistrib;
}
unsigned long selectComponent(bool selectCompA[],double wFactor,MixtureGD &inputM){
  unsigned long nbInputDistrib=inputM.getDistribCount(); 
  for (unsigned long idx=0;idx<nbInputDistrib;idx++) selectCompA[idx]=true;
  unsigned long nbOutputDistrib=nbInputDistrib;
  for (unsigned long idx=0;idx<nbInputDistrib;idx++)
    if ((selectCompA[idx])&&(inputM.weight(idx)<wFactor)){
      selectCompA[idx]=false;
      nbOutputDistrib--;
    } 
  return nbOutputDistrib;
}
unsigned long selectComponent(bool selectCompA[],unsigned long nbTop,MixtureGD &inputM){
  unsigned long nbInputDistrib=inputM.getDistribCount(); 
  for (unsigned long idx=0;idx<nbInputDistrib;idx++) selectCompA[idx]=false;
  TabWeight tabW(inputM);
  for (unsigned long i=0; i<nbTop;i++){
    selectCompA[tabW.getDistrib(i)]=true;
  }
  return nbTop;
}
double reduceModel(bool selectCompA[],MixtureGD &inputM,MixtureGD &outputM){
  unsigned long nbInputDistrib=inputM.getDistribCount(); 
  unsigned long outputIdx=0;
  double totW=0;
  for (unsigned idx=0;idx<nbInputDistrib;idx++){
    if (selectCompA[idx]){
      DistribGD &s=inputM.getDistrib(idx);
      DistribGD &d=outputM.getDistrib(outputIdx);
      copyMixture(s,d);
      outputM.weight(outputIdx)=inputM.weight(idx);
      totW+=outputM.weight(outputIdx);
      outputIdx++;
    }
  }
  return totW;
}

void normalizeWeights(MixtureGD &outputM){
  double totW=0;
  for (unsigned long idx=0;idx<outputM.getDistribCount();idx++)
    totW+=outputM.weight(idx);
  for (unsigned long idx=0;idx<outputM.getDistribCount();idx++)
    outputM.weight(idx)/=totW;
}

// TRAINING AND ADAPTING FUNCTIONS
//*******************************

// fuseModels. this function builds a GMM using two statistic estimations, represented by two models and the number of frames
// the life of this function should be very short, as a +operator on the stat accumulator will replace it
void fuseModels(const MixtureGD & model1,unsigned long nbFrameM1,const MixtureGD & model2,unsigned long nbFrameM2,MixtureGD &result){
  unsigned long vectSize=model1.getVectSize();  
  if (vectSize!=model2.getVectSize()) throw Exception("Feature vector size should be the same" , __FILE__, __LINE__);
  unsigned long nbDistrib=model1.getDistribCount();
  if (nbDistrib!=model2.getDistribCount()) throw Exception("Number of components should be the same" , __FILE__, __LINE__);
  for (unsigned long idx=0;idx<nbDistrib;idx++){
    DistribGD & d1=model1.getDistrib(idx);
    DistribGD & d2=model2.getDistrib(idx);
    DistribGD & dr=result.getDistrib(idx);
    double totOcc1=nbFrameM1*model1.weight(idx);
    double totOcc2=nbFrameM2*model2.weight(idx);
    for (unsigned long c=0;c<vectSize;c++){
      dr.setMean(((d1.getMean(c)*totOcc1)+(d2.getMean(c)*totOcc2))/(totOcc1+totOcc2),c);
      dr.setCov(((d1.getCov(c)*totOcc1)+(d2.getCov(c)*totOcc2))/(totOcc1+totOcc2),c);// TO BE CHANGED !!!!!!!
    }
    dr.computeAll();
    result.weight(idx)=((model1.weight(idx)*nbFrameM1)+(model2.weight(idx)*nbFrameM2))/(nbFrameM1+nbFrameM2);
  }
}

// gaussian/mixture fusion using less likelihood loss criterion
void gaussianFusion(const DistribGD &g1,double w1,const DistribGD & g2,double w2, DistribGD &res,double &w)
{
 unsigned long vectSize = g1.getVectSize();
 double a1=w1/(w1+w2);
 double a2=1.0-a1;
 for (unsigned long k=0;k<vectSize;k++)
 {
  double d=g1.getMean(k)-g2.getMean(k);
  res.setCov(a1*g1.getCov(k)+a2*g2.getCov(k)+a1*a2*d*d,k);
  res.setMean((a1*g1.getMean(k))+(a2*g2.getMean(k)),k);
 }
 res.computeAll();
 w=w1+w2;
}
void mixtureFusion(const MixtureGD &mixt,DistribGD &res,double &wres){
  unsigned long distribCount = mixt.getDistribCount();
  
  DistribGD& tmp=mixt.getDistrib(0);
  double wtmp=mixt.weight(0);
  // used only for monocomponent mixture
  res=tmp; 
  wres=wtmp;
  for (unsigned long i=1;i<distribCount;i++){
   gaussianFusion(mixt.getDistrib(i),mixt.weight(i),res,wtmp,res,wres);
   wtmp=wres;
  }	
}
// normalizeMixture() normmalizes  the mixture to fit
// the data distriobution
// Usually used with a mean=0, cov=1 distribution
void normalizeMixture(MixtureGD &mixt,const DoubleVector &meanSignal,const DoubleVector &covSignal,bool zeroOne,
		      unsigned long nbIt, bool meanOnly,Config &config){
  unsigned long vectSize = mixt.getVectSize();
  unsigned long distribCount = mixt.getDistribCount();
  double wtmp;
  MixtureServer ms(config);
  DistribGD &tmp=ms.createDistribGD();
  for (unsigned long it=0;it<nbIt;it++){ // NbIt used only for mean only option
    mixtureFusion(mixt,tmp,wtmp); // Compute in tmp the mono gaussian ML corresponding to the initial mixture
    for (unsigned long c=0;c<distribCount;c++){ // normalize the components
      DistribGD &d=mixt.getDistrib(c);
      for (unsigned long i=0;i<vectSize;i++){
		double newMean=d.getMean(i)-tmp.getMean(i);
		newMean/=sqrt(tmp.getCov(i));
		if (!zeroOne) {// the target is not N(m=0,std=1){
	  		newMean*=covSignal[i];
	  		newMean+=meanSignal[i];
		}
		d.setMean(newMean,i);
		if (!meanOnly){
	  		double newCov=d.getCov(i)/tmp.getCov(i);
	  		if (zeroOne==false) newCov*=covSignal[i]; // the target is not N(m=0,std=1)
	  		d.setCov(newCov,i);
		}
      }
      d.computeAll();
    }
  }
}
void normalizeMixture(MixtureGD &mixt,const DoubleVector &meanSignal,const DoubleVector &covSignal,bool zeroOne,TrainCfg &trainCfg,Config & config){
  normalizeMixture(mixt,meanSignal,covSignal,zeroOne,trainCfg.getNormalizeModelNbIt(),trainCfg.getNormalizeModelMeanOnly(),config);
}
void normalizeMixture(MixtureGD &mixt,const DoubleVector &meanSignal,const DoubleVector covSignal,bool zeroOne,Config &config){
  TrainCfg trainCfg(config);
  normalizeMixture(mixt,meanSignal,covSignal,zeroOne,trainCfg,config);
}

// for compatibility reasons only, same function but the target is always D(mean=0,std=1)
void normalizeMixture(MixtureGD &mixt,TrainCfg &trainCfg,Config &config){
  DoubleVector fake1,fake2;
  normalizeMixture(mixt,fake1,fake2,true,trainCfg,config);
}
void normalizeMixture(MixtureGD &mixt,Config &config){
  TrainCfg trainCfg(config);
  DoubleVector fake1,fake2;
  normalizeMixture(mixt,fake1,fake2,true,trainCfg,config);
}
void normalizeMixture(MixtureGD &mixt,MAPCfg &mapCfg,Config &config){
  DoubleVector fake1,fake2;
  normalizeMixture(mixt,fake1,fake2,true,mapCfg.getNormalizeModelNbIt(),mapCfg.getNormalizeModelMeanOnly(),config);
}
double likelihoodLoss(const DistribGD &g1,double w1,const DistribGD & g2,double w2)
{
  unsigned long vectSize = g1.getVectSize();
  double a1=w1/(w1+w2);
  double a2=1.0-a1;
  double d1=1.0;
  double d2=1.0;
  for (unsigned long k=0;k<vectSize;k++){
    double d=g1.getMean(k)-g2.getMean(k);
    double var=a1*g1.getCov(k)+a2*g2.getCov(k)+a1*a2*d*d;
    d1*=var/g1.getCov(k);
    d2*=var/g2.getCov(k);	
  }
  return (0.5*(w1*log(d1)+w2*log(d2)));
}

// Bayesian Adaptation functions (MAP)
//Direct Mean Only interpolation: constant apriori for each component
void computeMAPConst(const MixtureGD& initModel, MixtureGD &client, MAPCfg &cfg)
{
  unsigned long vectSize = initModel.getVectSize();
  unsigned long distribCount = initModel.getDistribCount();
  Config config;
  MixtureServer ms(config);
  MixtureGD & tmp=ms.duplicateMixture(initModel,DUPL_DISTRIB);
  if (cfg.getMeanAdapt()){
    double alpha=cfg.getMeanAlpha();
    for ( unsigned long indC=0; indC < distribCount; indC++) {
      DistribGD& t = tmp.getDistrib(indC);            // A  priori data for a component (from initModel model)
      DistribGD& c  = client.getDistrib(indC);        // Estimate statistics for the component
      for (unsigned long coef=0;coef<vectSize;coef++){
	double res=(alpha*t.getMean(coef)) +((1-alpha)*c.getMean(coef));
	t.setMean(res, coef);
      }
    }
  }
  if (cfg.getVarAdapt()){
    //    double alpha=cfg.varAlpha();
    // TODO
  } 
  if (cfg.getWeightAdapt()){
    //double alpha=cfg.weightAlpha();
    // TODO
  } 
  client=tmp;	    
}
/************************************************************
* Calcul MAP avec une constante
* 
* Origine : Corinne Fredouille
************************************************************/
void computeMAPConst2(const MixtureGD& initModel, MixtureGD &client,MAPCfg &cfg )
{
  unsigned long vectSize = initModel.getVectSize();
  unsigned long distribCount = initModel.getDistribCount();
  Config config;
  MixtureServer ms(config);
  MixtureGD & tmp=ms.duplicateMixture(initModel,DUPL_DISTRIB);
  //MixtureGD tmp(Id,vectSize,distribCount);
  //tmp=initModel;
  if (cfg.getMeanAdapt()){
    double alpha=cfg.getMeanAlpha();
    for ( unsigned long indC=0; indC < distribCount; indC++) {
      DistribGD& t = tmp.getDistrib(indC); // A  priori data for a component (from initModel model)
      DistribGD& c  = client.getDistrib(indC);  // Estimate statistics for the component
      for (unsigned long coef=0;coef<vectSize;coef++){
	double res=((alpha*tmp.weight(indC)*t.getMean(coef)) +
		    ((1-alpha)*client.weight(indC)*c.getMean(coef)))/(tmp.weight(indC)*alpha+client.weight(indC)*(1-alpha));
	t.setMean(res, coef);
      }
    } 
  }
  if (cfg.getVarAdapt()){
    //	  double alpha=cfg.varAlpha();
    // TODO
  }   
  if (cfg.getWeightAdapt()){
    //double alpha=cfg.weightAlpha();
    // TODO
  }
  client=tmp;	    
}

void copyMean(MixtureGD & mixtS, MixtureGD & mixtD){
  for ( unsigned long indC=0; indC <mixtS.getDistribCount(); indC++) {
    DistribGD& s = mixtS.getDistrib(indC);                 
    DistribGD& d = mixtD.getDistrib(indC);             
    for (unsigned long coef=0;coef<mixtS.getVectSize();coef++)
      d.setMean(s.getMean(coef), coef);;
    }
  mixtD.computeAll(); 
}
void copyVar(MixtureGD & mixtS, MixtureGD & mixtD){
  for ( unsigned long indC=0; indC <mixtS.getDistribCount(); indC++) {
    DistribGD& s = mixtS.getDistrib(indC);                 
    DistribGD& d = mixtD.getDistrib(indC);             
    for (unsigned long coef=0;coef<mixtS.getVectSize();coef++)
      d.setCov(s.getCov(coef), coef);;
    }
  mixtD.computeAll(); 
}
void copyWeight(MixtureGD & mixtS, MixtureGD & mixtD){
  for ( unsigned long indC=0; indC <mixtS.getDistribCount(); indC++) 
    mixtD.weight(indC)=mixtS.weight(indC);
}

// Classical LIMSI or MIT mean only MAP adaptation
void computeMAPOccDep(const MixtureGD& initModel, MixtureGD &client, MAPCfg &cfg,double frameCount)
{
  unsigned long vectSize = initModel.getVectSize();
  unsigned long distribCount = initModel.getDistribCount();
  if (verbose) cfg.showConfig();
  // TODO : suppress the need of a temporary mixture server
  Config config;
  MixtureServer ms(config);
  MixtureGD & tmp=ms.duplicateMixture(initModel,DUPL_DISTRIB);
  //
  if ((cfg.getMeanAdapt()) || (cfg.getVarAdapt()))
    for ( unsigned long indC=0; indC < distribCount; indC++) {
      DistribGD& t = tmp.getDistrib(indC);                  // A  priori data for a component (from initModel model)
      DistribGD& c = client.getDistrib(indC);               // Estimate statistics for the component
      DistribGD& w = initModel.getDistrib(indC);                // The initModel
      double alpha= (client.weight(indC)*frameCount);       // occupation for the component
      if (cfg.getMeanAdapt()){
	double alphaMean=alpha/(alpha+cfg.getMeanReg());				          // weight for the component
	for (unsigned long coef=0;coef<vectSize;coef++){
	  double res=((1-alphaMean)*w.getMean(coef)) +(alphaMean*c.getMean(coef));
	  t.setMean(res, coef);
	}
      }
      if (cfg.getVarAdapt()){
	double alphaVar=alpha/(alpha+cfg.getVarReg());	
	for (unsigned long coef=0;coef<vectSize;coef++){
	  double res=((1-alphaVar)*w.getCov(coef)) +(alphaVar*c.getCov(coef))
	    +(((1-alphaVar)*alphaVar)*pow(w.getMean(coef)-c.getMean(coef),2));
	  t.setCov(res, coef);
	}  
      }
    } 
  if (cfg.getWeightAdapt()){
    double sumWeight=0.0;
    for ( unsigned long indC=0; indC < distribCount; indC++) {
      double alpha= (client.weight(indC)*frameCount); 
      double alphaWeight=alpha/(alpha+cfg.getWeightReg());	
      tmp.weight(indC)=alphaWeight*client.weight(indC)+((1-alphaWeight)*initModel.weight(indC));
      sumWeight+=tmp.weight(indC);
    }
    for ( unsigned long indC=0; indC < distribCount; indC++) tmp.weight(indC)/=sumWeight;
  }
  tmp.computeAll();
  client=tmp;	  
}
// Classical LIMSI or MIT mean only MAP adaptation but based on ML estimate of the training data
void computeModelBasedMAPOccDep(const MixtureGD& initModel, MixtureGD &client, MAPCfg &cfg,double frameCount)
{
  unsigned long vectSize = initModel.getVectSize();
  unsigned long distribCount = initModel.getDistribCount();
  if (verbose) {
    cout << "modelBased MAP Adaptation"<<endl;
    cfg.showConfig();
  }
  // TODO : suppress the need of a temporary mixture server
  Config config;
  MixtureServer ms(config);
  MixtureGD & tmp=ms.duplicateMixture(initModel,DUPL_DISTRIB);

  //
  if ((cfg.getMeanAdapt()) || (cfg.getVarAdapt()))
    for ( unsigned long indC=0; indC < distribCount; indC++) {
      DistribGD& t = tmp.getDistrib(indC);                  // A  priori data for a component (from initModel model)
      DistribGD& c = client.getDistrib(indC);               // Estimate statistics for the component (Full ML in this case)
      DistribGD& w = initModel.getDistrib(indC);            // The initModel
      double alpha= (client.weight(indC)*frameCount);       // occupation for the component
      if (cfg.getMeanAdapt()){
	double alphaMean=alpha/(alpha+cfg.getMeanReg());				          // weight for the component
	for (unsigned long coef=0;coef<vectSize;coef++){
	  double res=((1-alphaMean)*w.getMean(coef)) +(alphaMean*c.getMean(coef));
	  t.setMean(res, coef);
	}
      }
      if (cfg.getVarAdapt()){
	double alphaVar=alpha/(alpha+cfg.getVarReg());	
	for (unsigned long coef=0;coef<vectSize;coef++){
	  double res=((1-alphaVar)*w.getCov(coef)) +(alphaVar*c.getCov(coef))
	    +(((1-alphaVar)*alphaVar)*pow(w.getMean(coef)-c.getMean(coef),2));
	  t.setCov(res, coef);
	}  
      }
    } 
  if (cfg.getWeightAdapt()){
    double sumWeight=0.0;
    for ( unsigned long indC=0; indC < distribCount; indC++) {
      double alpha= (client.weight(indC)*frameCount); 
      double alphaWeight=alpha/(alpha+cfg.getWeightReg());	
      tmp.weight(indC)=alphaWeight*client.weight(indC)+((1-alphaWeight)*initModel.weight(indC));
      sumWeight+=tmp.weight(indC);
    }
    for ( unsigned long indC=0; indC < distribCount; indC++) tmp.weight(indC)/=sumWeight;
  }
  tmp.computeAll();
  client=tmp;	  
}
// Main MAP functions
void computeMAP(MixtureServer &ms,const MixtureGD& initModel,MixtureGD &client,unsigned long frameCount,Config &config)
{
  MAPCfg cfg(config);                          // Get all the needed parameters in the config
  computeMAP(ms,initModel,client,frameCount,cfg);
}
void computeMAP(MixtureServer &ms,const MixtureGD& initModel,MixtureGD &client,unsigned long frameCount,MAPCfg &cfg)
{
  if (debug) cfg.showConfig(cout);
  if (cfg.getMethod()=="MAPConst")
    computeMAPConst(initModel,client,cfg);
  else if (cfg.getMethod()=="MAPOccDep")
    computeMAPOccDep(initModel,client,cfg,frameCount);
  else if (cfg.getMethod()=="MAPConst2")
    computeMAPConst2(initModel,client,cfg);
  else if (cfg.getMethod()=="MAPModelBased")
    computeModelBasedMAPOccDep(initModel,client,cfg,frameCount);
 else cerr <<"mapAlgo["<<cfg.getMethod()<<"] unknown, No adaptation will be perform"<<endl; 
}
//-------------------------------------------------------------------------
double setItParameter(double begin, double end, int nbIt, int it){
  if (nbIt<2) return begin;
  double itVal=((begin-end)/((double)nbIt-1));
  return begin-(itVal*it);
}
//-------------------------------------------------------------------------
// varianceControl for GMM models (florring and ceiling)
void varianceControl(MixtureGD& model,double flooring,double ceiling,const DoubleVector &covSignal)
{
  unsigned long vectSize     = model.getVectSize();
  unsigned long distribCount = model.getDistribCount();
  unsigned long cptFlooring  = 0;
  unsigned long cptCeiling   = 0;

  for (unsigned long c=0; c<distribCount; c++){
    DistribGD& d = model.getDistrib(c); 
    for (unsigned long v=0; v<vectSize; v++){
      double cov = d.getCov(v); 
      if (cov <= (flooring*covSignal[v])) { cov = flooring*covSignal[v]; cptFlooring++; }
      if (cov >= (ceiling*covSignal[v]))  { cov = ceiling*covSignal[v];  cptCeiling++; }
      d.setCov(cov, v);
    } 
    d.computeAll();  
  }
  if (verbose) cout << " total variance flooring = " << cptFlooring
		    << " ceiling = " << cptCeiling << endl;
  
}
 
// **********************************************************************
// Mixture Initialization stuff
//***********************************************************************
// cov and mean initialization
unsigned long computeMeanCov(Config &config,FeatureServer **fsTab,SegCluster ** segTab,unsigned long nbStream,DoubleVector &mean,DoubleVector &cov){
  initialize01(fsTab[0]->getVectSize(),mean,cov); // Set the good size to mean and cov
  FrameAccGD frameAcc;
  for (unsigned long stream=0;stream<nbStream;stream++)
    accumulateStatFrame(frameAcc,*fsTab[stream],*segTab[stream],config);
  mean=frameAcc.getMeanVect(); 
  cov=frameAcc.getCovVect();
  unsigned long nbFrame= frameAcc.getCount();
  return nbFrame;
}
unsigned long computeMeanCov(Config &config,FeatureServer &fs,SegCluster &seg,DoubleVector &mean,DoubleVector &cov){
  initialize01(fs.getVectSize(),mean,cov); // Set the good size to mean and cov
  FrameAccGD frameAcc;
  accumulateStatFrame(frameAcc,fs,seg,config);
  mean=frameAcc.getMeanVect(); 
  cov=frameAcc.getCovVect();
  unsigned long nbFrame= frameAcc.getCount();
  return nbFrame;
}
void initialize01(unsigned long vectSize,DoubleVector &mean,DoubleVector &cov){
  for(unsigned long i=0;i<vectSize;i++){
    mean.addValue(0);
    cov.addValue(1);
  }
}
// Mixture initialisation, based on random picking of frames
MixtureGD &mixtureInit(MixtureServer &ms,FeatureServer &fs,MixtureGD &world,
	               SegCluster &selectedSegments,const DoubleVector &globalCov, Config& config,TrainCfg & trainCfg)
{
  unsigned long vectSize = world.getVectSize();	
  unsigned long minimumLength=3;
  unsigned long maximumLength=7;
  if (config.existsParam("baggedMinimalLength")) minimumLength=config.getParam("baggedMinimalLength").toLong();
  if (config.existsParam("baggedMaximalLength")) maximumLength=config.getParam("baggedMaximalLength").toLong();
  if (verbose) cout <<"baggedMinimalLength["<< minimumLength<<"] MaximalLength["<<maximumLength<<"]"<<endl;
  unsigned long distribCount = world.getDistribCount(); 
  FrameAcc** frameAcc=new FrameAcc* [distribCount];             // Initialise the frame accumulator tab
  for (unsigned long indg=0;indg<distribCount;indg++) {         // for the component indg
    frameAcc[indg]= new FrameAccGD();                           // Create the accumulators
    frameAcc[indg]->reset();                                    // Reset it
  } ;
  double baggedProba=trainCfg.getBaggedFrameProbabilityInit()/distribCount;
  unsigned long nbBaggedIt=1;
  if (baggedProba>1){
    nbBaggedIt=(unsigned long) baggedProba+1;
    baggedProba/=nbBaggedIt;
  }
  for (unsigned long baggedIt=0;baggedIt<nbBaggedIt;baggedIt++){
    for(unsigned long indg=0;indg<distribCount;indg++){
      SegServer segServer;                                        // Create a local segment server 
      SegCluster & baggedFramesCluster=segServer.createCluster(1,"",""); // Create the cluster for describing the selected frames
      srand((indg+1)*(baggedIt+1)); 
      baggedSegments(selectedSegments,baggedFramesCluster,
		     baggedProba,minimumLength,maximumLength);
      accumulateStatFrame(*frameAcc[indg],fs,baggedFramesCluster,config); // Accumulate picked frames for each component
    }
  }
  for (unsigned long indg=0;indg<distribCount;indg++)  {        // For each component
    DistribGD& d = world.getDistrib(indg);                      // Get it
    DoubleVector mean=frameAcc[indg]->getMeanVect();            // Get the mean
    for (unsigned long c=0;c<vectSize;c++){                     // copy it
      d.setCov(globalCov[c], c);
      d.setMean(mean[c], c);
    }
    d.computeAll();
  }  
  world.equalizeWeights();                                      // set weight = 1/distribCount for each distrib 
  if (verbose){
    cout << "Initialize model"<<endl;
    for (unsigned long indg=0;indg<distribCount;indg++) 
      cout <<"Nb Frame for mean["<<indg<<"] init =["<<frameAcc[indg]->getCount()<<"]"<<endl;
  }
  for (unsigned long indg=0;indg<distribCount;indg++)           // Free the memory
    delete frameAcc[indg];
  delete [] frameAcc;
  return world;
}


// Mixture initialisation, based on random picking of frames
// Multiple data frames
MixtureGD &mixtureInit(MixtureServer &ms,FeatureServer **fsTab,SegCluster **segTab,double *weightTab,unsigned long nbStream,MixtureGD &world,
	               const DoubleVector &globalCov, Config& config,TrainCfg & trainCfg)
{
  if (debug|| verbose) cout << "begin of mixtureInit"<<endl;	
  unsigned long minimumLength=3;
  unsigned long maximumLength=7;
  if (config.existsParam("baggedMinimalLength")) minimumLength=config.getParam("baggedMinimalLength").toLong();
  if (config.existsParam("baggedMaximalLength")) maximumLength=config.getParam("baggedMaximalLength").toLong();
  if (verbose) cout <<"baggedMinimalLength["<< minimumLength<<"] MaximalLength["<<maximumLength<<"]"<<endl;
  unsigned long vectSize = world.getVectSize();	
  unsigned long distribCount = world.getDistribCount(); 
  FrameAcc** frameAcc=new FrameAcc* [distribCount];             // Initialise the frame accumulator tab
  for (unsigned long indg=0;indg<distribCount;indg++) {         // for the component indg
    frameAcc[indg]= new FrameAccGD();                           // Create the accumulators
    frameAcc[indg]->reset();                                    // Reset it
  }
  unsigned long nbTotalFrame=0;
  for (unsigned long stream=0;stream<nbStream;stream++)
    nbTotalFrame+=(unsigned long) ((double)totalFrame(*segTab[stream])*weightTab[stream]);
  double nbFrameToSelect=50;
  if (config.existsParam("nbFrameToSelect")) nbFrameToSelect=config.getParam("nbFrameToSelect").toLong();
  
  if (verbose) cout <<"mixtureInit, total frames with stream weighting["<<nbTotalFrame<<"] total frame to select by component["<<nbFrameToSelect<<"]"<<endl;                                                        
  // Init the bagged probablity vector : one probability by frame
  DoubleVector baggedProbaA(nbStream,nbStream);
  ULongVector baggedItA(nbStream,nbStream);
  for (unsigned long stream=0;stream<nbStream;stream++){         // For each input stream
      baggedProbaA[stream]=(nbFrameToSelect*weightTab[stream])/(double)totalFrame(*segTab[stream]);
      baggedItA[stream]=1;
      double baggedTmp=baggedProbaA[stream];
      while (baggedTmp>1){
	  	baggedItA[stream]++;
	  	baggedTmp/=baggedProbaA[stream]/baggedItA[stream];
      	}
      baggedProbaA[stream]=baggedTmp;
  }
  if (debug || verbose) cout << "bagged proba init ok"<<endl;
  /*
  for (unsigned long stream=0;stream<nbStream;stream++){         // For each input stream
      for (unsigned long baggedIt=0;baggedIt<baggedItA[stream];baggedIt++){
	  	if (debug || verbose) cout <<"Init Stream["<<stream<<"] Bagged by comp="<<baggedProbaA[stream]<<" It["<<baggedIt<<"]"<<endl;
	  	// Initialise the array of bagged segment clusters, one by component
	  	RefVector <SegCluster> baggedA(distribCount); // Will be deleted automatically at the end of each loop it
	  	SegServer segServer; // Create a local segment server, one for all the bagged clusters
	  	for (unsigned long indG=0;indG<distribCount;indG++)
	    	  baggedA.addObject(segServer.createCluster(indG,"",""),indG); // Create the cluster for describing the selected frames
	  	// Compute the bagging - random selection, for all components
	  	srand(((stream+1)*100)+(baggedIt+1)); // 100 to be sure to not have 2 equal init
	  	baggedSegments(*segTab[stream],baggedA,baggedProbaA[stream],minimumLength,maximumLength);
	  	// Accumulate the statistics
	  	for(unsigned long indG=0;indG<distribCount;indG++){ // For each component
	  		if (debug || verbose) cout << "AccumulateStat stream["<<stream<<"] component["<<indG<<"]"<<endl;
	      	accumulateStatFrame(*frameAcc[indG],*fsTab[stream],baggedA[indG],config); // Accumulate picked frames for the component and the input stream
	  	}
      }
  }*/
  for (unsigned long stream=0;stream<nbStream;stream++){         // For each input stream
      for (unsigned long baggedIt=0;baggedIt<baggedItA[stream];baggedIt++){
	  	if (debug || verbose) cout <<"Init Stream["<<stream<<"] Bagged by comp="<<baggedProbaA[stream]<<" It["<<baggedIt<<"]"<<endl;
	  	// Initialise the array of bagged segment cluster, who will contain all the segments for all the streams
	  	SegServer segServer; // Create a local segment server, one for all the bagged clusters
	    SegCluster & baggedSeg=segServer.createCluster(1,"",""); // Create the cluster for describing the selected frames
	  	// Compute the bagging - random selection, for all components
	  	srand(((stream+1)*100)+(baggedIt+1)); // 100 to be sure to not have 2 equal init
	  	baggedSegments(*segTab[stream],baggedSeg,distribCount,baggedProbaA[stream],minimumLength,maximumLength);
	  	// Accumulate the statistics
	  	Seg* seg;                                                     // reset the reader at the begin of the input stream
        baggedSeg.rewind();      
  		while((seg=baggedSeg.getSeg())!=NULL)                  // For each of the selected segments
    		accumulateStatFrame(*frameAcc[seg->labelCode()],*fsTab[stream],seg,config);
      }
  }
  for (unsigned long indg=0;indg<distribCount;indg++)  {        // For each component
    DistribGD& d = world.getDistrib(indg);                      // Get it
    DoubleVector mean=frameAcc[indg]->getMeanVect();            // Get the mean
    for (unsigned long c=0;c<vectSize;c++){                     // copy it
      d.setCov(globalCov[c], c);
      d.setMean(mean[c], c);
    }
    d.computeAll();
  }  
  
  world.equalizeWeights();                                      // set weight = 1/distribCount for each distrib 
  if (verbose){
    cout << "Initialize model"<<endl;
    for (unsigned long indg=0;indg<distribCount;indg++) 
      cout <<"Nb Frame for mean["<<indg<<"] init =["<<frameAcc[indg]->getCount()<<"]"<<endl;
  }
  for (unsigned long indg=0;indg<distribCount;indg++)           // Free the memory
    delete frameAcc[indg];
  delete [] frameAcc;
  return world;
}
MixtureGD &mixtureInit(MixtureServer &ms,FeatureServer &fs,MixtureGD &world,
	               SegCluster &selectedSegments,const DoubleVector &globalCov, Config& config){
  TrainCfg trainCfg(config);
  return mixtureInit(ms,fs,world,selectedSegments,globalCov,config,trainCfg);
}
MixtureGD &mixtureInit(MixtureServer &ms,FeatureServer **fsTab, SegCluster **segTab,double *weightTab,unsigned long nbStream,MixtureGD &world,
		       const DoubleVector &globalCov, Config& config){
  TrainCfg trainCfg(config);
  return mixtureInit(ms,fsTab,segTab,weightTab,nbStream,world,globalCov,config,trainCfg);
}



//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//Compute the MLLR transformation for means of a input Model and store it in an output Model

// remplacement de toutes les matrices carrées par des Matrices rectangulaires
// Dimension de W (vectsize , vectsize+1)
// Dimension de Z (vectsize, vectsize+1)
// Dimension de G (vectsize +1 , vectsize +1) (en fait on a vectisze Matrices G : une par distribution)

Matrix <double> computeMLLR (MixtureGD &inM,MixtureGD& outM,unsigned long frameCount, Config &config){
	StatServer ss(config);


if(verbose)cout<<"MLLR begins"<<endl;
	unsigned long distribCount=inM.getDistribCount();	        // number of gaussians in the Mixture	
	unsigned long vectSize= inM.getVectSize();				// dimension of the gaussians
	Matrix <double> W (vectSize, vectSize+1); 				// MMLR transformation matrix 
	Matrix <double> Z(vectSize, vectSize+1);				// matrix  (one unused colunn)
	W.setAllValues(0);
	Z.setAllValues(0);
	
	// Define and initialize the G matrices
	DoubleSquareMatrix **G;
	G=new DoubleSquareMatrix*[vectSize];
	for (unsigned i=0;i<vectSize;i++) {
		G[i]= new DoubleSquareMatrix();
		G[i]->setSize(vectSize+1);
		G[i]->setAllValues(0);
		}
	if(verboseLevel>=1 )cout<<"Initialization complete"<<endl<<"Begining the Statistic accumulation"<<endl;

	for (unsigned long j=0; j<distribCount;j++){
		DistribGD & distri = inM.getDistrib(j);
		DistribGD & distriMl = outM.getDistrib(j);
		double occ=outM.weight(j)*frameCount;
		DoubleVector mean(vectSize+1,1);
		DoubleVector &cov=distri.getCovVect();

		//Initialization of the mean doubleVector (offset and other values)
		mean[0]=1;
		mean.addValue(distri.getMeanVect());

		// Compute the Z matrix
		for (unsigned long p=0; p<vectSize;p++){
			double occFrame=distriMl.getMean(p)*occ;
			for (unsigned long q=0; q<vectSize+1;q++) 
				Z(p,q)+=occFrame*mean[q]/cov[p];
		}
		//Compute the G matrix
		for (unsigned long p=0; p<vectSize;p++){
			for (unsigned long q=0; q<vectSize+1; q++){	
				for (unsigned long r=0; r<vectSize+1; r++) (*G[p])(q,r) += occ*mean[q]*mean[r]/cov[p];
			}
		}
	}	
	
	if(verboseLevel >=1)cout<<"G and Z matrix calculated"<<endl;
	
	//Compute the W matrix using the G inverse and Z
	for (unsigned long l=0; l<vectSize;l++) {
		DoubleSquareMatrix G_inv(vectSize+1);
		(*G[l]).invert(G_inv);
		for (unsigned long c=0; c<vectSize+1;c++) 
			for( unsigned long k=0;k<vectSize+1;k++) 
				W(l,c) += G_inv(c,k)*Z(l,k);
		
	}
	if(verboseLevel >= 1)cout<<"W complete"<<endl<<"Computing new means"<<endl;
	
	//Compute the mean vector of outModel
	for (unsigned long j=0; j<distribCount;j++){
		DoubleVector meanOut(vectSize, vectSize);
		DistribGD & distri = inM.getDistrib(j);
		DoubleVector &mean= distri.getMeanVect();
		for (unsigned long i=0;i<vectSize;i++) meanOut[i]=W(i,0);	//initialize the mean vector with the offset
		for (unsigned long i=0;i<vectSize;i++) {
			for (unsigned long k=0;k<vectSize;k++) meanOut[i]+=W(i,k+1)*mean[k];
			}
		outM.getDistrib(j).setMeanVect(meanOut);	//update the mean vector of the out Mixture
	}
	
	if(verbose) cout<<"MLLR finished, New Means computed"<<endl;
		
	copyVar(inM, outM);
	copyWeight(inM,outM);

return W;
}


// The main function for estimate a client model by bayesian adaptattion of a aprioriModel model
// Using EM and MAP 
void adaptModel(Config& config,StatServer &ss,MixtureServer &ms,FeatureServer &fs,SegCluster& selectedSegments,
		MixtureGD &aprioriModel,MixtureGD &clientMixture,MAPCfg &mapCfg){ 
  if (verbose) mapCfg.showConfig();
  // if (verboseLevel>1) cout << "Mean LLK Init = " << meanLikelihood(ss,fs,clientMixture,selectedSegments,config)<< endl;    
  for (unsigned long trainIt=0; trainIt<mapCfg.getNbTrainIt(); trainIt++){    // Begin the initial adaptation loop (with bagged frames)
    MixtureStat &emAcc=ss.createAndStoreMixtureStat(clientMixture);             // Create a statistic accumulator using the curent model
    SegServer segServer;                                                      // Create a local segment server 
    SegCluster & baggedFramesCluster=segServer.createCluster(1,"","");        // Create the cluster for describing the selected frames
    baggedSegments(selectedSegments,baggedFramesCluster,mapCfg.getBaggedFrameProbability());
    emAcc.resetEM();
    srand(trainIt);
    double llkPreviousIt=accumulateStatEM(ss,fs,emAcc,baggedFramesCluster,config);// Accumulate the EM statistics
    clientMixture = emAcc.getEM();                                             // Get the EM estimate	 
    unsigned long frameCount=(unsigned long) emAcc.getEMFeatureCount();
    llkPreviousIt=llkPreviousIt/(double) frameCount;
    if (verbose)cout <<"ML (partial) estimate it["<<trainIt<<"] (take care, it corresponds to the previous it,0 means init likelihood) = "
		     <<llkPreviousIt<<endl;
    if(mapCfg.getMethod()=="MLLR") {
	    //DoubleSquareMatrix W=computeMLLR(aprioriModel,clientMixture,frameCount,config);			//MLLR adaptation
	    Matrix <double> W=computeMLLR(aprioriModel,clientMixture,frameCount,config);			//MLLR adaptation
            //outputDSM(W,clientMixture,config) ;
	    String MLLR_matrix = "MLLR_matrix.mat";
	    W.save(MLLR_matrix, config);
     } 
    else {
                computeMAP(ms,aprioriModel,clientMixture,frameCount,config);               // Bayesian Adaptation client=MAP(aprioriModel,client)
	}
    if (mapCfg.getNormalizeModel()) normalizeMixture(clientMixture,mapCfg,config);    // Normalize/fit the model if needed
    ss.deleteMixtureStat(emAcc);  
    if (verboseLevel>2) cout << "Likelihood on all frames= "<< meanLikelihood(ss,fs,clientMixture, selectedSegments,config) << endl;        
    }    
  if (verboseLevel>1) cout << "Final likelihood on all frames= "<< meanLikelihood(ss,fs,clientMixture, selectedSegments,config) << endl;                               
  if (debug) cout <<"adaptModel nb distrib:"<<ms.getDistribCount() <<"nb mixt:"<<ms.getMixtureCount()<<endl;
}

void adaptModel(Config& config,StatServer &ss,MixtureServer &ms,FeatureServer &fs,SegCluster& selectedSegments,
		MixtureGD &aprioriModel,MixtureGD &clientMixture){ 
  MAPCfg mapCfg(config);
  adaptModel(config,ss,ms,fs,selectedSegments,aprioriModel,clientMixture,mapCfg);
}


//------------------------------------------------------------------------
// ** New training algo based on a true EM/ML estimate of the training data before to apply MAP
void modelBasedadaptModel(Config& config,StatServer &ss,MixtureServer &ms,FeatureServer &fs,SegCluster& selectedSegments,
			  MixtureGD &aprioriModel,MixtureGD &clientMixture, MixtureGD &initModel){
  MAPCfg mapCfg(config);
  if (verbose){
    cout <<"Model adaptation based on true EM/ML estimate of training data"<<endl;
    mapCfg.showConfig();
  }
  MixtureServer msTmp(config);
  MixtureGD & data=msTmp.duplicateMixture(initModel,DUPL_DISTRIB);
  
  unsigned long totalFrameCount=totalFrame(selectedSegments);
  //if (verboseLevel>1) cout << "Mean LLK Init = " << meanLikelihood(ss,fs,data,selectedSegments,config)<< endl;    
  for (unsigned long emIt=0; emIt<mapCfg.getNbEmIt(); emIt++){                 // begin the true EM/ML estimate of the adpatation data 
    MixtureStat &emAcc=ss.createAndStoreMixtureStat(data);                     // Create a statistic accumulator using the curent model
    SegServer segServer;                                                       // Create a local segment server
    SegCluster & baggedFramesCluster=segServer.createCluster(1,"","");         // Create the cluster for describing the selected frames
    baggedSegments(selectedSegments,baggedFramesCluster,mapCfg.getBaggedFrameProbability());
    emAcc.resetEM();
    double llkPreviousIt=accumulateStatEM(ss,fs,emAcc,baggedFramesCluster,config);// Accumulate the EM statistics
    data = emAcc.getEM();                                                         // Get the EM estimate	 
    unsigned long frameCount=(unsigned long) emAcc.getEMFeatureCount();
    llkPreviousIt=llkPreviousIt/(double) frameCount;
    if (verbose)cout <<"ML (partial) estimate it["<<emIt<<"] (take care, it corresponds to the previous it,0 means init likelihood) = "
		     <<llkPreviousIt<<endl;
    ss.deleteMixtureStat(emAcc);
  }
  unsigned long modelNbComp =aprioriModel.getDistribCount();
  // Begin the estimation of the statistic using the EM/ML model of the adaptation data
  unsigned long vectSize=fs.getVectSize();                               // Complete log likelihood of the adaptation data given the apriori model
  for  (unsigned long idxModel=0;idxModel<modelNbComp;idxModel++){ // Initialize the client mixture
    DistribGD &c=clientMixture.getDistrib(idxModel);
    for (unsigned long idxC=0;idxC<vectSize;idxC++){
      c.setMean(0.00,idxC);
      c.setCov(0.00,idxC);
    }
  }
  DoubleVector apProbaTot(modelNbComp,modelNbComp); 
  apProbaTot.setAllValues(0.0);
  for (unsigned long idxData=0;idxData<data.getDistribCount();idxData++){
    if (debug) cout <<"Distrib Data["<<idxData<<"]"<<endl;
    DistribGD &d=data.getDistrib(idxData);
    double totLk=0.0;   // Likelihood of the current data component given the apriori model
    DoubleVector apProba(modelNbComp,modelNbComp); 
    apProba.setAllValues(0.0);
    for  (unsigned long idxModel=0;idxModel<modelNbComp;idxModel++){
      if (debug) cout <<"Distrib A Priori model["<<idxModel<<"]"<<endl;
      DistribGD &m=aprioriModel.getDistrib(idxModel);
      apProba[idxModel]=aprioriModel.weight(idxModel)*likelihoodGD(d,m);
      totLk+=apProba[idxModel]; 
    }
    for  (unsigned long idxModel=0;idxModel<modelNbComp;idxModel++){
      DistribGD &c=clientMixture.getDistrib(idxModel);
      apProba[idxModel]/=totLk;
      for (unsigned long idxC=0;idxC<vectSize;idxC++){
	c.setMean(c.getMean(idxC)+(d.getMean(idxC)*apProba[idxModel]*data.weight(idxData)),idxC);
	c.setCov(c.getCov(idxC)+((d.getMean(idxC)*d.getMean(idxC))*apProba[idxModel]*data.weight(idxData)),idxC);
      }
      apProbaTot[idxModel]+=apProba[idxModel]*data.weight(idxData);
    }
  }
  for (unsigned long idxModel=0;idxModel<modelNbComp;idxModel++){
    DistribGD &c=clientMixture.getDistrib(idxModel);
    for (unsigned long idxC=0;idxC<vectSize;idxC++){
      c.setMean(c.getMean(idxC)/apProbaTot[idxModel],idxC);
      //	c.setCov(c.getCov(idxC),idxC);
    }
  }
  for (unsigned long idxModel=0;idxModel<modelNbComp;idxModel++){
    DistribGD &c=clientMixture.getDistrib(idxModel);
    c.computeAll();
  }
  //
  computeMAP(ms,aprioriModel,clientMixture,totalFrameCount,config);                      // Bayesian Adaptation client=MAP(aprioriModel,client)
  if (mapCfg.getNormalizeModel()) normalizeMixture(clientMixture,mapCfg,config);    // Normalize/fit the model if needed 
  if (verboseLevel>2) cout << "Likelihood on all frames= "<< meanLikelihood(ss,fs,clientMixture, selectedSegments,config) << endl;       
  if (verboseLevel>1) cout << "Final likelihood on all frames= "<< meanLikelihood(ss,fs,clientMixture, selectedSegments,config) << endl;  if (debug) cout <<"adaptModel nb distrib:"<<ms.getDistribCount() <<"nb mixt:"<<ms.getMixtureCount()<<endl;  
}
//-------------------------------------------------------------------------
void trainModel(Config& config,StatServer &ss,FeatureServer &fs,SegCluster& selectedSegments,
			 DoubleVector & globalMean,DoubleVector &globalCov,MixtureGD &world,TrainCfg &trainCfg){ 
  float varianceFlooring = (float)trainCfg.getInitVarFloor();
  float varianceCeiling  = (float)trainCfg.getInitVarCeil();   
  if (verbose) trainCfg.showConfig();
  if (verboseLevel>1) cout << "Train It : initial model: ll world = " <<  meanLikelihood(ss,fs,world,selectedSegments,config)<< endl;
  try{   	
    for (unsigned long trainIt=0; trainIt<trainCfg.getNbTrainIt(); trainIt++){                           // Begin the initial EM iteration loop 
      MixtureStat& emAcc=ss.createAndStoreMixtureStat(world);             // Create a statistic accumulator using the init model
      // Set the variance control parameters
      varianceFlooring = (float)setItParameter(trainCfg.getInitVarFloor(),trainCfg.getFinalVarFloor(),trainCfg.getNbTrainIt(),trainIt);
      varianceCeiling  = (float)setItParameter(trainCfg.getInitVarCeil(),trainCfg.getFinalVarCeil(), trainCfg.getNbTrainIt(), trainIt);
      if (verbose) cout << "Train it["<<trainIt<<"], Variance floor["<<varianceFlooring
			<<"] ceiling["<<varianceCeiling<<"]"<<endl; 
      SegServer segServer;                                        // Create a local segment server 
      SegCluster & baggedFramesCluster=segServer.createCluster(1,"",""); // Create the cluster for describing the selected frames
      srand(trainIt);    
      baggedSegments(selectedSegments,baggedFramesCluster,trainCfg.getBaggedFrameProbability());
      emAcc.resetEM();                                                                                    // EM stuff
      double llkPreviousIt=accumulateStatEM(ss,fs,emAcc,baggedFramesCluster,config);                        // Compute EM statistics
      llkPreviousIt=llkPreviousIt/(double) emAcc.getEMFeatureCount();
      world = emAcc.getEM();                                                                             // Get the EM estimate
      varianceControl(world,varianceFlooring,varianceCeiling,globalCov);
      if (trainCfg.getNormalizeModel()) normalizeMixture(world,config);	
      if (verbose) cout << "Partial Train it["<<trainIt<<"] (take care, it corresponds to the previous it,0 means init likelihood) ="<<llkPreviousIt<<" Nb Frames="<<emAcc.getEMFeatureCount()<<endl;
      if (verboseLevel>2) cout << "Train it["<<trainIt <<"] Complete LLK=[" << meanLikelihood(ss,fs,world,selectedSegments,config)<< "]"<<endl;       
      ss.deleteMixtureStat(emAcc);
    }
    if (verboseLevel>1) cout << " Final ll world = " << meanLikelihood(ss,fs,world,selectedSegments,config)<< endl;
  }
  catch (Exception& e){
    cout << e.toString() << endl;
  }
}

//-------------------------------------------------------------------------
// stream version, should become the only one !!!
void trainModelStream(Config& config,MixtureServer &ms,StatServer &ss,FeatureServer **fsTab,SegCluster** segTab,
		double *weightTab,unsigned long nbStream,
		DoubleVector & globalMean,DoubleVector &globalCov,MixtureGD* &world,TrainCfg &trainCfg){ 
  unsigned long initialDistribCount=world->getDistribCount();
  double varianceFlooring = trainCfg.getInitVarFloor();
  double varianceCeiling  = trainCfg.getInitVarCeil(); 
  unsigned long initRand=0; //initRand allows to run several TrainWorld using different set of (random) data
  if (config.existsParam("initRand")) initRand=config.getParam("initRand").toLong();  
  unsigned long minimumLength=3;
  unsigned long maximumLength=7;
  if (config.existsParam("baggedMinimalLength")) minimumLength=config.getParam("baggedMinimalLength").toLong();
  if (config.existsParam("baggedMaximalLength")) maximumLength=config.getParam("baggedMaximalLength").toLong();
  if (verbose) trainCfg.showConfig();
  if (verbose) cout <<"baggedMinimalLength["<< minimumLength<<"] MaximalLength["<<maximumLength<<"]"<<endl;
  if (verboseLevel>1) cout <<"LLK Initial model= " <<  meanLikelihood(ss,fsTab,segTab,nbStream,*world,config)<< endl;
  try{   	
    for (unsigned long trainIt=0; trainIt<trainCfg.getNbTrainIt(); trainIt++){             // Begin the initial EM iteration loop;
      MixtureStat& emAcc=ss.createAndStoreMixtureStat(*world);                             // Create a statistic accumulator using the init model
      // Set the variance control parameters
      varianceFlooring = setItParameter(trainCfg.getInitVarFloor(),trainCfg.getFinalVarFloor(),trainCfg.getNbTrainIt(),trainIt);
      varianceCeiling  = setItParameter(trainCfg.getInitVarCeil(),trainCfg.getFinalVarCeil(), trainCfg.getNbTrainIt(), trainIt);
      if (verbose) cout << "Train it["<<trainIt<<"], Variance floor["<<varianceFlooring
			<<"] ceiling["<<varianceCeiling<<"]"<<endl; 
      unsigned long nbTotalFrame=0;
      for (unsigned long stream=0;stream<nbStream;stream++) nbTotalFrame+=(unsigned long) ((double)totalFrame(*segTab[stream])*weightTab[stream]);
      double nbFrameToSelect=trainCfg.getBaggedFrameProbability()*nbTotalFrame;
      if (verbose) cout <<"total frames with stream weighting["<<nbTotalFrame<<"] total frame to select ["<<(unsigned long) nbFrameToSelect<<"]"<<endl;
      emAcc.resetEM();  
      double llkPreviousIt=0;
      for(unsigned long stream=0;stream<nbStream;stream++){
	unsigned long nbBaggedIt=1;
	double baggedProba=(nbFrameToSelect*weightTab[stream])/(double)totalFrame(*segTab[stream]);
	if (baggedProba>1){
	  nbBaggedIt=(unsigned long) baggedProba+1;
	  baggedProba/=nbBaggedIt;
	}
	for (unsigned long baggedIt=0;baggedIt<nbBaggedIt;baggedIt++){
	  if (debug || verboseLevel) cout <<"Stream["<<stream<<"] Bagged="<<baggedProba<<endl;
	  SegServer segServer;                                        // Create a local segment server 
	  SegCluster & baggedFramesCluster=segServer.createCluster(1,"",""); // Create the cluster for describing the selected frames
      srand(((trainIt+1+initRand)*200)+(((stream+1)*20)+(baggedIt+1)));   
	  baggedSegments(*segTab[stream],baggedFramesCluster,baggedProba,minimumLength,maximumLength);
	  llkPreviousIt+=accumulateStatEM(ss,*fsTab[stream],emAcc,baggedFramesCluster,config);                        // Compute EM statistics
	}
      }
      llkPreviousIt=llkPreviousIt/(double) emAcc.getEMFeatureCount();
      (*world) = emAcc.getEM(); // Get the EM estimate
      varianceControl(*world,varianceFlooring,varianceCeiling,globalCov);
      if (trainCfg.getComponentReduction()){
	// Tableau dynamique instanciÃ© sans pointeur (ancien code)
	//bool selectCompA[world->getDistribCount()] TODO netoyage memoire
	bool * selectCompA = new bool[world->getDistribCount()];

	double diff=(initialDistribCount-trainCfg.getTargetDistribCount())/(double) trainCfg.getNbTrainIt();
	unsigned long nbTop=initialDistribCount-(unsigned long)((double)(trainIt+1)*diff);
	if (trainIt== (trainCfg.getNbTrainIt()-1)) nbTop=trainCfg.getTargetDistribCount();
	if (nbTop<world->getDistribCount())
	  {
	    if (verbose) cout <<"target number of d:"<<nbTop<<endl;
	    unsigned long nbOutputDistrib=selectComponent(selectCompA,nbTop,*world);
	    if (verbose) cout <<" suppress distrib nb initial:"<<world->getDistribCount()<<" nb final :"<<nbOutputDistrib;
	    MixtureGD &outputM=ms.createMixtureGD(nbOutputDistrib);
	    double totWeights=reduceModel(selectCompA,*world,outputM);
	    if (verbose) cout <<" Total weights["<<totWeights<<"] normalized to 1"<<endl;
	    normalizeWeights(outputM);
	    world=&outputM;
	  }
      }
      if (trainCfg.getNormalizeModel()) normalizeMixture(*world,trainCfg,config);	
 
      if (verbose) cout << "Partial Train it["<<trainIt<<"] (take care, it corresponds to the previous it,0 means init likelihood) ="<<llkPreviousIt<<" Nb Frames="<<(unsigned long) emAcc.getEMFeatureCount()<<endl;
       if (verboseLevel>2) cout << "Train it["<<trainIt <<"] Complete LLK=[" << meanLikelihood(ss,fsTab,segTab,nbStream,*world,config)<< "]"<<endl;    
      ss.deleteMixtureStat(emAcc);

    }
    if (verboseLevel>1) cout << " Final ll world = " << meanLikelihood(ss,fsTab,segTab,nbStream,*world,config)<< endl;
  }
  catch (Exception& e){
    cout << e.toString() << endl;
  }
}
void trainModel(Config& config,StatServer &ss,FeatureServer &fs,SegCluster& selectedSegments,
			 DoubleVector & globalMean,DoubleVector &globalCov,MixtureGD &world){
  TrainCfg trainCfg(config);
  trainModel(config,ss,fs,selectedSegments,globalMean,globalCov,world,trainCfg) ;
}
void trainModel(Config& config,MixtureServer &ms,StatServer &ss,FeatureServer **fsTab,SegCluster** segTab,double *weightTab,unsigned long nbStream,
		DoubleVector & globalMean,DoubleVector &globalCov,MixtureGD* &world){ 
  TrainCfg trainCfg(config);
  if (verbose) trainCfg.showConfig();
  trainModelStream(config,ms,ss,fsTab,segTab,weightTab,nbStream,globalMean,globalCov,world,trainCfg);
}
#endif //!defined(ALIZE_TrainTools_cpp)

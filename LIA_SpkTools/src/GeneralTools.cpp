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

#if !defined(ALIZE_GeneralTools_cpp)
#define ALIZE_GeneralTools_cpp

#include <iostream>
#include <fstream>  // pour outFile
#include <cstdio>   // pour printf()
#include <cassert> // pour le debug pratique
#include <cmath>

#include "liatools.h"

//----------------------- ComputeTest stuff --------------------------------------

// --------------------------------------------------------------------------------
unsigned long TabClientLine::nbClientLine(){
  return nbModelsLine; 
}
// --------------------------------------------------------------------------------
String& TabClientLine::getClientName(unsigned long nClient){
  //First element of line is the segment
  return pline->getElement(nClient+1);
}
// --------------------------------------------------------------------------------
MixtureGD& TabClientLine::getClientModel(unsigned long nClient){
  return *tabModel[nClient];   
} 

// --------------------------------------------------------------------------------
// Should be added - one step initialisation with a complete mixture server - TODO
/*unsigned long initialize(MixtureServer * m,Config &config, unsigned long maxModel=CST_MAX_MODEL_LINE)
  { 
  return 0; // Number of client models pre-loaded
  }
*/
// --------------------------------------------------------------------------------
unsigned long TabClientLine::loadLine(XLine* linep,String label,bool useDefaultClientModel,bool byUserRep){
  pline=linep; //Take care; the XLine is not stored but only a pointor on it !
  String * pModelName;  
  nbModelsLine=0;  
  bool useLabel=(label!=""); 
  pline->getElement(1);
  while (((pModelName=pline->getElement()) != NULL) && (nbModelsLine<nbMaxModelLine)){ 
    String modelName;
    if (useLabel)    
      modelName=*pModelName+"_"+label;                            // Use ID+label mode for multipl e models by speaker
    else  
      modelName=*pModelName;                                      // Normal mode
    int indexModel=ms->getMixtureIndex(modelName);
    if (indexModel==-1){                                          // the model is not already loaded in the mixture server
      String file;
      if (byUserRep) file=*pModelName+"/"+modelName;
      else file=modelName;
      try{  
	MixtureGD & m=ms->loadMixtureGD(file);
	tabModel[nbModelsLine]=&m; 
	nbModelsLine++;                                          // A model is added to the line
      }
      catch (Exception& e){ 
	cout << "WARNING, model ["<<modelName <<"] not found "<< endl;
 
	if ((label!="")&&(useDefaultClientModel)){
	  cout << "Trying to use speaker general model["<<*pModelName<<"]"<<endl;
	  MixtureGD & m=ms->loadMixtureGD(*pModelName);
	  ms->setMixtureId(m,modelName); 
	  tabModel[nbModelsLine]=&m;
	  if (verbose) cout << "model loaded"<<endl;
	  nbModelsLine++;                                          // A model is added to the line
	}
      }
    }
    else {
      tabModel[nbModelsLine]=&(ms->getMixtureGD(indexModel));     // the model is already in the memory and is added to the line
      nbModelsLine++;
    }
  }
  if (nbModelsLine>=nbMaxModelLine) cerr << "TabClientLine::loadLine() - nb model line > nbMaxModelLine)" << endl;	
  return nbModelsLine;
}

// --------------------------------------------------------------------------------
TabClientLine::TabClientLine(MixtureServer & m, Config &config, unsigned long maxModel)
{
  nbMaxModelLine=maxModel;
  ms=&m;
  conf=&config;
  tabModel=new MixtureGD*[nbMaxModelLine];
  pline=NULL;
}
//---------------------------------------------------------------------------------	
TabClientLine::~TabClientLine() {delete []tabModel;}
//---------------------------------------------------------------------------------
//----------------------- end of TabClientLine definition -------------------------


//---------------------------------------------------------------------------------
//-----------------------  TabHisto definition ----------------------------------

  Histo & TabHisto::getHistoFromVect(unsigned long n) {
	Object & tmp=_tabHisto.getObject(n);
	return (Histo&)tmp;
  }

void TabHisto::accumulateValueInTab(const String &id,double score){
	long n=_id.getIndex(id); // was unsigned long N.S. 23/09/05
	if (n==-1) {
		if (debug) cerr << id << " unknown :";
		_id.addElement(id);
		_nb++;
		n=_id.getIndex(id);
		if (debug) cerr << id << " now index "<< n << endl;
	}
	getHistoFromVect(n).accumulateValue(score);
}

void TabHisto::computeHistoInTab(const String &id){
  return getHisto(id).computeHisto();
}
unsigned long TabHisto::getIndex(const String &id){
  long n=_id.getIndex(id);
  if (n==-1)throw Exception("out of array", __FILE__, __LINE__);
  else return n;
}
String& TabHisto::getId(unsigned long n){
  if (n>_nb)  throw Exception("out of array" , __FILE__, __LINE__);
  else return _id.getElement(n);
}
Histo & TabHisto::getHisto(const String & id)
{
	unsigned long n=getIndex(id);
	if (n>_nb)  throw Exception("out of array" , __FILE__, __LINE__);
	else {unsigned long n=getIndex(id);
	return getHistoFromVect(n);}
}
//---------------------------------------------------------------------------------
//----------------------- end of TabHisto definition ----------------------------


//-----------------------------------------------------------------------------------
// Compute the entropy from an Histo
double computeEntropy(Histo & hist)
{
	// integral is done by summing with a step (weighted by log(x) function)
	double prob=0.0;
	double entropy=0.0;
	double bound=hist.lowerBound(0);
	double step=0.001;
	while (bound < hist.higherBound(hist.size()-1)) {
		prob=hist(bound); // get density
		if (prob<1e-20) {entropy+=0.0;} // take care of lim x->0 x*log(x) (=0)
		else {entropy+=-prob*step*log(prob);} //get the area of the bin to get P(x), indeed P(x)=sum_x(p(x))
		bound+=step;
		//cerr << "prob: "<<prob<<" entropy: "<<entropy<<endl;
	}
return entropy;
}

//-----------------------------------------------------------------------------------
// Compute mean from an Histo
double computeMean(Histo & hist)
{
	// integral is done by summing with a step (weighted by log(x) function)
	double prob=0.0;
	double mean=0.0;
	double bound=hist.lowerBound(0);
	double step=0.001;
	while (bound < hist.higherBound(hist.size()-1)) {
		prob=step*hist(bound);
		mean+=bound*prob; //get the area of the bin to get P(x), indeed P(x)=sum_x(p(x))
		bound+=step;
		//cerr << "prob: "<<prob<<" mean: "<<mean<<endl;
	}
return mean;
}


// --------------------------------------------------------
// Decision function ...
long setDecision(double LLRClient, double decisionThreshold)
{
  if (LLRClient>=decisionThreshold) return 1; else return  0;
}

//---------------------------------------------------------------------------------
//-----------------------  ScoreAccum definition ----------------------------------

void ScoreAccum::addAndAccumulate(const String &id,double score, unsigned long nbFrames){
  double value=score*(double)nbFrames;
  long n=_id.getIndex(id);
  if (n==-1){
    _score.addValue(value);
    _nbFrame.addValue(nbFrames);
    _id.addElement(id);
    _nb++;
  }
  else{
    _score[n]+=value;
    _nbFrame[n]+=nbFrames;
  }
}
double ScoreAccum::getScore(const String &id){
  return getScore(getIndex(id));
}
unsigned long ScoreAccum::getIndex(const String &id){
  long n=_id.getIndex(id);
  if (n==-1)throw Exception("out of array", __FILE__, __LINE__);
  else return n;
}
double  ScoreAccum::getScore(unsigned long n){
  if (_nbFrame[n]==0)throw Exception("No score accumulated", __FILE__, __LINE__); // Could also make and exception index out of bound
  double ret=_score[n]; // Could return an exception index out of bound
  ret/=_nbFrame[n];
  return ret;
}
String& ScoreAccum::getId(unsigned long n){
  if (n>_nb)  throw Exception("out of array" , __FILE__, __LINE__);
  else return _id.getElement(n);
}
//---------------------------------------------------------------------------------
//----------------------- end of ScoreAccum definition ----------------------------

//
// Component selection functions
int _compF(const void * op1, const void *op2){
  if (((TabWeightElem*)op1)->weight>((TabWeightElem*)op2)->weight)
    return -1; else return 1;
}
void TabWeight::init(const MixtureGD &model,unsigned long topDistribs){
  _size=model.getDistribCount();
  _tab=new TabWeightElem[_size];
  _sortByWeight(model);
  _nbTop=topDistribs;
}
void TabWeight::init(const MixtureGD &model,double threshold){
  _size=model.getDistribCount();
  _tab=new TabWeightElem[_size];
  _sortByWeight(model);
  _nbTopDyn(threshold);
}
TabWeight::TabWeight(const MixtureGD &model){
  init(model,model.getDistribCount());
}
TabWeight::TabWeight(const MixtureGD &model,unsigned long topDistribs){
  init(model,topDistribs);
}
TabWeight::TabWeight(const MixtureGD &model,double threshold){
  init(model,threshold);
}

// Random picking of frames (bagging) functions, based on segment/cluster processing
// It is independent of the segment length
// Try to decrease the number of segments for fasting the world model training
// Used mainly in TrainTools.cpp 
// Author: JFB
// Just a function for selecting or not randomly a frame
bool baggedFrame(double baggedFrameProbability){
  // return (drand48()< baggedFrameProbability);
  double res= ((double)rand()/ (double)RAND_MAX);
  //cout << " baggedProba= " << baggedFrameProbability<<" res=" << res<<endl; //TODO revenir a drand48
return (res < baggedFrameProbability);
}

unsigned long correctedLength(unsigned long length,unsigned long &minimumLength,unsigned long&maximumLength){
  if (length<minimumLength) length=minimumLength;
  if (length>maximumLength) length=maximumLength;
  return length;
}

void baggedSegmentsConstraint(SegCluster &selectedSegments,SegCluster &baggedFrameSegment,double baggedProbability,
		    unsigned long minimumLength,unsigned long maximumLength){
	do{
		baggedSegments(selectedSegments,baggedFrameSegment,baggedProbability,minimumLength,maximumLength);
	}while(totalFrame(baggedFrameSegment) == 0);
}
// Works on a set of bagged clusters - only one reading of the
// segments and multiple selections, one by bagged cluster
void baggedSegments(SegCluster &selectedSegments,SegCluster &baggedSeg,unsigned long nbBagged,double & baggedProbability,
		    unsigned long minimumLength,unsigned long maximumLength){  
 if (debug) cout << "begin of baggedSegments !!!"<<endl;
  Seg* seg;                                                     // reset the reader at the begin of the input stream
  selectedSegments.rewind();      
  seg=selectedSegments.getSeg();
  bool end=(seg==NULL);
  unsigned long beginSeg=0,lengthSeg=0;
  if (!end){
    beginSeg=seg->begin();
    lengthSeg=seg->length();
  }
  while(!end){
    if (debug) cout << "bagged, current input seg ["<<beginSeg<<","<<lengthSeg<<"]"<<endl;
    unsigned long  verifyLength=correctedLength(lengthSeg,minimumLength,maximumLength);
    bool moveSeg=true;
    unsigned long length=0;
    if (lengthSeg<=verifyLength){
      moveSeg=true;
      length=lengthSeg;
      if (debug) cout <<"change seg"<<endl;
    }
    else{
      moveSeg=false;
      length=verifyLength;
    }
    // for all cluster in baggedA
    SegServer &segServerOutput=baggedSeg.getServer();
    if (length>0){
		for (unsigned long idx=0;idx<nbBagged;idx++) // For each component
	    	if(baggedFrame(baggedProbability)){
			Seg &newSeg=segServerOutput.createSeg(beginSeg,length,idx,seg->string(),seg->sourceName());       
			baggedSeg.add(newSeg);
			if (debug) cout << "bagged - Adding in bagged["<<idx<<"] the seg ["<<seg->sourceName()<<"]"<<newSeg.begin()<<" "<<newSeg.length()<<endl;   
	    	}
    }
    if (moveSeg){
	seg=selectedSegments.getSeg();
      end=(seg==NULL);
      if (!end){
	beginSeg=seg->begin();
	lengthSeg=seg->length();
      }
    }
    else{
      lengthSeg-=length;
      beginSeg+=length;
    }
  } 
  
  if ((debug) || (verboseLevel>3)){
    cout <<"Bagged segments"<<endl;
	showCluster(baggedSeg);
    }
  if (verbose){
    unsigned long total=totalFrame(selectedSegments);
	unsigned long selected=totalFrame(baggedSeg);
	double percent=(double)selected*100/(double) total;
	cout <<"Bagged segments, Initial frames["<<total<<"] Selected frames["<<selected<<"] % selected["<<percent<<"]"<<endl;
  }
}
/*
void baggedSegments(SegCluster &selectedSegments,RefVector<SegCluster> &baggedA,double & baggedProbability,
		    unsigned long minimumLength,unsigned long maximumLength){  
  Seg* seg;                                                     // reset the reader at the begin of the input stream
  selectedSegments.rewind();      
  seg=selectedSegments.getSeg();
  bool end=(seg==NULL);
  unsigned long beginSeg=0,lengthSeg=0;
  if (!end){
    beginSeg=seg->begin();
    lengthSeg=seg->length();
  }
  while(!end){
    if (debug) cout << "bagged, current input seg ["<<beginSeg<<","<<lengthSeg<<"]"<<endl;
    unsigned long  verifyLength=correctedLength(lengthSeg,minimumLength,maximumLength);
    bool moveSeg=true;
    unsigned long length=0;
    if (lengthSeg<=verifyLength){
      moveSeg=true;
      length=lengthSeg;
      if (debug) cout <<"change seg"<<endl;
    }
    else{
      moveSeg=false;
      length=verifyLength;
    }
    // for all cluster in baggedA
    if (length>0)
	for (unsigned long idx=0;idx<baggedA.size();idx++) // For each component
	    if(baggedFrame(baggedProbability)){
		SegServer &segServerOutput=baggedA[idx].getServer();
		Seg &newSeg=segServerOutput.createSeg(beginSeg,length,0,seg->string(),seg->sourceName());       
		baggedA[idx].add(newSeg);
		if (debug) cout << "bagged - Adding in bagged["<<idx<<"] the seg ["<<seg->sourceName()<<"]"<<newSeg.begin()<<" "<<newSeg.length()<<endl;   
	    }
    if (moveSeg){
	seg=selectedSegments.getSeg();
      end=(seg==NULL);
      if (!end){
	beginSeg=seg->begin();
	lengthSeg=seg->length();
      }
    }
    else{
      lengthSeg-=length;
      beginSeg+=length;
    }
  } 
  if ((debug) || (verboseLevel>3)){
    cout <<"Bagged segments"<<endl;
    for (unsigned long idx=0;idx<baggedA.size();idx++){
	cout << "Bagged cluster["<<idx<<"]"<<endl;
	showCluster(baggedA[idx]);
    }
  }
  if (verbose){
    unsigned long total=totalFrame(selectedSegments);
    for (unsigned long idx=0;idx<baggedA.size();idx++){
	unsigned long selected=totalFrame(baggedA[idx]);
	double percent=(double)selected*100/(double) total;
	cout <<"Bagged segments["<<idx<<"] Initial frames["<<total<<"] Selected frames["<<selected<<"] % selected["<<percent<<"]"<<endl;
    }
  }
}*/
void baggedSegments(SegCluster &selectedSegments,SegCluster &baggedFrameSegment,double baggedProbability,
		    unsigned long minimumLength,unsigned long maximumLength){  
  SegServer &segServerOutput=baggedFrameSegment.getServer();
  Seg* seg;                                                     // reset the reader at the begin of the input stream
  selectedSegments.rewind();      
  seg=selectedSegments.getSeg();
  bool end=(seg==NULL);
  unsigned long beginSeg=0,lengthSeg=0;
  if (!end){
    beginSeg=seg->begin();
    lengthSeg=seg->length();
  }
  while(!end){
    if (debug) cout << "bagged, current input seg ["<<beginSeg<<","<<lengthSeg<<"]"<<endl;
    unsigned long  verifyLength=correctedLength(lengthSeg,minimumLength,maximumLength);
    double segBaggedProbability=baggedProbability;
    bool moveSeg=true;
    unsigned long length=0;
    if (lengthSeg<=verifyLength){
      moveSeg=true;
      length=lengthSeg;
      if (debug) cout <<"change seg"<<endl;
    }
    else{
      moveSeg=false;
      length=verifyLength;
    }
    if ((length>0) &&(baggedFrame(segBaggedProbability))){
      Seg &newSeg=segServerOutput.createSeg(beginSeg,length,0,seg->string(),seg->sourceName());       
      baggedFrameSegment.add(newSeg);
      if (debug) cout << "bagged - Adding the seg ["<<seg->sourceName()<<"]"<<newSeg.begin()<<" "<<newSeg.length()<<endl;   
    }
    if (moveSeg){
      seg=selectedSegments.getSeg();
      end=(seg==NULL);
      if (!end){
	beginSeg=seg->begin();
	lengthSeg=seg->length();
      }
    }
    else{
      lengthSeg-=length;
      beginSeg+=length;
    }
  } 
  if ((debug) || (verboseLevel>3)){
    cout <<"Bagged segments"<<endl;
    showCluster(baggedFrameSegment);
  }
  if (verbose){
    unsigned long total=totalFrame(selectedSegments);
    unsigned long selected=totalFrame(baggedFrameSegment);
    double percent=(double)selected*100/(double) total;
    cout <<"Bagged segments, Initial frames["<<total<<"] Selected frames["<<selected<<"] % selected["<<percent<<"]"<<endl;
  }
}
// Same but returns both selected and unselected clusters
// Take care, both clusters should be created in the same server - NOT TESTED
void baggedSegments(SegCluster &selectedSegments,SegCluster &baggedSelected,SegCluster &baggedUnselected,double baggedProbability,
		    unsigned long minimumLength,unsigned long maximumLength){  
  SegServer &segServerOutput=baggedSelected.getServer();
  Seg* seg;                                                     // reset the reader at the begin of the input stream
  selectedSegments.rewind();      
  seg=selectedSegments.getSeg();
  bool end=(seg==NULL);
  unsigned long beginSeg=0,lengthSeg=0;
  if (!end){
    beginSeg=seg->begin();
    lengthSeg=seg->length();
  }
  while(!end){
    if (debug) cout << "bagged, current input seg ["<<beginSeg<<","<<lengthSeg<<"]"<<endl;
    unsigned long  verifyLength=correctedLength(lengthSeg,minimumLength,maximumLength);
    double segBaggedProbability=baggedProbability;
    bool moveSeg=true;
    unsigned long length=0;
    if (lengthSeg<=verifyLength){
      moveSeg=true;
      length=lengthSeg;
      if (debug) cout <<"change seg"<<endl;
    }
    else{
      moveSeg=false;
      length=verifyLength;
    }
    if (length>0){
    	Seg &newSeg=segServerOutput.createSeg(beginSeg,length,0,seg->string(),seg->sourceName());       
    	if(baggedFrame(segBaggedProbability)){
       	 baggedSelected	.add(newSeg);
      	 if (debug) cout << "baggedSelected - Adding the seg ["<<seg->sourceName()<<"]"<<newSeg.begin()<<" "<<newSeg.length()<<endl;   
    	}
    	else{
    	  baggedUnselected.add(newSeg);
      	  if (debug) cout << "baggedUnselected - Adding the seg ["<<seg->sourceName()<<"]"<<newSeg.begin()<<" "<<newSeg.length()<<endl;   
    	 }
    }
    if (moveSeg){
      seg=selectedSegments.getSeg();
      end=(seg==NULL);
      if (!end){
	beginSeg=seg->begin();
	lengthSeg=seg->length();
      }
    }
    else{
      lengthSeg-=length;
      beginSeg+=length;
    }
  } 
  if ((debug) || (verboseLevel>3)){
    cout <<"BaggedSelected"<<endl;
    showCluster(baggedSelected);
    cout <<"BaggedUnselected"<<endl;
    showCluster(baggedUnselected);
  }
  if (verboseLevel>1){
    unsigned long total=totalFrame(selectedSegments);
    unsigned long selected=totalFrame(baggedSelected);
    double percent=(double)selected*100/(double) total;
    cout <<"Bagged segments, Initial frames["<<total<<"] Selected frames["<<selected<<"] % selected["<<percent<<"]"<<endl;
  }
}

//-------------------------------------------------------------------------
//-- Compute the mean log likelihood for the Selected frames and a given model
double meanLikelihood(StatServer &ss,FeatureServer &fs,MixtureGD &model,unsigned long idxBeginFrame,unsigned long nbFrames,Config &config)  {
  MixtureStat &llkAcc=ss.createAndStoreMixtureStat(model);
  llkAcc.resetLLK();
  accumulateStatLLK(ss,fs,llkAcc,idxBeginFrame,nbFrames,config); 
  double llk=llkAcc.getMeanLLK();
  ss.deleteMixtureStat(llkAcc);
  return llk;
}
// one a cluster
double meanLikelihood(StatServer &ss,FeatureServer &fs,MixtureGD &model,SegCluster &selectedSegments,
		      Config &config){
  MixtureStat &llkAcc=ss.createAndStoreMixtureStat(model);
  llkAcc.resetLLK();
  accumulateStatLLK(ss,fs,llkAcc,selectedSegments,config);
  double llk=llkAcc.getMeanLLK();
  ss.deleteMixtureStat(llkAcc);
  return llk;
}
// on a set of input streams
double meanLikelihood(StatServer &ss,FeatureServer **fsTab,SegCluster **segTab,unsigned long nbStream,MixtureGD &model,Config &config){
  MixtureStat &llkAcc=ss.createAndStoreMixtureStat(model);
  llkAcc.resetLLK();
  for (unsigned long stream=0;stream<nbStream;stream++)
    accumulateStatLLK(ss,*fsTab[stream],llkAcc,*segTab[stream],config);
  double llk=llkAcc.getMeanLLK();
  ss.deleteMixtureStat(llkAcc);
  return llk;
}

//A.P.
double meanLikelihood(StatServer & ss, ObjectRefVector & FeatServ,
  ObjectRefVector & ClusterSeg, MixtureGD & model, DoubleVector & decision,
  Config & config)
{
  MixtureStat & llkAcc = ss.createAndStoreMixtureStat(model);
  llkAcc.resetLLK();
  for (unsigned long nbFs = 0; nbFs < FeatServ.size(); nbFs++)
    accumulateStatLLK(ss,
      (static_cast < FeatureServer & >(FeatServ.getObject(nbFs))), llkAcc,
      (static_cast < SegCluster & >(ClusterSeg.getObject(nbFs))),
      decision[nbFs], config);
  double llk = llkAcc.getMeanLLK();
  ss.deleteMixtureStat(llkAcc);
  return llk;
}


//-------------------------------------------------------------------------
//-- Compute the mean and cov of selected the data using segmental mode
void globalMeanCov (FeatureServer &fs,SegCluster &selectedSegments,FrameAcc & globalFrameAcc,Config &config)  {
  globalFrameAcc.reset();  
  accumulateStatFrame(globalFrameAcc,fs,selectedSegments,config);  
}
// On a complete feature stream
void globalMeanCov (FeatureServer &fs,FrameAcc & globalFrameAcc,Config &config)  {
  globalFrameAcc.reset();  
  accumulateStatFrame(globalFrameAcc,fs,0,fs.getFeatureCount(),config);  
}

// ----------------------------------------------------------------------------------------------------------
// Feature Warping giving a source(tab of histo, one by coeff) and a target distribution 
// for a segment and cluster (segment is the minimum time unit to perform this)
void computeWarp(Histo *histoT,Histo &destH,FeatureServer & fs,unsigned long begin, unsigned long length,Config &config) {  
  unsigned long vectsize=fs.getVectSize();                                       // Get the vect size (number of coeff)
  Feature f;
  fs.seekFeature(begin);
  for (unsigned long idxFrame=0;idxFrame<length;idxFrame++){                     // for all the features of the segment
    fs.readFeature(f,0);  
    // Get the feature;
    for (unsigned int i = 0; i < vectsize; i++){    // For each coeff
      f[i]=warping(f[i],histoT[i],destH);           // Apply the warping function
    }
    fs.writeFeature(f);
  }
}
// on a segment
void computeWarp(Histo *histoT,Histo &destH, FeatureServer & fs, Seg* seg,Config & config){
  unsigned long begin=seg->begin()+fs.getFirstFeatureIndexOfASource(seg->sourceName()); // Idx of the first frame of the current file in the feature server
  computeWarp(histoT,destH,fs,begin,seg->length(),config);      
}
// on a cluster
void computeWarp(Histo *histoT,Histo &destH, FeatureServer & fs, SegCluster & selectedSegments, Config & config){
  Seg *seg;                                                             // current selectd segment
  selectedSegments.rewind();                                            // reset the reader at the begin of the input stream
  while((seg=selectedSegments.getSeg())!=NULL)                          // For each of the selected segments
    computeWarp(histoT,destH,fs,seg,config);                            // Normalize the features
}

// ----------------------------------------------------------------------------------------------------------
// Feature Mean subtraction and Cov reduction for a segment and cluster (segment is considerred to be the minimum time unit to perform this).
void computeZeroOne(const DoubleVector &featureMean,const DoubleVector &featureStd,FeatureServer & fs,unsigned long begin, unsigned long length,Config &config) {
  unsigned long vectsize=fs.getVectSize();                                       // Get the vect size (number of coeff)
  Feature f;
  fs.seekFeature(begin);
  for (unsigned long idxFrame=0;idxFrame<length;idxFrame++){                     // for all the features of the segment
    fs.readFeature(f,0);  
    // Get the feature;
    for (unsigned int i = 0; i < vectsize; i++)       {                  // For each coeff
      f[i]=(f[i]-featureMean[i])/featureStd[i];    // Apply the 0 mean 1 cov normalisation
    }
    fs.writeFeature(f);
  }
}
void computeZeroOne(FrameAccGD &frameAccu,FeatureServer & fs,unsigned long begin, unsigned long length,Config &config) {  
  const DoubleVector & featureMean = frameAccu.getMeanVect();                     // Get the mean vector
  const DoubleVector & featureStd = frameAccu.getStdVect();                       // Get the std vector (sqrt(cov))
  computeZeroOne(featureMean,featureStd,fs,begin,length,config);
}
// on a segment
void computeZeroOne(FrameAccGD & frameAccu, FeatureServer & fs, Seg* seg, Config & config){
  const DoubleVector & featureMean = frameAccu.getMeanVect();                     // Get the mean vector
  const DoubleVector & featureStd = frameAccu.getStdVect();                       // Get the std vector (sqrt(cov))
  unsigned long begin=seg->begin()+fs.getFirstFeatureIndexOfASource(seg->sourceName());          // Idx of the first frame of the current file in the feature server
  computeZeroOne(featureMean,featureStd,fs,begin,seg->length(),config);                      // Normalize the feature to fit 0 mean, 1 cov
}
void computeZeroOne(const DoubleVector &featureMean,const DoubleVector &featureStd, FeatureServer & fs, Seg* seg, Config & config){
  unsigned long begin=seg->begin()+fs.getFirstFeatureIndexOfASource(seg->sourceName());          // Idx of the first frame of the current file in the feature server
  computeZeroOne(featureMean,featureStd,fs,begin,seg->length(),config);                      // Normalize the feature to fit 0 mean, 1 cov
}
// on a cluster
void computeZeroOne(FrameAccGD & frameAccu, FeatureServer & fs, SegCluster & selectedSegments, Config & config){
  const DoubleVector & featureMean = frameAccu.getMeanVect();                     // Get the mean vector
  const DoubleVector & featureStd = frameAccu.getStdVect();                       // Get the std vector (sqrt(cov))
  computeZeroOne(featureMean,featureStd,fs,selectedSegments,config); 
}
void computeZeroOne(const DoubleVector &featureMean,const DoubleVector &featureStd, FeatureServer & fs, SegCluster & selectedSegments, Config & config){
  if (verbose) cout << "(General Tools) Compute CMS on Feature Server" << endl;
  Seg *seg;                                                                         // current selectd segment
  selectedSegments.rewind();                                                      // reset the reader at the begin of the input stream
  while((seg=selectedSegments.getSeg())!=NULL) {                 // For each of the selected segments
    computeZeroOne(featureMean,featureStd,fs,seg,config);}                      // Normalize the feature to fit 0 mean, 1 cov
}

void cms(String & featureFileName,FeatureServer &fs,Config &config) {
        unsigned long begin=fs.getFirstFeatureIndexOfASource(featureFileName);
	fs.seekFeature(begin);
	SegServer segmentsServer;
	LabelServer labelServer;
	initializeClusters(featureFileName,segmentsServer,labelServer,config);
	verifyClusterFile(segmentsServer,fs,config);
	unsigned long codeSelectedFrame=labelServer.getLabelIndexByString(config.getParam("labelSelectedFrames"));
	SegCluster& selectedSegments=segmentsServer.getCluster(codeSelectedFrame);  
	selectedSegments.rewind();  
	RealVector <double> mean,cov;
	FrameAccGD frameAccu;
	frameAccu.reset();   
	accumulateStatFrame(frameAccu,fs, selectedSegments, config);
	mean = frameAccu.getMeanVect();     // Get the mean vector
	cov  = frameAccu.getStdVect();      // Get the std vector
	computeZeroOne(mean,cov,fs, selectedSegments, config);	
}

// Feature writing in an output stream w - could be used for multiple segmen,ts from multiple files to one file
void outputFeatureFile(Config &config, FeatureServer &fs, Feature & f, FeatureFileWriter &w) { 
  fs.readFeature(f);
  w.writeFeature(f); 
}
// on a part of a file
void outputFeatureFile(Config &config, FeatureServer &fs, unsigned long begin,unsigned long length, FeatureFileWriter &w) {
  Feature f;
  fs.seekFeature(begin);
  for (unsigned long idxFrame=0;idxFrame<length;idxFrame++){
    outputFeatureFile(config, fs, f,w);
  }  
}
// on a segment
void outputFeatureFile(Config &config, FeatureServer &fs, Seg * seg, FeatureFileWriter &w) {
  unsigned long begin=seg->begin()+fs.getFirstFeatureIndexOfASource(seg->sourceName());          // Idx of the first frame of the current file in the feature server
  if (verbose) cout <<"(GeneralTools) Writing ["<<seg->sourceName()<<"]"<<" begin:"<<begin<<" length:"<<seg->length()<<endl;
  outputFeatureFile(config,fs,begin,seg->length(),w);
}

// on a cluster
void outputFeatureFile(Config &config, FeatureServer &fs, SegCluster & selectedSegments, FeatureFileWriter &w) {
  Seg *seg;                                                                         // current selectd segment
  selectedSegments.rewind();                                                        // reset the reader at the begin of the input stream
  while((seg=selectedSegments.getSeg())!=NULL){
    outputFeatureFile(config, fs, seg, w); 
  }
}

// Feature Mapping related functions
unsigned long getBestGaussian(Mixture & M, Feature & f) {
  double vrais=0.0;
  double _v=0.0;
  unsigned long idx=0;
  unsigned long model_size=M.getDistribCount();
  for (unsigned long g=0;g<model_size;g++){
    _v=M.getDistrib(g).computeLK(f)*M.weight(g);
    if(_v>vrais) {
      idx=g;
      vrais=_v;
    }
  }
  return idx;
}

void mapDataToDistrib(double & data, const double meanData, const double covData, const double meanMap, const double covMap) {
  //cout<< data << " mu,std: "<<meanData<<","<<covData<<" "<<meanMap<<","<<covMap << endl;
  data=sqrt(covMap/covData)*(data-meanData)+meanMap;
}

void featureMapping(MixtureServer & ms, Feature & f,Config &config) {
  unsigned long vectsize=f.getVectSize();                                       // Get the vect size (number of coeff)
  unsigned long idCD=getBestGaussian(ms.getMixture(1),f);						// Finding the best gaussian in the sub-model	
  DistribGD & DCD=ms.getMixtureGD(1).getDistrib(idCD); // get the gaussian in both models
  DistribGD & DCI=ms.getMixtureGD(0).getDistrib(idCD);
  for (unsigned int i = 0; i < vectsize; i++) {					// map
    if (debug && (i==0 || i==1)) cout<<"C:"<<idCD<<"["<<i<<"] CDm,v:"<<DCD.getMean(i)<<","<<DCD.getCov(i)<<" CIm,v: "<<DCI.getMean(i)<<","<<DCI.getCov(i)<< endl;	
    mapDataToDistrib(f[i],DCD.getMean(i),DCD.getCov(i),DCI.getMean(i),DCI.getCov(i));}
}

// on a segment
void featureMapping(MixtureServer & ms, FeatureServer & fs,Seg * seg,Config &config) {
  unsigned long begin=seg->begin()+fs.getFirstFeatureIndexOfASource(seg->sourceName()); // Idx of the first frame of the current file in the feature server
  fs.seekFeature(begin);
  Feature f;
  for (unsigned long idxFrame=0;idxFrame<seg->length();idxFrame++){                          // for all the features of the segment
    fs.readFeature(f,0); 		
    featureMapping(ms,f,config);
    fs.writeFeature(f);
  }
}

// on a cluster
void featureMapping(MixtureServer & ms, FeatureServer & fs,SegCluster & selectedSegments,Config &config) {
  Seg *seg;                                                                         // current selectd segment
  selectedSegments.rewind();                                                      // reset the reader at the begin of the input stream
  while((seg=selectedSegments.getSeg())!=NULL){                                   // For each of the selected segments
    featureMapping(ms,fs,seg,config);
  }	
}

// Model based distance
// Authors: JF Bonastre  - Driss Matrouf
//
double likelihoodGD(const DistribGD& d,const DistribGD &m){
  double partial=0;
  double *dMean=d.getMeanVect().getArray();
  double *mMean=m.getMeanVect().getArray();
  double *dCov=d.getCovVect().getArray();
  double *mCov=m.getCovVect().getArray();
  for (unsigned long idxC=0;idxC<d.getVectSize();idxC++){
    double meanDiff=dMean[idxC]-mMean[idxC];
    partial+=(dCov[idxC]+(meanDiff*meanDiff))/mCov[idxC];
  }
  return m.getCst()*exp(-0.5*partial);
}
double likelihoodGD(const MixtureGD& data,const MixtureGD &model,TabWeight &tabWD,TabWeight &tabWM){
  // could work with a subset of components in the data model
  double result=0;
  TabWeightElem *dataArray=tabWD.getArray();
  TabWeightElem *modelArray=tabWM.getArray();

  for (unsigned long idxDataG=0;idxDataG<tabWD.getNbTop();idxDataG++){
    if (debug) cout <<"Distrib Data["<<tabWD.getDistrib(idxDataG)<<"] weight["<<tabWD.getWeight(idxDataG)<<" "<<dataArray[idxDataG].weight<<"]"<<endl;
    DistribGD &d=*(dataArray[idxDataG].distribP);
    double lkCompData=0.0;
    for  (unsigned long idxModelG=0;idxModelG<tabWM.getNbTop();idxModelG++){
      if (debug) cout <<"Distrib Model["<<tabWM.getDistrib(idxModelG)<<"] weight["<<tabWM.getWeight(idxModelG)<<" "<<modelArray[idxModelG].weight<<"]"<<endl;
      DistribGD &m=*(modelArray[idxModelG].distribP);
      double tmp=modelArray[idxModelG].weight*likelihoodGD(d,m);
      lkCompData+=tmp; 
    }
    result+=dataArray[idxDataG].weight*log(lkCompData);
  }
  return result;
}
double likelihoodGD(const MixtureGD& data,const MixtureGD &model,TabWeight &tabWD){
  TabWeight tabWM(model);
  return likelihoodGD(data,model,tabWD,tabWM);
}
double likelihoodGD(const MixtureGD& data,const MixtureGD &model){
  TabWeight tabWM(model);
  TabWeight tabWD(data);
  return likelihoodGD(data,model,tabWD,tabWM);
}


// Information on the quantity of data available for each client
// Outputs a list with the selected files for a defined quantity of data
int ExtractTargetDataInfo(Config& config)
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

#endif //!defined(ALIZE_GeneralTools_cpp)

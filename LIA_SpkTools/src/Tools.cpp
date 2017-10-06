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

#if !defined(ALIZE_Tools_cpp)
#define ALIZE_Tools_cpp

#include "Tools.h"


/*************************************************************
* learnModelFromSplit: building of a model by EM and split initialisation
**************************************************************/
MixtureGD& learnModelFromSplit(Config& config,SegCluster& clusterSeg,
		   StatServer& ss,FeatureServer &fs,MixtureServer& ms, unsigned long maxDistribCount, unsigned long deleteMixture, unsigned long nbIt){
  

TrainCfg trainCfg(config);
trainCfg.setNbTrainIt(3);

FrameAccGD globalFrameAcc;                                                    	//compute global Mean and cov
globalMeanCov (fs,clusterSeg,globalFrameAcc,config);                             // Compute the global mean and covariance
DoubleVector globalMean=globalFrameAcc.getMeanVect(); 
DoubleVector globalCov=globalFrameAcc.getCovVect();

MixtureGD& world=mixtureInitBySplit(config,ms,ss,fs,clusterSeg,globalMean,globalCov,maxDistribCount,trainCfg, deleteMixture);                   // Init the model
if(verbose) cout << "Init Model" << endl;

trainCfg.setBaggedFrameProbability(1.0);
trainCfg.setNbTrainIt(nbIt);
trainModel(config,ss,fs,clusterSeg,globalMean,globalCov,world, trainCfg);            // Train it by EM/ML


return world;

}//learnModelFromSplit


/*************************************************************
* learnModel: building of a model by EM and model already initialised
**************************************************************/
void learnModel(Config& config,SegCluster& clusterSeg,
		   StatServer& ss,FeatureServer &fs,MixtureServer& ms, MixtureGD& mInit){
  
TrainCfg trainCfg(config);

FrameAccGD globalFrameAcc;                                                    	//compute global Mean and cov
globalMeanCov (fs,clusterSeg,globalFrameAcc,config);                             // Compute the global mean and covariance
DoubleVector globalMean=globalFrameAcc.getMeanVect(); 
DoubleVector globalCov=globalFrameAcc.getCovVect();

trainCfg.setBaggedFrameProbability(1.0);
trainCfg.setNbTrainIt(5);
trainModel(config,ss,fs,clusterSeg,globalMean,globalCov,mInit, trainCfg);            // Train it by EM/ML


}//learnModel



/*************************************************************
* learnModelFromBaggedFrame: building of a model by EM and bagged frame initialisation
**************************************************************/
void learnModelFromBaggedFrame(Config& config,SegCluster& clusterSeg,MixtureGD& world,
		   StatServer& ss,FeatureServer &fs,MixtureServer& ms){
  
FrameAccGD globalFrameAcc;                                                    //compute global Mean and cov
globalMeanCov (fs,clusterSeg,globalFrameAcc,config);                             // Compute the global mean and covariance
DoubleVector globalMean=globalFrameAcc.getMeanVect(); 
DoubleVector globalCov=globalFrameAcc.getCovVect();
double factor=300.0;

TrainCfg trainCfg(config);

// mixture init
if(totalFrame(clusterSeg) < 6000){
	factor=30.0;	
}

trainCfg.setInitVarFlooring(0.5);
trainCfg.setFinalVarFlooring(0.1);

double baggedFrameProbabilityInit = (double)(world.getDistribCount()*factor/(double)totalFrame(clusterSeg));
cout << "Training Model: " << world.getDistribCount() << " " << totalFrame(clusterSeg) << " " << baggedFrameProbabilityInit << endl;
trainCfg.setBaggedFrameProbabilityInit(baggedFrameProbabilityInit);
trainCfg.setBaggedFrameProbability(0.1);
trainCfg.setNbTrainIt(2);
mixtureInit(ms,fs,world,clusterSeg,globalCov,config,trainCfg);                   // Init the model


if(verbose) cout<<" Training EM for world"<<endl;	
trainCfg.setBaggedFrameProbability(0.1);
trainCfg.setNbTrainIt(10);
trainModel(config,ss,fs,clusterSeg,globalMean,globalCov,world, trainCfg);            // Train it by EM/ML
trainCfg.setBaggedFrameProbability(1.0);
trainCfg.setNbTrainIt(2);
trainModel(config,ss,fs,clusterSeg,globalMean,globalCov,world, trainCfg);            // Train it by EM/ML

}//learnModelFromBaggedFrame


/************************************************************
* modifyMAPConfig: Modify MAP adaptation parameters according to the context
*************************************************************/
void modifyMAPConfig(Config& config,MAPCfg& cfg, String context){

if(context=="selection"){
	cfg.setMethod(config.getParam("MAPAlgoSelection"));	
	cfg.setMeanAdapt(config.getParam("meanAdaptSelection").toBool());
	if(cfg.getMeanAdapt()){
		cfg.setMeanReg(config.getParam("MAPRegFactorMeanSelection").toDouble());
		cfg.setMeanAlpha(config.getParam("MAPAlphaMeanSelection").toDouble());
	}
	cfg.setWeightAdapt(config.getParam("weightAdaptSelection").toBool());
	if(cfg.getWeightAdapt()){
		cfg.setWeightReg(config.getParam("MAPRegFactorWeightSelection").toDouble());
		cfg.setWeightAlpha(config.getParam("MAPAlphaWeightSelection").toDouble());
	}
	cfg.setVarAdapt(config.getParam("varAdaptSelection").toBool());
	if(cfg.getVarAdapt()){
		cfg.setVarReg(config.getParam("MAPRegFactorVarSelection").toDouble());
		cfg.setVarAlpha(config.getParam("MAPAlphaVarSelection").toDouble());
	}
}
else if(context=="normal"){
	cfg.setMethod(config.getParam("MAPAlgo"));	
	cfg.setMeanAdapt(config.getParam("meanAdapt").toBool());
	if(cfg.getMeanAdapt()){
		cfg.setMeanReg(config.getParam("MAPRegFactorMean").toDouble());
		cfg.setMeanAlpha(config.getParam("MAPAlphaMean").toDouble());
	}		
	cfg.setWeightAdapt(config.getParam("weightAdapt").toBool());
	if(cfg.getWeightAdapt()){
		cfg.setWeightReg(config.getParam("MAPRegFactorWeight").toDouble());
		cfg.setWeightAlpha(config.getParam("MAPAlphaWeight").toDouble());
	}
	cfg.setVarAdapt(config.getParam("varAdapt").toBool());
	if(cfg.getVarAdapt()){
		cfg.setVarReg(config.getParam("MAPRegFactorVar").toDouble());
		cfg.setVarAlpha(config.getParam("MAPAlphaVar").toDouble());
	}
}
}//modifyMAPConfig

/************************************************************
* findAlreadyUsed: verify if the selected segment is suitable or not
*************************************************************/
bool findAlreadyUsed(SegCluster& clusterUsedSeg,unsigned long start,unsigned long length,String fileLabel){

Seg* segment;
unsigned long startc,lengthc;
String fileLabelc;

clusterUsedSeg.rewind();
while((segment=clusterUsedSeg.getSeg())!=NULL){
	startc=segment->begin();	
	lengthc=segment->length();
	fileLabelc=segment->sourceName();
	
	if(fileLabelc==fileLabel){
		if((start>=startc)&&(start<lengthc+startc-1))return true;	
	}
		
}//while segments
return false;
}//findAlreadyUsed

/*************************************************************
* CreateInitHMM: building of an HMM with 1 state, GMM initialized by model adaptation
**************************************************************/
void CreateInitHMM(Config& config,MAPCfg &cfg,SegCluster& clusterSeg,MixtureGD& world,
		   hmm& hmmv,StatServer& ss,FeatureServer &fs,MixtureServer& ms){
  
hmmv.reset();
hmmv.LoadState(world,clusterSeg.string());
MixtureGD& m0=(MixtureGD&)hmmv.getDensity(0);
if(verbose) cout<<"Adaptation (MAP) du modele L0"<<endl;

adaptModel(config,ss,ms,fs,clusterSeg,world,m0,cfg);        
}//CreateInitHMM


/*************************************************************
* CreateInitHMM: building of an HMM with 1 state, GMM initialized by EM
**************************************************************/
void CreateInitHMM(Config& config,SegCluster& clusterSeg, hmm& hmmv, StatServer& ss,FeatureServer &fs,MixtureServer& ms){
  
hmmv.reset();
unsigned long maxDistribCount=config.getParam("mixtureDistribCount").toLong();
MixtureGD& world=learnModelFromSplit(config,clusterSeg,ss,fs,ms,maxDistribCount,0,20);
hmmv.setDensity(world, 0);
}

/*******************************************************
* findNameClass: Extract the acoustic class from the speaker label
********************************************************/
void findNameClass(String& name,String& nameClass){
char temp2[300];
int i;

int ind=name.find("_");
const char *temp1=name.c_str();

for(i=0;i<ind;i++)
	temp2[i]=temp1[i];
temp2[i]=0;
	
nameClass=temp2;	
}//findNameClass


/************************************************************
* isNewSpeaker: verify wether the speaker label looks like L0x (new speaker)
*************************************************************/
bool isNewSpeaker(String nameCluster)
{
	const char *temp=nameCluster.c_str();
	long ind=nameCluster.find("_L");
	
	return (temp[ind+2]=='a');
}

// A function for initializing an hmm with several states
void InitHMM(Config& config,SegServer& segmentsServer,MixtureGD& world,hmm& actualHMM)
{
  if(verbose) cout<<">> HMM Initialisation <<"<<endl;
  actualHMM.reset();
  for(unsigned long i=0;i<segmentsServer.getClusterCount();i++){
    SegCluster& cluster=segmentsServer.getCluster(i);
    actualHMM.LoadState(world,cluster.string());
  }
  if(verbose) cout<<" HMM has "<<actualHMM.getNbState()<<" states"<<endl;
  String transitionMethod=config.getParam("transitionMethod"); // Compute the transitions
  double gamma=config.getParam("gammaTransition").toDouble();
  computeTransitions(config,actualHMM,transitionMethod,gamma,segmentsServer);
}//InitHMM


// A function for initializing an hmm with several states
void InitHMM(Config& config,SegServer& segmentsServer,MixtureGD& speech,MixtureGD& nonspeech, hmm& actualHMM)
{
  if(verbose) cout<<">> HMM Initialisation <<"<<endl;
  actualHMM.reset();
  for(unsigned long i=0;i<segmentsServer.getClusterCount();i++){
    SegCluster& cluster=segmentsServer.getCluster(i);
    if(cluster.string() == "speech") actualHMM.LoadState(speech,cluster.string());
    else if(cluster.string() == "nonspeech") actualHMM.LoadState(nonspeech,cluster.string());
  }
  if(verbose) cout<<" HMM has "<<actualHMM.getNbState()<<" states"<<endl;
  String transitionMethod=config.getParam("transitionMethod"); // Compute the transitions
  double gamma=config.getParam("gammaTransition").toDouble();
  computeTransitions(config,actualHMM,transitionMethod,gamma,segmentsServer);
}//InitHMM


// A function for initializing an hmm with several states
void InitHMMOverlap(Config& config,SegServer& segmentsServer,MixtureGD& world,hmm& actualHMM)
{
  if(verbose) cout<<">> HMM Initialisation <<"<<endl;
  actualHMM.reset();
  for(unsigned long i=0;i<segmentsServer.getClusterCount();i++){
    SegCluster& cluster=segmentsServer.getCluster(i);
    actualHMM.LoadState(world,cluster.string());
  }
  actualHMM.LoadState(world,"speech_L99");
  if(verbose) cout<<" HMM has "<<actualHMM.getNbState()<<" states"<<endl;
  String transitionMethod=config.getParam("transitionMethod"); // Compute the transitions
  double gamma=config.getParam("gammaTransition").toDouble();
  computeTransitions(config,actualHMM,transitionMethod,gamma,segmentsServer);
}//InitHMM

// 

/****************************************************************
* computeTransitions: HMM transition computation: 
* either equiprobable or with a constant pii and pij ratio 
*****************************************************************/
void computeTransitions(Config& config,hmm& actualHMM,String transitionMethod,double gamma,SegServer& Server){

if(transitionMethod=="Equiprob"){
	if((gamma<0)||(gamma>1)){
		cout<<"Gamma is not in interval 0 1! "<<gamma<<endl;
		exit(-1);
	}
	for(unsigned long i=0;i<actualHMM.getNbState();i++){
		for(unsigned long j=0;j<actualHMM.getNbState();j++){
			if(i==j){
				actualHMM.setTransition(gamma,i,j);
				if ((verbose) && (verboseLevel == 2)) cout<<gamma<<" ";
			}
			else{
				actualHMM.setTransition((1.0-gamma)/(double)(actualHMM.getNbState()-1.0),i,j);
				if ((verbose) && (verboseLevel == 2)) cout<<(1.0-gamma)/(double)(actualHMM.getNbState()-1.0)<<" ";
			}
		}//for j
		if ((verbose) && (verboseLevel == 2)) cout<<endl;
	}//for i
}//if equiprob
if(transitionMethod=="RapConst"){
	double transij = ((double)1.0 / (gamma + (double)(actualHMM.getNbState()-1.0)));
	double transii = (gamma * transij);
	for(unsigned long i=0;i<actualHMM.getNbState();i++){
		for(unsigned long j=0;j<actualHMM.getNbState();j++){
			if(i==j){
				actualHMM.setTransition(transii,i,j);
				if ((verbose) && (verboseLevel == 2)) cout<<transii<<" ";
			}
			else{
				actualHMM.setTransition(transij,i,j);
				if (verbose) cout<<transij<<" ";
			}
		}//for j
	if ((verbose) && (verboseLevel == 2)) cout<<endl;
	}//for i
}//if RapConst


if(transitionMethod=="Unity"){
	for(unsigned long i=0;i<actualHMM.getNbState();i++){
		for(unsigned long j=0;j<actualHMM.getNbState();j++){
			actualHMM.setTransition(1.0,i,j);
			if ((verbose) && (verboseLevel == 2)) cout<<1.0<<" ";
		}//for j
	if ((verbose) && (verboseLevel == 2)) cout<<endl;
	}//for i
}//if RapConst

// (used for the re-segmentation!) the idea here is that the transitions are dependent on the speaker time duration computed during the segmentation phase 
if(transitionMethod=="SpeakerTime"){
	unsigned long *speakerTime,timeTotal=0;
	double transij,*speakerPercent; 
	speakerTime=new unsigned long[Server.getClusterCount()];
	speakerPercent=new double[Server.getClusterCount()];
	for(unsigned long icluster=0;icluster<Server.getClusterCount();icluster++){
		//cluster extraction 
		SegCluster& cluster=Server.getCluster(icluster);
		speakerTime[icluster]=totalFrame(cluster);
		timeTotal+=speakerTime[icluster];
	}//for icluster
	for(unsigned long icluster=0;icluster<Server.getClusterCount();icluster++){
		speakerPercent[icluster]=(double)speakerTime[icluster]/(double)timeTotal;
	}//for icluster
	for(unsigned long i=0;i<actualHMM.getNbState();i++){
		for(unsigned long j=0;j<actualHMM.getNbState();j++){
			if(i==j){
				actualHMM.setTransition(gamma,i,j);
				if (verbose) cout<<gamma<<" ";
			}
			else{
				transij=(1-gamma)*(speakerPercent[j]+(speakerPercent[i]/Server.getClusterCount()));
				actualHMM.setTransition(transij,i,j);
				if (verbose) cout<<transij<<" ";
			}
		}//for j
		if (verbose) cout<<endl;
	}//for i
}//SpeakerTime
}//computeTransitions

void computeTransitions(Config& config,hmm& actualHMM,SegServer& Server){ 

String transitionMethod=config.getParam("transitionMethod");
double gamma=config.getParam("gammaTransition").toDouble();
computeTransitions(config,actualHMM,transitionMethod,gamma,Server);
}


// comparison between two segmentations according to different criteria
bool isComparableAndVerifSpeaker(Config& config,SegServer& segTemp,SegServer& actualSeg,DoubleVector& actualViterbiProb,DoubleVector& previousViterbiProb){

/*if(isValidLengthOfLastSpeaker(config,actualSeg)){
	return true; 
}*/

return (isComparable(config,segTemp,actualSeg,actualViterbiProb,previousViterbiProb));
}//isComparableAndVerif


// comparison between two segmentations according to different criteria
bool isComparable(Config& config,SegServer& segTemp,SegServer& actualSeg,DoubleVector& actualViterbiProb,DoubleVector& previousViterbiProb){

double epsilon=config.getParam("decodingEpsilon").toDouble();

if(verbose) cout<<">> Comparison between 2 segmentations <<"<<endl;

if(segTemp.getClusterCount()!=actualSeg.getClusterCount()){
	cout<<"ERROR : Problem with number of clusters! ("<<segTemp.getClusterCount()<<" != "<<actualSeg.getClusterCount()<<endl;
	exit(-1);
}

for(unsigned long icluster=0;icluster<segTemp.getClusterCount();icluster++){
	if(isDifferentSegmentation(config, segTemp.getCluster(icluster), actualSeg.getCluster(icluster))){
		if(!isLessEpsilon(config,actualViterbiProb,previousViterbiProb,epsilon))
						return false;
	}
}
if(debug) cout << "isComparable => true !!!"<< endl;
return true;
}//isComparable

// comparison between two segmentations according to the boundaries
bool isDifferentSegmentation(Config& config, SegCluster& clusterA, SegCluster& clusterB){
       
Seg *segment1,*segment2;
 	
if(debug) cout << "Same number of segments ?" << clusterA.getCount() << " vs " << clusterB.getCount() << endl;
if(clusterA.getCount()!=clusterB.getCount()){ 
	if(debug) cout << "Different Segmentations = True" << endl;
 	return true;
}

clusterA.rewind();clusterB.rewind();
while(((segment1=clusterA.getSeg())!=NULL)){
	segment2=clusterB.getSeg();
	if(debug){
		cout << "Segment start :" << segment1->begin() << " vs " << segment2->begin() << endl;
		cout << "Segment length:" << segment1->length() << " vs " << segment2->length() << endl;
	}
	if((segment1->begin()!=segment2->begin())||(segment1->length()!=segment2->length())){	
		if(debug) cout << "Different segmentation = True" << endl;
		return true;
	}
}
if(debug) cout << "Identical segmentation => false" << endl;

return false;
} //isDifferentSegmentation

// comparison between two viterbi paths according to likelihoods
bool isLessEpsilon(Config& config,DoubleVector& actualViterbiProb,DoubleVector& previousViterbiProb,double& epsilon){
double diff;

for(unsigned long iseg=0;iseg<actualViterbiProb.size();iseg++){
	diff=fabs(actualViterbiProb[iseg]-previousViterbiProb[iseg]);
	if(debug) cout<<"Diff :"<<diff<<" ? Epsilon : "<<epsilon<<endl;
	if(diff>epsilon){
		for(unsigned long jseg=0;jseg<actualViterbiProb.size();jseg++)
			previousViterbiProb[jseg]=actualViterbiProb[jseg];
		return false;
	}
}
return true;

}//isLessEpsilon

// Verify the validity of the last speaker according to the corresponding duration
bool isValidLengthOfLastSpeaker(Config& config,SegServer& actualSeg){
	unsigned long speakerMinTime=config.getParam("speakerMinTime").toLong();
	unsigned long lastSpeakerTime=totalFrame(actualSeg.getCluster(actualSeg.getClusterCount()-1));
	
	return (lastSpeakerTime <= speakerMinTime);
}


// Compute the likelihood betwwen a segmentation and an HMM
double ComputePath(Config& config,hmm& actualHMM,DoubleVector transitionsa,SegServer& actualSeg,XList& listseg,StatServer& ss){
String file;
int *validity,dr,imem=0,etatPrec=-1;
Seg *segmin;
unsigned long startMin=0;
FeatureServer *fs;
double llr=0;
XLine *linep;

if((verbose) && (verboseLevel == 2)){
	cout<<"Compute LLR of segmentation"<<endl;
 	cout<<"Segmentation to be compared"<<endl;
	displayAllSegments(config, actualSeg);
}       

listseg.rewind();
while ((linep=listseg.getLine()) != NULL){
	file = linep->getElement(0,false);
	fs=new FeatureServer(config, file);
	
	validity=new int[actualSeg.getSegCount()];
	unsigned long i;
	for(i=0;i<actualSeg.getSegCount();i++) validity[i]=0;
	do{
		dr=0;
		for(i=0;i<actualSeg.getSegCount();i++){
			Seg& segment=actualSeg.getSeg(i);
			if((segment.sourceName()==file)&&(validity[i]!=1)){
				if(dr==0){
					startMin=segment.begin();	
					segmin=&segment;
					imem=i;dr=1;
				}else
					if(segment.begin()<startMin){
						startMin=segment.begin();	
						segmin=&segment;	
						imem=i;
					}
			}else validity[i]=1;
		}
	
		validity[imem]=1;
		i=0;		
		while((i < actualHMM.getNbState()) && (actualHMM.getStateName(i)!=segmin->string()))i++;
		if(i == actualHMM.getNbState()){
			cout << i << ": problem with segment label: " << segmin->string() << " " << segmin->begin()<< " " <<segmin->begin()+segmin->length() << " Does not correspond to any hmmState name" << endl;
			for(unsigned long j=0; j<actualHMM.getNbState(); j++){
		 		cout << actualHMM.getStateName(j) << endl;
			 }		
		}
	       else{
			MixtureGD& m=(MixtureGD&)actualHMM.getDensity(i);
			if(etatPrec==-1) etatPrec=i;
			Feature f;
			for(unsigned long j=segmin->begin();j<=(segmin->begin()+segmin->length()-1);j++){
				fs->seekFeature(j);
				fs->readFeature(f);
				llr+=log(transitionsa[etatPrec+i*actualHMM.getNbState()])+ss.computeLLK(m,f);
				etatPrec=i;	
			}
	      	} // fin else
		dr=0;
		for(unsigned long j=0;j<actualSeg.getSegCount();j++)
			if(validity[j]==0) dr=1;
	}while(dr==1);	
	
	delete []validity;delete fs;
}//while linep

return (llr);

}//ComputePath

// Compute the likelihood betwwen a segmentation and an HMM
double ComputePath(Config& config,hmm& actualHMM,DoubleVector transitionsa,SegServer& actualSeg,FeatureServer &fs,StatServer& ss, MixtureGD& world){

int *validity,dr,imem=0,etatPrec=-1;
Seg *segmin;
unsigned long startMin=0;
double llr=0;

if((verbose) && (verboseLevel == 2)){
	cout<<"Compute LLR of segmentation"<<endl;
 	cout<<"Segmentation to be compared"<<endl;
	displayAllSegments(config, actualSeg);
}       

validity=new int[actualSeg.getSegCount()];
unsigned long i;
for(i=0;i<actualSeg.getSegCount();i++) validity[i]=0;
do{
	dr=0;
	for(i=0;i<actualSeg.getSegCount();i++){
		Seg& segment=actualSeg.getSeg(i);
		if(validity[i]!=1){
			if(dr==0){
				startMin=segment.begin();	
				segmin=&segment;
				imem=i;dr=1;
			}else
				if(segment.begin()<startMin){
					startMin=segment.begin();	
					segmin=&segment;	
					imem=i;
				}
				}else validity[i]=1;
	}
	validity[imem]=1;
	i=0;		
	while((i < actualHMM.getNbState()) && (actualHMM.getStateName(i)!=segmin->string()))i++;
	if(i == actualHMM.getNbState()){
		cout << i << ": problem with segment label: " << segmin->string() << " " << segmin->begin()<< " " <<segmin->begin()+segmin->length() << " Does not correspond to any hmmState name" << endl;
		for(unsigned long j=0; j<actualHMM.getNbState(); j++){
	 		cout << actualHMM.getStateName(j) << endl;
		 }		
	}
	else{
		MixtureGD& m=(MixtureGD&)actualHMM.getDensity(i);
		if(etatPrec==-1) etatPrec=i;
		Feature f;
		for(unsigned long j=segmin->begin();j<=(segmin->begin()+segmin->length()-1);j++){
			fs.seekFeature(j);
			fs.readFeature(f);
			llr+=log(transitionsa[etatPrec+i*actualHMM.getNbState()])+(ss.computeLLK(m,f)-ss.computeLLK(world,f));
			etatPrec=i;	
		}
	} // fin else
	dr=0;
	for(unsigned long j=0;j<actualSeg.getSegCount();j++)
		if(validity[j]==0) dr=1;
}while(dr==1);	
	
delete []validity;

return (llr);

}//ComputePath


// Compute the likelihood betwwen a state of an HMM and associated set of segments
double computeLLKMeanOnOneCluster(Config& config,hmm& actualHMM,SegServer& actualSeg,unsigned long icluster,StatServer& ss,FeatureServer &fs){

SegCluster& clusterT=actualSeg.getCluster(icluster);
MixtureGD& m=(MixtureGD&)actualHMM.getDensity(icluster);
return meanLikelihood(ss,fs,m,clusterT,config); 
}//ComputeLLROneCluster


void saveSegmentation(Config& config,SegServer& Resultat,FeatureServer &fs,String& outputFilesPath, int sorted){

double frameLength=config.getParam("frameLength").toDouble();
String outputFileName,label;
int *validity, dr;
unsigned long startMin,imem=0;
Seg *segmin;
FILE *f;
double start,stop;
String saveSegmentationExtension=config.getParam("saveSegmentationExtension");
  	
try{

if(sorted==0){
	for(unsigned long fi=0;fi<fs.getSourceCount();fi++){
		const String &file=fs.getNameOfASource(fi);
		outputFileName=outputFilesPath+file+saveSegmentationExtension;
		f=fopen(outputFileName.c_str(),"wt");
		if(verbose) cout<<"Save into the file "<<file<<saveSegmentationExtension<< " the following segmentation: " << endl;
		if(verbose) displayAllSegments(config, Resultat);
     
		validity=new int[Resultat.getSegCount()];
		unsigned long i;
		for(i=0;i<Resultat.getSegCount();i++) validity[i]=0;
		do{
			dr=0;
			for(i=0;i<Resultat.getSegCount();i++){
				Seg& segment=Resultat.getSeg(i);
				if((segment.sourceName()==file)&&(validity[i]!=1)){
					if(dr==0){
						startMin=segment.begin();	
						segmin=&segment;
						imem=i;dr=1;
					}else
					if(segment.begin()<=startMin){
						startMin=segment.begin();	
						segmin=&segment;	
						imem=i;
					}
				}else validity[i]=1;
			}
			
			validity[imem]=1;
			i=0;
			label=segmin->string();
			start=(double)segmin->begin();
			start*=frameLength;
			stop=(double)(segmin->begin()+segmin->length()-1);
			stop*=frameLength;
			fprintf(f,"%lf %lf %s \n",start,stop,label.c_str());
			printf("%lf %lf %s \n",start,stop,label.c_str());
			dr=0;
			for(unsigned long j=0;j<Resultat.getSegCount();j++)
				if(validity[j]==0) dr=1;
		}while(dr==1);	
		delete []validity;
		fclose(f);
	}//for
}else if(sorted==1){
	for(unsigned long fi=0;fi<fs.getSourceCount();fi++){
		const String &file=fs.getNameOfASource(fi);
		outputFileName=outputFilesPath+file+saveSegmentationExtension;
		f=fopen(outputFileName.c_str(),"wt");
		if(verbose) cout<<"Save into the file "<<file<<saveSegmentationExtension<< " the following segmentation: " << endl;
		if(verbose) displayAllSegments(config, Resultat);
     
		for(unsigned long i=0;i<Resultat.getSegCount();i++){
			Seg& segment=Resultat.getSeg(i);
			label=segment.string();
			start=(double)segment.begin();
			start*=frameLength;
			stop=(double)(segment.begin()+segment.length()-1);
			stop*=frameLength;
			fprintf(f,"%lf %lf %s \n",start,stop,label.c_str());
		}
		fclose(f);
	}//for
	
}	
}//try
catch (Exception& e){ 
	cout << e.toString().c_str() << endl;
}
}//saveSegmentation


void saveAcousticSegmentation(Config& config,SegServer& Resultat,FeatureServer &fs,String& outputFilesPath, int sorted){

double frameLength=config.getParam("frameLength").toDouble();
long limitClusterDuration=0;
long cumulClusterDuration=0;
long indexClusterDuration=0;

if(config.existsParam("limitClusterDuration")) limitClusterDuration=config.getParam("limitClusterDuration").toLong();

String outputFileName,label;
int *validity, dr;
unsigned long startMin,imem=0;
Seg *segmin;
FILE *f;
double start,stop;
String saveSegmentationExtension=config.getParam("saveSegmentationExtension");
  	
try{

if(sorted==0){
	for(unsigned long fi=0;fi<fs.getSourceCount();fi++){
		const String &file=fs.getNameOfASource(fi);
		outputFileName=outputFilesPath+file+saveSegmentationExtension;
		f=fopen(outputFileName.c_str(),"wt");
		if(verbose) cout<<"Save into the file "<<file<<saveSegmentationExtension<< " the following segmentation: " << endl;
		if(verbose) displayAllSegments(config, Resultat);
     
		validity=new int[Resultat.getSegCount()];
		unsigned long i;
		for(i=0;i<Resultat.getSegCount();i++) validity[i]=0;
		do{
			dr=0;
			for(i=0;i<Resultat.getSegCount();i++){
				Seg& segment=Resultat.getSeg(i);
				if((segment.sourceName()==file)&&(validity[i]!=1)){
					if(dr==0){
						startMin=segment.begin();	
						segmin=&segment;
						imem=i;dr=1;
					}else
					if(segment.begin()<=startMin){
						startMin=segment.begin();	
						segmin=&segment;	
						imem=i;
					}
				}else validity[i]=1;
			}
			
			validity[imem]=1;
			i=0;
			label=segmin->string();
			start=(double)segmin->begin();
			start*=frameLength;
			stop=(double)(segmin->begin()+segmin->length()-1);
			stop*=frameLength;
			String L=label;
			if(limitClusterDuration!=0){
				// Reset of counters
				if(cumulClusterDuration > limitClusterDuration){
					indexClusterDuration++;
					cumulClusterDuration=0;
				}
				// to reinit at each step
				L=label+String::valueOf(indexClusterDuration);
			}
			fprintf(f,"%lf %lf %s \n",start,stop,L.c_str());
			printf("%lf %lf %s \n",start,stop,L.c_str());
			dr=0;
			cumulClusterDuration+=segmin->length();
			for(unsigned long j=0;j<Resultat.getSegCount();j++)
				if(validity[j]==0) dr=1;
		}while(dr==1);	
		delete []validity;
		fclose(f);
	}//for
}else if(sorted==1){
	for(unsigned long fi=0;fi<fs.getSourceCount();fi++){
		const String &file=fs.getNameOfASource(fi);
		outputFileName=outputFilesPath+file+saveSegmentationExtension;
		f=fopen(outputFileName.c_str(),"wt");
		if(verbose) cout<<"Save into the file "<<file<<saveSegmentationExtension<< " the following segmentation: " << endl;
		if(verbose) displayAllSegments(config, Resultat);
     
		for(unsigned long i=0;i<Resultat.getSegCount();i++){
			Seg& segment=Resultat.getSeg(i);
			label=segment.string();
			start=(double)segment.begin();
			start*=frameLength;
			stop=(double)(segment.begin()+segment.length()-1);
			stop*=frameLength;
			fprintf(f,"%lf %lf %s \n",start,stop,label.c_str());
		}
		fclose(f);
	}//for
	
}	
}//try
catch (Exception& e){ 
	cout << e.toString().c_str() << endl;
}
}//saveAcousticSegmentation


// clean the last speaker cluster for a best segment selection
void CleanSpeaker(Config& config,SegCluster& clusterToClean,SegCluster& clusterL0, unsigned long startSeg,LabelServer& labelServer){

Seg *segment;

unsigned long start, length;
clusterToClean.rewind();
while((segment=clusterToClean.getSeg())!=NULL){
	start=segment->begin();	
	length=segment->length();
	if(!(((startSeg+1)>=start)&&(startSeg<length+start-1))){//not the segment to find
		moveSegmentFromOneClusterToAnother(labelServer, segment, clusterToClean, clusterL0);
		clusterToClean.rewind();
	}
}

}//CleanSpeaker

// look for empty speaker (without any segments and delete them
void NoDataSpeakerVerification(Config& config,hmm& actualHMM,SegServer& actualSeg){
int indCorr=0;
unsigned long speakerTime;
bool deletedSpeakerFlag=false;

for(unsigned long icluster=0;icluster<actualSeg.getClusterCount();icluster++){
	SegCluster& cluster=actualSeg.getCluster(icluster);
	speakerTime=totalFrame(cluster);
	
	if(speakerTime==0){//delete speaker
		if((verbose) && (verboseLevel==2)) cout<<" Because has 0 segments the speaker "<<cluster.string()<<" is deleted"<<endl;
		actualHMM.deleteState(icluster-indCorr);
		indCorr++;
		//transition computation
		String transitionMethod=config.getParam("transitionMethod");
		double gamma=config.getParam("gammaTransition").toDouble();
		computeTransitions(config,actualHMM,transitionMethod,gamma,actualSeg);
		
		
		//delete all corresponding cluster
		while(cluster.getCount()!=0){
			actualSeg.remove(cluster.get(0));		
		}
		cluster.setString("ToDelete");
		deletedSpeakerFlag=true;
	}

}//for icluster

if(deletedSpeakerFlag){
	//delete clusters
	bool flagDelete=false;
	do{
		flagDelete=false;
		unsigned long icluster=0;
		while((icluster<actualSeg.getClusterCount())&&(!flagDelete)){
			SegCluster& cluster=actualSeg.getCluster(icluster);
			if(cluster.string()=="ToDelete") {
				actualSeg.remove(cluster);
				flagDelete=true;
			}
			icluster++;
		}//while icluster
	}while(flagDelete);
}//if deletedSpeakerFlag

}//NoDataSpeakerVerification


//-- Accumulate the log likelihood for the selected frames and a given model
void accumulateStatLLKExt(StatServer &ss,FeatureServer &fs,MixtureStat &llkAcc, unsigned long idxBeginFrame,unsigned long nbFrames,
		       Config &config)  
{

	float * tabDouble;
	tabDouble = new float[nbFrames];
	double total=0;
	fs.seekFeature(idxBeginFrame);                                           // go to the frame in the buffer (and load it if needed)
	for (unsigned long n=0;n<nbFrames;n++){
		Feature f;
		if(fs.readFeature(f)==false)
			cout<<"No more features"<<endl;
		double val=llkAcc.computeLLK(f);
		tabDouble[n]=(float)val;
		total+=val;
		if (debug)
			cout << "likelihood Frame["<<idxBeginFrame+n<<"]="<<val<<endl;
	}
	qsort(tabDouble,nbFrames,sizeof(float),compare);
	int dix_pourcent = (int)(nbFrames*0.1);
	double t_min = tabDouble[nbFrames-dix_pourcent];
	double t_max = tabDouble[dix_pourcent];
	if(debug) {
		cout<<"T_min : "<<t_min<<endl;
		cout<<"T_max : "<<t_max<<endl;
	}
	fs.seekFeature(idxBeginFrame);                                           // go to the frame in the buffer (and load it if needed)
	for(unsigned long n=0;n<nbFrames;n++){
		Feature f;
		//fs.init(config);
		// cout << "n: " << idxBeginFrame+n ;
		if(fs.readFeature(f)==false)
			cout<<"No more features"<<endl;
		double to = llkAcc.computeLLK(f);
		if(to > t_min && to > t_max)
		{
			double toto=llkAcc.computeAndAccumulateLLK(f);
			if(debug) cout << "n: " << idxBeginFrame+n << " => likelihood Frame["<<"]="<<toto<<endl;
		}
	}

}

//-------------------------------------------------------------------------
//-- Compute the mean log likelihood for the Selected frames and a given model
// Required for selection on best frames !

int compare (const void  * a, const void * b){
  	if(*(float*)a - *(float*)b > 0)
		return -1;
	if(*(float*)a - *(float*)b == 0)
		return 0;
	return 1;
}

double meanLikelihoodExt(StatServer &ss,FeatureServer &fs,MixtureGD &model,unsigned long idxBeginFrame,unsigned long nbFrames,Config &config)  {
  double llk=-200.0;
  MixtureStat &llkAcc=ss.createAndStoreMixtureStat(model);
  llkAcc.resetLLK();
  accumulateStatLLKExt(ss,fs,llkAcc,idxBeginFrame,nbFrames,config); 
  if(llkAcc.getAccumulatedLLKFeatureCount() != 0){
  	llk=llkAcc.getMeanLLK();
  }
  ss.deleteMixtureStat(llkAcc);
  return llk;
}



void copySeg(SegCluster &cOld, SegCluster& cNew, SegServer &segTemp, LabelServer& labelServer,String l){
	Seg *segment;
	cOld.rewind();
	while((segment=cOld.getSeg()) != NULL){
		cNew.add(segTemp.createSeg(segment->begin(),segment->length(),0,l,segment->sourceName()));
	}
}


void createOneClusterPerSegment(SegServer& segmentServer, SegServer& segTemp){

for(unsigned long i=0;i<segmentServer.getSegCount();i++){
	Seg &segment = segmentServer.getSeg(i);
	SegCluster& cluster=segTemp.createCluster();
	cluster.add(segTemp.duplicateSeg(segment));
}
}

void createOneClusterPerFixedSegment(SegServer& segmentServer, SegServer& segTemp){

unsigned long fixed=200;
for(unsigned long i=0;i<segmentServer.getSegCount();i++){
	Seg &segment = segmentServer.getSeg(i);
	unsigned long start=segment.begin();
	while(start+fixed <= (segment.begin()+segment.length()-1)){
		SegCluster& cluster=segTemp.createCluster();
		cluster.add(segTemp.createSeg(start,fixed,0,segment.string(),segment.sourceName()));
		cout << start << " " << start+fixed-1 << endl;
		start+=fixed;
	}
	if((start+fixed > (segment.begin()+segment.length()-1)) && (segment.begin()+segment.length()-start>100)){
		SegCluster& cluster=segTemp.createCluster();
		cluster.add(segTemp.createSeg(start,segment.begin()+segment.length()-start,0,segment.string(),segment.sourceName()));
		cout << start << " " << start+(segment.begin()+segment.length()-start)-1 << endl;
	}
}	
}


//-------------------------------------------------------------------------
void viterbiDecoding(Config& config,hmm& actualHMM,SegCluster& cluster,
		 SegServer& segTemp,StatServer& ss,FeatureServer & fs,LabelServer& labelServer,DoubleVector& actualViterbiProb){

if (verbose) cout<<">> Decoding with HMM at  "<<actualHMM.getNbState()<<" states <<"<<endl;

Seg* segment;

ViterbiAccum& va=createAndInitializeViterbiAccum(ss,actualHMM);        // CREATE AND INITIALIZE the Viterbi accumulator
initializeCluster(segTemp,actualHMM,labelServer);                      // Initialize the segmentation (empty with one cluster by state/speaker)
if(cluster.getCount()==0){
	cout<<"ERROR: Not found Segments to decode !"<<endl;
	exit(-1);
}

cluster.rewind();
unsigned long iprob=0;
while((segment=cluster.getSeg())!=NULL){                               // Loop on each segment comming from the macro-class acoustic segmentation
//	cout << "Segment to decode => " << segment->string() << " " << segment->begin() << endl;
	va.reset();	                                                         // Reset the viterbi accumulator
	accumulateStatViterbi(fs,va,segment,config);                             // Accumulate the stats for a segment
	const ULongVector& path=va.getPath();                                // Get the path;
/*	cout << "Path: " ;
	for(unsigned long i=0; i<path.size(); i++) 
		cout << path[i] << " " ;
	cout << endl;*/
	double llp=va.getLlp();   
	                                           // Get the log likelihood of the path
	actualViterbiProb[iprob]=llp;
	if((verbose) && (verboseLevel==2)) cout << "Path probability: " << llp << endl;
	iprob++;
	copyPathInCluster(config, segTemp,path,actualHMM,labelServer,segment->begin(),segment->sourceName()); // Put the viterbi path (for the current macro segment)  in the cluster
}//while segments	
}//Viterbi


//-------------------------------------------------------------------------
MixtureGD &mixtureInitBySplit(Config& config, MixtureServer &ms, StatServer &ss, FeatureServer &fs, 
		SegCluster& selectedSegments, DoubleVector & globalMean,DoubleVector &globalCov,  unsigned long maxDistribCount, TrainCfg &trainCfg, unsigned long deleteDistrib){ 
  float varianceFlooring = trainCfg.getInitVarFloor();
  float varianceCeiling  = trainCfg.getInitVarCeil();   

 	 /* prepare the first gaussian */
	  MixtureGD &world=ms.createMixtureGD(0);
	  DistribGD &m=ms.createDistribGD();
	  unsigned long vectSize=fs.getVectSize();
	  if (verbose) trainCfg.showConfig(cout);
  	
	  for(unsigned long i=0; i<vectSize; i++){
		m.setMean(globalMean[i], i);
		m.setCov(globalCov[i], i);
	  }
	  m.computeAll();
	  ms.addDistribToMixture(world, m, 1.0);
  
	  if((verbose) && (verboseLevel > 1)){
		  cout << "Initialization of the first Gaussian: " << world.toString() << endl;
		  for(unsigned long d=0; d<world.getDistribCount(); d++){
	  		cout << world.getDistrib(d).toString() << endl;
	  	}
	  }

	  while((2*world.getDistribCount()) <= maxDistribCount){
  		/* Each gaussian will be split */
		if(verbose) cout << "Gaussian Nb: " << world.getDistribCount() << endl;
		unsigned long distribCount=world.getDistribCount();
		for(unsigned long d=0; d<distribCount; d++){
			/* Gaussian to split */
			DistribGD &m1=world.getDistrib(d);
			/* New gaussian */
			DistribGD &m2=ms.duplicateDistrib(m1);
		
			for(unsigned long i=0; i<vectSize; i++){
				/* first gaussian mean is fixed to mean(gaussian to split)+squared(cov(gaussian to split))
				   second gaussian mean is fixed to mean(gaussian to split)-squared(cov(gaussian to split)) */
				double mean=m1.getMean(i);
				m1.setMean(mean+sqrt(m1.getCov(i)), i);
				m2.setMean(mean-sqrt(m1.getCov(i)), i);
			}
			double w1=world.weight(d);
			world.weight(d)=w1/2;
			ms.addDistribToMixture(world, m2, w1/2);				
		}
	        world.computeAll();
	    
	    	if (verboseLevel>1) cout << "Train It : initial model: ll world = " <<  meanLikelihood(ss,fs,world,selectedSegments,config)<< endl;
	    	for (unsigned long trainIt=0; trainIt<trainCfg.getNbTrainIt(); trainIt++){                           // Begin the initial EM iteration loop 
			MixtureStat& emAcc=ss.createAndStoreMixtureStat(world);             // Create a statistic accumulator using the init model
			// Set the variance control parameters
		        varianceFlooring = setItParameter(trainCfg.getInitVarFloor(),trainCfg.getFinalVarFloor(),trainCfg.getNbTrainIt(),trainIt);
			varianceCeiling  = setItParameter(trainCfg.getInitVarCeil(),trainCfg.getFinalVarCeil(), trainCfg.getNbTrainIt(), trainIt);
	   	        if (verboseLevel>2) cout << "Train it["<<trainIt<<"], Variance floor["<<varianceFlooring
			<<"] ceiling["<<varianceCeiling<<"]"<<endl; 
	
			emAcc.resetEM();                                                                                    // EM stuff
			double llkPreviousIt=accumulateStatEM(ss,fs,emAcc,selectedSegments,config);                        // Compute EM statistics
			llkPreviousIt=llkPreviousIt/(double) emAcc.getEMFeatureCount();
			world = emAcc.getEM();                                                                             // Get the EM estimate
			varianceControl(world,varianceFlooring,varianceCeiling,globalCov);
			if (trainCfg.getNormalizeModel()) normalizeMixture(world,config);	
			if (verbose) cout << "Partial Train it["<<trainIt<<"] (take care, it corresponds to the previous it,0 means init likelihood) ="<<llkPreviousIt<<" Nb Frames="<<emAcc.getEMFeatureCount()<<endl;
			if (verboseLevel>2) cout << "Train it["<<trainIt <<"] Complete LLK=[" << meanLikelihood(ss,fs,world,selectedSegments,config)<< "]"<<endl;       
			ss.deleteMixtureStat(emAcc);
		}
		if (verboseLevel>1) cout << " Final ll world = " << meanLikelihood(ss,fs,world,selectedSegments,config)<< endl;
		if(verbose && verboseLevel > 2){
			cout << world.toString() << endl;
		        for(unsigned long d=0; d<world.getDistribCount(); d++){
	  			cout << world.getDistrib(d).toString() << endl;
  			}
		}
  	  }
	  while(world.getDistribCount() < maxDistribCount){
		if(verbose) cout << "Gaussian Nb (unitary split): " << world.getDistribCount() << endl;
		if(verbose) cout << world.toString() << endl;
	  	double maxWeight=0.0;
		unsigned long indD=0;
		unsigned long distribCount=world.getDistribCount();
		
		for(unsigned long d=0; d<distribCount; d++){
			if(world.weight(d) > maxWeight){
				maxWeight=world.weight(d);
				indD=d;
			}
		}
		DistribGD &m1=world.getDistrib(indD);
		/* New gaussian */
		DistribGD &m2=ms.duplicateDistrib(m1);
		
		for(unsigned long i=0; i<vectSize; i++){
			/* first gaussian mean is fixed to mean(gaussian to split)+squared(cov(gaussian to split))
			   second gaussian mean is fixed to mean(gaussian to split)-squared(cov(gaussian to split)) */
			double mean=m1.getMean(i);
			m1.setMean(mean+sqrt(m1.getCov(i)), i);
			m2.setMean(mean-sqrt(m1.getCov(i)), i);
		}
		double w1=world.weight(indD);
		world.weight(indD)=w1/2;
		ms.addDistribToMixture(world, m2, w1/2);				
		
	        world.computeAll();
	    
	    	if (verboseLevel>1) cout << "Train It : initial model: ll world = " <<  meanLikelihood(ss,fs,world,selectedSegments,config)<< endl;
	    	for (unsigned long trainIt=0; trainIt<trainCfg.getNbTrainIt(); trainIt++){                           // Begin the initial EM iteration loop 
			MixtureStat& emAcc=ss.createAndStoreMixtureStat(world);             // Create a statistic accumulator using the init model
			// Set the variance control parameters
		        varianceFlooring = setItParameter(trainCfg.getInitVarFloor(),trainCfg.getFinalVarFloor(),trainCfg.getNbTrainIt(),trainIt);
			varianceCeiling  = setItParameter(trainCfg.getInitVarCeil(),trainCfg.getFinalVarCeil(), trainCfg.getNbTrainIt(), trainIt);
	   	        if (verboseLevel>2) cout << "Train it["<<trainIt<<"], Variance floor["<<varianceFlooring
			<<"] ceiling["<<varianceCeiling<<"]"<<endl; 
	
			emAcc.resetEM();                                                                                    // EM stuff
			double llkPreviousIt=accumulateStatEM(ss,fs,emAcc,selectedSegments,config);                        // Compute EM statistics
			llkPreviousIt=llkPreviousIt/(double) emAcc.getEMFeatureCount();
			world = emAcc.getEM();                                                                             // Get the EM estimate
			varianceControl(world,varianceFlooring,varianceCeiling,globalCov);
			if (trainCfg.getNormalizeModel()) normalizeMixture(world,config);	
			if (verbose) cout << "Partial Train it["<<trainIt<<"] (take care, it corresponds to the previous it,0 means init likelihood) ="<<llkPreviousIt<<" Nb Frames="<<emAcc.getEMFeatureCount()<<endl;
			if (verboseLevel>2) cout << "Train it["<<trainIt <<"] Complete LLK=[" << meanLikelihood(ss,fs,world,selectedSegments,config)<< "]"<<endl;       
			ss.deleteMixtureStat(emAcc);
		}
		if (verboseLevel>1) cout << " Final ll world = " << meanLikelihood(ss,fs,world,selectedSegments,config)<< endl;
		if(verbose && verboseLevel > 2){
			cout << world.toString() << endl;
		        for(unsigned long d=0; d<world.getDistribCount(); d++){
	  			cout << world.getDistrib(d).toString() << endl;
  			}
		}
	  }
	  if(verbose) cout << "Final Gaussian Nb: " << world.getDistribCount() << endl;
	  
	  if(deleteDistrib==1){
        	MixtureGD &reducedWorld=ms.createMixtureGD(0);

		for(unsigned long d=0; d<world.getDistribCount(); d++){
			double w=world.weight(d);
			if(w*totalFrame(selectedSegments) > 10){
				DistribGD &m1=world.getDistrib(d);
				/* New gaussian */
				DistribGD &m2=ms.duplicateDistrib(m1);
				ms.addDistribToMixture(reducedWorld, m2, w);
			}
		}
		normalizeWeights(reducedWorld);
		reducedWorld.computeAll();
		if(verbose) cout << "World gaussian number goes from "<< world.getDistribCount() << " to " << reducedWorld.getDistribCount() << endl;
		return reducedWorld;
			  	
	  }
	  else 
	  	return world;
}


//-------------------------------------------------------------------------
void copyPathInCluster(Config& config,SegServer &segTemp,const ULongVector &path,hmm &cHmm,LabelServer &labelServer,unsigned long start, String sourceName){

unsigned long startV,stopV,indiceLoc;                                      // Begin/end of the current segment
indiceLoc=path[0]; startV=start;                                                 // did we begin a segment ?
for(unsigned long i=1;i<path.size();i++){
	if((path[i]!=indiceLoc)){                                 // End of a segment
		stopV=(i-1)+start;
		Label l;
		l.setString(cHmm.getStateName(indiceLoc));
		unsigned long cl=labelServer.addLabel(l);	
		SegCluster& clusterV=segTemp.getCluster(indiceLoc);
		clusterV.add(segTemp.createSeg(startV,stopV-startV+1,cl,cHmm.getStateName(indiceLoc),sourceName));
	
		indiceLoc=path[i];
		startV=i+start;
	} //if		
} //for
// Deal with the last segment if needed
stopV=(path.size()-1)+start;
Label l;
l.setString(cHmm.getStateName(indiceLoc));
unsigned long cl=labelServer.addLabel(l);	
SegCluster& clusterV=segTemp.getCluster(indiceLoc);
clusterV.add(segTemp.createSeg(startV,stopV-startV+1,cl,cHmm.getStateName(indiceLoc),sourceName));
}//copyPath


/**********************************************************
* createWorld: 
***********************************************************/
MixtureGD& createWorld(Config& config,SegCluster& cluster,StatServer& ss,FeatureServer &fs,MixtureServer& ms){
	if((config.existsParam("loadWorldFromExternalFile")) && (config.getParam("loadWorldFromExternalFile").toBool()) == true){
		String worldFile=config.getParam("worldModel");
		MixtureGD &world=ms.loadMixtureGD(worldFile);
		return world;
	}
	else{
		unsigned long maxDistribCount=config.getParam("mixtureDistribCount").toLong();
		MixtureGD &world=learnModelFromSplit(config,cluster,ss,fs,ms,maxDistribCount,0,20);
		return world;
	}
}


MixtureGD& createWorld(Config& config,SegCluster& cluster,StatServer& ss,FeatureServer &fs,MixtureServer& ms, unsigned long maxDistribCount){
	if((config.existsParam("loadWorldFromExternalFile")) && (config.getParam("loadWorldFromExternalFile").toBool()) == true){
		String worldFile=config.getParam("worldModel");
		MixtureGD &world=ms.loadMixtureGD(worldFile);
		return world;
	}
	else{
		MixtureGD &world=learnModelFromSplit(config,cluster,ss,fs,ms,maxDistribCount,0,20);
		return world;
	}
}


/****************************************
* segAdaptation: Adapts all the speaker/gmm models of a EHMM depending on the current segmentation
****************************************/

void segAdaptation(Config& config,MAPCfg &mapCfg, hmm& actualHMM,MixtureGD& world,SegServer& previousSeg,SegServer& actualSeg,
		   StatServer& ss,FeatureServer &fs,MixtureServer& ms){
  
NoDataSpeakerVerification(config,actualHMM,actualSeg);

for(unsigned long icluster=0;icluster<actualSeg.getClusterCount();icluster++){ // For each cluster
    
	SegCluster& clusterA=actualSeg.getCluster(icluster);         // Get the current cluster of segments 
	
	if (verbose) cout<<" Adapt segments for speaker "<<clusterA.string()<<endl;
	if(previousSeg.getClusterCount() != 0){
		SegCluster& clusterP=previousSeg.getCluster(icluster);
		if(isDifferentSegmentation(config, clusterP, clusterA)){
			actualHMM.setDensity(world, icluster);		
			if ((verbose) && (verboseLevel==2)) cout<<" altered segmentation for this speaker => adaptation "<<icluster<<endl;
			MixtureGD& m=(MixtureGD&)actualHMM.getDensity(icluster);
			adaptModel(config,ss,ms,fs,clusterA,world,m, mapCfg);
      		}	
	      	else{
      			if ((verbose) && (verboseLevel==2)) cout<<" identical segmentation for this speaker => no adaptation "<<icluster<<endl;		 
      		}
	}
	else{// 1st  adaptation/decodage step 
		actualHMM.setDensity(world, icluster);
		MixtureGD& m=(MixtureGD&)actualHMM.getDensity(icluster);
      
		adaptModel(config,ss,ms,fs,clusterA,world,m, mapCfg);
	}
}//for icluster
}//segAdaptation


/****************************************
* segEM: Train all the speaker/gmm models of a EHMM depending on the current segmentation based on the EM algorithm
****************************************/

void segEM(Config& config,hmm& actualHMM,SegServer& previousSeg,SegServer& actualSeg,
		   StatServer& ss,FeatureServer &fs,MixtureServer& ms,unsigned long nbDistrib){
  
NoDataSpeakerVerification(config,actualHMM,actualSeg);

for(unsigned long icluster=0;icluster<actualSeg.getClusterCount();icluster++){ // For each cluster
    
	SegCluster& clusterA=actualSeg.getCluster(icluster);         // Get the current cluster of segments 

	if(icluster==actualSeg.getClusterCount()-1) nbDistrib /= 2 ;

	if (verbose) cout<<" Adapt segments for speaker "<<clusterA.string()<< " NbDistrib: " << nbDistrib << endl;
	if(previousSeg.getClusterCount() != 0){
		SegCluster& clusterP=previousSeg.getCluster(icluster);
		if(isDifferentSegmentation(config, clusterP, clusterA)){
			if ((verbose) && (verboseLevel==2)) cout<<" altered segmentation for this speaker => adaptation "<<icluster<<endl;
			MixtureGD& m=learnModelFromSplit(config,clusterA,ss,fs,ms,nbDistrib,0,5);
			actualHMM.setDensity(m, icluster);
      		}	
	      	else{
      			if ((verbose) && (verboseLevel==2)) cout<<" identical segmentation for this speaker => no adaptation "<<icluster<<endl;		 
      		}
	}
	else{// 1st  adaptation/decodage step 
		MixtureGD& m=learnModelFromSplit(config,clusterA,ss,fs,ms,nbDistrib,0,5);
		actualHMM.setDensity(m, icluster);
	
	}
}//for icluster

}//segEM


#endif

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

#include<iostream>
#include<fstream>
#include<cstdio>
#include<cassert>
#include<cmath>
#include "TrainTools.h"
#include "SegTools.h"

using namespace alize;
using namespace std;



/**********************************************************
* clrCrit: Application of CLR criterion for clustering
***********************************************************/
double clrCrit(Config& config, SegCluster& c1, SegCluster &c2, StatServer& ss, FeatureServer& fs, MixtureGD& m1, MixtureGD& m2, MixtureGD& world){

if(verbose) cout << "Likelihood computation" << endl;	
MixtureStat &m1c2Acc=ss.createAndStoreMixtureStat(m1);                    // Create a statistic accumulator
MixtureStat &m2c1Acc=ss.createAndStoreMixtureStat(m2);                    // Crbeeate a statistic accumulator
MixtureStat &mWc1Acc=ss.createAndStoreMixtureStat(world);                    // Create a statistic accumulator
MixtureStat &mWc2Acc=ss.createAndStoreMixtureStat(world);                    // Create a statistic accumulator

accumulateStatLLK(ss,fs,m1c2Acc,c2,config);
accumulateStatLLK(ss,fs,mWc2Acc,c2,config);
accumulateStatLLK(ss,fs,m2c1Acc,c1,config);
accumulateStatLLK(ss,fs,mWc1Acc,c1,config);

double m1c2LLK = m1c2Acc.getMeanLLK();
double m2c1LLK = m2c1Acc.getMeanLLK();
double mWc1LLK = mWc1Acc.getMeanLLK();
double mWc2LLK = mWc2Acc.getMeanLLK();

if ((verbose) && (verboseLevel == 2)){
	cout << "LLKm2c1 ("<<totalFrame(c1)<<") => " << m2c1LLK << endl;
	cout << "LLKWc1 ("<<totalFrame(c1)<<") => " << mWc1LLK << endl;
	cout << "LLKm1c2 ("<<totalFrame(c2)<<") => " << m1c2LLK << endl;
	cout << "LLKWc2 ("<<totalFrame(c2)<<") => " << mWc2LLK << endl;	
}

double clr = (m1c2LLK - mWc2LLK) + (m2c1LLK - mWc1LLK);
if (verbose) cout << "CLR ("<<totalFrame(c1)<<"/"<<totalFrame(c2)<< ") => " << clr << endl;
return clr;
}

/**********************************************************
* glrCrit: Application of GLR criterion for clustering
***********************************************************/
double gllrCrit(Config& config, SegCluster& c1, SegCluster &c2, SegCluster &c12, StatServer& ss, FeatureServer& fs, MixtureGD& m1, MixtureGD& m2, MixtureGD& m12){

MixtureStat &m1Acc=ss.createAndStoreMixtureStat(m1);                    // Create a statistic accumulator
MixtureStat &m2Acc=ss.createAndStoreMixtureStat(m2);                    // Create a statistic accumulator
MixtureStat &m12Acc=ss.createAndStoreMixtureStat(m12);                    // Create a statistic accumulator

accumulateStatLLK(ss,fs,m1Acc,c1,config);
accumulateStatLLK(ss,fs,m2Acc,c2,config);
accumulateStatLLK(ss,fs,m12Acc,c12,config);

double m1LLK = m1Acc.getAccumulatedLLK();
double m2LLK = m2Acc.getAccumulatedLLK();
double m12LLK = m12Acc.getAccumulatedLLK();
if ((verbose) && (verboseLevel == 2)){
	cout << "LLK1 => " << m1LLK << endl;
	cout << "LLK2 => " << m2LLK << endl;
	cout << "LLK12 => " << m12LLK << endl;
}
double gllr = m12LLK - (m1LLK + m2LLK);
if (verbose) cout << "GLLR ("<<totalFrame(c1)<<"/"<<totalFrame(c2)<< ") => " << gllr << endl;
return gllr;
}

/**********************************************************
* bicCrit: Application of BIC criterion for clustering
***********************************************************/
double bicCrit(Config& config, SegCluster& c1, SegCluster &c2, SegCluster &c12, StatServer& ss, FeatureServer& fs, MixtureGD& m1, MixtureGD& m2, MixtureGD& m12){

double gllr = gllrCrit(config, c1, c2, c12, ss, fs, m1, m2, m12);
double lambda = 1.0;
double P = 0.5 * ((2 * fs.getVectSize()+1)*m1.getDistribCount()) * log((double)totalFrame(c1)+totalFrame(c2));
double bic = (-1 * gllr) - lambda * P;
if (verbose) cout << "BIC ("<<totalFrame(c1)<<"/"<<totalFrame(c2)<< ") => " << bic << " P => " << P << endl;

return bic;	
}

/**********************************************************
* bicCrit: Application of Delta BIC criterion for clustering (IDIAP)
***********************************************************/
double deltabicCrit(Config& config, SegCluster& c1, SegCluster &c2, SegCluster &c12, StatServer& ss, FeatureServer& fs, MixtureGD& m1, MixtureGD& m2, MixtureGD& m12){

double deltabic = gllrCrit(config, c1, c2, c12, ss, fs, m1, m2, m12);
if (verbose) cout << "Delta BIC (=GLLR) ("<<totalFrame(c1)<<"/"<<totalFrame(c2)<< ") => " << deltabic << endl;

return deltabic;	
}

/**********************************************************
* clusteringCriterionByAdapt: Application of different criterion for clustering (models are trained by adaptation)
***********************************************************/
double clusteringCriterionByAdapt(Config& config, Seg *segment1, Seg *segment2, StatServer& ss, FeatureServer& fs, MixtureGD& world, String crit){

if(verbose){
	cout << "Computation between: " << segment1->begin() << " " << endSeg(segment1); 
	cout << " and " << segment2->begin() << " " << endSeg(segment2) << endl; 
}

SegServer segTemp;
MixtureServer ms(config); 	
MAPCfg mapCfg(config);
mapCfg.setMethod("MAPOccDep");
mapCfg.setMeanReg(16);
mapCfg.setBaggedFrameProbability(1.0);	

if(verbose) cout << "Mixture building m1 and m2" << endl; 
SegCluster& c1=segTemp.createCluster();
c1.add(segTemp.createSeg(segment1->begin(),segment1->length(),0,"null",segment1->sourceName()));
SegCluster& c2=segTemp.createCluster();
c2.add(segTemp.createSeg(segment2->begin(),segment2->length(),0,"null",segment2->sourceName()));

MixtureGD &m1=ms.duplicateMixture(world,DUPL_DISTRIB); 	
MixtureGD &m2=ms.duplicateMixture(world,DUPL_DISTRIB); 	
adaptModel(config,ss,ms,fs,c1,world,m1, mapCfg);
adaptModel(config,ss,ms,fs,c2,world,m2, mapCfg);

if((crit=="GLR") || (crit=="BIC")){
	if(verbose) cout << "Mixture building m12" << endl; 
	SegCluster& c12=segTemp.createCluster();
	c12.add(segTemp.createSeg(segment1->begin(),segment1->length(),0,"null",segment1->sourceName()));
	c12.add(segTemp.createSeg(segment2->begin(),segment2->length(),0,"null",segment2->sourceName()));

	MixtureGD &m12=ms.duplicateMixture(world,DUPL_DISTRIB); 	
	adaptModel(config,ss,ms,fs,c12,world,m12, mapCfg);	

	if(crit=="GLR"){
		return gllrCrit(config, c1, c2, c12, ss, fs, m1, m2, m12);
	}
	if(crit=="BIC"){
		return bicCrit(config, c1, c2, c12, ss, fs, m1, m2, m12);
	}
}

else if(crit == "CLR"){
	return clrCrit(config, c1, c2, ss, fs, m1, m2, world);
}
else{
	cout << "ERROR: unknown clustering criterion !!!!" << endl;
	return -1;
}
return -1;

}


/**********************************************************
* clusteringCriterion: Application of different criterion for clustering (models are trained by EM/ML - Initialization by external model)
***********************************************************/
double clusteringCriterion(Config& config, Seg *segment1, Seg *segment2, StatServer& ss, FeatureServer& fs, MixtureGD& world, String crit){

	
if(verbose){
	cout << "Computation between: " << segment1->begin() << " " << endSeg(segment1); 
	cout << " and " << segment2->begin() << " " << endSeg(segment2) << endl; 
}

if(verbose) cout << "Mixture building m1 and m2" << endl;
TrainCfg trainCfg(config);
trainCfg.setBaggedFrameProbability(0.8);
trainCfg.setNbTrainIt(10);

SegServer segTemp;
SegCluster& c1=segTemp.createCluster();
c1.add(segTemp.createSeg(segment1->begin(),segment1->length(),0,"null",segment1->sourceName()));
SegCluster& c2=segTemp.createCluster();
c2.add(segTemp.createSeg(segment2->begin(),segment2->length(),0,"null",segment2->sourceName()));

MixtureServer ms(config); 	
MixtureGD &m1=ms.duplicateMixture(world,DUPL_DISTRIB); 	
MixtureGD &m2=ms.duplicateMixture(world,DUPL_DISTRIB); 	

FrameAccGD globalFrameAcc1;                                                    //compute global Mean and cov
globalMeanCov (fs,c1,globalFrameAcc1,config);                             // Compute the global mean and covariance
DoubleVector globalMean1=globalFrameAcc1.getMeanVect(); 
DoubleVector globalCov1=globalFrameAcc1.getCovVect();
FrameAccGD globalFrameAcc2;                                                    //compute global Mean and cov
globalMeanCov (fs,c2,globalFrameAcc2,config);                             // Compute the global mean and covariance
DoubleVector globalMean2=globalFrameAcc2.getMeanVect(); 
DoubleVector globalCov2=globalFrameAcc2.getCovVect();

if(verbose) cout<<" Training EM for speaker m1"<<endl;	
trainModel(config,ss,fs,c1,globalMean1,globalCov1,m1, trainCfg);            // Train it by EM/ML
if(verbose) cout<<" Training EM for speaker m2"<<endl;	
trainModel(config,ss,fs,c2,globalMean2,globalCov2,m2, trainCfg);            // Train it by EM/ML

if((crit=="GLR") || (crit=="BIC")){
	if(verbose) cout << "Mixture building m12" << endl;
	SegCluster& c12=segTemp.createCluster();
	c12.add(segTemp.createSeg(segment1->begin(),segment1->length(),0,"null",segment1->sourceName()));
	c12.add(segTemp.createSeg(segment2->begin(),segment2->length(),0,"null",segment2->sourceName()));
	MixtureGD &m12=ms.duplicateMixture(world,DUPL_DISTRIB); 	
	FrameAccGD globalFrameAcc12;                                                    //compute global Mean and cov
	globalMeanCov (fs,c12,globalFrameAcc12,config);                             // Compute the global mean and covariance
	DoubleVector globalMean12=globalFrameAcc12.getMeanVect(); 
	DoubleVector globalCov12=globalFrameAcc12.getCovVect();
	if(verbose) cout<<" Training EM for speaker m12"<<endl;	
	trainModel(config,ss,fs,c12,globalMean12,globalCov12,m12, trainCfg);            // Train it by EM/ML

	if(crit=="GLR"){
		return gllrCrit(config, c1, c2, c12, ss, fs, m1, m2, m12);
	}
	if(crit=="BIC"){
		return bicCrit(config, c1, c2, c12, ss, fs, m1, m2, m12);
	}	
}

else if(crit == "CLR"){
	return clrCrit(config, c1, c2, ss, fs, m1, m2, world);
}
else{
	cout << "ERROR: unknown clustering criterion !!!!" << endl;
	return -1;
}
return -1;
}


/**********************************************************
* mergeCluster: merge two clusters
***********************************************************/
SegCluster& mergeCluster(SegCluster& c1, SegCluster& c2, SegServer& segTemp, String merge="NULL"){

SegCluster& c12=segTemp.createCluster();
c12.setString(merge);
Seg *segment;
c1.rewind();
while((segment=c1.getSeg()) != NULL){
	c12.add(segTemp.createSeg(segment->begin(),segment->length(),0,merge,segment->sourceName()));
}
c2.rewind();
while((segment=c2.getSeg()) != NULL){
	c12.add(segTemp.createSeg(segment->begin(),segment->length(),0,merge,segment->sourceName()));
}
return c12;
}

/**********************************************************
* clusteringCriterionWithoutWorldInit: Application of different criterion for clustering (models are trained by EM/ML - Initialization by external model)
***********************************************************/
double clusteringCriterionWithoutWorldInit(Config& config, SegCluster& c1, SegCluster& c2, StatServer& ss, 
		FeatureServer& fs,MixtureGD& world,String crit){
	
double criterion=0.0;
SegServer segTemp;
MixtureServer ms(config); 

if(verbose) cout << "Mixture building m1 and m2" << endl;
TrainCfg trainCfg(config);
trainCfg.setInitVarFlooring(0.0);
trainCfg.setFinalVarFlooring(0.0);


FrameAccGD globalFrameAcc1;                                                    //compute global Mean and cov
globalMeanCov (fs,c1,globalFrameAcc1,config);                             // Compute the global mean and covariance
DoubleVector globalMean1=globalFrameAcc1.getMeanVect(); 
DoubleVector globalCov1=globalFrameAcc1.getCovVect();
FrameAccGD globalFrameAcc2;                                                    //compute global Mean and cov
globalMeanCov (fs,c2,globalFrameAcc2,config);                             // Compute the global mean and covariance
DoubleVector globalMean2=globalFrameAcc2.getMeanVect(); 
DoubleVector globalCov2=globalFrameAcc2.getCovVect();


if(crit == "DELTABIC"){
	if (verbose) cout <<"m1 and m2 model init from scratch (hola !)"<<endl;
	trainCfg.setBaggedFrameProbabilityInit(0.4);

	int gausNb1=5;
	int gausNb2=5;
	
	/*int gausNb1 = totalFrame(c1)/800;
	int gausNb2 = totalFrame(c2)/800;
	
	if(gausNb1 == 0) gausNb1=1;
	if(gausNb2 == 0) gausNb2=1;*/
	
	MixtureGD &m1=ms.createMixtureGD(gausNb1);
	MixtureGD &m2=ms.createMixtureGD(gausNb2);
	trainCfg.setNbTrainIt(1);

	mixtureInit(ms,fs,m1,c1,globalCov1,config,trainCfg);                   // Init the model
	mixtureInit(ms,fs,m2,c2,globalCov2,config,trainCfg);                   // Init the model

	trainCfg.setNbTrainIt(5);
	trainCfg.setBaggedFrameProbability(1.0);
	if(verbose) cout<<" Training EM for speaker m1 (nbGauss:"<<gausNb1<<")"<<endl;	
	trainModel(config,ss,fs,c1,globalMean1,globalCov1,m1, trainCfg);            // Train it by EM/ML
	if(verbose) cout<<" Training EM for speaker m2 (nbGauss:"<<gausNb2<<")"<<endl;	
	trainModel(config,ss,fs,c2,globalMean2,globalCov2,m2, trainCfg);            // Train it by EM/ML

	if(verbose) cout << "Mixture building m12" << endl;
	SegCluster& c12=mergeCluster(c1, c2, segTemp);
	
	int gausNb12=gausNb1+gausNb2;
	MixtureGD &m12=ms.createMixtureGD(gausNb12);

	FrameAccGD globalFrameAcc12;                                                    //compute global Mean and cov
	globalMeanCov (fs,c12,globalFrameAcc12,config);                             // Compute the global mean and covariance
	DoubleVector globalMean12=globalFrameAcc12.getMeanVect(); 
	DoubleVector globalCov12=globalFrameAcc12.getCovVect();

	if (verbose) cout <<"m12 model init from scratch"<<endl;
	trainCfg.setBaggedFrameProbabilityInit(0.4);
	trainCfg.setBaggedFrameProbability(0.4);
	trainCfg.setNbTrainIt(1);
	mixtureInit(ms,fs,m12,c12,globalCov12,config,trainCfg);                   // Init the model

	if(verbose) cout<<" Training EM for speaker m12"<<endl;	
	trainCfg.setNbTrainIt(5);
	trainCfg.setBaggedFrameProbability(1.0);
	trainModel(config,ss,fs,c12,globalMean12,globalCov12,m12, trainCfg);            // Train it by EM/ML

	criterion = deltabicCrit(config, c1, c2, c12, ss, fs, m1, m2, m12);
	ms.deleteMixture(m1);
	ms.deleteMixture(m2);
	ms.deleteMixture(m12);
	ms.deleteUnusedDistribs();

	return criterion;
}
else{
	if (verbose) cout <<"m1 and m2 model init from scratch"<<endl;
	trainCfg.setBaggedFrameProbabilityInit(0.1);
	trainCfg.setBaggedFrameProbability(0.4);
	MixtureGD &m1=ms.createMixtureGD();
	MixtureGD &m2=ms.createMixtureGD();

	mixtureInit(ms,fs,m1,c1,globalCov1,config,trainCfg);                   // Init the model
	mixtureInit(ms,fs,m2,c2,globalCov2,config,trainCfg);                   // Init the model

	trainCfg.setNbTrainIt(10);
	trainCfg.setBaggedFrameProbability(1.0);
	if(verbose) cout<<" Training EM for speaker m1"<<endl;	
	trainModel(config,ss,fs,c1,globalMean1,globalCov1,m1, trainCfg);            // Train it by EM/ML
	if(verbose) cout<<" Training EM for speaker m2"<<endl;	
	trainModel(config,ss,fs,c2,globalMean2,globalCov2,m2, trainCfg);            // Train it by EM/ML


	if((crit=="GLR") || (crit=="BIC")){
		if(verbose) cout << "Mixture building m12" << endl;
		SegCluster& c12=mergeCluster(c1, c2, segTemp);

		FrameAccGD globalFrameAcc12;                                                    //compute global Mean and cov
		globalMeanCov (fs,c12,globalFrameAcc12,config);                             // Compute the global mean and covariance
		DoubleVector globalMean12=globalFrameAcc12.getMeanVect(); 
		DoubleVector globalCov12=globalFrameAcc12.getCovVect();

		if (verbose) cout <<"m12 model init from scratch"<<endl;	
		MixtureGD &m12=ms.createMixtureGD();
		trainCfg.setBaggedFrameProbabilityInit(0.1);
		trainCfg.setBaggedFrameProbability(0.4);
		trainCfg.setNbTrainIt(2);
		mixtureInit(ms,fs,m12,c12,globalCov12,config,trainCfg);                   // Init the model
		if(verbose) cout<<" Training EM for speaker m12"<<endl;	
		trainCfg.setNbTrainIt(10);
		trainCfg.setBaggedFrameProbability(1.0);
		trainModel(config,ss,fs,c12,globalMean12,globalCov12,m12, trainCfg);            // Train it by EM/ML

		if(crit=="GLR"){
			criterion = gllrCrit(config, c1, c2, c12, ss, fs, m1, m2, m12);
			ms.deleteMixture(m1);
			ms.deleteMixture(m2);
			ms.deleteMixture(m12);
			ms.deleteUnusedDistribs();
			return criterion;
		}
		if(crit=="BIC"){
			criterion = bicCrit(config, c1, c2, c12, ss, fs, m1, m2, m12);
			ms.deleteMixture(m1);
			ms.deleteMixture(m2);
			ms.deleteMixture(m12);
			ms.deleteUnusedDistribs();
			return criterion;
		}		
	}

	else if(crit == "CLR"){
		criterion = clrCrit(config, c1, c2, ss, fs, m1, m2, world);
		ms.deleteMixture(m1);
		ms.deleteMixture(m2);
		ms.deleteUnusedDistribs();
		return criterion;
	}
	else{
		cout << "ERROR: unknown clustering criterion !!!!" << endl;
		return -1;
	}
}

return -1;
}

/**********************************************************
* clusteringCriterionWithoutWorldInit: Application of different criterion for clustering (models are trained by EM/ML - Initialization by external model)
***********************************************************/
double clusteringCriterionWithoutWorldInit(Config& config, SegCluster& c1, MixtureGD& m1, SegCluster& c2, MixtureGD& m2, StatServer& ss, 
		FeatureServer& fs,MixtureGD& world,String crit){
	
SegServer segTemp;

TrainCfg trainCfg(config);
trainCfg.setBaggedFrameProbability(0.8);
trainCfg.setNbTrainIt(10);
trainCfg.setInitVarFlooring(0.0);
trainCfg.setFinalVarFlooring(0.0);

if(crit == "DELTABIC"){

	SegCluster& c12=mergeCluster(c1, c2, segTemp);
	
	MixtureServer msc12(config); 
	MixtureGD &m12=msc12.createMixtureGD((c1.getCount()+c2.getCount())*5);
	//MixtureGD &m12=msc12.createMixtureGD((c1.getCount()+c2.getCount()));
	
	FrameAccGD globalFrameAcc12;                                                    //compute global Mean and cov
	globalMeanCov (fs,c12,globalFrameAcc12,config);                             // Compute the global mean and covariance
	DoubleVector globalMean12=globalFrameAcc12.getMeanVect(); 
	DoubleVector globalCov12=globalFrameAcc12.getCovVect();

	if(totalFrame(c12) < 150){
		trainCfg.setBaggedFrameProbabilityInit(1.0);
		trainCfg.setBaggedFrameProbability(1.0);
	}else{
		trainCfg.setBaggedFrameProbabilityInit(0.4);
		trainCfg.setBaggedFrameProbability(0.4);	
	}
	trainCfg.setNbTrainIt(2);
	mixtureInit(msc12,fs,m12,c12,globalCov12,config,trainCfg);                   // Init the model

	trainCfg.setNbTrainIt(5);
	trainCfg.setBaggedFrameProbability(1.0);
	trainModel(config,ss,fs,c12,globalMean12,globalCov12,m12, trainCfg);            // Train it by EM/ML

	return deltabicCrit(config, c1, c2, c12, ss, fs, m1, m2, m12);

}
else{
	if((crit=="GLR") || (crit=="BIC")){
		if(verbose) cout << "Mixture building m12" << endl;
		MixtureServer ms(config); 	
		SegCluster& c12=mergeCluster(c1, c2, segTemp);

		FrameAccGD globalFrameAcc12;                                                    //compute global Mean and cov
		globalMeanCov (fs,c12,globalFrameAcc12,config);                             // Compute the global mean and covariance
		DoubleVector globalMean12=globalFrameAcc12.getMeanVect(); 
		DoubleVector globalCov12=globalFrameAcc12.getCovVect();

		if (verbose) cout <<"m12 model init from scratch"<<endl;	
		MixtureGD &m12=ms.createMixtureGD(c1.getCount());
		trainCfg.setBaggedFrameProbabilityInit(0.4);
		trainCfg.setBaggedFrameProbability(0.4);
		trainCfg.setNbTrainIt(2);
		mixtureInit(ms,fs,m12,c12,globalCov12,config,trainCfg);                   // Init the model
		if(verbose) cout<<" Training EM for speaker m12"<<endl;	
		trainCfg.setNbTrainIt(5);
		trainCfg.setBaggedFrameProbability(1.0);
		trainModel(config,ss,fs,c12,globalMean12,globalCov12,m12, trainCfg);            // Train it by EM/ML

		if(crit=="GLR"){
			return gllrCrit(config, c1, c2, c12, ss, fs, m1, m2, m12);
		}
		if(crit=="BIC"){
			return bicCrit(config, c1, c2, c12, ss, fs, m1, m2, m12);
		}		
	}

	else if(crit == "CLR"){
		return clrCrit(config, c1, c2, ss, fs, m1, m2, world);
	}
	else{
		cout << "ERROR: unknown clustering criterion !!!!" << endl;
		return -1;
	}
}
return -1;
}

/**********************************************************
* clusteringCriterion: Application of different criterion for clustering (models are trained by EM/ML - Initialization by external model)
***********************************************************/
double clusteringCriterionWithoutWorldInit(Config& config, Seg *segment1, Seg *segment2, StatServer& ss, FeatureServer& fs,MixtureGD& world,String crit){

SegServer segTemp;	
if(verbose){
	cout << "Computation between: " << segment1->begin() << " " << endSeg(segment1); 
	cout << " and " << segment2->begin() << " " << endSeg(segment2) << endl; 
}

SegCluster& c1=segTemp.createCluster();
c1.add(segTemp.createSeg(segment1->begin(),segment1->length(),0,"null",segment1->sourceName()));
SegCluster& c2=segTemp.createCluster();
c2.add(segTemp.createSeg(segment2->begin(),segment2->length(),0,"null",segment2->sourceName()));

return clusteringCriterionWithoutWorldInit(config, c1, c2, ss, fs,world,crit);

}


/**********************************************************
* isSimilarSegment: compare two segments with BIC, GLR, ..., criterions
***********************************************************/
bool isSimilarSegment(Config& config, Seg *segment1, Seg *segment2, StatServer& ss, FeatureServer& fs, MixtureGD& world, String crit){
	// if deltabic positive => two segments are similar
	// if deltabic negative => two segments are different
	double threshold=0.0;

	//double critValue=clusteringCriterionByAdapt(config, segment1, segment2, ss, fs, world, crit);
	double critValue=clusteringCriterionWithoutWorldInit(config, segment1, segment2, ss, fs, world, crit);
	if(config.existsParam("clusteringCritThresh"))
		threshold=config.getParam("clusteringCritThresh").toDouble();
	if((crit == "BIC") || (crit == "CLR") || (crit == "DELTABIC")) 
		if(critValue > threshold) return true;
		//else return false;
	else if(crit == "GLR")
		if(critValue < threshold) return true;
		//else return false;
return false; // avoid warning NS (check that!)
}



/**********************************************************
* cohortMaxLikelihood: compute the maximum likelihood over a set of cluster for a given segment
***********************************************************/
double cohortMaxLikelihood(StatServer& ss,FeatureServer &fs,hmm& actualHMM,unsigned long except,unsigned long begin,unsigned long longSelection,Config& config){
DoubleVector llr;

cout << "nb: " << actualHMM.getNbState() << endl;
for(unsigned long i=0; i<actualHMM.getNbState(); i++){	
	if(i!=except){
		MixtureGD& m=(MixtureGD&)actualHMM.getDensity(i);
		double mean = meanLikelihood(ss,fs,m,begin,longSelection,config);
		cout << "mean: " << mean << endl;
		llr.addValue(mean);	
	}
}
double max = llr[llr.getIndexOfLargestValue()];
if((verbose) && (verboseLevel == 2))
	cout << "Max of cohort: " << max << endl;
return max;
}


/**********************************************************
* bestFittingSegment: search the best segments of a cluster (in terms of likelihood) and return it
***********************************************************/
Seg *bestFittingSegment(Config& config, SegCluster& cluster, MixtureGD& m, StatServer& ss, FeatureServer& fs){

DoubleVector llr,starts;                                        // The LLR value and the startindex of each of the potential trial
cluster.rewind();
Seg *segment;	
while((segment=cluster.getSeg())!=NULL){
		String fileLabel=segment->sourceName();
		unsigned long begin=segment->begin()+fs.getFirstFeatureIndexOfASource(fileLabel);//?????????????
		
		double llrTmp = meanLikelihood(ss,fs,m,begin,segment->length(),config);
		if(verbose) cout << segment->begin() << " " << endSeg(segment) << " => " << llrTmp << endl;
		llr.addValue(llrTmp);
}//while segment

if(llr.size() == 0) return NULL;

bool goon=true;
while(goon){
	unsigned long ind = llr.getIndexOfLargestValue();
	if((verbose) && (verboseLevel == 2)) cout << "Best llk " << ind << " => " << llr[ind] << endl;
	if(llr[ind] == -200)
		return NULL;

	cluster.rewind();
	for (unsigned long i=0; i<=ind; i++){
		segment=cluster.getSeg();
	}
	if(segment->length() > 600)
		goon = false;
	else{
		llr[ind]=-200;
	}
}	
cluster.rewind();

return segment;
}

/**********************************************************
* bestFittingSegment: search the best segments of a cluster (in terms of likelihood) and return it
***********************************************************/
Seg *bestFittingSegment(Config& config, SegCluster& cluster, MixtureGD& m, StatServer& ss, FeatureServer& fs, MixtureGD& world){

DoubleVector llr,starts;                                        // The LLR value and the startindex of each of the potential trial
cluster.rewind();
Seg *segment;	
while((segment=cluster.getSeg())!=NULL){
		String fileLabel=segment->sourceName();
		unsigned long begin=segment->begin()+fs.getFirstFeatureIndexOfASource(fileLabel);//?????????????
		
		// Which normalization ?
		double normalizedFactor=meanLikelihood(ss,fs,world,begin,segment->length(),config);

		double llrTmp = meanLikelihood(ss,fs,m,begin,segment->length(),config) - normalizedFactor;
		if(verbose) cout << segment->begin() << " " << endSeg(segment) << " => " << llrTmp << endl;
		llr.addValue(llrTmp);
}//while segment

if(llr.size() == 0) return NULL;

bool goon=true;
while(goon){
	unsigned long ind = llr.getIndexOfLargestValue();
	if((verbose) && (verboseLevel == 2)) cout << "Best llk " << ind << " => " << llr[ind] << endl;
	if(llr[ind] == -200)
		return NULL;

	cluster.rewind();
	for (unsigned long i=0; i<=ind; i++){
		segment=cluster.getSeg();
	}
	if(segment->length() > 600)
		goon = false;
	else{
		llr[ind]=-200;
	}
}	
cluster.rewind();

return segment;
}

/**********************************************************
* bestFittingSegment: search the best segments of a cluster (in terms of likelihood) and return it
***********************************************************/
Seg *bestFittingSegment(Config& config, SegCluster& cluster, MixtureGD& m, StatServer& ss, FeatureServer& fs, hmm& actualHMM,unsigned long except){

DoubleVector llr,starts;                                        // The LLR value and the startindex of each of the potential trial
cluster.rewind();
Seg *segment;	
while((segment=cluster.getSeg())!=NULL){
		String fileLabel=segment->sourceName();
		unsigned long begin=segment->begin()+fs.getFirstFeatureIndexOfASource(fileLabel);//?????????????
		
		// Which normalization ?
		double normalizedFactor=cohortMaxLikelihood(ss,fs,actualHMM,except,segment->begin(),segment->length(),config);
			
		double llrTmp = meanLikelihood(ss,fs,m,begin,segment->length(),config) - normalizedFactor;
		if(verbose) cout << segment->begin() << " " << endSeg(segment) << " => " << llrTmp << endl;
		llr.addValue(llrTmp);
}//while segment

if(llr.size() == 0) return NULL;

bool goon=true;
while(goon){
	unsigned long ind = llr.getIndexOfLargestValue();
	if((verbose) && (verboseLevel == 2)) cout << "Best llk " << ind << " => " << llr[ind] << endl;
	if(llr[ind] == -200)
		return NULL;

	cluster.rewind();
	for (unsigned long i=0; i<=ind; i++){
		segment=cluster.getSeg();
	}
	if(segment->length() > 600)
		goon = false;
	else{
		llr[ind]=-200;
	}
}	
cluster.rewind();

return segment;
}


/**********************************************************
* bestFittingCLuster: search the best cluster for a segment (in terms of normalized likelihood) and return it
***********************************************************/
unsigned long bestFittingCluster(Config& config, hmm& actualHMM, SegServer& actualSeg, Seg *segment, StatServer& ss, FeatureServer& fs, unsigned long exceptInd=200){

DoubleVector llr,allInd;                                        // The LLR value and the startindex of each of the potential trial

for(unsigned long icluster=0; icluster<actualSeg.getClusterCount(); icluster++){
	if(icluster != exceptInd){
		MixtureGD& m=(MixtureGD&)actualHMM.getDensity(icluster);
		double llrTmp = meanLikelihood(ss,fs,m,segment->begin(),segment->length(),config);
		if(verbose) cout << "meanLikelihood of cluster: " << icluster << " => " << llrTmp << endl;
		llr.addValue(llrTmp);
		allInd.addValue((double)icluster);
	}
}
unsigned long ind = (unsigned long)(allInd[llr.getIndexOfLargestValue()]);
if(verbose) cout << "Best cluster: " << ind << endl;

return ind;
}

/**********************************************************
* intraCluster: evaluate purity intra cluster
***********************************************************/
void intraCluster(Config& config, hmm& actualHMM, SegServer& actualSeg, StatServer& ss, FeatureServer& fs, MixtureGD& world, String crit){
	for(unsigned long icluster=0;icluster<actualSeg.getClusterCount();icluster++){ // For each cluster
		SegCluster& clusterA=actualSeg.getCluster(icluster);         // Get the current cluster of segments 
		if (verbose) cout<<" Purify cluster for speaker "<<clusterA.string()<<endl;
		clusterA.rewind();
		Seg *segment1 = bestFittingSegment(config, clusterA, (MixtureGD&)actualHMM.getDensity(icluster), ss, fs,world);
		if ((verbose) && (verboseLevel == 2)) cout<<"Best segment "<<segment1->begin()<< " " << endSeg(segment1) << endl;
		Seg *segment2;	
		while((segment2=clusterA.getSeg())!=NULL)                               // Loop on each segment comming from the macro-class acoustic segmentation
			isSimilarSegment(config, segment1, segment2, ss, fs, world, crit);
	}//for icluster
}


/**********************************************************
* interCluster: Evaluate purity inter cluster
***********************************************************/
void interCluster(Config& config, hmm& actualHMM, SegServer& actualSeg, StatServer& ss, FeatureServer& fs, MixtureGD& world, String crit){

	for(unsigned long icluster=0;icluster<actualSeg.getClusterCount();icluster++){ // For each cluster
	    
		SegCluster& clusterA=actualSeg.getCluster(icluster);         // Get the current cluster of segments 
		if (verbose) cout<<" Evaluation of purity inter cluster for speaker "<<clusterA.string()<<endl;
		clusterA.rewind();

		Seg *segment1 = bestFittingSegment(config, clusterA, (MixtureGD&)actualHMM.getDensity(icluster), ss, fs,world);
		if ((verbose) && (verboseLevel == 2)) cout<<"Best segment "<<segment1->begin()<< " " << endSeg(segment1) << endl;
		Seg *segment2;	
		
		for(unsigned long other=0; other<actualSeg.getClusterCount(); other++){
			if(other != icluster){
				SegCluster& clusterB=actualSeg.getCluster(other);         // Get the current cluster of segments 
				clusterB.rewind();
				if (verbose) cout<<" Against speaker "<<clusterB.string()<<endl;
				
				while((segment2=clusterB.getSeg())!=NULL)                               // Loop on each segment comming from the macro-class acoustic segmentation
					isSimilarSegment(config, segment1, segment2, ss, fs, world, crit);
			}//if other
		}
	}//for icluster
}


/**********************************************************
* InitOneGausMixture: Create a mixture with only one Gaussian and initialized to the global covariance and mean
***********************************************************/
MixtureGD & InitOneGausMixture(Config& config, SegCluster& c, FeatureServer& fs, MixtureServer& ms){

unsigned long vectSize=fs.getVectSize();
MixtureGD &m=ms.createMixtureGD(1);
DistribGD &d=m.getDistrib(0);
	
FrameAccGD globalFrameAcc;                                                    //compute global Mean and cov
globalMeanCov (fs,c,globalFrameAcc,config);                             // Compute the global mean and covariance
DoubleVector globalMean=globalFrameAcc.getMeanVect(); 
DoubleVector globalCov=globalFrameAcc.getCovVect();
	
for(unsigned long i=0; i<vectSize; i++){
	d.setMean(globalMean[i], i);
	d.setCov(globalCov[i], i);
}
d.computeAll();	
	
return m;
}


/**********************************************************
* clusteringCriterionWithoutWorldInitOneGaus: Application of different criterion for clustering (models are trained by EM/ML - Initialization by external model)
***********************************************************/
double clusteringCriterionWithoutWorldInitOneGaus(Config& config, SegCluster& c1, SegCluster& c2, StatServer& ss, 
		FeatureServer& fs,String crit){
	
double criterion=0.0;
SegServer segTemp;
MixtureServer ms(config); 

TrainCfg trainCfg(config);

unsigned long vectSize=fs.getVectSize();

SegCluster& c12=mergeCluster(c1, c2, segTemp);


if(crit == "DELTABIC"){

	MixtureGD& m1=InitOneGausMixture(config, c1, fs, ms);
	MixtureGD& m2=InitOneGausMixture(config, c2, fs, ms);
	MixtureGD& m12=InitOneGausMixture(config, c12, fs, ms);

	
	/* Split of the unique Gaussian + 2 iterations of EM */
        DistribGD &d12=m12.getDistrib(0);
	DistribGD &d12Bis=ms.duplicateDistrib(d12);

	for(unsigned long i=0; i<vectSize; i++){
		/* first gaussian mean is fixed to mean(gaussian to split)+squared(cov(gaussian to split))
		second gaussian mean is fixed to mean(gaussian to split)-squared(cov(gaussian to split)) */
		double mean=d12.getMean(i);
		d12.setMean(mean+sqrt(d12.getCov(i)), i);
		d12Bis.setMean(mean-sqrt(d12.getCov(i)), i);
	}
	
	m12.weight(0)=0.5;
	ms.addDistribToMixture(m12, d12Bis, 0.5);
	
	trainCfg.setNbTrainIt(2);
	for (unsigned long trainIt=0; trainIt<trainCfg.getNbTrainIt(); trainIt++){                           // Begin the initial EM iteration loop 
		MixtureStat& emAcc=ss.createAndStoreMixtureStat(m12);             // Create a statistic accumulator using the init model
		// Set the variance control parameters

		emAcc.resetEM();                                                                                    // EM stuff
		double llkPreviousIt=accumulateStatEM(ss,fs,emAcc,c12,config);                        // Compute EM statistics
		m12 = emAcc.getEM();                                                                             // Get the EM estimate
		//varianceControl(world,varianceFlooring,varianceCeiling,globalCov12);
		ss.deleteMixtureStat(emAcc);
	}

	criterion = deltabicCrit(config, c1, c2, c12, ss, fs, m1, m2, m12);
	ms.deleteMixture(m1);
	ms.deleteMixture(m2);
	ms.deleteMixture(m12);
	ms.deleteUnusedDistribs();

	return criterion;
}
else{
	MixtureGD& m1=InitOneGausMixture(config, c1, fs, ms);
	MixtureGD& m2=InitOneGausMixture(config, c2, fs, ms);

	if((crit=="GLR") || (crit=="BIC") || (crit=="DGLR")){
		MixtureGD& m12=InitOneGausMixture(config, c12, fs, ms);

		if(crit=="GLR"){
			criterion = gllrCrit(config, c1, c2, c12, ss, fs, m1, m2, m12);
			ms.deleteMixture(m1);
			ms.deleteMixture(m2);
			ms.deleteMixture(m12);
			ms.deleteUnusedDistribs();
			return criterion;
		}
		if(crit=="DGLR"){
			criterion = gllrCrit(config, c1, c2, c12, ss, fs, m1, m2, m12);
			ms.deleteMixture(m1);
			ms.deleteMixture(m2);
			ms.deleteMixture(m12);
			ms.deleteUnusedDistribs();
			return (-1.0*criterion);
		}
		if(crit=="BIC"){
			criterion = bicCrit(config, c1, c2, c12, ss, fs, m1, m2, m12);
			ms.deleteMixture(m1);
			ms.deleteMixture(m2);
			ms.deleteMixture(m12);
			ms.deleteUnusedDistribs();
			return criterion;
		}		
	}

	else{
		cout << "ERROR: unknown clustering criterion !!!!" << endl;
		return -1;
	}
}

return -1;
}

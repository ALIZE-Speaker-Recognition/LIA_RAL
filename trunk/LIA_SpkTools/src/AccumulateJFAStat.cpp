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

#if !defined(ALIZE_JFAAcc_cpp)
#define ALIZE_JFAAcc_cpp

#if defined(_WIN32)
  #include <cfloat> // for _isnan()
  #define ISNAN(x) _isnan(x)
  #define ISINF(x) (!_finite(x))
#elif defined(linux) || defined(__linux) || defined(__CYGWIN__) || defined(__APPLE__)
  #define ISNAN(x) isnan(x)
  #define ISINF(x) isinf(x)
#else
  #error "Unsupported OS\n"
#endif

#include<AccumulateJFAStat.h>
#include<iostream>
#include<fstream>
#include<cstdio>
#include<cassert>
#include<cmath>
#include "SuperVectors.h"
#include "RealVector.h"
#include <liatools.h>
#include <limits>
#ifdef THREAD
#include <pthread.h>
#endif


using namespace alize;
using namespace std;


//***************************************************************//
//			Constructeurs et Destructeurs				 //
//***************************************************************//


//-----------------------------------------------------------------------------------------------------------------------------------------------------------
JFAAcc::JFAAcc(String & featFilename,Config & config)
	:_ms(config),_ss(config){ // constructor for a single file
		XList jfaNdx;

		if(featFilename.endsWith(".ndx")){
			jfaNdx.load(featFilename,config);
		}
		else{

			if(verboseLevel >1)cout<<"Init JFAAcc with only one file"<<endl;
			XLine& tmpLine = XLine::create();
			tmpLine.addElement(featFilename);

			//XList jfaNdx;
			jfaNdx.addLine()=tmpLine;
		}
		_init(jfaNdx,config);
}





//-----------------------------------------------------------------------------------------------------------------------------------------------------------
JFAAcc::JFAAcc(String & featFilename,Config & config, String task)
	:_ms(config),_ss(config){ // constructor for a single file
		XList jfaNdx;

		if(featFilename.endsWith(".ndx")){
			jfaNdx.load(featFilename,config);
		}
		else{

			if(verboseLevel >1)cout<<"Init JFAAcc with only one file"<<endl;
			XLine& tmpLine = XLine::create();
			tmpLine.addElement(featFilename);

			//XList jfaNdx;
			jfaNdx.addLine()=tmpLine;
		}
		_init(jfaNdx,config,task);
}


//-----------------------------------------------------------------------------------------------------------------------------------------------------------
JFAAcc::JFAAcc(XList & ndx,Config & config)
	:_ms(config),_ss(config){ // constructor
	_init(ndx,config);
}


//-----------------------------------------------------------------------------------------------------------------------------------------------------------
JFAAcc::JFAAcc(XList & ndx,Config & config,String task)
	:_ms(config),_ss(config){ // constructor
	_init(ndx,config,task);
}


//-----------------------------------------------------------------------------------------------------------------------------------------------------------
JFAAcc::~JFAAcc(){
	_vEvT.deleteAllObjects();
	_uEuT.deleteAllObjects();
	_vuEvuT.deleteAllObjects();
	_l_spk_inv.deleteAllObjects();
	_l_sess_inv.deleteAllObjects();
	_l_yx_inv.deleteAllObjects();
	_Aev.deleteAllObjects();
	_Aec.deleteAllObjects();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
String JFAAcc::getClassName() const{	return "JFAAcc";}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::_init(XList &ndx, Config &config){

	///Convert the NDX file
	_fileList=ndx;
	_ndxTable=JFATranslate(ndx);

	///Load the UBM
	MixtureGD& UBM = _ms.loadMixtureGD(config.getParam("inputWorldFilename"));

	_vectSize = UBM.getVectSize();
	_n_distrib = UBM.getDistribCount();
	_svSize = _vectSize*_n_distrib;

	_rankEV=1;
	if(config.existsParam("eigenVoiceNumber"))
		_rankEV=config.getParam("eigenVoiceNumber").toULong();
	_rankEC=1;
	if(config.existsParam("eigenChannelNumber"))
		_rankEC=config.getParam("eigenChannelNumber").toULong();
	
	///Read NDX file
	_n_speakers=_fileList.getLineCount();
	_n_sessions=_fileList.getAllElements().getElementCount();
	
	_ubm_means.setSize(_svSize);
	_ubm_invvar.setSize(_svSize);

	///Create UBM supervectors
	for(unsigned long i=0;i<_n_distrib;i++){
		DistribGD & dis=UBM.getDistrib(i);
		DoubleVector & c=dis.getCovInvVect();
		DoubleVector & m=dis.getMeanVect();
		for(unsigned long j=0;j<_vectSize;j++){
			_ubm_means[i*_vectSize+j]=m[j];
			_ubm_invvar[i*_vectSize+j]=c[j];
		}
	}

	///Create and initialise statistics acumulators
	_matN= Matrix<double>(_n_speakers, _n_distrib);
	_N_h= Matrix<double>(_n_sessions, _n_distrib);
	_F_X= Matrix<double>(_n_speakers, _n_distrib*_vectSize);
	_F_X_h= Matrix<double>(_n_sessions, _n_distrib*_vectSize);
	_matN.setAllValues(0.0);
	_N_h.setAllValues(0.0);
	_F_X.setAllValues(0.0);
	_F_X_h.setAllValues(0.0);
	
	///Create and initialise statistics acumulators
	_cN= Matrix<double>(_n_speakers, _n_distrib);
	_cN_h= Matrix<double>(_n_sessions, _n_distrib);
	_cF_X= Matrix<double>(_n_speakers, _n_distrib*_vectSize);
	_cF_X_h= Matrix<double>(_n_sessions, _n_distrib*_vectSize);
	_cN.setAllValues(0.0);
	_cN_h.setAllValues(0.0);
	_cF_X.setAllValues(0.0);
	_cF_X_h.setAllValues(0.0);
	
	///Initialise Matrices
	_D.setSize(_svSize);
	_D.setAllValues(0.0);
	
	_V.setDimensions(_rankEV, _svSize);
	_V.setAllValues(0.0);
	
	_matU.setDimensions(_rankEC, _svSize);
	_matU.setAllValues(0.0);
	
	_VU.setDimensions(_rankEV + _rankEC, _svSize);
	_VU.setAllValues(0.0);
	
	_Z.setDimensions(_n_speakers, _svSize);
	_Z.setAllValues(0.0);
	
	_Y.setDimensions(_n_speakers, _rankEV);
	_Y.setAllValues(0.0);
	
	_matX.setDimensions(_n_sessions, _rankEC);
	_matX.setAllValues(0.0);
	
	_YX.setDimensions(_n_speakers, _rankEV + _rankEC);
	_YX.setAllValues(0.0);
	
	///Create the accumulators for L matrices computation
	for(unsigned long s =0; s<_n_distrib; s++){
		_vEvT.addObject(*new DoubleSquareMatrix(_rankEV));
		_vEvT[s].setAllValues(0.0);
	}
	for(unsigned long s =0; s<_n_distrib; s++){
		_uEuT.addObject(*new DoubleSquareMatrix(_rankEC));
		_uEuT[s].setAllValues(0.0);
	}

	for(unsigned long s =0; s<_n_distrib; s++){
		_vuEvuT.addObject(*new DoubleSquareMatrix(_rankEV + _rankEC));
		_vuEvuT[s].setAllValues(0.0);
	}

	///Create accumulators to compute and inverse L matrices for Eigenvoices
	for(unsigned long spk =0; spk<_n_speakers; spk++){
		_l_spk_inv.addObject(*new DoubleSquareMatrix(_rankEV));
		_l_spk_inv[spk].setAllValues(0.0);
	}
	
	///Create accumulators to compute and inverse L matrices for EigenChannel
	for(unsigned long sess =0; sess<_n_sessions; sess++){
		_l_sess_inv.addObject(*new DoubleSquareMatrix(_rankEC));
		_l_sess_inv[sess].setAllValues(0.0);
	}

	///Create accumulators to compute and inverse L matrices for training phase
	for(unsigned long spk =0; spk<_n_speakers; spk++){
		_l_yx_inv.addObject(*new DoubleSquareMatrix(_rankEV + _rankEC));
		_l_yx_inv[spk].setAllValues(0.0);
	}
	
	///Create accumulators to compute the EigenVoices matrix
	for(unsigned long spk =0; spk<_n_distrib; spk++){
		_Aev.addObject(*new DoubleSquareMatrix(_rankEV));
		_Aev[spk].setAllValues(0.0);
	}
	_Cev.setDimensions(_rankEV,_svSize);
	_Cev.setAllValues(0.0);

	///Create accumulators to compute the EigenChannel matrix
	for(unsigned long spk =0; spk<_n_distrib; spk++){
		_Aec.addObject(*new DoubleSquareMatrix(_rankEC));
		_Aec[spk].setAllValues(0.0);
	}
	_Cec.setDimensions(_rankEC,_svSize);
	_Cec.setAllValues(0.0);
	
}


//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::_init(XList &ndx, Config &config, String task){

	//task can take values:
	// Accumulate
	// EigenVoice
	// EigenChannelLFA
	// EigenChannelJFA
	// EstimateD
	// TrainTarget
	// ComputeTest

	///Convert the NDX file
	_fileList=ndx;
	_ndxTable=JFATranslate(ndx); 
	///Load the UBM
	MixtureGD& UBM = _ms.loadMixtureGD(config.getParam("inputWorldFilename"));

	_vectSize = UBM.getVectSize();
	_n_distrib = UBM.getDistribCount();
	_svSize = _vectSize*_n_distrib;

	_rankEV=1;
	if(config.existsParam("eigenVoiceNumber"))
		_rankEV=config.getParam("eigenVoiceNumber").toULong();
	_rankEC=1;
	if(config.existsParam("eigenChannelNumber"))
		_rankEC=config.getParam("eigenChannelNumber").toULong();
	
	///Read NDX file
	_n_speakers=_fileList.getLineCount();
	_n_sessions=_fileList.getAllElements().getElementCount();
	
	_ubm_means.setSize(_svSize);
	_ubm_invvar.setSize(_svSize);

	///Create UBM supervectors
	for(unsigned long i=0;i<_n_distrib;i++){
		DistribGD & dis=UBM.getDistrib(i);
		DoubleVector & c=dis.getCovInvVect();
		DoubleVector & m=dis.getMeanVect();
		for(unsigned long j=0;j<_vectSize;j++){
			_ubm_means[i*_vectSize+j]=m[j];
			_ubm_invvar[i*_vectSize+j]=c[j];
		}
	}

	///Create and initialise statistics acumulators
	_matN= Matrix<double>(_n_speakers, _n_distrib);
	_N_h= Matrix<double>(_n_sessions, _n_distrib);
	_F_X= Matrix<double>(_n_speakers, _n_distrib*_vectSize);
	_F_X_h= Matrix<double>(_n_sessions, _n_distrib*_vectSize);
	_matN.setAllValues(0.0);
	_N_h.setAllValues(0.0);
	_F_X.setAllValues(0.0);
	_F_X_h.setAllValues(0.0);
	
	///Create and initialise statistics acumulators
	_cN= Matrix<double>(_n_speakers, _n_distrib);
	_cN_h= Matrix<double>(_n_sessions, _n_distrib);
	_cF_X= Matrix<double>(_n_speakers, _n_distrib*_vectSize);
	_cF_X_h= Matrix<double>(_n_sessions, _n_distrib*_vectSize);
	_cN.setAllValues(0.0);
	_cN_h.setAllValues(0.0);
	_cF_X.setAllValues(0.0);
	_cF_X_h.setAllValues(0.0);
	
	///Initialise Matrices
	_D.setSize(_svSize);
	_D.setAllValues(0.0);
	
	_V.setDimensions(_rankEV, _svSize);
	_V.setAllValues(0.0);
	
	_matU.setDimensions(_rankEC, _svSize);
	_matU.setAllValues(0.0);

	if((task == "TrainTarget")||(task == "ComputeTest")){
		_VU.setDimensions(_rankEV + _rankEC, _svSize);
		_VU.setAllValues(0.0);
	}
	else{
		_VU.setDimensions(1,1);
		_VU.setAllValues(0.0);
	}

	_Z.setDimensions(_n_speakers, _svSize);
	_Z.setAllValues(0.0);

	_Y.setDimensions(_n_speakers, _rankEV);
	_Y.setAllValues(0.0);
	
	_matX.setDimensions(_n_sessions, _rankEC);
	_matX.setAllValues(0.0);


	if((task == "TrainTarget")||(task == "ComputeTest")){
		_YX.setDimensions(_n_speakers, _rankEV + _rankEC);
		_YX.setAllValues(0.0);
	}
	else{
		_YX.setDimensions(1, _rankEV + _rankEC);
		_YX.setAllValues(0.0);
	}
	
	///Create the accumulators for L matrices computation
	for(unsigned long s =0; s<_n_distrib; s++){
		_vEvT.addObject(*new DoubleSquareMatrix(_rankEV));
		_vEvT[s].setAllValues(0.0);
	}

	if((task == "EigenChannelJFA")||(task == "EigenChannelLFA")||(task == "EstimateD")||(task == "TrainTarget")||(task == "ComputeTest")){
		for(unsigned long s =0; s<_n_distrib; s++){
			_uEuT.addObject(*new DoubleSquareMatrix(_rankEC));
			_uEuT[s].setAllValues(0.0);
		}
	}
	else{
		for(unsigned long s =0; s<1; s++){
			_uEuT.addObject(*new DoubleSquareMatrix(1));
			_uEuT[s].setAllValues(0.0);
		}
	}

	if((task == "TrainTarget")||(task == "ComputeTest")){
		for(unsigned long s =0; s<_n_distrib; s++){
			_vuEvuT.addObject(*new DoubleSquareMatrix(_rankEV + _rankEC));
			_vuEvuT[s].setAllValues(0.0);
		}
	}
	else{
		for(unsigned long s =0; s<1; s++){
			_vuEvuT.addObject(*new DoubleSquareMatrix(1));
			_vuEvuT[s].setAllValues(0.0);
		}
	}

	///Create accumulators to compute and inverse L matrices for Eigenvoices
	for(unsigned long spk =0; spk<_n_speakers; spk++){
		_l_spk_inv.addObject(*new DoubleSquareMatrix(_rankEV));
		_l_spk_inv[spk].setAllValues(0.0);
	}
	
	//Create accumulators to compute and inverse L matrices for EigenChannel
	if((task == "EigenChannelJFA")||(task == "EigenChannelLFA")||(task == "EstimateD")||(task == "TrainTarget")||(task == "ComputeTest")){
		for(unsigned long sess =0; sess<_n_sessions; sess++){
			_l_sess_inv.addObject(*new DoubleSquareMatrix(_rankEC));
			_l_sess_inv[sess].setAllValues(0.0);
		}
	}
	else{
		for(unsigned long sess =0; sess<1; sess++){
			_l_sess_inv.addObject(*new DoubleSquareMatrix(1));
			_l_sess_inv[sess].setAllValues(0.0);
		}
	}

	///Create accumulators to compute and inverse L matrices for training phase
	if((task == "TrainTarget")||(task == "ComputeTest")){
		for(unsigned long spk =0; spk<_n_speakers; spk++){
			_l_yx_inv.addObject(*new DoubleSquareMatrix(_rankEV + _rankEC));
			_l_yx_inv[spk].setAllValues(0.0);
		}
	}
	else{
		for(unsigned long spk =0; spk<1; spk++){
			_l_yx_inv.addObject(*new DoubleSquareMatrix(1));
			_l_yx_inv[spk].setAllValues(0.0);
		}
	}
	
	///Create accumulators to compute the EigenVoices matrix
	for(unsigned long spk =0; spk<_n_distrib; spk++){
		_Aev.addObject(*new DoubleSquareMatrix(_rankEV));
		_Aev[spk].setAllValues(0.0);
	}
	_Cev.setDimensions(_rankEV,_svSize);
	_Cev.setAllValues(0.0);


	///Create accumulators to compute the EigenChannel matrix
	if((task == "EigenChannelJFA")||(task == "EigenChannelLFA")||(task == "EstimateD")||(task == "TrainTarget")||(task == "ComputeTest")){
		for(unsigned long spk =0; spk<_n_distrib; spk++){
			_Aec.addObject(*new DoubleSquareMatrix(_rankEC));
			_Aec[spk].setAllValues(0.0);
		}
		_Cec.setDimensions(_rankEC,_svSize);
		_Cec.setAllValues(0.0);
	}
	else{
		for(unsigned long spk =0; spk<1; spk++){
			_Aec.addObject(*new DoubleSquareMatrix(1));
			_Aec[spk].setAllValues(0.0);
		}
		_Cec.setDimensions(1,1);
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::computeAndAccumulateJFAStat(Config& config){
	
	//STILL TO IMPLEMENT A THREADED VERSION OF THIS FUNCTION
	
//	#ifdef THREAD          
//	if (config.existsParam("numThread") && config.getParam("numThread").toULong() >0){
//		computeAndAccumulateJFAStatThreaded(config);				//accumulate stats
//	}
//	else	computeAndAccumulateJFAStatUnThreaded(config); 			//unthreaded version
//	#else
		computeAndAccumulateJFAStatUnThreaded(config);			//accumute stats
//	#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::computeAndAccumulateJFAStatUnThreaded(Config& config){

	MixtureGD& UBM = _ms.getMixtureGD(0);
	MixtureGDStat &acc=_ss.createAndStoreMixtureStat(UBM);

	///Create and initialise the feature server
	FeatureServer fs;
	XLine allFiles=_fileList.getAllElements();

	fs.init(config, allFiles);

	///Create and initialise feature clusters
	SegServer segmentsServer;
	LabelServer labelServer;
	initializeClusters(allFiles,segmentsServer,labelServer,config);

	verifyClusterFile(segmentsServer,fs,config);
	unsigned long codeSelectedFrame=labelServer.getLabelIndexByString(config.getParam("labelSelectedFrames"));

	SegCluster& selectedSegments=segmentsServer.getCluster(codeSelectedFrame);
	acc.resetOcc();
	Seg *seg; 
	selectedSegments.rewind();

	///Fast access to vector and matrices array
	double *n_h, *n, *f_x_h, *f_x;	
	n_h=_N_h.getArray(); n=_matN.getArray(); f_x_h=_F_X_h.getArray();f_x=_F_X.getArray();
	
	String currentSource="";unsigned long loc=0;unsigned long session=0;
	while((seg=selectedSegments.getSeg())!=NULL){

		unsigned long begin=seg->begin()+fs.getFirstFeatureIndexOfASource(seg->sourceName()); 				/// Idx of the first frame of the current file in the feature server

		JFATranslate ndxTable=JFATranslate(_fileList);
		
		if (currentSource!=seg->sourceName()) {
			currentSource=seg->sourceName();
			loc=ndxTable.locNb(currentSource);
			session=ndxTable.sessionNb(currentSource);
			if (verboseLevel >= 1)cout << "Processing file["<<currentSource<<"]"<< endl;	
		}

		fs.seekFeature(begin);
		Feature f;

		for (unsigned long idxFrame=0;idxFrame<seg->length();idxFrame++){
			fs.readFeature(f);
			acc.computeAndAccumulateOcc(f);
			DoubleVector aPost=acc.getOccVect();
			double *ff=f.getDataVector();

			for(unsigned long k=0;k<_n_distrib;k++) {

				n_h[session*_n_distrib+k]+=aPost[k];
				n[loc*_n_distrib+k]   +=aPost[k];

				for (unsigned long i=0;i<_vectSize;i++) {
					f_x_h[session*_svSize+(k*_vectSize+i)]+=aPost[k]*ff[i];
					f_x[loc*_svSize+(k*_vectSize+i)]   +=aPost[k]*ff[i];
				}
			}
		}
	}	
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::computeAndAccumulateJFAStat(FeatureServer &fs,Config & config){
	SegServer segmentsServer;
	LabelServer labelServer;
	initializeClusters(_fileList,segmentsServer,labelServer,config);
	verifyClusterFile(segmentsServer,fs,config);
	unsigned long codeSelectedFrame=labelServer.getLabelIndexByString(config.getParam("labelSelectedFrames"));
	SegCluster& selectedSegments=segmentsServer.getCluster(codeSelectedFrame); 
	this->computeAndAccumulateJFAStat(selectedSegments,fs,config);
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::computeAndAccumulateJFAStat(SegCluster &selectedSegments,FeatureServer &fs,Config & config){

	MixtureGD& UBM = _ms.getMixtureGD(0);
	MixtureGDStat &acc=_ss.createAndStoreMixtureStat(UBM);

	///Fast access to vector and matrices array
	double *n_h, *n, *f_x_h, *f_x;	
	n_h=_N_h.getArray(); n=_matN.getArray(); f_x_h=_F_X_h.getArray();f_x=_F_X.getArray();
	
	///Compute Occupations and Statistics	
	acc.resetOcc();
	Seg *seg; 
	selectedSegments.rewind();
	
	String currentSource="";unsigned long loc=0;unsigned long session=0;
	while((seg=selectedSegments.getSeg())!=NULL){

		unsigned long begin=seg->begin()+fs.getFirstFeatureIndexOfASource(seg->sourceName()); 				/// Idx of the first frame of the current file in the feature server

		JFATranslate ndxTable=JFATranslate(_fileList);
		
		if (currentSource!=seg->sourceName()) {
			currentSource=seg->sourceName();
			loc=ndxTable.locNb(currentSource);
			session=ndxTable.sessionNb(currentSource);
			if (verboseLevel >= 1)cout << "Processing file["<<currentSource<<"]"<< endl;	
		}

		fs.seekFeature(begin);
		Feature f;

		for (unsigned long idxFrame=0;idxFrame<seg->length();idxFrame++){

			fs.readFeature(f);
			acc.computeAndAccumulateOcc(f);
			DoubleVector aPost=acc.getOccVect();
			double *ff=f.getDataVector();

			for(unsigned long k=0;k<_n_distrib;k++) {

				n_h[session*_n_distrib+k]+=aPost[k];
				n[loc*_n_distrib+k]   +=aPost[k];

				for (unsigned long i=0;i<_vectSize;i++) {
					f_x_h[session*_svSize+(k*_vectSize+i)]+=aPost[k]*ff[i];
					f_x[loc*_svSize+(k*_vectSize+i)]   +=aPost[k]*ff[i];
				}
			}
		}
	}	
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::resetAcc(){
	_F_X.setAllValues(0.0);	
	_F_X_h.setAllValues(0.0);	
	_N_h.setAllValues(0.0);
	_matN.setAllValues(0.0);	
	if (verboseLevel >= 1) cout << "# JFA Accumulators reset" << endl;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::resetTmpAcc(){
	_Cev.setAllValues(0.0);
	_Cec.setAllValues(0.0);
	///Reinitialise accumulators for L matrices computation
	for(unsigned long s =0; s<_n_distrib; s++){
		_vEvT[s].setAllValues(0.0);
		_uEuT[s].setAllValues(0.0);

		_Aev[s].setAllValues(0.0);
		_Aec[s].setAllValues(0.0);
	}

	///Reinitialise accumulators for L matrices computation and inversion
	for(unsigned long spk =0; spk<_n_speakers; spk++){
		_l_spk_inv[spk].setAllValues(0.0);
	}

	for(unsigned long session =0; session<_n_speakers; session++){
		_l_sess_inv[session].setAllValues(0.0);
	}
}


//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::resetTmpAcc(String task){
	
	_Cev.setAllValues(0.0);
	_Cec.setAllValues(0.0);
	///Reinitialise accumulators for L matrices computation
	for(unsigned long s =0; s<_n_distrib; s++){
		_vEvT[s].setAllValues(0.0);
		_Aev[s].setAllValues(0.0);
	}

	///Reinitialise accumulators for L matrices computation
	unsigned long tmpDistribNb;
	if((task == "EigenChannelJFA")||(task == "EigenChannelLFA")||(task == "EstimateD")||(task == "TrainTarget")||(task == "ComputeTest")){
		tmpDistribNb = _n_distrib;
	}
	else{
		tmpDistribNb = 1;
	}

	for(unsigned long s =0; s<tmpDistribNb; s++){
		_uEuT[s].setAllValues(0.0);
		_Aec[s].setAllValues(0.0);
	}

	///Reinitialise accumulators for L matrices computation and inversion
	for(unsigned long spk =0; spk<_n_speakers; spk++){
		_l_spk_inv[spk].setAllValues(0.0);
	}

	unsigned long tmpSpkNb;
	if((task == "EigenChannelJFA")||(task == "EigenChannelLFA")||(task == "EstimateD")||(task == "TrainTarget")||(task == "ComputeTest")){
		tmpSpkNb = _n_speakers;
	}
	else{
		tmpSpkNb = 1;
	}
	for(unsigned long session =0; session<tmpSpkNb; session++){
		_l_sess_inv[session].setAllValues(0.0);
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//void JFAAcc::resetTmpAccIV(){
//	
//	_Cev.setAllValues(0.0);
//	_Cec.setAllValues(0.0);
//	///Reinitialise accumulators for L matrices computation
//	for(unsigned long s =0; s<_n_distrib; s++){
//		_vEvT[s].setAllValues(0.0);
//		//_uEuT[s].setAllValues(0.0);
//		
//		_Aev[s].setAllValues(0.0);
//		//_Aec[s].setAllValues(0.0);
//	}
//
//	///Reinitialise accumulators for L matrices computation and inversion
//	for(unsigned long spk =0; spk<_n_speakers; spk++){
//		_l_spk_inv[spk].setAllValues(0.0);
//	}
//
////	for(unsigned long session =0; session<_n_speakers; session++){
////		_l_sess_inv[session].setAllValues(0.0);
////	}
//}


//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::loadEV(const String& evFilename, Config& config){	///load an EigenVoices Matrix
	String filename = config.getParam("matrixFilesPath") + evFilename +  config.getParam("saveMatrixFilesExtension");
	_V.load (filename, config);

	//Transpose the matrix if the number of line is higher than the numbe of columns
	if(_V.cols()<_V.rows()){
		cout<<"Load EV : number of lines ( "<<_V.rows() <<" ) higher than the number of columns ( "<<_V.cols()<<" ("<<endl;
		_V.transpose();
	}
	_rankEV=_V.rows();

	if(_rankEV != config.getParam("eigenVoiceNumber").toULong()){
		throw Exception("Incorrect dimension of EigenVoice Matrix",__FILE__,__LINE__);
	}

	cout << "(AccumulateJFAStat) Init EV matrix from "<< filename <<"  for EigenVoices Matrix: "<<", rank: ["<<_V.rows() << "] sv size: [" << _V.cols() <<"]"<<endl;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::loadEV(Matrix<double> & V, Config& config){
	_rankEV=V.rows();
	_V=V;

	//Transpose the matrix if the number of line is higher than the numbe of columns
	if(_V.cols()<_V.rows()){
		cout<<"Load EV : number of lines higher than the number of columns"<<endl;
		_V.transpose();
	}

	unsigned long rankEV = 1;
	if(config.existsParam("eigenVoiceNumber"))
		rankEV = config.getParam("eigenVoiceNumber").toULong();

	if(_rankEV != config.getParam("eigenVoiceNumber").toULong()){
		throw Exception("Incorrect dimension of EigenVoice Matrix",__FILE__,__LINE__);
	}

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::loadEC(const String& ecFilename, Config& config){	///load an EigenChannel Matrix
	String filename = config.getParam("matrixFilesPath") + ecFilename +  config.getParam("saveMatrixFilesExtension");
	_matU.load (filename, config);
	_rankEC=_matU.rows();

	if(_rankEC != config.getParam("eigenChannelNumber").toULong()){
		throw Exception("Incorrect dimension of EigenChannel Matrix",__FILE__,__LINE__);
	}

	cout << "(AccumulateJFAStat) Init EC matrix from "<< filename <<"  for EigenChannel Matrix: "<<", rank: ["<<_matU.rows() << "] sv size: [" << _matU.cols() <<"]"<<endl;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::loadEC(Matrix<double> & U, Config& config){
	_rankEC=U.rows();
	_matU=U;

	unsigned long rankEC = 1;
	if(config.existsParam("eigenChannelNumber"))
		rankEC = config.getParam("eigenChannelNumber").toULong();

	if(_rankEC != rankEC){
		throw Exception("Incorrect dimension of EigenChannel Matrix",__FILE__,__LINE__);
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::loadD(const String& dFilename, Config& config){	///load a D Matrix
	String filename = config.getParam("matrixFilesPath") + dFilename +  config.getParam("saveMatrixFilesExtension");
	
	Matrix<double> D(filename, config);
	double *d = D.getArray();
	
	if( (D.rows() != 1) || ( D.cols() != _svSize ) ){
		throw Exception("Incorrect dimension of D Matrix",__FILE__,__LINE__);
	}
	
	else{
		for(unsigned long i=0; i<_svSize; i++){
			_D[i] = d[i];
		}
	}
	cout << "(AccumulateJFAStat) Init D from "<< filename <<"  for D Matrix: "<<endl;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::loadD(DoubleVector & D){
	_D=D;

	if(_D.size() != _svSize ){
		throw Exception("Incorrect dimension of D Matrix",__FILE__,__LINE__);
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::loadN(Config& config){
	String filname=config.getParam("matrixFilesPath") + config.getParam("nullOrderStatSpeaker") + config.getParam("loadMatrixFilesExtension");
	_matN.load (filname, config);

	if((_matN.rows() != _fileList.getLineCount()) || (_matN.cols() != _n_distrib)){
		throw Exception("Incorrect dimension of N Matrix",__FILE__,__LINE__);
	}
	cout<<"(AccumulateJFAStat) ---------load statistics N [ "<<_matN.rows()<<" ] [ "<<_matN.cols()<<" ]"<<endl;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::loadN_h(Config& config){
	String filname=config.getParam("matrixFilesPath") + config.getParam("nullOrderStatSession") + config.getParam("loadMatrixFilesExtension");
	_N_h.load (filname, config);


	if((_N_h.rows() != _fileList.getAllElements().getElementCount()) || (_N_h.cols() != _n_distrib)){
		throw Exception("Incorrect dimension of N_h Matrix",__FILE__,__LINE__);
	}
	cout<<"(AccumulateJFAStat) ---------load statistics N_h [ "<<_N_h.rows()<<" ] [ "<<_N_h.cols()<<" ]"<<endl;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::loadF_X(Config& config){
	String filname=config.getParam("matrixFilesPath") + config.getParam("firstOrderStatSpeaker") + config.getParam("loadMatrixFilesExtension");
	_F_X.load (filname, config);

	if((_F_X.rows() != _fileList.getLineCount()) || (_F_X.cols() != _svSize)){
		throw Exception("Incorrect dimension of F_X Matrix",__FILE__,__LINE__);
	}
	cout<<"(AccumulateJFAStat) ---------load statistics F_X [ "<<_F_X.rows()<<" ] [ "<<_F_X.cols()<<" ]"<<endl;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::loadF_X_h(Config& config){
	String filname=config.getParam("matrixFilesPath") + config.getParam("firstOrderStatSession") + config.getParam("loadMatrixFilesExtension");
	_F_X_h.load (filname, config);

	if((_F_X_h.rows() != _fileList.getAllElements().getElementCount()) || (_F_X_h.cols() != _svSize)){
		throw Exception("Incorrect dimension of F_X_h Matrix",__FILE__,__LINE__);
	}
	cout<<"(AccumulateJFAStat) ---------load statistics F_X_h [ "<<_F_X_h.rows()<<" ] [ "<<_F_X_h.cols()<<" ]"<<endl;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::initEV(Config& config){		///random initialisation of the EigenVoices Matrix
		_rankEV=config.getParam("eigenVoiceNumber").toULong();
		_V.setDimensions(_rankEV,_svSize);

		//Two type of random initializations
		//	normal  : random with a normal law by using a Box-Muller Generator
		//	uniform : random with a uniform distribution
try{
		String randomLaw = "normal";
		if(config.existsParam("randomInitLaw"))	randomLaw = config.getParam("randomInitLaw");

		//Initialize the matrix by generating random values following a uniform law
		if(config.getParam("randomInitLaw") == "uniform"){

			srand48(_svSize*_rankEV);
			_V.randomInit();

			double norm=0;
			for(unsigned long k=0; k<_svSize; k++){
				norm += _ubm_invvar[k];
			}
			norm = norm/_svSize;
			for(unsigned long i=0; i<_V.rows(); i++){
				for(unsigned long j=0; j<_V.cols(); j++){
					_V(i,j) = _V(i,j)*norm;
				}
			}
		}

		//Initialize the matrix by generating random values following a normal law whose mean is 0 and variance is 
		//th norm of the Covariance matrix from the UBM multiplied by 0.001
		else if (config.getParam("randomInitLaw") == "normal"){
			double norm=0;
			for(unsigned long k=0; k<_svSize; k++){
				norm += _ubm_invvar[k];
			}

			boxMullerGeneratorInit();
			for(unsigned long i=0; i<_V.rows(); i++){
				for(unsigned long j=0; j<_V.cols(); j++){
					double val = boxMullerGenerator(0.0, 1.0);
					while((ISNAN(val)) || (ISINF(val))){
						val = boxMullerGenerator(0.0, 1.0);
					}
					_V(i,j) = val * norm * 0.001;
				}
			}
		}
		else{
			throw Exception("Selected random initialization law does not exist",__FILE__,__LINE__);
		}
		if (verboseLevel >=1) cout << "(AccumulateJFAStat) Random Init for EigenVoices Matrix with "<<randomLaw<<" law "<<", rank: ["<<_V.rows() << "] sv size: [" << _V.cols() <<"]"<<endl;
}
catch (Exception& e){ 
	cout << e.toString().c_str() << endl;
}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::initEC(Config& config){		///random initialisation of the EigenChannel Matrix
		_rankEC=config.getParam("eigenChannelNumber").toULong();
		_matU.setDimensions(_rankEC,_svSize);

		//Two type of random initializations
		//	normal  : random with a normal law by using a Box-Muller Generator
		//	uniform : random with a uniform distribution
try{
	String randomLaw = "normal";
	if(config.existsParam("randomInitLaw"))	randomLaw = config.getParam("randomInitLaw");

	//Initialize the matrix by generating random values following a uniform law
	if(config.getParam("randomInitLaw") == "uniform"){
		srand48(_svSize*_rankEC);
		_matU.randomInit();
	}

	//Initialize the matrix by generating random values following a normal law whose mean is 0 and variance is 
	//th norm of the Covariance matrix from the UBM multiplied by 0.001
	else if (config.getParam("randomInitLaw") == "normal"){
		double norm=0;
		for(unsigned long k=0; k<_svSize; k++){
			norm += _ubm_invvar[k];
		}

		boxMullerGeneratorInit();
		for(unsigned long i=0; i<_matU.rows(); i++){
			for(unsigned long j=0; j<_matU.cols(); j++){
				double val = boxMullerGenerator(0.0, 1.0);
				while((ISNAN(val)) || (ISINF(val))){
					val = boxMullerGenerator(0.0, 1.0);
				}
				_matU(i,j) = val * norm * 0.001;
			}
		}
	}
	else{
		throw Exception("Slected random initialization law does not exist",__FILE__,__LINE__);
	}
	if (verboseLevel >=1) cout << "(AccumulateJFAStat) Random Init for EigenChannel Matrix with "<<randomLaw<<" law "<<", rank: ["<<_matU.rows() << "] sv size: [" << _matU.cols() <<"]"<<endl;
}
catch (Exception& e){ 
	cout << e.toString().c_str() << endl;
}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::initD(Config& config){			///initialisation of the D Matrix
		Matrix<double> D;
		D.setDimensions(1,_svSize);

try{

	//Choose the initialisation type for D
	String initDType = "MAP";
	if(config.existsParam("initDType"))	initDType = config.getParam("initDType");

	//Random initalisation using a uniform law
	if(initDType == "uniform"){
		srand48(_svSize);
		D.randomInit();
		for(unsigned long i=0; i<_svSize; i++){
			_D[i] = D(1,i);
		}
		cout<< "(AccumulateJFAStat) Init D Matrix with uniform random law :  sv size: [" << _D.size() <<"]"<<endl;
	}
	else if(initDType == "normal"){
		double norm=0;
		for(unsigned long k=0; k<_svSize; k++){
			norm += _ubm_invvar[k];
		}

		boxMullerGeneratorInit();
		for(unsigned long i=0; i<_svSize; i++){
			double val = boxMullerGenerator(0.0, 1.0);
			while((ISNAN(val)) || (ISINF(val))){
				val = boxMullerGenerator(0.0, 1.0);
			}
			_D[i] = val * norm * 0.001;
		}
		cout<< "(AccumulateJFAStat) Init D Matrix with Box-Muller random generator :  sv size: [" << _D.size() <<"]"<<endl;
	}
	else if(initDType == "MAP"){
		for(unsigned long i=0; i<_svSize; i++){
			_D[i] = sqrt(1.0/(_ubm_invvar[i]*config.getParam("regulationFactor").toDouble()));
		}
		cout<< "(AccumulateJFAStat) Init D Matrix with MAP :  sv size: [" << _D.size() <<"]"<<endl;
	}
	else{
		throw Exception("Selected D matrix initialization parameter does not exist",__FILE__,__LINE__);
	}
}
catch (Exception& e){ 
	cout << e.toString().c_str() << endl;
}	
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::initVU(){
	_VU.concatCols(_V,_matU);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::saveV(const String& filename, Config& config){
	String vName = config.getParam("matrixFilesPath") + filename +  config.getParam("saveMatrixFilesExtension");
	_V.save(vName, config);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::saveU(const String& filename, Config& config){
	String uName = config.getParam("matrixFilesPath") + filename +  config.getParam("saveMatrixFilesExtension");
	_matU.save(uName, config);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::saveD(const String& filename, Config& config){
	String dName = config.getParam("matrixFilesPath") + filename +  config.getParam("saveMatrixFilesExtension");
	Matrix<double> D;
	D.setDimensions(1,_svSize);
	for(unsigned long i=0; i<_svSize; i++){
		D(0,i) = _D[i];
	}
	D.save(dName, config);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::estimateVEVT(Config &config){
	
	#ifdef THREAD          
	if (config.existsParam("numThread") && config.getParam("numThread").toULong() >0)	estimateVEVTThreaded(config.getParam("numThread").toULong());
	else estimateVEVTUnThreaded();
	#else
	estimateVEVTUnThreaded();
	#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::estimateVEVTUnThreaded(){	
        if (verboseLevel >= 1) cout << "(AccumulateJFAStat) compute V * Sigma-1 * V "<<endl;

	for(unsigned long d=0; d<_n_distrib; d++){

		Matrix<double>ssV= _V.crop(0,d*_vectSize, _rankEV,_vectSize);

		///Compute  _vEvT matrices
		double *vEvT, *v, *E;
		_vEvT[d].setAllValues(0.0);		//ajout du 2011/07/06 pour etre comme la partie threadee
		vEvT= _vEvT[d].getArray();
		v= ssV.getArray();
		E= _ubm_invvar.getArray();

		for(unsigned long i = 0; i<_rankEV; i++){
			for(unsigned long j = 0; j<=i; j++){
				for(unsigned long k=0; k<_vectSize;k++){
					vEvT[i*_rankEV+j] += v[i*_vectSize+k] * E[d*_vectSize+k] * v[j*_vectSize+k];
				}
			}
		}
		
		for(unsigned long i=0;i<_rankEV;i++){
			for(unsigned long j=i+1;j<_rankEV;j++) {
				vEvT[i*_rankEV+j] = vEvT[j*_rankEV+i];
			}
		}
	}
}

#ifdef THREAD
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Data strucutre of thread
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
struct VEVthread_data{

	double *V;
	double *ubm_invvar;
	unsigned long disBottom;
	unsigned long disUp;	
	unsigned long rankEV;
	unsigned long svSize;
	unsigned long vectSize;
	RefVector <DoubleSquareMatrix>* vevT;
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Thread Routine
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void *VEVthread(void *threadarg) {
	struct VEVthread_data *my_data;
	my_data = (struct VEVthread_data *) threadarg;
	
	double *v = my_data->V;
	unsigned long disBottom = my_data->disBottom;	
	unsigned long disUp = my_data->disUp;
	double *E=my_data->ubm_invvar;
	unsigned long _rankEV=my_data->rankEV;
	unsigned long _svSize=my_data->svSize;
	unsigned long _vectSize=my_data->vectSize;

	for (unsigned long d=disBottom;d <disUp;d++){
		
		DoubleSquareMatrix &vEvT=(*(my_data->vevT))[d];
		vEvT.setAllValues(0.0);
		double *vevt = vEvT.getArray();
		
		for(unsigned long i = 0; i<_rankEV; i++){
			for(unsigned long j = 0; j<=i; j++){
				for(unsigned long k=0; k<_vectSize;k++){
					vevt[i*_rankEV+j] += v[(i*_svSize)+(d*_vectSize)+k] * E[d*_vectSize+k] * v[(j*_svSize)+(d*_vectSize)+k];
				}
			}
		}

		for(unsigned long i=0;i<_rankEV;i++){
			for(unsigned long j=i+1;j<_rankEV;j++) {
				vevt[i*_rankEV+j]=vevt[j*_rankEV+i];
			}
		}
	}

	pthread_exit((void*) 0);
	return (void*)0 ;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::estimateVEVTThreaded(unsigned long NUM_THREADS){
	
	if (verboseLevel >=1) cout << "(AccumulateJFAStat) Compute vEvT Threaded"<<endl;
	if (NUM_THREADS==0) throw Exception("Num threads can not be 0",__FILE__,__LINE__);
	
	double *V=_V.getArray();
	double *ubm_invvar=_ubm_invvar.getArray();
	
	int rc, status;
	if (NUM_THREADS > _n_distrib) NUM_THREADS=_n_distrib;
	
	struct VEVthread_data *thread_data_array = new VEVthread_data[NUM_THREADS];
	pthread_t *threads = new pthread_t[NUM_THREADS];

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	unsigned long offset=_n_distrib/NUM_THREADS;
	
	unsigned long disBottom = 0;
	unsigned long disUp = 0;
	unsigned long re=_n_distrib - NUM_THREADS*offset;
	
	//Create threads : one per distribution as a maximum
	for(unsigned long t=0; t<NUM_THREADS; t++){
	
		disUp = disBottom + offset;
		if(t<re) disUp +=1;

		thread_data_array[t].V = V;
		thread_data_array[t].ubm_invvar=ubm_invvar;
		thread_data_array[t].disBottom=disBottom;
		thread_data_array[t].disUp=disUp;
		thread_data_array[t].rankEV=_rankEV;
		thread_data_array[t].svSize=_svSize;
		thread_data_array[t].vectSize=_vectSize;
		thread_data_array[t].vevT=&(_vEvT);

		if (verboseLevel > 1) cout<<"(AccumulateJFAStat) Creating thread n["<< t<< "] for distributions["<<disBottom<<"-->"<<disUp<<"]"<<endl;
		rc = pthread_create(&threads[t], &attr, VEVthread, (void *)&thread_data_array[t]);
		if (rc) throw Exception("ERROR; return code from pthread_create() is ",__FILE__,rc);

		disBottom = disUp;
	}

	pthread_attr_destroy(&attr);

	for(unsigned long t=0; t<NUM_THREADS; t++) {
		rc = pthread_join(threads[t], (void **)&status);
		if (rc)  throw Exception("ERROR; return code from pthread_join() is ",__FILE__,rc);
		if (verboseLevel >1) cout <<"(AccumulateJFAStat) Completed join with thread ["<<t<<"] status["<<status<<"]"<<endl;
	}


	free(thread_data_array);
	free(threads);

	if (verboseLevel >= 1) cout << "(AccumulateJFAStat) Done " << endl;
}
#endif

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::estimateUEUT(Config &config){
	#ifdef THREAD          
	if (config.existsParam("numThread") && config.getParam("numThread").toULong() >0)	estimateUEUTThreaded(config.getParam("numThread").toULong());
	else estimateUEUTUnThreaded();
	#else
	estimateUEUTUnThreaded();
	#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::estimateUEUTUnThreaded(){	
        if (verboseLevel >= 1) cout << "(AccumulateJFAStat) compute U * Sigma-1 * U "<<endl;

	for(unsigned long d=0; d<_n_distrib; d++){

		Matrix<double>ssU= _matU.crop(0,d*_vectSize, _rankEC,_vectSize);

		///Compute _uEuT matrices
		double *uEuT, *u, *E;
		uEuT= _uEuT[d].getArray();
		u= ssU.getArray();
		E= _ubm_invvar.getArray();
		
		for(unsigned long i = 0; i<_rankEC; i++){
			for(unsigned long j = 0; j<=i; j++){
				for(unsigned long k=0; k<_vectSize;k++){
					uEuT[i*_rankEC+j] += u[i*_vectSize+k] * E[d*_vectSize+k] * u[j*_vectSize+k];
				}
			}
		}
		for(unsigned long i=0;i<_rankEC;i++){
			for(unsigned long j=i+1;j<_rankEC;j++) {
				uEuT[i*_rankEC+j] = uEuT[j*_rankEC+i];
			}
		}
	}
}

#ifdef THREAD
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Data strucutre of thread
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
struct UEUthread_data{

	double *U;
	double *ubm_invvar;
	unsigned long disBottom;
	unsigned long disUp;	
	unsigned long rankEC;
	unsigned long svSize;
	unsigned long vectSize;
	RefVector <DoubleSquareMatrix>* ueuT;
	unsigned long nT;
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Thread Routine
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void *UEUthread(void *threadarg) {
	struct UEUthread_data *my_data;
	my_data = (struct UEUthread_data *) threadarg;
	
	double *u = my_data->U;
	unsigned long disBottom = my_data->disBottom;	
	unsigned long disUp = my_data->disUp;
	double *E=my_data->ubm_invvar;
	unsigned long _rankEC=my_data->rankEC;
	unsigned long _svSize=my_data->svSize;
	unsigned long _vectSize=my_data->vectSize;

	for (unsigned long d=disBottom;d <disUp;d++){
		
		DoubleSquareMatrix &uEuT=(*(my_data->ueuT))[d];
		uEuT.setAllValues(0.0);
		double *ueut = uEuT.getArray();
		
		for(unsigned long i = 0; i<_rankEC; i++){
			for(unsigned long j = 0; j<=i; j++){
				for(unsigned long k=0; k<_vectSize;k++){
					ueut[i*_rankEC+j] += u[(i*_svSize)+(d*_vectSize)+k] * E[d*_vectSize+k] * u[(j*_svSize)+(d*_vectSize)+k];
				}
			}
		}
		for(unsigned long i=0;i<_rankEC;i++){
			for(unsigned long j=i+1;j<_rankEC;j++) {
				ueut[i*_rankEC+j] = ueut[j*_rankEC+i];
			}
		}
	}
	pthread_exit((void*) 0);
	return (void*)0 ;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::estimateUEUTThreaded(unsigned long NUM_THREADS){
	
	if (verboseLevel >=1) cout << "(AccumulateJFAStat) Compute uEuT Threaded"<<endl;
	if (NUM_THREADS==0) throw Exception("Num threads can not be 0",__FILE__,__LINE__);
	
	double *U=_matU.getArray();
	double *ubm_invvar=_ubm_invvar.getArray();
	
	int rc, status;
	if (NUM_THREADS > _n_distrib) NUM_THREADS=_n_distrib;
	
	//struct UEUthread_data thread_data_array[NUM_THREADS];
	struct UEUthread_data *thread_data_array = new UEUthread_data[NUM_THREADS];

	//pthread_t threads[NUM_THREADS];
	 pthread_t *threads = new pthread_t[NUM_THREADS];

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	unsigned long offset=_n_distrib/NUM_THREADS;
	
	unsigned long disBottom = 0;
	unsigned long disUp=0;
	unsigned long re=_n_distrib - NUM_THREADS*offset;

	//Create threads
	for(unsigned long t=0; t<NUM_THREADS; t++){
	
		disUp = disBottom +offset;
		if(t<re) disUp +=1;
		
		thread_data_array[t].U=U;
		thread_data_array[t].ubm_invvar=ubm_invvar;
		thread_data_array[t].disBottom=disBottom;
		thread_data_array[t].disUp=disUp;
		thread_data_array[t].rankEC=_rankEC;
		thread_data_array[t].svSize=_svSize;
		thread_data_array[t].vectSize=_vectSize;
		thread_data_array[t].ueuT=&(_uEuT);
		thread_data_array[t].nT=t;

		if (verboseLevel >1) cout<<"(AccumulateJFAStat) Creating thread n["<< t<< "] for distributions["<<disBottom<<"-->"<<disUp<<"]"<<endl;
		rc = pthread_create(&threads[t], &attr, UEUthread, (void *)&thread_data_array[t]);
		if (rc) throw Exception("ERROR; return code from pthread_create() is ",__FILE__,rc);
		
		disBottom = disUp;
	}
	
	pthread_attr_destroy(&attr);
	for(unsigned long t=0; t<NUM_THREADS; t++) {
		rc = pthread_join(threads[t], (void **)&status);
		if (rc)  throw Exception("ERROR; return code from pthread_join() is ",__FILE__,rc);
		if (verboseLevel >1) cout <<"(AccumulateJFAStat) Completed join with thread ["<<t<<"] status["<<status<<"]"<<endl;
	}

	free(thread_data_array);
	free(threads);

	if (verboseLevel >= 1) cout << "(AccumulateJFAStat) Done " << endl;
}
#endif

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::estimateVUEVUT(Config &config){
	
	#ifdef THREAD          
	if (config.existsParam("numThread") && config.getParam("numThread").toULong() >0)	estimateVUEVUTThreaded(config.getParam("numThread").toULong());
	else estimateVUEVUTUnThreaded();
	#else
	estimateVUEVUTUnThreaded();
	#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::estimateVUEVUTUnThreaded(){	
        if (verboseLevel >= 1) cout << "(AccumulateJFAStat) compute VU * Sigma-1 * VU "<<endl;

	for(unsigned long d=0; d<_n_distrib; d++){
		Matrix<double>ssVU= _VU.crop(0,d*_vectSize, _rankEV + _rankEC,_vectSize);

		///Compute _uEuT matrices
		double *vuEvuT, *vu, *E;
		vuEvuT= _vuEvuT[d].getArray();
		vu= ssVU.getArray();
		E= _ubm_invvar.getArray();

		for(unsigned long i = 0; i<_rankEV + _rankEC; i++){
			for(unsigned long j = 0; j<=i; j++){
				for(unsigned long k=0; k<_vectSize;k++)
					vuEvuT[i*(_rankEV + _rankEC)+j] += vu[i*_vectSize+k] * E[d*_vectSize+k] * vu[j*_vectSize+k];
			}
		}
		for(unsigned long i=0;i<_rankEV + _rankEC;i++){
			for(unsigned long j=i+1;j<_rankEV + _rankEC;j++) {
				vuEvuT[i*(_rankEV + _rankEC)+j] = vuEvuT[j*(_rankEV + _rankEC)+i];
			}
		}
	}
}

#ifdef THREAD
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Data strucutre of thread
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
struct VUEVUthread_data{

	double *VU;
	double *ubm_invvar;
	unsigned long disBottom;
	unsigned long disUp;	
	unsigned long rankECV;
	unsigned long svSize;
	unsigned long vectSize;
	RefVector <DoubleSquareMatrix>* vuevuT;
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Thread Routine
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void *VUEVUthread(void *threadarg) {
	struct VUEVUthread_data *my_data;
	my_data = (struct VUEVUthread_data *) threadarg;
	
	double *vu = my_data->VU;
	unsigned long disBottom = my_data->disBottom;	
	unsigned long disUp = my_data->disUp;
	double *E=my_data->ubm_invvar;
	unsigned long _rankECV=my_data->rankECV;
	unsigned long _svSize=my_data->svSize;
	unsigned long _vectSize=my_data->vectSize;

	for (unsigned long d=disBottom;d <disUp;d++){
	
		DoubleSquareMatrix &vuEvuT=(*(my_data->vuevuT))[d];
		vuEvuT.setAllValues(0.0);
		double *vuevut = vuEvuT.getArray();

		for(unsigned long i = 0; i<_rankECV; i++){
			for(unsigned long j = 0; j<=i; j++){
				for(unsigned long k=0; k<_vectSize;k++){
					vuevut[i*_rankECV+j] += vu[(i*_svSize)+(d*_vectSize)+k] * E[d*_vectSize+k] * vu[(j*_svSize)+(d*_vectSize)+k];
				}
			}
		}
		for(unsigned long i=0;i<_rankECV;i++){
			for(unsigned long j=i+1;j<_rankECV;j++) {
				vuevut[i*_rankECV+j]=vuevut[j*_rankECV+i];
			}
		}

	}
	pthread_exit((void*) 0);
	return (void*)0 ;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::estimateVUEVUTThreaded(unsigned long NUM_THREADS){
	
	if (verboseLevel >=1) cout << "(AccumulateJFAStat) Compute vuEvuT Threaded"<<endl;
	if (NUM_THREADS==0) throw Exception("Num threads can not be 0",__FILE__,__LINE__);
	
	double *VU=_VU.getArray();
	double *ubm_invvar=_ubm_invvar.getArray();
	
	int rc, status;
	if (NUM_THREADS > _n_distrib) NUM_THREADS=_n_distrib;
	
	struct VUEVUthread_data *thread_data_array = new VUEVUthread_data[NUM_THREADS];
	pthread_t *threads = new pthread_t[NUM_THREADS];

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	unsigned long offset=_n_distrib/NUM_THREADS;
	
	unsigned long disBottom = 0;
	unsigned long disUp=0;
	unsigned long re=_n_distrib - NUM_THREADS*offset;

	//Create threads
	for(unsigned long t=0; t<NUM_THREADS; t++){
	
		disUp = disBottom +offset;
		if(t<re) disUp +=1;
		
		thread_data_array[t].VU=VU;
		thread_data_array[t].ubm_invvar=ubm_invvar;
		thread_data_array[t].disBottom=disBottom;
		thread_data_array[t].disUp=disUp;
		thread_data_array[t].rankECV=_rankEC+_rankEV;
		thread_data_array[t].svSize=_svSize;
		thread_data_array[t].vectSize=_vectSize;
		thread_data_array[t].vuevuT=&(_vuEvuT);

		if (verboseLevel >1) cout<<"(AccumulateJFAStat) Creating thread n["<< t<< "] for distributions["<<disBottom<<"-->"<<disUp<<"]"<<endl;
		rc = pthread_create(&threads[t], &attr, VUEVUthread, (void *)&thread_data_array[t]);
		if (rc) throw Exception("ERROR; return code from pthread_create() is ",__FILE__,rc);
		
		disBottom = disUp;
	}
	
	pthread_attr_destroy(&attr);
	for(unsigned long t=0; t<NUM_THREADS; t++) {
		rc = pthread_join(threads[t], (void **)&status);
		if (rc)  throw Exception("ERROR; return code from pthread_join() is ",__FILE__,rc);
		if (verboseLevel >1) cout <<"(AccumulateJFAStat) Completed join with thread ["<<t<<"] status["<<status<<"]"<<endl;
	}

	free(thread_data_array);
	free(threads);

	if (verboseLevel >= 1) cout << "(AccumulateJFAStat) Done " << endl;
}
#endif

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void  JFAAcc::getVY(DoubleVector &vy, String &file){		///Compute VY for the speaker corresponding to the given file given file
	vy.setAllValues(0.0);	
	
	double *v, *y; y=_Y.getArray(); v=_V.getArray();

	unsigned long idx=_ndxTable.locNb(file);
	for (unsigned long i=0;i<_svSize;i++){
		for (unsigned long j=0;j<_rankEV;j++){
			vy[i] += v[j*_svSize + i] * y[idx*_rankEV + j ];
		}
	}	
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void  JFAAcc::getVY(DoubleVector &vy, unsigned long spk) {		///Compute VY for speaker spk
	vy.setAllValues(0.0);
	
	double *v, *y; y=_Y.getArray(); v=_V.getArray();
	for (unsigned long i=0;i<_svSize;i++){
		for (unsigned long j=0;j<_rankEV;j++){
			vy[i] += v[j*_svSize + i] * y[spk*_rankEV + j ];
		}
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void  JFAAcc::getVUYX(DoubleVector &vuyx, String &file){
	vuyx.setAllValues(0.0);
	
	double *vu, *yx; vu=_VU.getArray(); yx=_YX.getArray();
	
	unsigned long idx=_ndxTable.locNb(file);
	for (unsigned long i=0;i<_svSize;i++){
		for (unsigned long j=0;j<_rankEC + _rankEV;j++) {
			vuyx[i] += vu[j*_svSize + i] * yx[idx*(_rankEC + _rankEV)+j];
		}
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void  JFAAcc::getVUYX(DoubleVector &vuyx, unsigned long spk) {
	vuyx.setAllValues(0.0);
	
	double *vu, *yx; vu=_VU.getArray(); yx=_YX.getArray();
	
	for (unsigned long i=0;i<_svSize;i++){
		for (unsigned long j=0;j<_rankEC + _rankEV;j++){
			vuyx[i] += vu[j*_svSize + i] * yx[spk*(_rankEC + _rankEV)+j];
		}
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void  JFAAcc::getDZ(DoubleVector &dz, unsigned long spk) {	///Compute DZ for speaker spk
	dz.setAllValues(0.0);	
	double *z=_Z.getArray();
	for (unsigned long i=0;i<_svSize;i++){
		dz[i] = _D[i] * z[spk*_svSize +i];
	}
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void  JFAAcc::getUX(DoubleVector &ux, unsigned long ses) {		///Compute UX for session ses
	ux.setAllValues(0.0);

	ux.setAllValues(0.0);
	double *u, *x;
	u=_matU.getArray(); x=_matX.getArray();

	for (unsigned long i=0;i<_svSize;i++){
		for (unsigned long j=0;j<_rankEC;j++){
			ux[i] = _matU(j,i)*_matX(ses,j);
		}
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::getUX(DoubleVector &ux, String& file){
	
	ux.setAllValues(0.0);	
	double *u, *x;
	u=_matU.getArray(); x=_matX.getArray();

	unsigned long idx=_ndxTable.sessionNb(file);
	
	for (unsigned long i=0;i<_svSize;i++){
		for (unsigned long j=0;j<_rankEC;j++) {
			ux[i]+=_matU(j,i)*_matX(idx,j);
		}
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::getVYplusDZ(DoubleVector & VYplusDZ, String & file){
	VYplusDZ.setAllValues(0.0);
	DoubleVector VY(_svSize,_svSize);
	getVY(VY,file);
	
	double *vy, *vyplusdz, *z;
	vy=VY.getArray(); vyplusdz=VYplusDZ.getArray(); z=_Z.getArray();
	
	unsigned long spk=_ndxTable.locNb(file);		
	for (unsigned long i=0;i<_svSize;i++){
		vyplusdz[i] = vy[i] + _D[i] * z[spk*_svSize+i];
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::getVYplusDZ(DoubleVector &VYplusDZ, unsigned long spk){
	VYplusDZ.setAllValues(0.0);
	DoubleVector VY(_svSize,_svSize);
	getVY(VY,spk);
	
	double *vy, *vyplusdz, *z;
	vy=VY.getArray(); vyplusdz=VYplusDZ.getArray(); z=_Z.getArray();
	
	for (unsigned long i=0;i<_svSize;i++){
		vyplusdz[i] = vy[i] + _D[i] * z[spk*_svSize+i];
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Compute supervector of client M_s_h=M+DZ_s
void JFAAcc::getMplusDZ(DoubleVector &Sp, String& file){
	Sp.setAllValues(0.0);
	
	unsigned long spk=_ndxTable.locNb(file);		
	for (unsigned long i=0;i<_svSize;i++)
		Sp[i]=_ubm_means[i]+_D[i]*_Z(spk,i);
	
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Compute supervector of client M_s_h=M+DZ_s
void JFAAcc::getMplusDZ(DoubleVector &Sp, unsigned long spk){
	Sp.setAllValues(0.0);
	
	for (unsigned long i=0;i<_svSize;i++)
		Sp[i]=_ubm_means[i]+_D[i]*_Z(spk,i);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Compute supervector of client M_s_h=M+ VY_s
void JFAAcc::getMplusVY(DoubleVector &Sp, String& file){
	Sp.setAllValues(0.0);
	DoubleVector vy(_svSize,_svSize);
	getVY(vy,file);		
	for (unsigned long i=0;i<_svSize;i++){
		Sp[i] = _ubm_means[i] + vy[i];
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Compute supervector of client M_s_h=M+ VY_s
void JFAAcc::getMplusVY(DoubleVector &Sp, unsigned long spk){
	Sp.setAllValues(0.0);
	DoubleVector vy(_svSize,_svSize);
	getVY(vy,spk);
	for (unsigned long i=0;i<_svSize;i++){
		Sp[i]=_ubm_means[i] + vy[i];
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Compute supervector of client M_s_h=M+ UX
void JFAAcc::getMplusUX(DoubleVector &Sp, String& file){
	Sp.setAllValues(0.0);
	DoubleVector ux(_svSize,_svSize);
	getUX(ux,file);		
	for (unsigned long i=0;i<_svSize;i++){
		Sp[i] = _ubm_means[i] + ux[i];
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Compute supervector of client M_s_h=M+ UX
void JFAAcc::getMplusUX(DoubleVector &Sp, unsigned long session){
	Sp.setAllValues(0.0);
	DoubleVector ux(_svSize,_svSize);
	getUX(ux,session);
	for (unsigned long i=0;i<_svSize;i++){
		Sp[i]=_ubm_means[i] + ux[i];
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Compute supervector of client M_s_h=M+ VY_s + DZ_s
void JFAAcc::getMplusVYplusDZ(DoubleVector &Sp, String& file){
	Sp.setAllValues(0.0);
	
	DoubleVector vy(_svSize,_svSize);
	getVY(vy,file);	
	unsigned long spk=_ndxTable.locNb(file);		
	for (unsigned long i=0;i<_svSize;i++){
		Sp[i] = _ubm_means[i] + vy[i] +_D[i]*_Z(spk,i);
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Compute supervector of client M_s_h=M+ VY_s + DZ_s
void JFAAcc::getMplusVYplusDZ(DoubleVector &Sp, unsigned long spk){
	Sp.setAllValues(0.0);
	DoubleVector vy(_svSize,_svSize);
	getVY(vy,spk);
	for (unsigned long i=0;i<_svSize;i++){
		Sp[i]=_ubm_means[i] + vy[i] +_D[i]*_Z(spk,i);
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Compute supervector of client M_s_h=M+ VUY_sX
void JFAAcc::getMplusVUYX(DoubleVector &Sp, String& file){
	Sp.setAllValues(0.0);
	DoubleVector vuyx(_svSize,_svSize);
	getVUYX(vuyx,file);		
	for (unsigned long i=0;i<_svSize;i++){
		Sp[i] = _ubm_means[i] + vuyx[i];
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Compute supervector of client M_s_h=M+ VUY_sX
void JFAAcc::getMplusVUYX(DoubleVector &Sp, unsigned long spk){
	Sp.setAllValues(0.0);
	DoubleVector vuyx(_svSize,_svSize);
	this->getVUYX(vuyx,spk);
	
	for (unsigned long i=0;i<_svSize;i++){
		Sp[i]=_ubm_means[i] + vuyx[i];
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::estimateAndInverseL_EV(Config& config){
	
	#ifdef THREAD          
	if (config.existsParam("numThread") && config.getParam("numThread").toULong() >0)	estimateAndInverseLThreaded_EV(config.getParam("numThread").toULong());
	else estimateAndInverseLUnThreaded_EV();
	#else
	estimateAndInverseLUnThreaded_EV();
	#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::estimateAndInverseLUnThreaded_EV(){
	
	if (verboseLevel >= 1) cout << "(AccumulateJFAStat) Compute and Inverse L Matrix for EigenVoices "<<endl;
	
	DoubleSquareMatrix L(_rankEV);

	for(unsigned long spk=0; spk<_n_speakers; spk++){	
		L.setAllValues(0.0);
		for(unsigned long i=0; i<_rankEV; i++){	L(i,i)=1.0;}
		
		double *l, *n;
		l=L.getArray();  n=_matN.getArray();

		for(unsigned long dis=0; dis<_n_distrib;dis++){
			double *vevt=_vEvT[dis].getArray();
			for(unsigned long i=0; i<_rankEV; i++){
				for(unsigned long j=0; j<=i; j++){
					l[i*_rankEV+j] =+ l[i*_rankEV+j] + vevt[i*_rankEV+j]*n[spk*_n_distrib+dis];
				}
			}
		}
		//As L is symetric, copy the second half of coefficients
		for(unsigned long i=0;i<_rankEV;i++){
			for(unsigned long j=i+1;j<_rankEV;j++) {
				l[i*_rankEV+j] = l[j*_rankEV+i];
			}
		}
		
		//Inverse L and stock it in _l_spk_inv[spk]
		L.invert(_l_spk_inv[spk]);
	}
}

#ifdef THREAD
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Data strucutre of thread
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
struct LVthread_data{

	double *N;
	unsigned long spkBottom;
	unsigned long spkUp;	
	unsigned long rankEV;
	unsigned long n_distrib;
	unsigned long nt;
	RefVector <DoubleSquareMatrix>* vevT;
	RefVector <DoubleSquareMatrix>* linv;
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Thread Routine
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void *LVthread(void *threadarg) {
	struct LVthread_data *my_data;
	my_data = (struct LVthread_data *) threadarg;

	double *n = my_data->N;
	unsigned long spkBottom = my_data->spkBottom;	
	unsigned long spkUp = my_data->spkUp;
	unsigned long _rankEV=my_data->rankEV;
	unsigned long _n_distrib=my_data->n_distrib;

	for(unsigned long spk=spkBottom; spk<spkUp; spk++){
	
		DoubleSquareMatrix L(_rankEV);	
		L.setAllValues(0.0);	
		double *l=L.getArray();
		for(unsigned long i=0; i<_rankEV; i++){	l[i*_rankEV+i]=1.0;}

		for(unsigned long d=0; d<_n_distrib;d++){
			DoubleSquareMatrix &vEvT=(*(my_data->vevT))[d];
			double *vevt = vEvT.getArray();
			for(unsigned long i=0; i<_rankEV; i++){
				for(unsigned long j=0; j<=i; j++){
					l[i*_rankEV+j] =+ l[i*_rankEV+j] + vevt[i*_rankEV+j]*n[spk*_n_distrib+d];
				}
			}
			
		}
		//As L is symetric, copy the second half of coefficients
		for(unsigned long i=0;i<_rankEV;i++){
			for(unsigned long j=i+1;j<_rankEV;j++) {
				l[i*_rankEV+j] = l[j*_rankEV+i];
			}
		}

		DoubleSquareMatrix &linv=(*(my_data->linv))[spk];		
		//Inverse L and stock it in _l_spk_inv[spk]
		L.invert(linv);
	}
	pthread_exit((void*) 0);
	return (void*)0 ;
}



//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::estimateAndInverseLThreaded_EV(unsigned long NUM_THREADS){
	
	if (verboseLevel >= 1) cout << "(AccumulateJFAStat) Compute and Inverse L Matrix for EigenVoices Threaded"<<endl;
	if (NUM_THREADS==0) throw Exception("Num threads can not be 0",__FILE__,__LINE__);

	int rc, status;
	if (NUM_THREADS > _n_speakers) NUM_THREADS=_n_speakers;

	struct LVthread_data *thread_data_array = new LVthread_data[NUM_THREADS];
	pthread_t *threads = new pthread_t[NUM_THREADS];

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	unsigned long offset=_n_speakers/NUM_THREADS;
	
	double *N=_matN.getArray(); 
	
	unsigned long spkBottom = 0;
	unsigned long spkUp=0;
	unsigned long re=_n_speakers - NUM_THREADS*offset;
	
	//Create threads
	for(unsigned long t=0; t<NUM_THREADS; t++){
	
		spkUp = spkBottom +offset;
		if(t<re) spkUp +=1;

		thread_data_array[t].N=N;
		thread_data_array[t].spkBottom=spkBottom;
		thread_data_array[t].spkUp=spkUp;
		thread_data_array[t].rankEV=_rankEV;
		thread_data_array[t].n_distrib=_n_distrib;
		thread_data_array[t].vevT=&(_vEvT);
		thread_data_array[t].linv=&(_l_spk_inv);
		thread_data_array[t].nt=t;

		if (verboseLevel>1) cout<<"(AccumulateJFAStat) Creating thread n["<< t<< "] for speakers["<<spkBottom<<"-->"<<spkUp-1<<"]"<<endl;
		rc = pthread_create(&threads[t], &attr, LVthread, (void *)&thread_data_array[t]);
		if (rc) throw Exception("ERROR; return code from pthread_create() is ",__FILE__,rc);
		
		spkBottom = spkUp;
	}
	
	pthread_attr_destroy(&attr);
	for(unsigned long t=0; t<NUM_THREADS; t++) {
		rc = pthread_join(threads[t], (void **)&status);
		if (rc)  throw Exception("ERROR; return code from pthread_join() is ",__FILE__,rc);
		if (verboseLevel >1) cout <<"(AccumulateJFAStat) Completed join with thread ["<<t<<"] status["<<status<<"]"<<endl;
	}

	free(thread_data_array);
	free(threads);

	if (verboseLevel >= 1) cout << "(AccumulateJFAStat) Done " << endl;
	
}
#endif

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::estimateAndInverseL_EC(Config& config){
	#ifdef THREAD    
	if (config.existsParam("numThread") && config.getParam("numThread").toULong() >0)	estimateAndInverseLThreaded_EC(config.getParam("numThread").toULong());
	else estimateAndInverseLUnThreaded_EC();
	#else
	estimateAndInverseLUnThreaded_EC();
	#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::estimateAndInverseLUnThreaded_EC(){

	if (verboseLevel >= 1) cout << "(AccumulateJFAStat) Compute and Inverse L Matrix for EigenChannel "<<endl;

	DoubleSquareMatrix L(_rankEC);

	for(unsigned long sess=0; sess<_n_sessions; sess++){
		L.setAllValues(0.0);
		for(unsigned long i=0; i<_rankEC; i++){	L(i,i)=1.0;}

		double *l, *n_h;
		l=L.getArray();  n_h=_N_h.getArray();

		for(unsigned long dis=0; dis<_n_distrib;dis++){
			double *ueut=_uEuT[dis].getArray();
			
			for(unsigned long i=0; i<_rankEC; i++){
				for(unsigned long j=0; j<=i; j++){
					l[i*_rankEC+j] =+ l[i*_rankEC+j] + ueut[i*_rankEC+j]*n_h[sess*_n_distrib+dis];
				}
			}
		}
		//Copy the second half of the symetric matrix
		for(unsigned long i=0;i<_rankEC;i++){
			for(unsigned long j=i+1;j<_rankEC;j++) {
				l[i*_rankEC+j] = l[j*_rankEC+i];
			}
		}
	
		//Inverse L and stock in  _l_sess_inv[spk]
		L.invert(_l_sess_inv[sess]);
	}
}

#ifdef THREAD
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Data strucutre of thread
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
struct LUthread_data{

	double *N_h;
	unsigned long sessBottom;
	unsigned long sessUp;	
	unsigned long rankEC;
	unsigned long n_distrib;
	RefVector <DoubleSquareMatrix>* ueuT;
	RefVector <DoubleSquareMatrix>* linv;
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Thread Routine
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void *LUthread(void *threadarg) {
	struct LUthread_data *my_data;
	my_data = (struct LUthread_data *) threadarg;

	double *n_h = my_data->N_h;
	unsigned long sessBottom = my_data->sessBottom;	
	unsigned long sessUp = my_data->sessUp;
	unsigned long _rankEC=my_data->rankEC;
	unsigned long _n_distrib=my_data->n_distrib;

	for(unsigned long session=sessBottom; session<sessUp; session++){
	
		DoubleSquareMatrix L(_rankEC);	
		L.setAllValues(0.0);	
		double *l=L.getArray();
		for(unsigned long i=0; i<_rankEC; i++){	l[i*_rankEC+i]=1.0;}

		for(unsigned long d=0; d<_n_distrib;d++){
			DoubleSquareMatrix &uEuT=(*(my_data->ueuT))[d];
			double *ueut = uEuT.getArray();
			
			for(unsigned long i=0; i<_rankEC; i++){
				for(unsigned long j=0; j<=i; j++){
					l[i*_rankEC+j] =+ l[i*_rankEC+j] + ueut[i*_rankEC+j]*n_h[session*_n_distrib+d];
				}
			}
			
		}
		//As L is symetric, copy the second half of coefficients
		for(unsigned long i=0;i<_rankEC;i++){
			for(unsigned long j=i+1;j<_rankEC;j++) {
				l[i*_rankEC+j] = l[j*_rankEC+i];
			}
		}
		
		DoubleSquareMatrix &linv=(*(my_data->linv))[session];		
		//Inverse L and stock it in _l_sess_inv[spk]
		L.invert(linv);
	}
	pthread_exit((void*) 0);
	return (void*)0 ;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::estimateAndInverseLThreaded_EC(unsigned long NUM_THREADS){
	
	if (verboseLevel >= 1) cout << "(AccumulateJFAStat) Compute and Inverse L Matrix for EigenChannel Threaded"<<endl;
	if (NUM_THREADS==0) throw Exception("Num threads can not be 0",__FILE__,__LINE__);

	int rc, status;
	if (NUM_THREADS > _n_distrib) NUM_THREADS=_n_distrib;

	struct LUthread_data *thread_data_array = new LUthread_data[NUM_THREADS];
	pthread_t *threads = new pthread_t[NUM_THREADS];

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	unsigned long offset=_n_sessions/NUM_THREADS;
	double *N_h=_N_h.getArray(); 
	
	unsigned long sessBottom = 0;
	unsigned long sessUp=0;
	unsigned long re=_n_sessions - NUM_THREADS*offset;
	
	//Create threads
	for(unsigned long t=0; t<NUM_THREADS; t++){
		sessUp = sessBottom +offset;
		if(t<re) sessUp +=1;

		thread_data_array[t].N_h=N_h;
		thread_data_array[t].sessBottom=sessBottom;
		thread_data_array[t].sessUp=sessUp;
		thread_data_array[t].rankEC=_rankEC;
		thread_data_array[t].n_distrib=_n_distrib;
		thread_data_array[t].ueuT=&(_uEuT);
		thread_data_array[t].linv=&(_l_sess_inv);

		if (verboseLevel >1) cout<<"(AccumulateJFAStat) Creating thread n["<< t<< "] for sessions ["<<sessBottom<<"-->"<<sessUp<<"]"<<endl;
		rc = pthread_create(&threads[t], &attr, LUthread, (void *)&thread_data_array[t]);
		if (rc) throw Exception("ERROR; return code from pthread_create() is ",__FILE__,rc);
		
		sessBottom = sessUp;
	}
	
	pthread_attr_destroy(&attr);
	for(unsigned long t=0; t<NUM_THREADS; t++) {
		rc = pthread_join(threads[t], (void **)&status);
		if (rc)  throw Exception("ERROR; return code from pthread_join() is ",__FILE__,rc);
		if (verboseLevel >1) cout <<"(AccumulateJFAStat) Completed join with thread ["<<t<<"] status["<<status<<"]"<<endl;
	}

	free(thread_data_array);
	free(threads);

	if (verboseLevel >= 1) cout << "(AccumulateJFAStat) Done " << endl;
}
#endif

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::estimateAndInverseL_VU(Config& config){

	#ifdef THREAD          
	if (config.existsParam("numThread") && config.getParam("numThread").toULong() >0)	estimateAndInverseLThreaded_VU(config.getParam("numThread").toULong());
	else estimateAndInverseLUnThreaded_VU();
	#else
	estimateAndInverseLUnThreaded_VU();
	#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::estimateAndInverseLUnThreaded_VU(){
	
	if (verboseLevel >= 1) cout << "(AccumulateJFAStat) Compute and Inverse L Matrix for training "<<endl;
	
	DoubleSquareMatrix L(_rankEV + _rankEC);

	for(unsigned long spk=0; spk<_n_speakers; spk++){	
		L.setAllValues(0.0);
		for(unsigned long i=0; i<_rankEV + _rankEC; i++){	L(i,i)=1.0;}

		double *l, *n;
		l=L.getArray();  n=_matN.getArray();

		for(unsigned long dis=0; dis<_n_distrib;dis++){
			double *vuevut=_vuEvuT[dis].getArray();
			
			for(unsigned long i=0; i<_rankEV + _rankEC; i++){
				for(unsigned long j=0; j<=i; j++){
					l[i*(_rankEV + _rankEC)+j] =+ l[i*(_rankEV + _rankEC)+j] + vuevut[i*(_rankEV + _rankEC)+j]*n[spk*_n_distrib+dis];
				}
			}
		}
		//Copy the second half of the symetric matrix
		for(unsigned long i=0;i<(_rankEV + _rankEC);i++){
			for(unsigned long j=i+1;j<(_rankEV + _rankEC);j++) {
				l[i*(_rankEV + _rankEC)+j] = l[j*(_rankEV + _rankEC)+i];
			}
		}

		//Inverse L and stock in _l_sess_inv[spk]
		L.invert(_l_yx_inv[spk]);
	}
}

#ifdef THREAD
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Data strucutre of thread
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
struct LVUthread_data{

	double *N;
	unsigned long spkBottom;
	unsigned long spkUp;	
	unsigned long rankEVC;
	unsigned long n_distrib;
	unsigned long nt;
	RefVector <DoubleSquareMatrix>* vuevuT;
	RefVector <DoubleSquareMatrix>* linv;
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Thread Routine
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void *LVUthread(void *threadarg) {
	struct LVUthread_data *my_data;
	my_data = (struct LVUthread_data *) threadarg;

	double *n = my_data->N;
	unsigned long spkBottom = my_data->spkBottom;	
	unsigned long spkUp = my_data->spkUp;
	unsigned long _rankEVC=my_data->rankEVC;
	unsigned long _n_distrib=my_data->n_distrib;

	for(unsigned long spk=spkBottom; spk<spkUp; spk++){
	
		DoubleSquareMatrix L(_rankEVC);
		L.setAllValues(0.0);	
		double *l=L.getArray();
		for(unsigned long i=0; i<_rankEVC; i++){	l[i*_rankEVC+i]=1.0;}

		for(unsigned long d=0; d<_n_distrib;d++){
			DoubleSquareMatrix &vuEvuT=(*(my_data->vuevuT))[d];
			double *vuevut = vuEvuT.getArray();
			
			for(unsigned long i=0; i<_rankEVC; i++){
				for(unsigned long j=0; j<=i; j++){
					l[i*_rankEVC+j] =+ l[i*_rankEVC+j] + vuevut[i*_rankEVC+j]*n[spk*_n_distrib+d];
				}
			}
			
		}
		//As L is symetric, copy the second half of coefficients
		for(unsigned long i=0;i<_rankEVC;i++){
			for(unsigned long j=i+1;j<_rankEVC;j++) {
				l[i*_rankEVC+j] = l[j*_rankEVC+i];
			}
		}

		DoubleSquareMatrix &linv=(*(my_data->linv))[spk];		
		//Inverse L and stock it in _l_yx_inv[spk]
		L.invert(linv);
	}
	pthread_exit((void*) 0);
	return (void*)0 ;
}



//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::estimateAndInverseLThreaded_VU(unsigned long NUM_THREADS){
	
	if (verboseLevel >= 1) cout << "(AccumulateJFAStat) Compute and Inverse L Matrix for Training Threaded"<<endl;
	if (NUM_THREADS==0) throw Exception("Num threads can not be 0",__FILE__,__LINE__);

	int rc, status;
	if (NUM_THREADS > _n_speakers) NUM_THREADS=_n_speakers;

	struct LVUthread_data *thread_data_array = new LVUthread_data[NUM_THREADS];
	pthread_t *threads = new pthread_t[NUM_THREADS];


	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	unsigned long offset=_n_speakers/NUM_THREADS;
	
	double *N=_N.getArray(); 
	
	unsigned long spkBottom = 0;
	unsigned long spkUp=0;
	unsigned long re=_n_speakers - NUM_THREADS*offset;
	
	//Create threads
	for(unsigned long t=0; t<NUM_THREADS; t++){
		spkUp = spkBottom +offset;
		if(t<re) spkUp +=1;

		thread_data_array[t].N=N;
		thread_data_array[t].spkBottom=spkBottom;
		thread_data_array[t].spkUp=spkUp;
		thread_data_array[t].rankEVC=_rankEV + _rankEC;
		thread_data_array[t].n_distrib=_n_distrib;
		thread_data_array[t].vuevuT=&(_vuEvuT);
		thread_data_array[t].linv=&(_l_yx_inv);

		if (verboseLevel >1) cout<<"(AccumulateJFAStat) Creating thread n["<< t<< "] for speakers["<<spkBottom<<"-->"<<spkUp-1<<"]"<<endl;
		rc = pthread_create(&threads[t], &attr, LVUthread, (void *)&thread_data_array[t]);
		if (rc) throw Exception("ERROR; return code from pthread_create() is ",__FILE__,rc);
		
		spkBottom = spkUp;
	}
	
	pthread_attr_destroy(&attr);
	for(unsigned long t=0; t<NUM_THREADS; t++) {
		rc = pthread_join(threads[t], (void **)&status);
		if (rc)  throw Exception("ERROR; return code from pthread_join() is ",__FILE__,rc);
		if (verboseLevel>1) cout <<"(AccumulateJFAStat) Completed join with thread ["<<t<<"] status["<<status<<"]"<<endl;
	}
	if (verboseLevel >= 1) cout << "(AccumulateJFAStat) Done " << endl;

	free(thread_data_array);
	free(threads);

}
#endif

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::estimateYandV(Config& config){
	#ifdef THREAD          
	if (config.existsParam("numThread") && config.getParam("numThread").toULong() >0)	estimateYandVThreaded(config.getParam("numThread").toULong());
	else estimateYandVUnThreaded();
	#else
	estimateYandVUnThreaded();
	#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::estimateYandVUnThreaded(){	//estimate Y for all speakers and V
	
	if (verboseLevel >=1) cout << "(AccumulateJFAStat) Compute Y  and V Estimate "<<endl;
	_Y.setAllValues(0.0);
	Matrix<double> AUX(_n_speakers,_rankEV);
	AUX.setAllValues(0.0);
	
	double *y, *v, *f_x, *aux, *invVar, *n, *c;
	y=_Y.getArray(); v=_V.getArray(); f_x=_F_X.getArray(); aux=AUX.getArray(); invVar=_ubm_invvar.getArray(); n=_matN.getArray(); c=_Cev.getArray();

	//For each speaker
	for(unsigned long spk=0; spk<_n_speakers; spk++){
		double *invl = _l_spk_inv[spk].getArray();

		for(unsigned long i=0;i<_rankEV;i++){
			for(unsigned long k=0;k<_svSize;k++) {
				aux[spk*_rankEV+i] += f_x[spk*_svSize+k] * invVar[k] * v[i*_svSize+k];
			}
		}
	
		//multiplication by invL
		for(unsigned long i=0; i<_rankEV;i++){
			for(unsigned long k=0; k<_rankEV; k++){
				y[spk*_rankEV+i] += aux[spk*_rankEV+k] * invl[i*_rankEV+k];
			}
		}

		for(unsigned long i=0;i<_rankEV;i++){
			for(unsigned long j=0;j<_rankEV;j++){
					invl[i*_rankEV+j] += y[spk*_rankEV+i]*y[spk*_rankEV+j];
			}
		}

		for(unsigned long dis = 0; dis<_n_distrib;dis++){
			for(unsigned long i=0;i<_rankEV;i++){
				for(unsigned long j = 0;j<_rankEV;j++){
					_Aev[dis](i,j) += invl[i*_rankEV+j] * n[spk*_n_distrib+dis];

				}
			}
		}

		for(unsigned long i=0;i<_rankEV;i++){
			for(unsigned long j=0;j<_svSize;j++){
					c[i*_svSize+j] += y[spk*_rankEV+i] * f_x[spk*_svSize+j];
			}
		}
	}
}

#ifdef THREAD
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Data strucutre of thread
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
struct estimateYandVthread_data{
	double *AUX;
	double *Y;
	double *V;
	double *F_X;
	double *ubm_invvar;
	unsigned long spkBottom;
	unsigned long spkUp;	
	unsigned long rankEV;
	unsigned long svSize;
	unsigned long n_distrib;
	unsigned long nt;
	RefVector <DoubleSquareMatrix>* linv;
};


struct computeC_data{
	double* tmpC;
	double *Y;
	double *F_X;
	unsigned long numThread;
	unsigned long spkBottom;
	unsigned long spkUp;	
	unsigned long rankEV;
	unsigned long svSize;
};


struct computeINVLandA_data{
	double *N;
	double *Y;
	unsigned long spk;
	unsigned long _n_distrib;
	unsigned long numThread;
	unsigned long rankBottom;
	unsigned long rankUp;	
	unsigned long rankEV;
	unsigned long svSize;
	RefVector <DoubleSquareMatrix>* A;
	RefVector <DoubleSquareMatrix>* linv;
};



//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Thread Routine
//-----------------------------------------------------------------------------------------------------------------------------------------------------------

void *estimateYandVthread(void *threadarg) {
	
	struct estimateYandVthread_data *my_data;
	my_data = (struct estimateYandVthread_data *) threadarg;

	double *v = my_data->V;
	double *y = my_data->Y;
	double *f_x = my_data->F_X;
	double *aux = my_data->AUX;
	double *ubm_invvar = my_data->ubm_invvar;
	unsigned long spkBottom = my_data->spkBottom;	
	unsigned long spkUp = my_data->spkUp;
	unsigned long _rankEV=my_data->rankEV;
	unsigned long _svSize=my_data->svSize;
	
	//For each session
	for(unsigned long spk= spkBottom; spk<spkUp; spk++){

		DoubleSquareMatrix &linv=(*(my_data->linv))[spk];
		double *invl = linv.getArray();

		for(unsigned long i=0;i<_rankEV;i++){
			for(unsigned long k=0;k<_svSize;k++) {
				aux[spk*_rankEV+i] += f_x[spk*_svSize+k] * ubm_invvar[k] * v[i*_svSize+k];
			}
		}
		
		//multiplication by invL
		for(unsigned long i=0; i<_rankEV;i++){
			for(unsigned long k=0; k<_rankEV; k++){
				y[spk*_rankEV+i] += aux[spk*_rankEV+k] * invl[i*_rankEV+k];
			}
		}
	}
	pthread_exit((void*) 0);
	return (void*)0 ;
}


void *computeC(void *threadarg){
	
	struct computeC_data *my_data;
	my_data = (struct computeC_data *) threadarg;

	unsigned long numThread = my_data->numThread;
	double *tmpc = my_data->tmpC;
	double *y = my_data->Y;
	double *f_x = my_data->F_X;
	unsigned long spkBottom = my_data->spkBottom;	
	unsigned long spkUp = my_data->spkUp;
	unsigned long _rankEV=my_data->rankEV;
	unsigned long _svSize=my_data->svSize;
	
	//For each session
	for(unsigned long spk= spkBottom; spk<spkUp; spk++){
		for(unsigned long i=0;i<_rankEV;i++){
			for(unsigned long j=0;j<_svSize;j++){
				tmpc[(numThread*_rankEV+i)*_svSize+j] += y[spk*_rankEV+i] * f_x[spk*_svSize+j];
			}
		}
	}
	pthread_exit((void*) 0);
	return (void*)0 ;
}



void *computeINVLandA(void *threadarg){

	struct computeINVLandA_data *my_data;
	my_data = (struct computeINVLandA_data *) threadarg;

	double *n = my_data->N;
	double *y = my_data->Y;
	unsigned long spk = my_data->spk;
	unsigned long _n_distrib = my_data->_n_distrib;
	unsigned long rankBottom = my_data->rankBottom;	
	unsigned long rankUp = my_data->rankUp;
	unsigned long _rankEV=my_data->rankEV;

	DoubleSquareMatrix &linv=(*(my_data->linv))[spk];
	double *invl = linv.getArray();

	for(unsigned long i= rankBottom; i<rankUp; i++){
		for(unsigned long j=0;j<_rankEV;j++){
			invl[i*_rankEV+j] += y[spk*_rankEV+i]*y[spk*_rankEV+j];
			for(unsigned long dis = 0; dis<_n_distrib;dis++){
				DoubleSquareMatrix &A=(*(my_data->A))[dis];
				double *a = A.getArray();
				a[i*_rankEV+j] += invl[i*_rankEV+j] * n[spk*_n_distrib+dis];
			}
		}
	}
	pthread_exit((void*) 0);
	return (void*)0 ;
}


//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::estimateYandVThreaded(unsigned long NUM_THREADS){

	if (verboseLevel >= 1) cout << "(AccumulateJFAStat) Estimate Y and V for each session Threaded"<<endl;
	if (NUM_THREADS==0) throw Exception("Num threads can not be 0",__FILE__,__LINE__);

	int rc, status;
	if (NUM_THREADS > _n_speakers) NUM_THREADS=_n_speakers;

	struct estimateYandVthread_data *thread_data_array = new estimateYandVthread_data[NUM_THREADS];
	pthread_t *threads = new pthread_t[NUM_THREADS];

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	unsigned long offset=_n_speakers/NUM_THREADS;

	_Y.setAllValues(0.0);
	Matrix<double> AUX(_n_speakers,_rankEV);
	AUX.setAllValues(0.0);
	
	double *y, *f_x, *n, *c;
	y=_Y.getArray(); f_x=_F_X.getArray(); n=_matN.getArray(); c=_Cev.getArray();
	
	double *V =_V.getArray();
	double *Y =_Y.getArray();
	double *F_X =_F_X.getArray();
	double *aux = AUX.getArray(); 
	double *ubm_invvar=_ubm_invvar.getArray();

	unsigned long spkBottom = 0;
	unsigned long spkUp=0;
	unsigned long re=_n_speakers - NUM_THREADS*offset;
	
	//Create threads
	for(unsigned long t=0; t<NUM_THREADS; t++){
		spkUp = spkBottom +offset;
		if(t<re) spkUp +=1;

		thread_data_array[t].V=V;
		thread_data_array[t].Y=Y;
		thread_data_array[t].F_X=F_X;
		thread_data_array[t].AUX=aux;
		thread_data_array[t].ubm_invvar=ubm_invvar;
		thread_data_array[t].spkBottom=spkBottom;
		thread_data_array[t].spkUp=spkUp;
		thread_data_array[t].rankEV=_rankEV;
		thread_data_array[t].svSize=_svSize;
		thread_data_array[t].linv=&(_l_spk_inv);

		if (verboseLevel >1) cout<<"(AccumulateJFAStat) Creating thread n [ "<< t<< " ] for speakers[ "<<spkBottom<<" --> "<<spkUp-1<<" ]"<<endl;
		rc = pthread_create(&threads[t], &attr, estimateYandVthread, (void *)&thread_data_array[t]);
		if (rc) throw Exception("ERROR; return code from pthread_create() is ",__FILE__,rc);
		
		spkBottom = spkUp;
	}
	
	pthread_attr_destroy(&attr);
	for(unsigned long t=0; t<NUM_THREADS; t++) {
		rc = pthread_join(threads[t], (void **)&status);
		if (rc)  throw Exception("ERROR; return code from pthread_join() is ",__FILE__,rc);
		if (verboseLevel>1) cout <<"(AccumulateJFAStat) Completed join with thread ["<<t<<"] status["<<status<<"]"<<endl;
	}
	
	free(thread_data_array);
	free(threads);
	

//*************************************************************************
	//Create a temporary matrix in by concatenating C matrices
	Matrix<double> _tmpC;
	_tmpC.setDimensions(_rankEV*NUM_THREADS,_svSize);

	//Initialize tmpC
	_tmpC.setAllValues(0.0);

	struct computeC_data *thread_data_array2 = new computeC_data[NUM_THREADS];
	pthread_t *threads2 = new pthread_t[NUM_THREADS];

	pthread_attr_t attr2;
	pthread_attr_init(&attr2);
	pthread_attr_setdetachstate(&attr2, PTHREAD_CREATE_JOINABLE);
	unsigned long offset2=_n_speakers/NUM_THREADS;
	
	double *tmpc;
	f_x=_F_X.getArray(); tmpc=_tmpC.getArray();
	
	spkBottom = 0;
	spkUp=0;
	re=_n_speakers - NUM_THREADS*offset2;

	//Create threads
	for(unsigned long t=0; t<NUM_THREADS; t++){
		spkUp = spkBottom +offset2;
		if(t<re) spkUp +=1;

		thread_data_array2[t].Y=Y;
		thread_data_array2[t].F_X = f_x;
		thread_data_array2[t].numThread = t;
		thread_data_array2[t].spkBottom = spkBottom;
		thread_data_array2[t].spkUp = spkUp;	
		thread_data_array2[t].rankEV = _rankEV;
		thread_data_array2[t].svSize = _svSize;
		thread_data_array2[t].tmpC = tmpc;

		if (verboseLevel >1) cout<<"(AccumulateJFAStat) Compute C creating thread n [ "<< t<< " ] for speakers[ "<<spkBottom<<" --> "<<spkUp-1<<" ]"<<endl;
			rc = pthread_create(&threads2[t], &attr2, computeC, (void *)&thread_data_array2[t]);
		if (rc) throw Exception("ERROR; return code from pthread_create() is ",__FILE__,rc);

		spkBottom = spkUp;
	}
	
	pthread_attr_destroy(&attr2);
	for(unsigned long t=0; t<NUM_THREADS; t++) {
		rc = pthread_join(threads2[t], (void **)&status);
		if (rc)  throw Exception("ERROR; return code from pthread_join() is ",__FILE__,rc);
		if (verboseLevel>1) cout <<"(AccumulateJFAStat) Completed join with thread ["<<t<<"] status["<<status<<"]"<<endl;
	}
	
	free(thread_data_array2);
	free(threads2);

	//Sum matrix C after multithreading
	for(unsigned long mt=0; mt<NUM_THREADS;mt++){
		for(unsigned long i=0; i<_rankEV;i++){
			for(unsigned long j=0; j<_svSize;j++){
				c[i*_svSize+j] += _tmpC(mt*_rankEV+i,j);
			}
		}
	}


///////////////////////////////////////////////////////////////////




	for(unsigned long spk=0; spk<_n_speakers; spk++){

		struct computeINVLandA_data *thread_data_array3 = new computeINVLandA_data[NUM_THREADS];
		pthread_t *threads3 = new pthread_t[NUM_THREADS];

		pthread_attr_t attr3;
		pthread_attr_init(&attr3);
		pthread_attr_setdetachstate(&attr3, PTHREAD_CREATE_JOINABLE);
		unsigned long offset3=_rankEV/NUM_THREADS;

		unsigned long rankBottom = 0;
		unsigned long rankUp=0;
		re=_rankEV - NUM_THREADS*offset3;

		//Create threads
		for(unsigned long t=0; t<NUM_THREADS; t++){
			rankUp = rankBottom +offset3;
			if(t<re) rankUp +=1;

			thread_data_array3[t].N = n;
			thread_data_array3[t].Y=Y;
			thread_data_array3[t].spk = spk;
			thread_data_array3[t]._n_distrib = _n_distrib;
			thread_data_array3[t].numThread = t;
			thread_data_array3[t].rankBottom = rankBottom;
			thread_data_array3[t].rankUp = rankUp;	
			thread_data_array3[t].rankEV = _rankEV;
			thread_data_array3[t].svSize = _svSize;
			thread_data_array3[t].A = &(_Aev);
			thread_data_array3[t].linv=&(_l_spk_inv);

			if (verboseLevel >2) cout<<"(AccumulateJFAStat) ComputeLinvandA creating thread n [ "<< t<< " ] for speakers[ "<<rankBottom<<" --> "<<rankUp-1<<" ]"<<endl;
				rc = pthread_create(&threads3[t], &attr3, computeINVLandA, (void *)&thread_data_array3[t]);
			if (rc) throw Exception("ERROR; return code from pthread_create() is ",__FILE__,rc);
			rankBottom = rankUp;
		}
	
		pthread_attr_destroy(&attr3);
		for(unsigned long t=0; t<NUM_THREADS; t++) {
			rc = pthread_join(threads3[t], (void **)&status);
			if (rc)  throw Exception("ERROR; return code from pthread_join() is ",__FILE__,rc);
			if (verboseLevel>2) cout <<"(AccumulateJFAStat) Completed join with thread ["<<t<<"] status["<<status<<"]"<<endl;
		}
	
		free(thread_data_array3);
		free(threads3);
	}

	if (verboseLevel >= 1) cout << "(AccumulateJFAStat) Done " << endl;
}
#endif

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::estimateY(Config& config){
	#ifdef THREAD          
	if (config.existsParam("numThread") && config.getParam("numThread").toULong() >0)	estimateYThreaded(config.getParam("numThread").toULong());
	else estimateYUnThreaded();
	#else
	estimateYUnThreaded();
	#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::estimateYUnThreaded(){	//estimate Y for all speakers
	
	if (verboseLevel>=1) cout << "(AccumulateJFAStat) Compute Y  Estimate "<<endl;
	_Y.setAllValues(0.0);
	Matrix<double> AUX(_n_speakers,_rankEV);
	AUX.setAllValues(0.0);
	
	double *y, *v, *f_x, *aux, *invVar;
	y=_Y.getArray(); v=_V.getArray(); f_x=_F_X.getArray(); aux=AUX.getArray(); invVar=_ubm_invvar.getArray();

	//For each speaker
	for(unsigned long spk=0; spk<_n_speakers; spk++){

		double *invl = _l_spk_inv[spk].getArray();

		for(unsigned long i=0;i<_rankEV;i++){
			for(unsigned long k=0;k<_svSize;k++) {
				aux[spk*_rankEV+i] += f_x[spk*_svSize+k] * invVar[k] * v[i*_svSize+k];
			}
		}
		
		//multiplication by invL
		for(unsigned long i=0; i<_rankEV;i++){
			for(unsigned long k=0; k<_rankEV; k++){
				y[spk*_rankEV+i] += aux[spk*_rankEV+k] * invl[i*_rankEV+k];
			}
		}
	}
}

#ifdef THREAD
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Data strucutre of thread
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
struct estimateYthread_data{
	double *AUX;
	double *Y;
	double *V;
	double *F_X;
	double *ubm_invvar;
	unsigned long spkBottom;
	unsigned long spkUp;	
	unsigned long rankEV;
	unsigned long svSize;
	unsigned long n_distrib;
	unsigned long nt;
	RefVector <DoubleSquareMatrix>* linv;
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Thread Routine
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void *estimateYthread(void *threadarg) {
	
	struct estimateYthread_data *my_data;
	my_data = (struct estimateYthread_data *) threadarg;

	double *v = my_data->V;
	double *y = my_data->Y;
	double *f_x = my_data->F_X;
	double *aux = my_data->AUX;
	double *ubm_invvar = my_data->ubm_invvar;
	unsigned long spkBottom = my_data->spkBottom;	
	unsigned long spkUp = my_data->spkUp;
	unsigned long _rankEV=my_data->rankEV;
	unsigned long _svSize=my_data->svSize;
	
	//For each session
	for(unsigned long spk= spkBottom; spk<spkUp; spk++){

		DoubleSquareMatrix &linv=(*(my_data->linv))[spk];
		double *invl = linv.getArray();

		for(unsigned long i=0;i<_rankEV;i++){
			for(unsigned long k=0;k<_svSize;k++) {
				aux[spk*_rankEV+i] += f_x[spk*_svSize+k] * ubm_invvar[k] * v[i*_svSize+k];
			}
		}
		
		//multiplication by invL
		for(unsigned long i=0; i<_rankEV;i++){
			for(unsigned long k=0; k<_rankEV; k++){
				y[spk*_rankEV+i] += aux[spk*_rankEV+k] * invl[i*_rankEV+k];
			}
		}
	}
	pthread_exit((void*) 0);
	return (void*)0 ;
}


//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::estimateYThreaded(unsigned long NUM_THREADS){

	if (verboseLevel >= 1) cout << "(AccumulateJFAStat) Estimate Y for each speaker Threaded"<<endl;
	if (NUM_THREADS==0) throw Exception("Num threads can not be 0",__FILE__,__LINE__);

	int rc, status;
	if (NUM_THREADS > _n_speakers) NUM_THREADS=_n_speakers;

	struct estimateYthread_data *thread_data_array = new estimateYthread_data[NUM_THREADS];
	pthread_t *threads = new pthread_t[NUM_THREADS];

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	unsigned long offset=_n_speakers/NUM_THREADS;

	_Y.setAllValues(0.0);
	Matrix<double> AUX(_n_speakers,_rankEV);
	AUX.setAllValues(0.0);
	
	double *y, *f_x, *c;
	y=_Y.getArray(); f_x=_F_X.getArray(); c=_Cev.getArray();
	
	double *V =_V.getArray();
	double *Y =_Y.getArray();
	double *F_X =_F_X.getArray();
	double *aux = AUX.getArray(); 
	double *ubm_invvar=_ubm_invvar.getArray();

	unsigned long spkBottom = 0;
	unsigned long spkUp=0;
	unsigned long re=_n_speakers - NUM_THREADS*offset;
	
	//Create threads
	for(unsigned long t=0; t<NUM_THREADS; t++){
		spkUp = spkBottom +offset;
		if(t<re) spkUp +=1;

		thread_data_array[t].V=V;
		thread_data_array[t].Y=Y;
		thread_data_array[t].F_X=F_X;
		thread_data_array[t].AUX=aux;
		thread_data_array[t].ubm_invvar=ubm_invvar;
		thread_data_array[t].spkBottom=spkBottom;
		thread_data_array[t].spkUp=spkUp;
		thread_data_array[t].rankEV=_rankEV;
		thread_data_array[t].svSize=_svSize;
		thread_data_array[t].linv=&(_l_spk_inv);

		if (verboseLevel>1) cout<<"(AccumulateJFAStat) Creating thread n [ "<< t<< " ] for speakers[ "<<spkBottom<<" --> "<<spkUp-1<<" ]"<<endl;
		rc = pthread_create(&threads[t], &attr, estimateYandVthread, (void *)&thread_data_array[t]);
		if (rc) throw Exception("ERROR; return code from pthread_create() is ",__FILE__,rc);
		
		spkBottom = spkUp;
	}
	
	pthread_attr_destroy(&attr);
	for(unsigned long t=0; t<NUM_THREADS; t++) {
		rc = pthread_join(threads[t], (void **)&status);
		if (rc)  throw Exception("ERROR; return code from pthread_join() is ",__FILE__,rc);
		if (verboseLevel>1) cout <<"(AccumulateJFAStat) Completed join with thread ["<<t<<"] status["<<status<<"]"<<endl;
	}

	free(thread_data_array);
	free(threads);

	if (verboseLevel >= 1) cout << "(AccumulateJFAStat) Done " << endl;
}
#endif

//----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::estimateXandU(Config &config){
	#ifdef THREAD          
	if (config.existsParam("numThread") && config.getParam("numThread").toULong() >0)	estimateXandUThreaded(config.getParam("numThread").toULong());
	else estimateXandUUnThreaded();
	#else
	estimateXandUUnThreaded();
	#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::estimateXandUUnThreaded(){	//estimate X for all sessions and U
	
	if (verboseLevel>=1) cout << "(AccumulateJFAStat) Compute X  andU Estimate "<<endl;
	_matX.setAllValues(0.0);
	Matrix<double> AUX(_n_sessions,_rankEC);
	AUX.setAllValues(0.0);
	
	double *x, *u, *f_x_h, *aux, *invVar, *n_h, *c;
	x=_matX.getArray(); u=_matU.getArray(); f_x_h=_F_X_h.getArray(); aux=AUX.getArray(); invVar=_ubm_invvar.getArray(); n_h=_N_h.getArray(); c=_Cec.getArray();

	//For each session
	for(unsigned long session=0; session<_n_sessions; session++){

		double *invl = _l_sess_inv[session].getArray();

		for(unsigned long i=0;i<_rankEC;i++){
			for(unsigned long k=0;k<_svSize;k++) {
				aux[session*_rankEC+i] += f_x_h[session*_svSize+k] * invVar[k] * u[i*_svSize+k];
			}
		}
		
		//multiplication by invL
		for(unsigned long i=0; i<_rankEC;i++){
			for(unsigned long k=0; k<_rankEC; k++){
				x[session*_rankEC+i] += aux[session*_rankEC+k] * invl[i*_rankEC+k];
			}
		}

		for(unsigned long i=0;i<_rankEC;i++){
			for(unsigned long j=0;j<_rankEC;j++){
					invl[i*_rankEC+j] += x[session*_rankEC+i]*x[session*_rankEC+j];
			}
		}

		for(unsigned long dis = 0; dis<_n_distrib;dis++){
			for(unsigned long i=0;i<_rankEC;i++){
				for(unsigned long j = 0;j<_rankEC;j++){
					_Aec[dis](i,j) += invl[i*_rankEC+j] * n_h[session*_n_distrib+dis];

				}
			}
		}

		for(unsigned long i=0;i<_rankEC;i++){
			for(unsigned long j=0;j<_svSize;j++){
					c[i*_svSize+j] += x[session*_rankEC+i] * f_x_h[session*_svSize+j];
			}
		}
	}
}

#ifdef THREAD
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Data strucutre of thread
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
struct estimateXandUthread_data{
	double *U;
	double *X;
	double *F_X_h;
	double *AUX;
	double *ubm_invvar;
	unsigned long sessBottom;
	unsigned long sessUp;	
	unsigned long rankEC;
	unsigned long svSize;
	unsigned long nt;
	RefVector <DoubleSquareMatrix>* linv;
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Thread Routine
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void *estimateXandUthread(void *threadarg) {
	
	struct estimateXandUthread_data *my_data;
	my_data = (struct estimateXandUthread_data *) threadarg;

	double *u = my_data->U;
	double *x = my_data->X;
	double *f_x_h = my_data->F_X_h;
	double *aux = my_data->AUX;
	double *ubm_invvar = my_data->ubm_invvar;
	unsigned long sessBottom = my_data->sessBottom;	
	unsigned long sessUp = my_data->sessUp;
	unsigned long _rankEC=my_data->rankEC;
	unsigned long _svSize=my_data->svSize;
	
	//For each session
	for(unsigned long session= sessBottom; session<sessUp; session++){

		DoubleSquareMatrix &linv=(*(my_data->linv))[session];
		double *invl = linv.getArray();

		for(unsigned long i=0;i<_rankEC;i++){
			for(unsigned long k=0;k<_svSize;k++) {
				aux[session*_rankEC+i] += f_x_h[session*_svSize+k] * ubm_invvar[k] * u[i*_svSize+k];
			}
		}
		
		//multiplication by invL
		for(unsigned long i=0; i<_rankEC;i++){
			for(unsigned long k=0; k<_rankEC; k++){
				x[session*_rankEC+i] += aux[session*_rankEC+k] * invl[i*_rankEC+k];
			}
		}
	}
	pthread_exit((void*) 0);
	return (void*)0 ;
}


//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::estimateXandUThreaded(unsigned long NUM_THREADS){

	if (verboseLevel >= 1) cout << "(AccumulateJFAStat) Estimate X and U for each session Threaded"<<endl;
	if (NUM_THREADS==0) throw Exception("Num threads can not be 0",__FILE__,__LINE__);

	int rc, status;
	if (NUM_THREADS > _n_sessions) NUM_THREADS=_n_sessions;

	struct estimateXandUthread_data *thread_data_array = new estimateXandUthread_data[NUM_THREADS];
	pthread_t *threads = new pthread_t[NUM_THREADS];

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	unsigned long offset=_n_sessions/NUM_THREADS;

	_matX.setAllValues(0.0);
	Matrix<double> AUX(_n_sessions,_rankEC);
	AUX.setAllValues(0.0);
	
	double *x, *f_x_h, *n_h, *c;
	x=_matX.getArray(); f_x_h=_F_X_h.getArray(); n_h=_N_h.getArray(); c=_Cec.getArray();
	
	double *U =_matU.getArray();
	double *X =_matX.getArray();
	double *F_X_h =_F_X_h.getArray();
	double *aux = AUX.getArray(); 
	double *ubm_invvar=_ubm_invvar.getArray();

	unsigned long sessBottom = 0;
	unsigned long sessUp=0;
	unsigned long re=_n_sessions - NUM_THREADS*offset;
	
	//Create threads
	for(unsigned long t=0; t<NUM_THREADS; t++){
		sessUp = sessBottom +offset;
		if(t<re) sessUp +=1;

		thread_data_array[t].U=U;
		thread_data_array[t].X=X;
		thread_data_array[t].F_X_h=F_X_h;
		thread_data_array[t].AUX=aux;
		thread_data_array[t].ubm_invvar=ubm_invvar;
		thread_data_array[t].sessBottom=sessBottom;
		thread_data_array[t].sessUp=sessUp;
		thread_data_array[t].rankEC=_rankEC;
		thread_data_array[t].svSize=_svSize;
		thread_data_array[t].linv=&(_l_sess_inv);

		if (verboseLevel>1) cout<<"(AccumulateJFAStat) Creating thread n["<< t<< "] for speakers["<<sessBottom<<"-->"<<sessUp-1<<"]"<<endl;
		rc = pthread_create(&threads[t], &attr, estimateXandUthread, (void *)&thread_data_array[t]);
		if (rc) throw Exception("ERROR; return code from pthread_create() is ",__FILE__,rc);
		
		sessBottom = sessUp;
	}
	
	pthread_attr_destroy(&attr);
	for(unsigned long t=0; t<NUM_THREADS; t++) {
		rc = pthread_join(threads[t], (void **)&status);
		if (rc)  throw Exception("ERROR; return code from pthread_join() is ",__FILE__,rc);
		if (verboseLevel>1) cout <<"(AccumulateJFAStat) Completed join with thread ["<<t<<"] status["<<status<<"]"<<endl;
	}

	free(thread_data_array);
	free(threads);
	
	
	//For each session
	for(unsigned long session=0; session<_n_sessions; session++){
		
		double *invl = _l_sess_inv[session].getArray();

		for(unsigned long i=0;i<_rankEC;i++){
			for(unsigned long j=0;j<_rankEC;j++){
					invl[i*_rankEC+j] += x[session*_rankEC+i]*x[session*_rankEC+j];
			}
		}

		for(unsigned long dis = 0; dis<_n_distrib;dis++){
			for(unsigned long i=0;i<_rankEC;i++){
				for(unsigned long j = 0;j<_rankEC;j++){
					_Aec[dis](i,j) += invl[i*_rankEC+j] * n_h[session*_n_distrib+dis];

				}
			}
		}

		for(unsigned long i=0;i<_rankEC;i++){
			for(unsigned long j=0;j<_svSize;j++){
					c[i*_svSize+j] += x[session*_rankEC+i] * f_x_h[session*_svSize+j];
			}
		}
	}

	if (verboseLevel >= 1) cout << "(AccumulateJFAStat) Done " << endl;
}
#endif


//----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::estimateX(Config &config){
	#ifdef THREAD    
	if (config.existsParam("numThread") && config.getParam("numThread").toULong() >0)	estimateXThreaded(config.getParam("numThread").toULong());
	else estimateXUnThreaded();
	#else
	estimateXUnThreaded();
	#endif
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::estimateXUnThreaded(){	//estimate X for all speakers

	if (verboseLevel>=1) cout << "(AccumulateJFAStat) Compute X Estimate "<<endl;
	_matX.setAllValues(0.0);
	Matrix<double> AUX(_n_sessions,_rankEC);
	AUX.setAllValues(0.0);

	double *x, *u, *f_x_h, *aux, *invVar;
	x=_matX.getArray(); u=_matU.getArray(); f_x_h=_F_X_h.getArray(); aux=AUX.getArray(); invVar=_ubm_invvar.getArray();

	//For each session
	for(unsigned long session=0; session<_n_sessions; session++){

		double *invl = _l_sess_inv[session].getArray();

		for(unsigned long i=0;i<_rankEC;i++){
			for(unsigned long k=0;k<_svSize;k++) {
				aux[session*_rankEC+i] += f_x_h[session*_svSize+k] * invVar[k] * u[i*_svSize+k];
			}
		}
	
		//multiplication by invL
		for(unsigned long i=0; i<_rankEC;i++){
			for(unsigned long k=0; k<_rankEC; k++){
				x[session*_rankEC+i] += aux[session*_rankEC+k] * invl[i*_rankEC+k];
			}
		}
	}
}

#ifdef THREAD
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Data strucutre of thread
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
struct estimateXthread_data{
	double *U;
	double *X;
	double *F_X_h;
	double *AUX;
	double *ubm_invvar;
	unsigned long sessBottom;
	unsigned long sessUp;	
	unsigned long rankEC;
	unsigned long svSize;
	unsigned long nt;
	RefVector <DoubleSquareMatrix>* linv;
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Thread Routine
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void *estimateXthread(void *threadarg) {
	
	struct estimateXthread_data *my_data;
	my_data = (struct estimateXthread_data *) threadarg;

	double *u = my_data->U;
	double *x = my_data->X;
	double *f_x_h = my_data->F_X_h;
	double *aux = my_data->AUX;
	double *ubm_invvar = my_data->ubm_invvar;
	unsigned long sessBottom = my_data->sessBottom;	
	unsigned long sessUp = my_data->sessUp;
	unsigned long _rankEC=my_data->rankEC;
	unsigned long _svSize=my_data->svSize;
	
	//For each session
	for(unsigned long session= sessBottom; session<sessUp; session++){

		DoubleSquareMatrix &linv=(*(my_data->linv))[session];
		double *invl = linv.getArray();

		for(unsigned long i=0;i<_rankEC;i++){
			for(unsigned long k=0;k<_svSize;k++) {
				aux[session*_rankEC+i] += f_x_h[session*_svSize+k] * ubm_invvar[k] * u[i*_svSize+k];
			}
		}
		
		//multiplication by invL
		for(unsigned long i=0; i<_rankEC;i++){
			for(unsigned long k=0; k<_rankEC; k++){
				x[session*_rankEC+i] += aux[session*_rankEC+k] * invl[i*_rankEC+k];
			}
		}
	}
	pthread_exit((void*) 0);
	return (void*)0 ;
}


//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::estimateXThreaded(unsigned long NUM_THREADS){

	if (verboseLevel >= 1) cout << "(AccumulateJFAStat) Estimate X for each session Threaded"<<endl;
	if (NUM_THREADS==0) throw Exception("Num threads can not be 0",__FILE__,__LINE__);

	int rc, status;
	if (NUM_THREADS > _n_sessions) NUM_THREADS=_n_sessions;

	struct estimateXthread_data *thread_data_array = new estimateXthread_data[NUM_THREADS];
	pthread_t *threads = new pthread_t[NUM_THREADS];

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	unsigned long offset=_n_sessions/NUM_THREADS;
	
	Matrix<double> AUX(_n_sessions,_rankEC);
	AUX.setAllValues(0.0);
	
	double *U =_matU.getArray();
	double *X =_matX.getArray();
	double *F_X_h =_F_X_h.getArray();
	double *aux = AUX.getArray(); 
	double *ubm_invvar=_ubm_invvar.getArray();

	unsigned long sessBottom = 0;
	unsigned long sessUp=0;
	unsigned long re=_n_sessions - NUM_THREADS*offset;
	
	//Create threads
	for(unsigned long t=0; t<NUM_THREADS; t++){
		sessUp = sessBottom +offset;
		if(t<re) sessUp +=1;

		thread_data_array[t].U=U;
		thread_data_array[t].X=X;
		thread_data_array[t].F_X_h=F_X_h;
		thread_data_array[t].AUX=aux;
		thread_data_array[t].ubm_invvar=ubm_invvar;
		thread_data_array[t].sessBottom=sessBottom;
		thread_data_array[t].sessUp=sessUp;
		thread_data_array[t].rankEC=_rankEC;
		thread_data_array[t].svSize=_svSize;
		thread_data_array[t].linv=&(_l_sess_inv);

		if (verboseLevel>1) cout<<"(AccumulateJFAStat) Creating thread n["<< t<< "] for speakers["<<sessBottom<<"-->"<<sessUp-1<<"]"<<endl;
		rc = pthread_create(&threads[t], &attr, estimateXthread, (void *)&thread_data_array[t]);
		if (rc) throw Exception("ERROR; return code from pthread_create() is ",__FILE__,rc);
		
		sessBottom = sessUp;
	}
	
	pthread_attr_destroy(&attr);
	for(unsigned long t=0; t<NUM_THREADS; t++) {
		rc = pthread_join(threads[t], (void **)&status);
		if (rc)  throw Exception("ERROR; return code from pthread_join() is ",__FILE__,rc);
		if (verboseLevel>1) cout <<"(AccumulateJFAStat) Completed join with thread ["<<t<<"] status["<<status<<"]"<<endl;
	}

	free(thread_data_array);
	free(threads);

	if (verboseLevel >= 1) cout << "(AccumulateJFAStat) Done " << endl;

}
#endif

//******************************************
// AJOUT POUR LE LFA VERSION DRISS
//******************************************
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::estimateU(){	//estimate X for all sessions and U
	
	if (verboseLevel>=1) cout << "(AccumulateJFAStat) Compute U Estimate "<<endl;
//	_X.setAllValues(0.0);
	Matrix<double> AUX(_n_sessions,_rankEC);
	AUX.setAllValues(0.0);
	
	double *x, *u, *f_x_h, /* *aux,*/ *invVar, *n_h, *c;
	x=_matX.getArray(); u=_matU.getArray(); f_x_h=_F_X_h.getArray();  invVar=_ubm_invvar.getArray(); n_h=_N_h.getArray(); c=_Cec.getArray();
	//aux=AUX.getArray();

	//For each session
	for(unsigned long session=0; session<_n_sessions; session++){

		double *invl = _l_sess_inv[session].getArray();

//		for(unsigned long i=0;i<_rankEC;i++){
//			for(unsigned long k=0;k<_svSize;k++) {
//				aux[session*_rankEC+i] += f_x_h[session*_svSize+k] * invVar[k] * u[i*_svSize+k];
//			}
//		}
		
//		//multiplication by invL
//		for(unsigned long i=0; i<_rankEC;i++){
//			for(unsigned long k=0; k<_rankEC; k++){
//				x[session*_rankEC+i] += aux[session*_rankEC+k] * invl[i*_rankEC+k];
//			}
//		}

		for(unsigned long i=0;i<_rankEC;i++){
			for(unsigned long j=0;j<_rankEC;j++){
					invl[i*_rankEC+j] += x[session*_rankEC+i]*x[session*_rankEC+j];
			}
		}

		for(unsigned long dis = 0; dis<_n_distrib;dis++){
			for(unsigned long i=0;i<_rankEC;i++){
				for(unsigned long j = 0;j<_rankEC;j++){
					_Aec[dis](i,j) += invl[i*_rankEC+j] * n_h[session*_n_distrib+dis];

				}
			}
		}

		for(unsigned long i=0;i<_rankEC;i++){
			for(unsigned long j=0;j<_svSize;j++){
					c[i*_svSize+j] += x[session*_rankEC+i] * f_x_h[session*_svSize+j];
			}
		}
	}
}
//******************************************
// FIN D'AJOUT POUR LE LFA VERSION DRISS
//******************************************

//----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::estimateZandD(){	//estimate Z and D
	
	if (verboseLevel>=1) cout << "(AccumulateJFAStat) Compute Z and D Estimate "<<endl;
	_Z.setAllValues(0.0);
	
	DoubleVector aux1(_svSize, _svSize), aux2(_svSize, _svSize);
	aux1.setAllValues(0.0); aux2.setAllValues(0.0);
	
	double *f_x, *n, *z;
	f_x=_F_X.getArray(); n=_matN.getArray(); z=_Z.getArray();

	for(unsigned long spk=0; spk<_n_speakers;spk++){
		
		DoubleVector L(_svSize, _svSize);
		L.setAllValues(0.0);
		
		for(unsigned long i=0; i<_n_distrib; i++){
			for(unsigned long j=0; j<_vectSize; j++){
				L[i*_vectSize+j] = 1 + _matN(spk,i) * _ubm_invvar[i*_vectSize+j] * _D[i*_vectSize+j] * _D[i*_vectSize+j];
				_Z(spk,i*_vectSize+j) = _F_X(spk,i*_vectSize+j) * _ubm_invvar[i*_vectSize+j] * _D[i*_vectSize+j] / L[i*_vectSize+j];
			}
		}

		for(unsigned long i=0;i<_n_distrib;i++){
			for(unsigned long j=0;j<_vectSize;j++){
				aux1[i * _vectSize +j] += ( 1 / L[i*_vectSize + j] + _Z(spk,i*_vectSize +j) * _Z(spk,i*_vectSize +j) ) * _matN(spk, i);
				aux2[i*_vectSize+j] += _Z(spk,i*_vectSize +j) * _F_X(spk,i*_vectSize+j);
				
			}
		}
	}
	
	for(unsigned long i=0;i<_svSize;i++){
		_D[i] = aux2[i] / aux1[i];
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::estimateYX(){	//estimate X and Y at the same time for all speakers

	if (verboseLevel>=1) cout << "(AccumulateJFAStat) Compute Y and X  Estimate "<<endl;

	_YX.setAllValues(0.0);
	Matrix<double> AUX(_n_speakers,_rankEV + _rankEC);
	AUX.setAllValues(0.0);
	
	double *yx, *vu, *f_x, *aux, *invVar;
	yx=_YX.getArray(); vu=_VU.getArray(); f_x=_F_X.getArray(); aux=AUX.getArray(); invVar=_ubm_invvar.getArray();

	//For each speaker
	for(unsigned long spk=0; spk<_n_speakers; spk++){

		double *invl = _l_yx_inv[spk].getArray();

		for(unsigned long i=0;i<_rankEV + _rankEC ;i++){
			for(unsigned long k=0;k<_svSize;k++) {
				aux[spk*(_rankEV + _rankEC)+i] += f_x[spk*_svSize+k] * invVar[k] * vu[i*_svSize+k];
			}
		}
		
		//multiplication by invL
		for(unsigned long i=0; i<_rankEV + _rankEC;i++){
			for(unsigned long k=0; k<_rankEV + _rankEC; k++){
				yx[spk*(_rankEV + _rankEC)+i] += aux[spk*(_rankEV + _rankEC)+k] * invl[i*(_rankEV + _rankEC)+k];
			}
		}
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::estimateZ(){	//estimate Z for all speakers during the training phase
	
	if (verboseLevel>=1) cout << "(AccumulateJFAStat) Compute Z Estimate during the training phase "<<endl;
	_Z.setAllValues(0.0);
	
	DoubleVector aux1(_svSize, _svSize), aux2(_svSize, _svSize);
	aux1.setAllValues(0.0); aux2.setAllValues(0.0);
	DoubleVector L(_svSize, _svSize);

	double *z, *n, *f_x;
	z=_Z.getArray(); n=_matN.getArray(); f_x=_F_X.getArray();
	
	for(unsigned long spk=0; spk<_n_speakers;spk++){

		L.setAllValues(0.0);
		
		for(unsigned long i=0; i<_n_distrib; i++){
			for(unsigned long j=0; j<_vectSize; j++){
				L[i*_vectSize+j] = 1 + _matN(spk,i) * _ubm_invvar[i*_vectSize+j] * _D[i*_vectSize+j] * _D[i*_vectSize+j];
				_Z(spk,i*_vectSize+j) = _F_X(spk,i*_vectSize+j) * _ubm_invvar[i*_vectSize+j] * _D[i*_vectSize+j] / L[i*_vectSize+j];
			}
		}
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::estimateZMAP(double tau){	//estimate Z for all speakers during the training phase
	
	if (verboseLevel>=1) cout << "(AccumulateJFAStat) Compute Z MAP estimate during the training phase "<<endl;
	_Z.setAllValues(0.0);
	
	DoubleVector aux1(_svSize, _svSize), aux2(_svSize, _svSize);
	aux1.setAllValues(0.0); aux2.setAllValues(0.0);

	double *z, *n, *f_x;
	z=_Z.getArray(); n=_matN.getArray(); f_x=_F_X.getArray();
	
	for(unsigned long spk=0; spk<_n_speakers;spk++){
		for(unsigned long i=0; i<_n_distrib; i++){
			for(unsigned long j=0; j<_vectSize; j++){
				_Z(spk,i*_vectSize+j) = (tau/(tau + _matN(spk,i))) * _D[i*_vectSize+j] * _ubm_invvar[i*_vectSize+j] * _F_X(spk,i*_vectSize+j);
			}
		}
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::updateVestimate(){

	Matrix<double> tmpC; tmpC.setDimensions(_rankEV, _svSize); tmpC.setAllValues(0.0);
	
	for(unsigned long d=0;d<_n_distrib;d++){
		DoubleSquareMatrix invA;
		invA.setSize(_rankEV);
		_Aev[d].invert(invA);
	
		double *tmpc, * inva, *cev;
		tmpc=tmpC.getArray(); inva=invA.getArray(); cev=_Cev.getArray();
	
		for(unsigned long i=0; i<_rankEV;i++){
			for(unsigned long j=0; j<_vectSize;j++){
				for(unsigned long k=0; k<_rankEV;k++){
					tmpC(i,d*_vectSize+j) += invA(i,k) *  _Cev(k,d*_vectSize+j);
				}
			}
		}
	}
	_Cev=tmpC;
	_V=_Cev;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::updateUestimate(){

	Matrix<double> tmpC; tmpC.setDimensions(_rankEC, _svSize); tmpC.setAllValues(0.0);
	for(unsigned long d=0;d<_n_distrib;d++){
		DoubleSquareMatrix invA;
		invA.setSize(_rankEC);
		_Aec[d].invert(invA);

		double *tmpc, * inva, *cec;
		tmpc=tmpC.getArray(); inva=invA.getArray(); cec=_Cec.getArray();
		
		for(unsigned long i=0; i<_rankEC;i++){
			for(unsigned long j=0; j<_vectSize;j++){
				for(unsigned long k=0; k<_rankEC;k++){
					tmpC(i,d*_vectSize+j) += invA(i,k) *  _Cec(k,d*_vectSize+j);
				}
			}
		}
	}
	_Cec=tmpC;
	_matU=_Cec;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
XList& JFAAcc::getXList(){
	return(_fileList);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned long JFAAcc::getNSpeakers(){
	return(_n_speakers);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned long JFAAcc::getNSessions(){
	return(_n_sessions);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned long JFAAcc::getNDistrib(){
	return(_n_distrib);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned long JFAAcc::getVectSize(){
	return(_vectSize);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned long JFAAcc::getSvSize(){
	return(_svSize);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned long JFAAcc::getRankEV(){
	return(_rankEV);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned long JFAAcc::getRankEC(){
	return(_rankEC);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
DoubleVector& JFAAcc::getUbmMeans(){
	return(_ubm_means);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
DoubleVector& JFAAcc::getUbmInvVar(){
	return(_ubm_invvar);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
Matrix<double> JFAAcc::getV(){
	return(_V);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
Matrix<double> JFAAcc::getU(){
	return(_matU);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
DoubleVector JFAAcc::getD(){
	return(_D);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
Matrix<double> JFAAcc::getVU(){
	return(_VU);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
Matrix<double> JFAAcc::getY(){
	return(_Y);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
Matrix<double> JFAAcc::getX(){
	return(_matX);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
Matrix<double> JFAAcc::getZ(){
	return(_Z);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
Matrix<double> JFAAcc::getYX(){
	return(_YX);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::saveX(String xFile,Config & config){
	_matX.save(xFile,config);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::saveY(String yFile,Config &config){
	_Y.save(yFile,config);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::saveZ(String xFile,Config &config){
	_Z.save(xFile,config);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
Matrix <double>& JFAAcc::getN() {
	return _matN;
}	
	
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
Matrix <double>& JFAAcc::getN_h() {
	return _N_h;
}
	
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
Matrix <double>& JFAAcc::getF_X() {
	return _F_X;
}	
	
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
Matrix <double>& JFAAcc::getF_X_h() {
	return _F_X_h;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::splitYX(){
	_Y= _YX.crop(0,0, _n_speakers, _rankEV);
	_matX= _YX.crop(0,_rankEV, _n_speakers, _rankEC);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::storeAccs(){ // save JFA state in temporary variables
		_cF_X=_F_X;
		_cF_X_h=_F_X_h;
		_cN_h=_N_h; 
		_cN=_matN;
		if (verboseLevel>=1) cout << "(AccumulateJFAStat) JFA Accs states stored" << endl;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::restoreAccs() {	// Restore Accumulators in temporary variables
	_matN=_cN;
	_N_h=_cN_h;	
	_F_X=_cF_X;
	_F_X_h=_cF_X_h;
	if (verboseLevel>=1) cout << "(AccumulateJFAStat) JFA Accs states restored" << endl;			
}
	
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::substractMplusDZ(Config & config){
	#ifdef THREAD          
	if (config.existsParam("numThread") && config.getParam("numThread").toULong() >0)	substractMplusDZThreaded(config.getParam("numThread").toULong());
	else substractMplusDZUnThreaded();
	#else
	substractMplusDZUnThreaded();
	#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::substractMplusDZUnThreaded(){
	
	if (verboseLevel >= 1) cout <<"(AccumulateJFAStat) Compute and Substract Speaker FA Stats M + DZ " << endl;	

	double *F_X = _F_X.getArray();
	
	for(unsigned long spk=0; spk<_n_speakers; spk++){
		//Compute M+DZ for the current speaker
		DoubleVector MplusDZ; MplusDZ.setSize(_svSize); MplusDZ.setAllValues(0.0);
		double *mplusdz=MplusDZ.getArray();
		this->getMplusDZ(MplusDZ,spk);

		double *n=_matN.getArray();
		
		//Substract M+DZ
		for(unsigned long i=0; i<_n_distrib; i++){
			for(unsigned long j = 0; j< _vectSize;j++){
				F_X[spk*_svSize+i*_vectSize+j] -= mplusdz[i*_vectSize+j]*n[spk*_n_distrib+i];
			}
		}
	}
}


#ifdef THREAD
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Data strucutre of thread
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
struct S_MplusDZthread_data{
	double *N;
	double *F_X;
	double *D;
	double *Z;
	double *ubm_means;
	unsigned long spkBottom;
	unsigned long spkUp;	
	unsigned long vectSize;
	unsigned long svSize;
	unsigned long n_distrib;
	unsigned long nt;
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Thread Routine
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void *S_MplusDZthread(void *threadarg) {
	struct S_MplusDZthread_data *my_data;
	my_data = (struct S_MplusDZthread_data *) threadarg;

	double *n = my_data->N;
	double *f_x = my_data->F_X;
	double *d = my_data->D;
	double *z = my_data->Z;
	double *ubm_means = my_data->ubm_means;
	unsigned long spkBottom = my_data->spkBottom;	
	unsigned long spkUp = my_data->spkUp;
	unsigned long _vectSize=my_data->vectSize;
	unsigned long _svSize=my_data->svSize;
	unsigned long _n_distrib=my_data->n_distrib;

	for(unsigned long spk=spkBottom; spk<spkUp; spk++){
		//Compute M+DZ for the current speaker
		DoubleVector MplusDZ; MplusDZ.setSize(_svSize); MplusDZ.setAllValues(0.0);
		for (unsigned long i=0;i<_svSize;i++){
			MplusDZ[i]=ubm_means[i]+d[i] * z[spk*_svSize+i];
		}

		//Substract M+DZ
		for(unsigned long i=0; i<_n_distrib; i++){
			for(unsigned long j = 0; j< _vectSize;j++){
				f_x[spk*_svSize+i*_vectSize+j] -= MplusDZ[i*_vectSize+j]*n[spk*_n_distrib+i];
			}
		}
	}
	pthread_exit((void*) 0);
	return (void*)0 ;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::substractMplusDZThreaded(unsigned long NUM_THREADS){
	
	if (verboseLevel >= 1) cout <<"(AccumulateJFAStat) Compute and Substract Speaker FA Stats M + DZ Threaded" << endl;
	if (NUM_THREADS==0) throw Exception("Num threads can not be 0",__FILE__,__LINE__);

	int rc, status;
	if (NUM_THREADS > _n_speakers) NUM_THREADS=_n_speakers;

	struct S_MplusDZthread_data *thread_data_array = new S_MplusDZthread_data[NUM_THREADS];
	pthread_t *threads = new pthread_t[NUM_THREADS];

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	unsigned long offset=_n_speakers/NUM_THREADS;
	
	double *N=_matN.getArray(); 
	double *F_X=_F_X.getArray();
	double *D=_D.getArray();
	double *Z=_Z.getArray();
	double *ubm_means=_ubm_means.getArray();
	unsigned long spkBottom = 0;
	unsigned long spkUp=0;
	unsigned long re=_n_speakers - NUM_THREADS*offset;
	
	//Create threads
	for(unsigned long t=0; t<NUM_THREADS; t++){
		spkUp = spkBottom +offset;
		if(t<re) spkUp +=1;

		thread_data_array[t].N=N;
		thread_data_array[t].F_X=F_X;
		thread_data_array[t].D=D;
		thread_data_array[t].Z=Z;
		thread_data_array[t].ubm_means=ubm_means;
		thread_data_array[t].spkBottom=spkBottom;
		thread_data_array[t].spkUp=spkUp;
		thread_data_array[t].vectSize=_vectSize;
		thread_data_array[t].svSize=_svSize;
		thread_data_array[t].n_distrib=_n_distrib;

		if (verboseLevel>1) cout<<"(AccumulateJFAStat) Creating thread n["<< t<< "] for speakers ["<<spkBottom<<"-->"<<spkUp-1<<"]"<<endl;
		rc = pthread_create(&threads[t], &attr, S_MplusDZthread, (void *)&thread_data_array[t]);
		if (rc) throw Exception("ERROR; return code from pthread_create() is ",__FILE__,rc);
		
		spkBottom = spkUp;
	}
	
	pthread_attr_destroy(&attr);
	for(unsigned long t=0; t<NUM_THREADS; t++) {
		rc = pthread_join(threads[t], (void **)&status);
		if (rc)  throw Exception("ERROR; return code from pthread_join() is ",__FILE__,rc);
		if (verboseLevel>1) cout <<"(AccumulateJFAStat) Completed join with thread ["<<t<<"] status["<<status<<"]"<<endl;
	}

	free(thread_data_array);
	free(threads);

	if (verboseLevel >= 1) cout << "(AccumulateJFAStat) Done " << endl;
}

#endif

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::substractMplusDZByChannel(){
	
	if (verboseLevel >= 1) cout <<"(AccumulateJFAStat) Compute and Substract Speaker FA Stats M + DZ " << endl;	

	double *F_X_h = _F_X_h.getArray();
	double *N_h = _N_h.getArray();

	XLine *pline; String *pFile; _fileList.rewind();

	unsigned long spk =0;
	unsigned long sess=0;

	while((pline=_fileList.getLine())!=NULL) {

		DoubleVector MplusDZ; MplusDZ.setSize(_svSize); MplusDZ.setAllValues(0.0);

		while((pFile=pline->getElement())!=NULL) {

		this->getMplusDZ(MplusDZ,spk);

			for(unsigned long k=0;k<_n_distrib;k++) 
				for (unsigned long i=0;i<_vectSize;i++) 
					F_X_h[sess*_svSize+(k*_vectSize+i)]-= N_h[sess*_n_distrib+k]*MplusDZ[i+k*_vectSize];
			sess++;
		}
		spk++;
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::substractMplusVY(Config &config){
	#ifdef THREAD          
	if (config.existsParam("numThread") && config.getParam("numThread").toULong() >0)	substractMplusVYThreaded(config.getParam("numThread").toULong());
	else substractMplusVYUnThreaded();
	#else
	substractMplusVYUnThreaded();
	#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::substractMplusVYUnThreaded(){
	
	if (verboseLevel >= 1) cout <<"(AccumulateJFAStat) Compute and Substract Speaker FA Stats M + VY" << endl;	

	double *f_x = _F_X.getArray();
	
	for(unsigned long spk=0; spk<_n_speakers; spk++){
		//compute M+VY for the current speaker
		DoubleVector MplusVY; MplusVY.setSize(_svSize); MplusVY.setAllValues(0.0);
		double *mplusvy=MplusVY.getArray();
		this->getMplusVY(MplusVY,spk);
		double *n=_matN.getArray();

		//substract M+VY
		for(unsigned long i=0; i<_n_distrib; i++){
			for(unsigned long j = 0; j< _vectSize;j++){
				f_x[spk*_svSize+i*_vectSize+j] -= mplusvy[i*_vectSize+j]*n[spk*_n_distrib+i];
			}
		}
	}
}

#ifdef THREAD
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Data strucutre of thread
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
struct S_MplusVYthread_data{
	double *N;
	double *F_X;
	double *V;
	double *Y;
	double *ubm_means;
	unsigned long spkBottom;
	unsigned long spkUp;	
	unsigned long vectSize;
	unsigned long svSize;
	unsigned long n_distrib;
	unsigned long rankEV;
	unsigned long nt;
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Thread Routine
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void *S_MplusVYthread(void *threadarg) {
	
	struct S_MplusVYthread_data *my_data;
	my_data = (struct S_MplusVYthread_data *) threadarg;

	double *n = my_data->N;
	double *f_x = my_data->F_X;
	double *v = my_data->V;
	double *y = my_data->Y;
	double *ubm_means = my_data->ubm_means;
	unsigned long spkBottom = my_data->spkBottom;	
	unsigned long spkUp = my_data->spkUp;
	unsigned long _vectSize=my_data->vectSize;
	unsigned long _svSize=my_data->svSize;
	unsigned long _n_distrib=my_data->n_distrib;
	unsigned long _rankEV=my_data->rankEV;

	DoubleVector MplusVY;MplusVY.setSize(_svSize);
	double* mplusvy=MplusVY.getArray();	

	for(unsigned long spk=spkBottom; spk<spkUp;spk++){

		//compute M+VY for the current speaker
		MplusVY.setAllValues(0.0);

		//Calcul de mplusvy
		for (unsigned long i=0;i<_svSize;i++){
			for (unsigned long j=0;j<_rankEV;j++){
				mplusvy[i] += v[j*_svSize + i] * y[spk*_rankEV + j ];
			}
			//Ajout de M
			mplusvy[i] += ubm_means[i];
		}
		//substract M+VY
		for(unsigned long i=0; i<_n_distrib; i++){
			for(unsigned long j = 0; j< _vectSize;j++){
				f_x[spk*_svSize+i*_vectSize+j] -= mplusvy[i*_vectSize+j]*n[spk*_n_distrib+i];
			}
		}
	}
	pthread_exit((void*) 0);
	return (void*)0 ;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::substractMplusVYThreaded(unsigned long NUM_THREADS){
	
	if (verboseLevel >= 1) cout <<"(AccumulateJFAStat) Compute and Substract Speaker FA Stats M + VY Threaded" << endl;
	if (NUM_THREADS==0) throw Exception("Num threads can not be 0",__FILE__,__LINE__);

	int rc, status;
	if (NUM_THREADS > _n_speakers) NUM_THREADS=_n_speakers;

	struct S_MplusVYthread_data *thread_data_array = new S_MplusVYthread_data[NUM_THREADS];
	pthread_t *threads = new pthread_t[NUM_THREADS];

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	unsigned long offset=_n_speakers/NUM_THREADS;
	
	double *N=_matN.getArray(); 
	double *F_X=_F_X.getArray();
	double *V=_V.getArray();
	double *Y=_Y.getArray();
	double *ubm_means = _ubm_means.getArray();
	unsigned long spkBottom = 0;
	unsigned long spkUp=0;
	unsigned long re=_n_speakers - NUM_THREADS*offset;
	
	//Create threads
	for(unsigned long t=0; t<NUM_THREADS; t++){
		spkUp = spkBottom +offset;
		if(t<re) spkUp +=1;

		thread_data_array[t].N=N;
		thread_data_array[t].F_X=F_X;
		thread_data_array[t].V=V;
		thread_data_array[t].Y=Y;
		thread_data_array[t].ubm_means=ubm_means;
		thread_data_array[t].spkBottom=spkBottom;
		thread_data_array[t].spkUp=spkUp;
		thread_data_array[t].vectSize=_vectSize;
		thread_data_array[t].svSize=_svSize;
		thread_data_array[t].n_distrib=_n_distrib;
		thread_data_array[t].rankEV=_rankEV;

		if (verboseLevel>1) cout<<"(AccumulateJFAStat) Creating thread n["<< t<< "] for speakers ["<<spkBottom<<"-->"<<spkUp-1<<"]"<<endl;
		rc = pthread_create(&threads[t], &attr, S_MplusVYthread, (void *)&thread_data_array[t]);
		if (rc) throw Exception("ERROR; return code from pthread_create() is ",__FILE__,rc);
		
		spkBottom = spkUp;
	}
	
	pthread_attr_destroy(&attr);
	for(unsigned long t=0; t<NUM_THREADS; t++) {
		rc = pthread_join(threads[t], (void **)&status);
		if (rc)  throw Exception("ERROR; return code from pthread_join() is ",__FILE__,rc);
		if (verboseLevel>1) cout <<"(AccumulateJFAStat) Completed join with thread ["<<t<<"] status["<<status<<"]"<<endl;
	}

	free(thread_data_array);
	free(threads);

	if (verboseLevel >= 1) cout << "(AccumulateJFAStat) Done " << endl;
}

#endif

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::substractUX(Config &config){
	#ifdef THREAD          
	if (config.existsParam("numThread") && config.getParam("numThread").toULong() >0)	substractUXThreaded(config.getParam("numThread").toULong());
	else substractUXUnThreaded();
	#else
	substractUXUnThreaded();
	#endif
}
	
//-----------------------------------------------------------------------------------------------------------------------------------------------------------	
void JFAAcc::substractUXUnThreaded(){

	if (verboseLevel >= 1) cout <<"(AccumulateJFAStat) Compute and Substract Speaker FA Stats UX " << endl;

	unsigned long spk=0;
	unsigned long session=0; 
	XLine *pline; String *pFile; 
	_fileList.rewind();

	DoubleVector UX;UX.setSize(_svSize); UX.setAllValues(0.0);
	double* ux=UX.getArray();	
	double *f_x = _F_X.getArray();
	double *n_h=_N_h.getArray(); 

	while((pline=_fileList.getLine())!=NULL) { 		
		while((pFile=pline->getElement())!=NULL) {	
			this->getUX(UX,*pFile);

			for(unsigned long k=0;k<_n_distrib;k++)
				for (unsigned long i=0;i<_vectSize;i++)
					f_x[spk*_svSize+k*_vectSize+i] -= n_h[session*_n_distrib+k]*ux[i+k*_vectSize];
			session++;
		}
		spk++;
	}
}

#ifdef THREAD
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Data strucutre of thread
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
struct S_UXthread_data{
	double *N_h;
	double *F_X;
	double *U;
	double *X;
	RefVector<DoubleVector>* sessIdx;
	unsigned long spkBottom;
	unsigned long spkUp;	
	unsigned long vectSize;
	unsigned long svSize;
	unsigned long n_distrib;
	unsigned long rankEC;
	unsigned long nt;
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Thread Routine
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void *S_UXthread(void *threadarg) {
	struct S_UXthread_data *my_data;
	my_data = (struct S_UXthread_data *) threadarg;

	double *n_h = my_data->N_h;
	double *f_x = my_data->F_X;
	double *u = my_data->U;
	double *x = my_data->X;
	unsigned long spkBottom = my_data->spkBottom;	
	unsigned long spkUp = my_data->spkUp;
	unsigned long _vectSize=my_data->vectSize;
	unsigned long _svSize=my_data->svSize;
	unsigned long _n_distrib=my_data->n_distrib;
	unsigned long _rankEC=my_data->rankEC;


	DoubleVector UX;UX.setSize(_svSize);
	double* ux=UX.getArray();	

	for(unsigned long spk=spkBottom; spk<spkUp;spk++){

		DoubleVector &sessIdx =(*(my_data->sessIdx))[spk];

		for(unsigned long sess=0; sess<sessIdx.size(); sess++){

			//getUX
			UX.setAllValues(0.0);
			unsigned long idx = (unsigned long)sessIdx[sess];

			for (unsigned long i=0;i<_svSize;i++){
				for (unsigned long j=0;j<_rankEC;j++) {
					ux[i]+=u[j*_svSize+i]*x[idx*_rankEC+j];
				}
			}
			//
			for(unsigned long k=0;k<_n_distrib;k++){
				for (unsigned long i=0;i<_vectSize;i++){
					f_x[spk*_svSize+k*_vectSize+i] -= n_h[idx*_n_distrib+k]*ux[i+k*_vectSize];
				}
			}
		}
	}
	pthread_exit((void*) 0);
	return (void*)0 ;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::substractUXThreaded(unsigned long NUM_THREADS){
	
	if (verboseLevel >= 1) cout <<"(AccumulateJFAStat) Compute and Substract Channel JFA Stats UX Threaded" << endl;
	if (NUM_THREADS==0) throw Exception("Num threads can not be 0",__FILE__,__LINE__);

	int rc, status;
	if (NUM_THREADS > _n_speakers) NUM_THREADS=_n_speakers;

	struct S_UXthread_data *thread_data_array = new S_UXthread_data[NUM_THREADS];
	pthread_t *threads = new pthread_t[NUM_THREADS];

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	unsigned long offset=_n_speakers/NUM_THREADS;
	
	double *N_h=_N_h.getArray(); 
	double *F_X=_F_X.getArray();
	double *U=_matU.getArray();
	double *X=_matX.getArray();

	RefVector<DoubleVector> sessIdx(_n_speakers);
	DoubleVector tmpSessIdx(0);

	for(unsigned long i=0; i<_n_speakers; i++){
		sessIdx.addObject(*new DoubleVector(tmpSessIdx));
	}

	_fileList.rewind(); XLine *pline; String *pfile;
	unsigned long tmpS=0;
	unsigned long tmpSpk =0;
	while((pline=_fileList.getLine()) != NULL){
		while((pfile=pline->getElement())!=NULL){
			sessIdx[tmpSpk].addValue(tmpS);
			tmpS++;
		}
		tmpSpk++;
	}

	unsigned long spkBottom = 0;
	unsigned long spkUp=0;
	unsigned long re=_n_speakers - NUM_THREADS*offset;
	
	//Create threads
	for(unsigned long t=0; t<NUM_THREADS; t++){
		spkUp = spkBottom +offset;
		if(t<re) spkUp +=1;

		thread_data_array[t].N_h=N_h;
		thread_data_array[t].F_X=F_X;
		thread_data_array[t].U=U;
		thread_data_array[t].X=X;
		thread_data_array[t].sessIdx=&(sessIdx);


		thread_data_array[t].spkBottom=spkBottom;
		thread_data_array[t].spkUp=spkUp;
		thread_data_array[t].vectSize=_vectSize;
		thread_data_array[t].svSize=_svSize;
		thread_data_array[t].n_distrib=_n_distrib;
		thread_data_array[t].rankEC=_rankEC;

		if (verboseLevel>1) cout<<"(AccumulateJFAStat) Creating thread n["<< t<< "] for speakers ["<<spkBottom<<"-->"<<spkUp-1<<"]"<<endl;
		rc = pthread_create(&threads[t], &attr, S_UXthread, (void *)&thread_data_array[t]);
		if (rc) throw Exception("ERROR; return code from pthread_create() is ",__FILE__,rc);
		
		spkBottom = spkUp;
	}


	pthread_attr_destroy(&attr);
	for(unsigned long t=0; t<NUM_THREADS; t++) {
		rc = pthread_join(threads[t], (void **)&status);
		if (rc)  throw Exception("ERROR; return code from pthread_join() is ",__FILE__,rc);
		if (verboseLevel>1) cout <<"(AccumulateJFAStat) Completed join with thread ["<<t<<"] status["<<status<<"]"<<endl;
	}

	free(thread_data_array);
	free(threads);

	sessIdx.deleteAllObjects();

	if (verboseLevel >= 1) cout << "(AccumulateJFAStat) Done " << endl;
}

#endif

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::substractMplusUX(){

	if (verboseLevel >= 1) cout <<"(AccumulateJFAStat) Compute and Substract Speaker FA Stats UX " << endl;

	unsigned long spk=0;
	unsigned long session=0; 
	XLine *pline; String *pFile; 
	_fileList.rewind();

	DoubleVector MplusUX;MplusUX.setSize(_svSize); MplusUX.setAllValues(0.0);
	double* mplusux=MplusUX.getArray();	
	double *f_x = _F_X.getArray();
	double *n_h=_N_h.getArray(); 

	while((pline=_fileList.getLine())!=NULL) { 		
		while((pFile=pline->getElement())!=NULL) {	
			this->getMplusUX(MplusUX,*pFile);			

			for(unsigned long k=0;k<_n_distrib;k++){
				for (unsigned long i=0;i<_vectSize;i++){
					f_x[spk*_svSize+k*_vectSize+i] -= n_h[session*_n_distrib+k] *  mplusux[i+k*_vectSize];
				}
			}
			session++;
		}
		spk++;
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::substractMplusVUYX(){
	
	if (verboseLevel >= 1) cout <<"(AccumulateJFAStat) Compute and Substract Speaker FA Stats M + VUYX" << endl;	

	double *f_x = _F_X.getArray();

	for(unsigned long spk=0; spk<_n_speakers; spk++){
		//compute M+VUYX for the current speaker
		DoubleVector MplusVUYX; MplusVUYX.setSize(_svSize); MplusVUYX.setAllValues(0.0);
		double *mplusvuyx=MplusVUYX.getArray();
		this->getMplusVUYX(MplusVUYX,spk);

		double *n=_matN.getArray();
		
		//substract M+VUYX
		for(unsigned long i=0; i<_n_distrib; i++){
			for(unsigned long j = 0; j< _vectSize;j++){
				f_x[spk*_svSize+i*_vectSize+j] -= mplusvuyx[i*_vectSize+j]*n[spk*_n_distrib+i];
			}
		}
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::substractMplusVYplusDZ(Config &config){
	#ifdef THREAD          
	if (config.existsParam("numThread") && config.getParam("numThread").toULong() >0)	substractMplusVYplusDZThreaded(config.getParam("numThread").toULong());
	else substractMplusVYplusDZUnThreaded();
	#else
	substractMplusVYplusDZUnThreaded();
	#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::substractMplusVYplusDZUnThreaded(){
	
	if (verboseLevel >= 1) cout <<"(AccumulateJFAStat) Compute and Substract Speaker FA Stats M + VY + DZ unthreaded" << endl;	

	DoubleVector MplusVYplusDZ(_svSize,_svSize);
	MplusVYplusDZ.setAllValues(0.0); 
	double *m_vy_dz=MplusVYplusDZ.getArray();  
	
	// Compute Occupations and Statistics
	unsigned long spk=0; 
	unsigned long session=0; 		
	XLine *pline; String *pFile; _fileList.rewind();

	double *n_h=_N_h.getArray(); 
	double *f_x_h=_F_X_h.getArray();
	while((pline=_fileList.getLine())!=NULL) {
		this->getMplusVYplusDZ(MplusVYplusDZ,spk);
		
		while((pFile=pline->getElement())!=NULL) {
			for(unsigned long k=0;k<_n_distrib;k++){
				for (unsigned long i=0;i<_vectSize;i++){ 
					f_x_h[session*_svSize+k*_vectSize+i] -= n_h[session*_n_distrib+k]*m_vy_dz[i+k*_vectSize];
				}
			}			
			session++;
		}
		spk++;
	}
}

#ifdef THREAD
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Data strucutre of thread
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
struct S_MplusVYplusDZthread_data{
	double *N_h;
	double *F_X_h;
	double *V;
	double *Y;
	double *D;
	double *Z;
	double *ubm_means;
	RefVector<DoubleVector>* sessIdx;

	unsigned long spkBottom;
	unsigned long spkUp;	
	unsigned long vectSize;
	unsigned long svSize;
	unsigned long n_distrib;
	unsigned long rankEV;
	unsigned long nt;
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Thread Routine
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void *S_MplusVYplusDZthread(void *threadarg) {
	struct S_MplusVYplusDZthread_data *my_data;
	my_data = (struct S_MplusVYplusDZthread_data *) threadarg;

	double *n_h = my_data->N_h;
	double *f_x_h = my_data->F_X_h;
	double *v = my_data->V;
	double *y = my_data->Y;
	double *d = my_data->D;
	double *z = my_data->Z;
	double *ubm_means = my_data->ubm_means;
	unsigned long spkBottom = my_data->spkBottom;	
	unsigned long spkUp = my_data->spkUp;
	unsigned long _vectSize=my_data->vectSize;
	unsigned long _svSize=my_data->svSize;
	unsigned long _n_distrib=my_data->n_distrib;
	unsigned long _rankEV=my_data->rankEV;

	DoubleVector Sp;Sp.setSize(_svSize);
	double* sp=Sp.getArray();
	
	for(unsigned long spk=spkBottom; spk<spkUp;spk++){


		DoubleVector &sessIdx =(*(my_data->sessIdx))[spk];

		for(unsigned long sess=0; sess<sessIdx.size(); sess++){

			unsigned long idx = (unsigned long)sessIdx[sess];

			Sp.setAllValues(0.0);
			
			//getMplusVYplusDZ
			for (unsigned long i=0;i<_svSize;i++){

				//getVY
				for (unsigned long j=0;j<_rankEV;j++){
					sp[i] += v[j*_svSize + i] * y[spk*_rankEV + j ];
				}

				//add DZ
				sp[i] +=ubm_means[i] + d[i]*z[spk*_svSize+i];
			}

			for(unsigned long k=0;k<_n_distrib;k++){
				for (unsigned long i=0;i<_vectSize;i++){ 
					f_x_h[idx*_svSize+k*_vectSize+i] -= n_h[idx*_n_distrib+k]*sp[i+k*_vectSize];
				}
			}
		}
	}
	pthread_exit((void*) 0);
	return (void*)0 ;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::substractMplusVYplusDZThreaded(unsigned long NUM_THREADS){
	
	if (verboseLevel >= 1) cout <<"(AccumulateJFAStat) Compute and Substract Speaker FA Stats M + VY + DZ Threaded" << endl;
	if (NUM_THREADS==0) throw Exception("Num threads can not be 0",__FILE__,__LINE__);

	int rc, status;
	if (NUM_THREADS > _n_speakers) NUM_THREADS=_n_speakers;

	struct S_MplusVYplusDZthread_data *thread_data_array = new S_MplusVYplusDZthread_data[NUM_THREADS];
	pthread_t *threads = new pthread_t[NUM_THREADS];

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	unsigned long offset=_n_speakers/NUM_THREADS;
	
	double *N_h=_N_h.getArray(); 
	double *F_X_h=_F_X_h.getArray();
	double *V=_V.getArray();
	double *Y=_Y.getArray();
	double *D=_D.getArray();
	double *Z=_Z.getArray();
	double *ubm_means=_ubm_means.getArray();

	RefVector<DoubleVector> sessIdx(_n_speakers);
	DoubleVector tmpSessIdx(0);

	for(unsigned long i=0; i<_n_speakers; i++){
		sessIdx.addObject(*new DoubleVector(tmpSessIdx));
	}

	_fileList.rewind(); XLine *pline; String *pfile;
	unsigned long tmpS=0;
	unsigned long tmpSpk =0;
	while((pline=_fileList.getLine()) != NULL){
		while((pfile=pline->getElement())!=NULL){
			sessIdx[tmpSpk].addValue(tmpS);
			tmpS++;
		}
		tmpSpk++;
	}

	unsigned long spkBottom = 0;
	unsigned long spkUp=0;
	unsigned long re=_n_speakers - NUM_THREADS*offset;
	
	//Create threads
	for(unsigned long t=0; t<NUM_THREADS; t++){
		spkUp = spkBottom +offset;
		if(t<re) spkUp +=1;

		thread_data_array[t].N_h=N_h;
		thread_data_array[t].F_X_h=F_X_h;
		thread_data_array[t].V=V;
		thread_data_array[t].Y=Y;
		thread_data_array[t].D=D;
		thread_data_array[t].Z=Z;
		thread_data_array[t].ubm_means=ubm_means;
		thread_data_array[t].sessIdx=&(sessIdx);

		thread_data_array[t].spkBottom=spkBottom;
		thread_data_array[t].spkUp=spkUp;
		thread_data_array[t].vectSize=_vectSize;
		thread_data_array[t].svSize=_svSize;
		thread_data_array[t].n_distrib=_n_distrib;
		thread_data_array[t].rankEV=_rankEV;

		if (verboseLevel>1) cout<<"(AccumulateJFAStat) Creating thread n["<< t<< "] for speakers [ "<<spkBottom<<" --> "<<spkUp-1<<" ]"<<endl;
		rc = pthread_create(&threads[t], &attr, S_MplusVYplusDZthread, (void *)&thread_data_array[t]);
		if (rc) throw Exception("ERROR; return code from pthread_create() is ",__FILE__,rc);
		
		spkBottom = spkUp;
	}
	
	pthread_attr_destroy(&attr);
	for(unsigned long t=0; t<NUM_THREADS; t++) {
		rc = pthread_join(threads[t], (void **)&status);
		if (rc)  throw Exception("ERROR; return code from pthread_join() is ",__FILE__,rc);
		if (verboseLevel>1) cout <<"(AccumulateJFAStat) Completed join with thread ["<<t<<"] status["<<status<<"]"<<endl;
	}

	free(thread_data_array);
	free(threads);

	sessIdx.deleteAllObjects();

	if (verboseLevel >= 1) cout << "(AccumulateJFAStat) Done " << endl;
}

#endif

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
///Fonction qui renvoie le modele M_s_h=M+Vy + Dz+Ux
void JFAAcc::getSpeakerModel(MixtureGD &mixture, String& file){
	
	DoubleVector Sp, MplusVYplusDZ, UX;
	Sp.setSize(_svSize); Sp.setAllValues(0.0);
	MplusVYplusDZ.setSize(_svSize);
	UX.setSize(_svSize);
	
	this->getMplusVYplusDZ(MplusVYplusDZ, file);
	this->getUX(UX, file);
	
	for(unsigned long i=0; i<_svSize; i++){
		Sp[i] = MplusVYplusDZ[i]+UX[i];
	}
	svToModel(Sp,mixture);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
/// Normalize features with a smooth mixture transformation o't=ot-sum(P(c|ot)Uc.x)
void JFAAcc::normalizeFeatures(SegCluster &selectedSegments,FeatureServer &fs,Config & config){
	
	if (verbose) cout << "(AccumulateJFAStat) Normalize Features" << endl;	
	//MixtureGD & clientMixture=_ms.getMixtureGD(1); 						// copy the UBM mixture
//	MixtureGD & world=_ms.getMixtureGD(0);

	MixtureGD & clientMixture=_ms.duplicateMixture (_ms.getMixtureGD(0), DUPL_DISTRIB);
	
	RealVector <double> ux; ux.setSize(_svSize); 						//vecteur utilise pour calculers Ux
	double *_ux=ux.getArray();
	
	bool topGauss = false;
	if(config.existsParam("topGauss"))	topGauss = config.getParam("topGauss").toBool();

	Seg *seg;          												// current selectd segment
	selectedSegments.rewind();
	String currentSource="";
	while((seg=selectedSegments.getSeg())!=NULL){
		unsigned long begin=seg->begin()+fs.getFirstFeatureIndexOfASource(seg->sourceName()); 
		if (currentSource!=seg->sourceName()) {
			currentSource=seg->sourceName();
			this->getUX(ux,currentSource);
			this->getSpeakerModel(clientMixture,currentSource);			//model contenant M_s_h=M+Vy + Dz+Ux
			if (verbose)cout << "Processing speaker["<<currentSource<<"]"<< endl;	
		}

		fs.seekFeature(begin);
		Feature f;

		if (!topGauss) {
			for (unsigned long idxFrame=0;idxFrame<seg->length();idxFrame++){
				fs.readFeature(f,0);
				double *ff=f.getDataVector();				
				double sum=0.0;
				DoubleVector P;
				P.setSize(_n_distrib);
				double *Prob=P.getArray();
				for(unsigned long k=0;k<_n_distrib;k++) {
					Prob[k]=clientMixture.weight(k)*clientMixture.getDistrib(k).computeLK(f);
					
//					Prob[k]=world.weight(k)* world.getDistrib(k).computeLK(f);
					sum+=Prob[k];
				}
				for(unsigned long k=0;k<_n_distrib;k++){
					Prob[k]/=sum;
				}

				//Substract the channel from the feature
				for(unsigned long k=0;k<_n_distrib;k++) {
					for (unsigned long i=0;i<_vectSize;i++) 
						ff[i]-= Prob[k]*_ux[k*_vectSize+i];				//Soustrait les statistiques du canal
					}
				fs.writeFeature(f);
			}	
		}
		else {
			throw Exception("no topgauss yet",__FILE__,__LINE__);
		}
	}
	
	_ms.deleteMixture(clientMixture);
	_ms.deleteUnusedDistribs();
	
};	

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::substractUXfromFeatures(FeatureServer &fs,Config &config){
	SegServer segmentsServer;
	LabelServer labelServer;
	initializeClusters(_fileList,segmentsServer,labelServer,config);
	verifyClusterFile(segmentsServer,fs,config);
	unsigned long codeSelectedFrame=labelServer.getLabelIndexByString(config.getParam("labelSelectedFrames"));
	SegCluster& selectedSegments=segmentsServer.getCluster(codeSelectedFrame);  
	normalizeFeatures(selectedSegments,fs,config);
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::orthonormalizeV(){
	// Gram-Schmidt algorithm
/*
[m,n] = size(B);
Q=zeros(m,n);
R=zeros(m,n);

%Orthonormalisation de A
for j=1:m
    v=B(j,:);
    for i=1:j-1
        R(i,j)=Q(i,:)*B(j,:)';
        v=v-R(i,j)*Q(i,:);
    end
    R(j,j)=norm(v);
    if ( R(j,j) == 0 )
        Q(j,:)=0;
    else
        Q(j,:)=v/R(j,j);
    end;
end*/

	Matrix<double> Q,R;
	Q.setDimensions(_rankEV, _svSize);
	Q.setAllValues(0.0);
	R.setDimensions(_rankEV, _svSize);
	R.setAllValues(0.0);

	for(unsigned long j=0;j<_rankEV;j++){

		DoubleVector v(_svSize,_svSize);
		for(unsigned long k=0; k<_svSize;k++){
			v[k] = _V(j,k);
		}
		for(unsigned long i=0;i<j;i++){
			//R(i,j)=Q(i,:)*B(j,:)';
			for(unsigned long k=0; k<_svSize;k++){
				R(i,j) += Q(i,k)*_V(j,k);
			}
			//v=v-R(i,j)*Q(i,:);
			for(unsigned long k=0; k<_svSize;k++){
				v[k] -= R(i,j)* Q(i,k);
			}
		}

		//R(j,j)=norm(v);
		double nV = 0;
		for(unsigned long k=0; k<_svSize;k++){
			nV += v[k]*v[k];
		}
		R(j,j) = sqrt(nV);

		//if ( R(j,j) == 0 )
		//    Q(j,:)=0;
		//else
		//    Q(j,:)=v/R(j,j);
		//end;
		if(R(j,j) ==0){
			for(unsigned long k=0; k<_svSize;k++){
				Q(j,k)=0;
			}
		}
		else{
			for(unsigned long k=0; k<_svSize;k++){
				Q(j,k)= v[k]/R(j,j);
			}
		}
	}

	//Copy the new matrix in _V
	for(unsigned long i=0;i<_rankEV;i++){
		for(unsigned long j=0;j<_svSize;j++){
			_V(i,j) = Q(i,j);
		}
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
///Save Accumulators on disk
void JFAAcc::saveAccs(Config &config) {

	String fxName, fxhName, nName, nhName;
	fxName = "F_X.mat"; fxhName = "F_X_h.mat"; nName = "N.mat"; nhName = "N_h.mat";
	if(config.existsParam("nullOrderStatSpeaker"))	nName = config.getParam("matrixFilesPath") + config.getParam("nullOrderStatSpeaker")+config.getParam("saveMatrixFilesExtension");
	if(config.existsParam("nullOrderStatSession"))	nhName = config.getParam("matrixFilesPath") + config.getParam("nullOrderStatSession")+config.getParam("saveMatrixFilesExtension");
	if(config.existsParam("firstOrderStatSpeaker")) fxName = config.getParam("matrixFilesPath") + config.getParam("firstOrderStatSpeaker")+config.getParam("saveMatrixFilesExtension");
	if(config.existsParam("firstOrderStatSession"))	fxhName = config.getParam("matrixFilesPath") + config.getParam("firstOrderStatSession")+config.getParam("saveMatrixFilesExtension");

	_F_X.save(fxName,config);
	_F_X_h.save(fxhName,config);
	_N_h.save(nhName,config); 
	_matN.save(nName,config);
	if (verboseLevel>=1) cout << "(AccumulateJFAStat) FA Accs states saved" << endl;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
double JFAAcc::getLLK(SegCluster &selectedSegments,MixtureGD &model,FeatureServer&fs,Config & config){

	if (verboseLevel >= 1) cout << "(FactorAnalysisStat) Compute Likelihood" << endl;
	double llk=0.0;

	MixtureGDStat &acc=_ss.createAndStoreMixtureStat(model);		
	Seg *seg;        
	selectedSegments.rewind();

	while((seg=selectedSegments.getSeg())!=NULL){                           	
		unsigned long begin=seg->begin()+fs.getFirstFeatureIndexOfASource(seg->sourceName()); 
		fs.seekFeature(begin);
		Feature f;
		for (unsigned long idxFrame=0;idxFrame<seg->length();idxFrame++){
			fs.readFeature(f); 
			acc.computeAndAccumulateLLK(f,1.0,TOP_DISTRIBS_NO_ACTION);
		}		
	}

	llk= acc.getMeanLLK();
	_ss.deleteMixtureStat(acc);

return llk;
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void JFAAcc::verifyEMLK(Config& config){

	XList ndx(_fileList);
	XLine *pline; String *pFile; ndx.rewind();

	double total=0.0;
	unsigned long maxLLKcomputed=1;
	maxLLKcomputed=config.getParam("computeLLK").toULong();	

	unsigned long cnt=0;
	while((pline=ndx.getLine())!=NULL && cnt < maxLLKcomputed) { 
		while((pFile=pline->getElement())!=NULL && cnt < maxLLKcomputed) {

			/// Compute JFA model
			MixtureServer ms(config);
			MixtureGD &model=ms.loadMixtureGD(config.getParam("inputWorldFilename"));
			this->getSpeakerModel(model,*pFile);
			
			/// Get LLK
			FeatureServer fs(config,*pFile);
			SegServer segmentsServer;
			LabelServer labelServer;
			initializeClusters(*pFile,segmentsServer,labelServer,config);
			verifyClusterFile(segmentsServer,fs,config);
			unsigned long codeSelectedFrame=labelServer.getLabelIndexByString(config.getParam("labelSelectedFrames"));
			SegCluster& selectedSegments=segmentsServer.getCluster(codeSelectedFrame);  
			double llk=this->getLLK(selectedSegments,model,fs,config); 
			if (verboseLevel>1) cout << "(EigenVoice) LLK["<<*pFile<<"]="<<llk<<endl;
			cnt++;
			total+=llk;
		}
	}
	if (verboseLevel >=1) cout << "*** (Verify LLK) Total LLK="<<total<<" ***"<<endl;
}


#endif

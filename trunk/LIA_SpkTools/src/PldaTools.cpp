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

#if !defined(ALIZE_PldaTools_cpp)
#define ALIZE_PldaTools_cpp

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

#include<PldaTools.h>
#include<iostream>
#include<fstream>
#include<cstdio>
#include<cassert>
#include<cmath>
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
PldaDev::PldaDev(String & ndxFilename,Config & config){

	// Create the XList and sort speakers by decreasing number of sessions
	_fileList.load(ndxFilename,config);
	_fileList.sortByElementNumber("descend");

	_n_speakers = _fileList.getLineCount();
	_session_per_speaker.setSize(_n_speakers);
	_n_sessions = 0;

	// From the first line, first element, get the vectSize
	XLine *linep;
	_fileList.getLine(0);
	linep=_fileList.getLine();
	String fileName = linep->getElement(0);

	String fName = config.getParam("vectorFilesPath") + "/" + fileName + config.getParam("vectorFilesExtension");
	Matrix<double> clientV(fName,config);
	_vectSize = clientV.cols();

	// Compute the total number of sessions
	_fileList.getLine(0);
	unsigned long c = 0;
	while ((linep=_fileList.getLine()) != NULL){
		_session_per_speaker[c] = linep->getElementCount();
		_n_sessions += linep->getElementCount();
		c++;
	}

	// Initialize data, mean vector and matrix
	_mean.setSize(_vectSize);
	_mean.setAllValues(0.0);

	_speaker_means.setDimensions(_vectSize,_n_speakers);
	_speaker_means.setAllValues(0.0);

	_class.setSize(_n_sessions);
	_style.setSize(_n_sessions);

	// Read vectors and fill the _data Matrix
	_data.setDimensions(_vectSize,_n_sessions);

	_fileList.getLine(0);
	unsigned long LCounter = 0;
	unsigned long SCounter = 0;
	while ((linep=_fileList.getLine()) != NULL){
		for(unsigned long s=0;s<linep->getElementCount();s++){

			// Read the vector from file
			fileName = linep->getElement(s);
			String fName = config.getParam("vectorFilesPath") + "/" + fileName + config.getParam("vectorFilesExtension");
			Matrix<double> tmpVect(fName,config);

			if((tmpVect.rows() !=1)||(tmpVect.cols()!=_vectSize))
				throw Exception("Incorrect dimension of vector to load",__FILE__,__LINE__);

			_class[SCounter]	= LCounter;
			_style[SCounter]	= 0;	// to modify for future implementation of Tied_FA and Tied-PLDA

			for(unsigned long k=0;k<_vectSize;k++){
				_data(k,SCounter) = tmpVect(0,k);
				_speaker_means(k,LCounter) += tmpVect(0,k);
				_mean[k] += tmpVect(0,k);
			}
			SCounter++;
		}
		LCounter++;
	}

	// Compute global and speaker means
	for(unsigned long k=0;k<_vectSize;k++){
		_mean[k] /= (double)_n_sessions;
		for(unsigned long s=0;s<_n_speakers;s++){
			_speaker_means(k,s) /= (double)_session_per_speaker[s];
		}
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
PldaDev::PldaDev(XList & fileList,Config & config){

	_fileList.sortByElementNumber("descend");

	_n_speakers = _fileList.getLineCount();
	_session_per_speaker.setSize(_n_speakers);
	_n_sessions = 0;

	// From the first line, first element, get the vectSize
	XLine *linep;
	_fileList.getLine(0);
	linep=_fileList.getLine();
	String fileName = linep->getElement(0);

	String fName = config.getParam("vectorFilesPath") + "/" + fileName + config.getParam("vectorFilesExtension");
	Matrix<double> clientV(fName,config);
	_vectSize = clientV.cols();

	// Compute the total number of sessions
	_fileList.getLine(0);
	unsigned long c = 0;
	while ((linep=_fileList.getLine()) != NULL){
		_session_per_speaker[c] = linep->getElementCount();
		_n_sessions += linep->getElementCount();
		c++;
	}

	// Initialize data, mean vector and matrix
	_mean.setSize(_vectSize);
	_mean.setAllValues(0.0);

	_speaker_means.setDimensions(_vectSize,_n_speakers);
	_speaker_means.setAllValues(0.0);

	_class.setSize(_n_sessions);
	_style.setSize(_n_sessions);

	// Read vectors and fill the _data Matrix
	_data.setDimensions(_vectSize,_n_sessions);

	_fileList.getLine(0);
	unsigned long LCounter = 0;
	unsigned long SCounter = 0;
	while ((linep=_fileList.getLine()) != NULL){
		for(unsigned long s=0;s<linep->getElementCount();s++){

			// Read the vector from file
			fileName = linep->getElement(s);
			String fName = config.getParam("vectorFilesPath") + "/" + fileName + config.getParam("vectorFilesExtension");
			Matrix<double> tmpVect(fName,config);

			// verifier la dimension de la matrice lue...
			if((tmpVect.rows() !=1)||(tmpVect.cols()!=_vectSize))
				throw Exception("Incorrect dimension of vector to load",__FILE__,__LINE__);

			_class[SCounter]	= LCounter;
			_style[SCounter]	= 0;	// to modify for future implementation of Tied_FA and Tied-PLDA

			for(unsigned long k=0;k<_vectSize;k++){
				_data(k,SCounter) = tmpVect(0,k);
				_speaker_means(k,LCounter) += tmpVect(0,k);
				_mean[k] += tmpVect(0,k);
			}
			SCounter++;
		}
		LCounter++;
	}

	// Compute global and speaker means
	for(unsigned long k=0;k<_vectSize;k++){
		_mean[k] /= (double)_n_sessions;
		for(unsigned long s=0;s<_n_speakers;s++){
			_speaker_means(k,s) /= (double)_session_per_speaker[s];
		}
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
PldaDev::~PldaDev(){}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
String PldaDev::getClassName() const{	return "PldaDev";}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaDev::computeAll(){

	// Reset all mean values
	_mean.setAllValues(0.0);
	_speaker_means.setAllValues(0.0);

	// From the first line, first element, get the vectSize
	XLine *linep;
	_fileList.getLine(0);
	unsigned long LCounter = 0;
	unsigned long SCounter = 0;
	while ((linep=_fileList.getLine()) != NULL){
		for(unsigned long s=0;s<linep->getElementCount();s++){
			for(unsigned long k=0;k<_vectSize;k++){
				_speaker_means(k,LCounter) += _data(k,SCounter);
				_mean[k] += _data(k,SCounter);
			}
			SCounter++;
		}
		LCounter++;
	}

	// Compute global and speaker means
	for(unsigned long k=0;k<_vectSize;k++){
		_mean[k] /= (double)_n_sessions;
		for(unsigned long s=0;s<_n_speakers;s++){
			_speaker_means(k,s) /= (double)_session_per_speaker[s];
		}
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned long PldaDev::getVectSize(){
	return _vectSize;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned long PldaDev::getSpeakerNumber(){
	return  _n_speakers;
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned long PldaDev::getSessionNumber(){
	return _n_sessions;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
Matrix<double> PldaDev::getData(){
	return _data;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
RealVector<double>& PldaDev::getMean(){
	return _mean;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaDev::getSpeakerMean(unsigned long spk,RealVector<double>& spkMean){
	spkMean.setSize(_vectSize);
	for(unsigned long i=0;i<_vectSize;i++)
		spkMean[i] = _speaker_means(i,spk);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaDev::lengthNorm(){

	// Compute the norm of all vectors in _data
	RealVector<double> vecNorm(0);
	vecNorm.setSize(_n_sessions);
	vecNorm.setAllValues(0.0);

	for(unsigned long s =0;s<_n_sessions;s++){
		double tmp = 0.0;
		for(unsigned long k=0;k<_vectSize;k++){
			tmp += _data(k,s)*_data(k,s);
		}
		vecNorm[s] = sqrt(tmp);
	}

	// Divide all vectors by their norm
	for(unsigned long s =0;s<_n_sessions;s++){
		for(unsigned long k=0;k<_vectSize;k++){
			_data(k,s) /=vecNorm[s];
		}
	}
	// Recompute the means
	this->computeAll();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaDev::center(RealVector<double> &mu){
	for(unsigned long s =0;s<_n_sessions;s++){
		for(unsigned long k=0;k<_vectSize;k++){
			_data(k,s) -= mu[k];
		}
	}
	// Recompute the means
	this->computeAll();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaDev::computeCovMat(DoubleSquareMatrix &Sigma, DoubleSquareMatrix &W, DoubleSquareMatrix &B, Config &config){

	#ifdef THREAD          
	if (config.existsParam("numThread") && config.getParam("numThread").toULong() > 0)	computeCovMatThreaded(Sigma,W,B,config.getParam("numThread").toULong());
	else computeCovMatUnThreaded(Sigma,W,B);
	#else
	computeCovMatUnThreaded(Sigma,W,B);
	#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaDev::computeCovMatUnThreaded(DoubleSquareMatrix &Sigma, DoubleSquareMatrix &W, DoubleSquareMatrix &B){

	// Initialize matrices
	Sigma.setSize(_vectSize);
	Sigma.setAllValues(0.0);
	W.setSize(_vectSize);
	W.setAllValues(0.0);
	B.setSize(_vectSize);
	B.setAllValues(0.0);

	for(unsigned long i=0;i<_vectSize;i++){
		for(unsigned long j=0;j<_vectSize;j++){
			for(unsigned long s=0;s<_n_sessions;s++){
				Sigma(i,j)	+= (_data(i,s)-_mean[i])*(_data(j,s)-_mean[j]);
				W(i,j)		+= (_data(i,s)-_speaker_means(i,_class[s]))*(_data(j,s)-_speaker_means(j,_class[s]));
			}

			for(unsigned long c=0;c<_n_speakers;c++){
				B(i,j)		+= _session_per_speaker[c]* (_speaker_means(i,c)-_mean[i]) * (_speaker_means(j,c)-_mean[j]);	
			}
		}
	}

	for(unsigned long i=0;i<_vectSize;i++){
		for(unsigned long j=0;j<_vectSize;j++){
			Sigma(i,j)	/= _n_sessions;
			W(i,j)		/= _n_sessions;
			B(i,j)		/= _n_sessions;
		}
	}
}

#ifdef THREAD
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Data structure of thread
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
struct Covthread_data{

	double *data;
	double *speaker_means;
	double *mean;
	unsigned long *_class;
	unsigned long *session_per_speaker;
	unsigned long n_speakers;
	unsigned long n_sessions;
	unsigned long vectSize;
	unsigned long startSession;
	unsigned long stopSession;
	unsigned long threadNb;
	RefVector <DoubleSquareMatrix> *sigma;
	RefVector <DoubleSquareMatrix> *w;
	RefVector <DoubleSquareMatrix> *b;

};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Thread Routine
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void *Covthread(void *threadarg){
	struct Covthread_data *my_data;
	my_data = (struct Covthread_data *) threadarg;

	double *data = my_data->data;
	double *speaker_means = my_data->speaker_means;
	double *mean = my_data->mean;
	unsigned long *_class = my_data->_class;
	unsigned long *session_per_speaker = my_data->session_per_speaker;
	unsigned long n_speakers = my_data->n_speakers;
	unsigned long n_sessions = my_data->n_sessions;	
	unsigned long vectSize = my_data->vectSize;	
	unsigned long startSession = my_data->startSession;
	unsigned long stopSession = my_data->stopSession;
	unsigned long threadNb = my_data->threadNb;

	DoubleSquareMatrix &_sigma=(*(my_data->sigma))[threadNb];
	DoubleSquareMatrix &_w=(*(my_data->w))[threadNb];
	DoubleSquareMatrix &_b=(*(my_data->b))[threadNb];

	double *sigma = _sigma.getArray();
	double *w = _w.getArray();
	double *b = _b.getArray();

	for(unsigned long i=0;i<vectSize;i++){
		for(unsigned long j=0;j<vectSize;j++){
			for(unsigned long s=startSession;s<stopSession;s++){
				_sigma(i,j)	+= (data[i*n_sessions+s]-mean[i])*(data[j*n_sessions+s]-mean[j]);
				_w(i,j)		+= (data[i*n_sessions+s]-speaker_means[i*n_speakers+_class[s]])*(data[j*n_sessions+s]-speaker_means[j*n_speakers+_class[s]]);
			}
			for(unsigned long c=_class[startSession];c<_class[stopSession-1]+1;c++){
				_b(i,j)		+= session_per_speaker[c]* (speaker_means[i*n_speakers+c]-mean[i]) * (speaker_means[j*n_speakers+c]-mean[j]);
			}
		}
	}
	
	pthread_exit((void*) 0);
	return (void*)0 ;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaDev::computeCovMatThreaded(DoubleSquareMatrix &Sigma, DoubleSquareMatrix &W, DoubleSquareMatrix &B, unsigned long NUM_THREADS){

	if (NUM_THREADS==0) throw Exception("Num threads can not be 0",__FILE__,__LINE__);

	Sigma.setSize(_vectSize); W.setSize(_vectSize); B.setSize(_vectSize);
	Sigma.setAllValues(0.0); W.setAllValues(0.0); B.setAllValues(0.0);

	// Compute the covariance matrices
	if(verboseLevel > 0) cout<<"	Compute Covariance matrices"<<endl;

	// split the list
	RealVector<unsigned long> startIndex;
	this->splitPerSpeaker(NUM_THREADS, startIndex);

	// create the 3 RefVector<DoubleSquareMatrix>
	RefVector<DoubleSquareMatrix> sigMatrices;
	for(unsigned long nt=0;nt<NUM_THREADS;nt++){
		sigMatrices.addObject(*new DoubleSquareMatrix(_vectSize));
		sigMatrices[nt].setAllValues(0.0);
	}

	RefVector<DoubleSquareMatrix> wMatrices;
	for(unsigned long nt=0;nt<NUM_THREADS;nt++){
		wMatrices.addObject(*new DoubleSquareMatrix(_vectSize));
		wMatrices[nt].setAllValues(0.0);
	}

	RefVector<DoubleSquareMatrix> bMatrices;
	for(unsigned long nt=0;nt<NUM_THREADS;nt++){
		bMatrices.addObject(*new DoubleSquareMatrix(_vectSize));
		bMatrices[nt].setAllValues(0.0);
	}

	// threads
	int rc, status;
	if (NUM_THREADS > _n_speakers) NUM_THREADS=_n_speakers;
	
	struct Covthread_data *thread_data_array = new Covthread_data[NUM_THREADS];
	pthread_t *threads = new pthread_t[NUM_THREADS];

	double *speaker_means, *data, *mean;
	unsigned long *Class, *session_per_speaker;
	speaker_means = _speaker_means.getArray(); data = _data.getArray(); mean = _mean.getArray(); Class=_class.getArray(); session_per_speaker = _session_per_speaker.getArray();

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	
	//Create threads : one per distribution as a maximum
	for(unsigned long t=0; t<NUM_THREADS; t++){

		thread_data_array[t].data = data;
		thread_data_array[t].speaker_means = speaker_means;
		thread_data_array[t].mean = mean;
		thread_data_array[t]._class = Class;
		thread_data_array[t].session_per_speaker = session_per_speaker;
		thread_data_array[t].n_speakers = _n_speakers;
		thread_data_array[t].n_sessions = _n_sessions;
		thread_data_array[t].vectSize = _vectSize;
		thread_data_array[t].startSession = startIndex[t];
		if(t<NUM_THREADS-1){
			thread_data_array[t].stopSession = startIndex[t+1];
		}
		else{
			thread_data_array[t].stopSession = _n_sessions;
		}

		thread_data_array[t].threadNb = t;
		thread_data_array[t].sigma=&(sigMatrices);
		thread_data_array[t].w=&(wMatrices);
		thread_data_array[t].b=&(bMatrices);

		if (verboseLevel > 0) cout<<"(computeCovMat) Creating thread n["<< t<< "] for sessions["<<startIndex[t]<<"-->"<<thread_data_array[t].stopSession<<"]"<<endl;
		rc = pthread_create(&threads[t], &attr, Covthread, (void *)&thread_data_array[t]);
		if (rc) throw Exception("ERROR; return code from pthread_create() is ",__FILE__,rc);

	}

	pthread_attr_destroy(&attr);

	for(unsigned long t=0; t<NUM_THREADS; t++) {
		rc = pthread_join(threads[t], (void **)&status);
		if (rc)  throw Exception("ERROR; return code from pthread_join() is ",__FILE__,rc);
		if (verboseLevel >1) cout <<"(computeCovMat) Completed join with thread ["<<t<<"] status["<<status<<"]"<<endl;
	}

	free(thread_data_array);
	free(threads);

	// join and sum the matrices
	for(unsigned long t=0;t<NUM_THREADS;t++){
		for(unsigned long i=0;i<_vectSize;i++){
			for(unsigned long j=0;j<_vectSize;j++){
				Sigma(i,j) += sigMatrices[t](i,j);
				W(i,j) += wMatrices[t](i,j);
				B(i,j) += bMatrices[t](i,j);
			}
		}
	}
	sigMatrices.deleteAllObjects();
	wMatrices.deleteAllObjects();
	bMatrices.deleteAllObjects();

	// Normalize 
	for(unsigned long i=0;i<_vectSize;i++){
		for(unsigned long j=0;j<_vectSize;j++){
			Sigma(i,j)	/= _n_sessions;
			W(i,j)		/= _n_sessions;
			B(i,j)		/= _n_sessions;
		}
	}
}
#endif

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaDev::computeWccnChol(DoubleSquareMatrix &WCCN, Config &config){

	#ifdef THREAD          
	if (config.existsParam("numThread") && config.getParam("numThread").toULong() > 0)	computeWccnCholThreaded(WCCN,config.getParam("numThread").toULong());
	else computeWccnCholUnThreaded(WCCN);
	#else
	computeWccnCholUnThreaded(WCCN);
	#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaDev::computeWccnCholUnThreaded(DoubleSquareMatrix &WCCN){

	WCCN.setSize(_vectSize);
	WCCN.setAllValues(0.0);

	DoubleSquareMatrix W(_vectSize);
	W.setAllValues(0.0);

	// Compute the covariance matrix per class
	if(verboseLevel > 0) cout<<"	Compute WCCN matrix"<<endl;
	unsigned long spk = 0;
	unsigned long s = 0;

	while(s < _n_sessions){
		DoubleSquareMatrix covSpk(_vectSize);
		covSpk.setAllValues(0.0);
		unsigned long sessionNumber = _session_per_speaker[spk];
		while((s<_n_sessions)&&(_class[s] == spk)){
			for(unsigned long i=0;i<_vectSize;i++){
				for(unsigned long j=0;j<_vectSize;j++){
					covSpk(i,j) += (_data(i,s)-_speaker_means(i,_class[s]))*(_data(j,s)-_speaker_means(j,_class[s]));
				}
			}
			s++;
		}
		for(unsigned long i=0;i<_vectSize;i++){
			for(unsigned long j=0;j<_vectSize;j++){
				W(i,j) += covSpk(i,j)/(double)sessionNumber;
			}
		}
		if(s<_n_sessions)	spk = _class[s];
	}

	// Normalize the WCCN Matrix by the number of speakers
	for(unsigned long i=0;i<_vectSize;i++){
		for(unsigned long j=0;j<_vectSize;j++){
			W(i,j) /= _n_speakers;
		}
	}

    // Invert the WCCN matrix
	if(verboseLevel > 0) cout<<"	Invert WCCN matrix"<<endl;
	DoubleSquareMatrix invW(_vectSize);
	W.invert(invW);

	// Choleski decompostion of WCCN inverse
	if(verboseLevel > 0) cout<<"	Cholesky decomposition of inv(WCCN)"<<endl;
	invW.upperCholesky(WCCN);

}

#ifdef THREAD
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Data structure of thread
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
struct WCCNthread_data{

	double *data;
	double *speaker_means;
	unsigned long *_class;
	unsigned long *session_per_speaker;
	unsigned long n_speakers;
	unsigned long n_sessions;
	unsigned long vectSize;
	unsigned long startSession;
	unsigned long stopSession;
	unsigned long threadNb;
	RefVector <DoubleSquareMatrix> *cov;
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Thread Routine
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void *WCCNthread(void *threadarg){
	struct WCCNthread_data *my_data;
	my_data = (struct WCCNthread_data *) threadarg;

	double *data = my_data->data;
	double *speaker_means = my_data->speaker_means;
	unsigned long *_class = my_data->_class;
	unsigned long *session_per_speaker = my_data->session_per_speaker;
	unsigned long n_speakers = my_data->n_speakers;
	unsigned long n_sessions = my_data->n_sessions;	
	unsigned long vectSize = my_data->vectSize;	
	unsigned long startSession = my_data->startSession;
	unsigned long stopSession = my_data->stopSession;
	unsigned long threadNb = my_data->threadNb;

	DoubleSquareMatrix &cov=(*(my_data->cov))[threadNb];
	double *c = cov.getArray();
	unsigned long spk = _class[startSession];
	unsigned long s = startSession;

	while(s < stopSession){
		DoubleSquareMatrix covSpk(vectSize);
		covSpk.setAllValues(0.0);
		unsigned long sessionNumber = session_per_speaker[spk];
		while((s<stopSession)&&(_class[s] == spk)){
			for(unsigned long i=0;i<vectSize;i++){
				for(unsigned long j=0;j<vectSize;j++){
					covSpk(i,j) += (data[i*n_sessions+s]-speaker_means[i*n_speakers+_class[s]])*(data[j*n_sessions+s]-speaker_means[j*n_speakers+_class[s]]);
				}
			}
			s++;
		}
		for(unsigned long i=0;i<vectSize;i++){
			for(unsigned long j=0;j<vectSize;j++){
				c[i*vectSize+j] += covSpk(i,j)/(double)sessionNumber;
			}
		}
		if(s<n_sessions)	spk = _class[s];
	}

	pthread_exit((void*) 0);
	return (void*)0 ;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaDev::computeWccnCholThreaded(DoubleSquareMatrix &WCCN, unsigned long NUM_THREADS){

	cerr<<"computeWccnCholThreaded"<<endl;
	if (NUM_THREADS==0) throw Exception("Num threads can not be 0",__FILE__,__LINE__);

	WCCN.setSize(_vectSize);
	WCCN.setAllValues(0.0);

	DoubleSquareMatrix W(_vectSize);
	W.setAllValues(0.0);

	// Compute the covariance matrix per class
	if(verboseLevel > 0) cout<<"	Compute WCCN matrix"<<endl;

	// split the list
	RealVector<unsigned long> startIndex;
	this->splitPerSpeaker(NUM_THREADS, startIndex);

	// create the RefVector<DoubleSquareMatrix>
	RefVector<DoubleSquareMatrix> covMatrices;
	for(unsigned long nt=0;nt<NUM_THREADS;nt++){
		covMatrices.addObject(*new DoubleSquareMatrix(_vectSize));
		covMatrices[nt].setAllValues(0.0);
	}

	// threads
	int rc, status;
	if (NUM_THREADS > _n_speakers) NUM_THREADS=_n_speakers;
	
	struct WCCNthread_data *thread_data_array = new WCCNthread_data[NUM_THREADS];
	pthread_t *threads = new pthread_t[NUM_THREADS];

	double *speaker_means, *data;
	unsigned long *Class, *session_per_speaker;
	speaker_means = _speaker_means.getArray(); data = _data.getArray(); Class=_class.getArray(); session_per_speaker = _session_per_speaker.getArray();

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	
	//Create threads : one per distribution as a maximum
	for(unsigned long t=0; t<NUM_THREADS; t++){

		thread_data_array[t].data = data;
		thread_data_array[t].speaker_means = speaker_means;
		thread_data_array[t]._class = Class;
		thread_data_array[t].session_per_speaker = session_per_speaker;
		thread_data_array[t].n_speakers = _n_speakers;
		thread_data_array[t].n_sessions = _n_sessions;
		thread_data_array[t].vectSize = _vectSize;
		thread_data_array[t].startSession = startIndex[t];
		if(t<NUM_THREADS-1){
			thread_data_array[t].stopSession = startIndex[t+1];
		}
		else{
			thread_data_array[t].stopSession = _n_sessions;
		}

		thread_data_array[t].threadNb = t;
		thread_data_array[t].cov=&(covMatrices);

		if (verboseLevel > 0) cout<<"(computeWccnChol) Creating thread n["<< t<< "] for sessions["<<startIndex[t]<<"-->"<<thread_data_array[t].stopSession<<"]"<<endl;
		rc = pthread_create(&threads[t], &attr, WCCNthread, (void *)&thread_data_array[t]);
		if (rc) throw Exception("ERROR; return code from pthread_create() is ",__FILE__,rc);

	}

	pthread_attr_destroy(&attr);

	for(unsigned long t=0; t<NUM_THREADS; t++) {
		rc = pthread_join(threads[t], (void **)&status);
		if (rc)  throw Exception("ERROR; return code from pthread_join() is ",__FILE__,rc);
		if (verboseLevel >1) cout <<"(computeWccnChol) Completed join with thread ["<<t<<"] status["<<status<<"]"<<endl;
	}

	free(thread_data_array);
	free(threads);

	// join and sum the matrices
	for(unsigned long t=0;t<NUM_THREADS;t++){
		for(unsigned long i=0;i<_vectSize;i++){
			for(unsigned long j=0;j<_vectSize;j++){
				W(i,j) += covMatrices[t](i,j);
			}
		}
	}
	covMatrices.deleteAllObjects();

	// Normalize the WCCN Matrix by the number of speakers
	for(unsigned long i=0;i<_vectSize;i++){
		for(unsigned long j=0;j<_vectSize;j++){
			W(i,j) /= _n_speakers;
		}
	}

    // Invert the WCCN matrix
	if(verboseLevel > 0) cout<<"	Invert WCCN matrix"<<endl;
	DoubleSquareMatrix invW(_vectSize);
	W.invert(invW);

	// Choleski decompostion of WCCN inverse
	if(verboseLevel > 0) cout<<"	Cholesky decomposition of inv(WCCN)"<<endl;
	invW.upperCholesky(WCCN);
}
#endif

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaDev::computeMahalanobis(DoubleSquareMatrix &M, Config& config){

	// Initialize matrices
	M.setSize(_vectSize);
	M.setAllValues(0.0);

	DoubleSquareMatrix W(_vectSize), B(_vectSize), Sigma(_vectSize);
	W.setAllValues(0.0); B.setAllValues(0.0); Sigma.setAllValues(0.0);

	this->computeCovMat(Sigma,W,B,config);

	W.invert(M);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaDev::splitPerSpeaker(unsigned long nbThread, RealVector<unsigned long> &startIndex){

	unsigned long spkPerThread = floor((double)_n_speakers/(double)nbThread);
	startIndex.setSize(nbThread);
	startIndex[0] = 0;
	
	unsigned long firstSpeakerNextList = spkPerThread;
	unsigned long spkCounter = 1;
	unsigned long currentSession = 0;

	while((currentSession<_n_sessions)&&(spkCounter<nbThread)){
		while(_class[currentSession]<spkCounter*spkPerThread){
			currentSession++;
		}
		startIndex[spkCounter] = currentSession;
		spkCounter++;
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaDev::computeScatterMat(DoubleSquareMatrix &SB, DoubleSquareMatrix &SW, Config &config){

	#ifdef THREAD          
	if (config.existsParam("numThread") && config.getParam("numThread").toULong() > 0)	computeScatterMatThreaded(SB,SW,config.getParam("numThread").toULong());
	else computeScatterMatUnThreaded(SB,SW);
	#else
	computeScatterMatUnThreaded(Sigma,W,B);
	#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaDev::computeScatterMatUnThreaded(DoubleSquareMatrix &SB, DoubleSquareMatrix &SW){

	// Initialize matrices
	SB.setSize(_vectSize);
	SB.setAllValues(0.0);
	SW.setSize(_vectSize);
	SW.setAllValues(0.0);

	unsigned long sessionCounter = 0;
	for(unsigned long c=0;c<_n_speakers;c++){

		DoubleSquareMatrix tmpSB;
		tmpSB.setSize(_vectSize);
		tmpSB.setAllValues(0.0);

		for(unsigned long i=0;i<_vectSize;i++){
			for(unsigned long j=0;j<_vectSize;j++){
				SB(i,j)		+= (_speaker_means(i,c)-_mean[i]) * (_speaker_means(j,c)-_mean[j]);	

				for(unsigned long s=0;s<_session_per_speaker[c];s++){
					tmpSB(i,j)		+= (_data(i,s)-_speaker_means(i,_class[s]))*(_data(j,s)-_speaker_means(j,_class[s]));
					sessionCounter++;
				}
			}
		}
		for(unsigned long i=0;i<_vectSize;i++){
			for(unsigned long j=0;j<_vectSize;j++){
				SW(i,j) = tmpSB(i,j) / _session_per_speaker[c];
			}
		}
	}
}

#ifdef THREAD
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Data structure of thread
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
struct Scatterthread_data{

	double *data;
	double *speaker_means;
	double *mean;
	unsigned long *_class;
	unsigned long *session_per_speaker;
	unsigned long n_speakers;
	unsigned long n_sessions;
	unsigned long vectSize;
	unsigned long startSession;
	unsigned long stopSession;
	unsigned long threadNb;
	RefVector <DoubleSquareMatrix> *sb;
	RefVector <DoubleSquareMatrix> *sw;

};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Thread Routine
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void *Scatterthread(void *threadarg){
	struct Scatterthread_data *my_data;
	my_data = (struct Scatterthread_data *) threadarg;

	double *data = my_data->data;
	double *speaker_means = my_data->speaker_means;
	double *mean = my_data->mean;
	unsigned long *_class = my_data->_class;
	unsigned long *session_per_speaker = my_data->session_per_speaker;
	unsigned long n_speakers = my_data->n_speakers;
	unsigned long n_sessions = my_data->n_sessions;	
	unsigned long vectSize = my_data->vectSize;	
	unsigned long startSession = my_data->startSession;
	unsigned long stopSession = my_data->stopSession;
	unsigned long threadNb = my_data->threadNb;

	DoubleSquareMatrix &_sb=(*(my_data->sb))[threadNb];
	DoubleSquareMatrix &_sw=(*(my_data->sw))[threadNb];

	double *sb = _sb.getArray();
	double *sw = _sw.getArray();

	unsigned long spk = _class[startSession];
	unsigned long s = startSession;

	DoubleSquareMatrix tmpSW;
	tmpSW.setSize(vectSize);
	tmpSW.setAllValues(0.0);

	unsigned long currentClass = _class[startSession];
	while(s<stopSession){
		while(currentClass == _class[s]){
			for(unsigned long i=0;i<vectSize;i++){
				for(unsigned long j=0;j<vectSize;j++){
					tmpSW(i,j)		+= (data[i*n_sessions+s]-speaker_means[i*n_speakers+_class[s]])*(data[j*n_sessions+s]-speaker_means[j*n_speakers+_class[s]]);
				}
			}
			s++;
		}
		for(unsigned long i=0;i<vectSize;i++){
			for(unsigned long j=0;j<vectSize;j++){
				sw[i*vectSize+j]		+= tmpSW(i,j) / session_per_speaker[currentClass];
				sb[i*vectSize+j]		+= (speaker_means[i*n_speakers+currentClass] - mean[i]) * (speaker_means[j*n_speakers+currentClass] - mean[j]);
			}
		}
		tmpSW.setAllValues(0.0);
		if(s<stopSession) currentClass = _class[s];
	}

	pthread_exit((void*) 0);
	return (void*)0 ;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaDev::computeScatterMatThreaded(DoubleSquareMatrix &SB, DoubleSquareMatrix &SW, unsigned long NUM_THREADS){

	if (NUM_THREADS==0) throw Exception("Num threads can not be 0",__FILE__,__LINE__);

	SB.setSize(_vectSize); SB.setAllValues(0.0);
	SW.setSize(_vectSize); SW.setAllValues(0.0);

	// Compute the covariance matrices
	if(verboseLevel > 0) cout<<"	Compute Scatter matrices"<<endl;

	// split the list
	RealVector<unsigned long> startIndex;
	this->splitPerSpeaker(NUM_THREADS, startIndex);

	// create the 2 RefVector<DoubleSquareMatrix>
	RefVector<DoubleSquareMatrix> sbMatrices;
	for(unsigned long nt=0;nt<NUM_THREADS;nt++){
		sbMatrices.addObject(*new DoubleSquareMatrix(_vectSize));
		sbMatrices[nt].setAllValues(0.0);
	}

	RefVector<DoubleSquareMatrix> swMatrices;
	for(unsigned long nt=0;nt<NUM_THREADS;nt++){
		swMatrices.addObject(*new DoubleSquareMatrix(_vectSize));
		swMatrices[nt].setAllValues(0.0);
	}

	// threads
	int rc, status;
	if (NUM_THREADS > _n_speakers) NUM_THREADS=_n_speakers;
	
	struct Scatterthread_data *thread_data_array = new Scatterthread_data[NUM_THREADS];
	pthread_t *threads = new pthread_t[NUM_THREADS];

	double *speaker_means, *data, *mean;
	unsigned long *Class, *session_per_speaker;
	speaker_means = _speaker_means.getArray(); data = _data.getArray(); mean = _mean.getArray(); Class=_class.getArray(); session_per_speaker = _session_per_speaker.getArray();

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	
	//Create threads : one per distribution as a maximum
	for(unsigned long t=0; t<NUM_THREADS; t++){

		thread_data_array[t].data = data;
		thread_data_array[t].speaker_means = speaker_means;
		thread_data_array[t].mean = mean;
		thread_data_array[t]._class = Class;
		thread_data_array[t].session_per_speaker = session_per_speaker;
		thread_data_array[t].n_speakers = _n_speakers;
		thread_data_array[t].n_sessions = _n_sessions;
		thread_data_array[t].vectSize = _vectSize;
		thread_data_array[t].startSession = startIndex[t];
		if(t<NUM_THREADS-1){
			thread_data_array[t].stopSession = startIndex[t+1];
		}
		else{
			thread_data_array[t].stopSession = _n_sessions;
		}

		thread_data_array[t].threadNb = t;
		thread_data_array[t].sb=&(sbMatrices);
		thread_data_array[t].sw=&(swMatrices);

		if (verboseLevel > 0) cout<<"(computeCovMat) Creating thread n["<< t<< "] for sessions["<<startIndex[t]<<"-->"<<thread_data_array[t].stopSession<<"]"<<endl;
		rc = pthread_create(&threads[t], &attr, Scatterthread, (void *)&thread_data_array[t]);
		if (rc) throw Exception("ERROR; return code from pthread_create() is ",__FILE__,rc);

	}

	pthread_attr_destroy(&attr);

	for(unsigned long t=0; t<NUM_THREADS; t++) {
		rc = pthread_join(threads[t], (void **)&status);
		if (rc)  throw Exception("ERROR; return code from pthread_join() is ",__FILE__,rc);
		if (verboseLevel >1) cout <<"(computeCovMat) Completed join with thread ["<<t<<"] status["<<status<<"]"<<endl;
	}

	free(thread_data_array);
	free(threads);

	// join and sum the matrices
	for(unsigned long t=0;t<NUM_THREADS;t++){
		for(unsigned long i=0;i<_vectSize;i++){
			for(unsigned long j=0;j<_vectSize;j++){
				SB(i,j) += sbMatrices[t](i,j);
				SW(i,j) += swMatrices[t](i,j);
			}
		}
	}
	sbMatrices.deleteAllObjects();
	swMatrices.deleteAllObjects();
}
#endif






























//***********************************************
//	PldaTest methods
//***********************************************

PldaTest::PldaTest(String & ndxFilename,Config &config){

	// Create the XList
	XList testList;
	testList.load(ndxFilename,config);

	_segLine.reset();
	_modelLine.reset();

	// Read the XList and fill the XLines when required
	XLine *linep;
	testList.getLine(0);
	String vect;

	while ((linep=testList.getLine()) != NULL){
		
		vect = linep->getElement(0);
		if(_segLine.getIndex(vect) == -1){
			_segLine.addElement(vect);
		}

		unsigned long e = 1;
		while(e<linep->getElementCount()){

			vect = linep->getElement(e);
			if(_modelLine.getIndex(vect) == -1){
				_modelLine.addElement(vect);
			}
			e++;
		}
	}
	_n_models	= _modelLine.getElementCount();
	_n_segments = _segLine.getElementCount();

	// Get _vectSize from the first model
	String tmp = _modelLine.getElement(0);
	String fName = config.getParam("testVectorFilesPath") + "/" + tmp + config.getParam("testVectorFilesExtension");
	Matrix<double> tmpVect(fName,config);
	_vectSize = tmpVect.cols();

	_models.setDimensions(_vectSize,_n_models);
	_segments.setDimensions(_vectSize,_n_segments);

	if(config.existsParam("verboseLevel") && (config.getParam("verboseLevel").toLong()>0)){
		cout<<"PldaTest found: "<<_n_models<<" models and "<<_n_segments<<" test segments"<<endl;
	}

	// Load model vectors
	String fileName;
	for(unsigned long m=0;m<_n_models;m++){

		// Read the vector from file
		fileName = _modelLine.getElement(m);
		String fName = config.getParam("testVectorFilesPath") + "/" + fileName + config.getParam("testVectorFilesExtension");
		Matrix<double> tmpVect(fName,config);

		for(unsigned long k=0;k<_vectSize;k++){
			_models(k,m) = tmpVect(0,k);
		}
	}

	// Load segment vectors
	for(unsigned long s=0;s<_n_segments;s++){
		// Read the vector from file
		fileName = _segLine.getElement(s);
		String fName = config.getParam("testVectorFilesPath") + "/" + fileName + config.getParam("testVectorFilesExtension");
		Matrix<double> tmpVect(fName,config);
		for(unsigned long k=0;k<_vectSize;k++){
			_segments(k,s) = tmpVect(0,k);
		}
	}

	// Initialize Score matrix with the correct dimensions
	_trials.setDimensions(_n_models,_n_segments);
	_trials.setAllValues(0);

	// Initialize Score matrix with the correct dimensions
	_scores.setDimensions(_n_models,_n_segments);
	_scores.setAllValues(0.0);

	//Initialize and fill _trials matrix
	_scores.setDimensions(_n_models,_n_segments);
	_scores.setAllValues(0);

	testList.getLine(0);
	while ((linep=testList.getLine()) != NULL){
	
		vect = linep->getElement(0);
		unsigned long seg = (unsigned long)_segLine.getIndex(vect);

		unsigned long e = 1;
		while(e<linep->getElementCount()){
			vect = linep->getElement(e);
			unsigned long model = (unsigned long)_modelLine.getIndex(vect);
			_trials(model,seg) = 1; 
			e++;
		}
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
PldaTest::PldaTest(XList &testList,Config & config){

	_segLine.reset();
	_modelLine.reset();

	// Read the XList and fill the XLines when required
	XLine *linep;
	testList.getLine(0);
	String vect;

	while ((linep=testList.getLine()) != NULL){
		
		vect = linep->getElement(0);
		if(_segLine.getIndex(vect) == -1){
			_segLine.addElement(vect);
		}

		unsigned long e = 1;
		while(e<linep->getElementCount()){

			vect = linep->getElement(e);
			if(_modelLine.getIndex(vect) == -1){
				_modelLine.addElement(vect);
			}
			e++;
		}
	}
	_n_models	= _modelLine.getElementCount();
	_n_segments = _segLine.getElementCount();

	// Get _vectSize from the first model
	String tmp = _modelLine.getElement(0);
	String fName = config.getParam("testVectorFilesPath") + "/" + tmp + config.getParam("testVectorFilesExtension");
	Matrix<double> tmpVect(fName,config);
	_vectSize = tmpVect.cols();

	_models.setDimensions(_vectSize,_n_models);
	_segments.setDimensions(_vectSize,_n_segments);

	if(config.existsParam("verboseLevel") && (config.getParam("verboseLevel").toLong()>0)){
		cout<<"PldaTest found: "<<_n_models<<" models and "<<_n_segments<<" test segments"<<endl;
	}

	// Load model vectors
	String fileName;
	for(unsigned long m=0;m<_n_models;m++){

		// Read the vector from file
		fileName = _modelLine.getElement(m);
		String fName = config.getParam("testVectorFilesPath") + "/" + fileName + config.getParam("testVectorFilesExtension");
		Matrix<double> tmpVect(fName,config);

		for(unsigned long k=0;k<_vectSize;k++){
			_models(k,m) = tmpVect(0,k);
		}
	}

	// Load segment vectors
	for(unsigned long s=0;s<_n_segments;s++){
		// Read the vector from file
		fileName = _segLine.getElement(s);
		String fName = config.getParam("testVectorFilesPath") + "/" + fileName + config.getParam("testVectorFilesExtension");
		Matrix<double> tmpVect(fName,config);
		for(unsigned long k=0;k<_vectSize;k++){
			_segments(k,s) = tmpVect(0,k);
		}
	}

	// Initialize Score matrix with the correct dimensions
	_trials.setDimensions(_n_models,_n_segments);
	_trials.setAllValues(0);

	// Initialize Score matrix with the correct dimensions
	_scores.setDimensions(_n_models,_n_segments);
	_scores.setAllValues(0.0);

	//Initialize and fill _trials matrix
	_scores.setDimensions(_n_models,_n_segments);
	_scores.setAllValues(0);

	testList.getLine(0);
	while ((linep=testList.getLine()) != NULL){
	
		vect = linep->getElement(0);
		unsigned long seg = (unsigned long)_segLine.getIndex(vect);

		unsigned long e = 1;
		while(e<linep->getElementCount()){
			vect = linep->getElement(e);
			unsigned long model = (unsigned long)_modelLine.getIndex(vect);
			_trials(model,seg) = 1; 
			e++;
		}
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
PldaTest::~PldaTest(){}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
String PldaTest::getClassName() const{	return "PldaTest";}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned long PldaTest::getVectSize(){
	return _vectSize;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned long PldaTest::getModelsNumber(){
	return _n_models;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned long PldaTest::getSegmentsNumber(){
	return _n_segments;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
Matrix<double> PldaTest::getModels(){
	return _models;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
Matrix<double> PldaTest::getSegments(){
	return _segments;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
Matrix<double> PldaTest::getScores(){
	return _scores;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
Matrix<unsigned long> PldaTest::getTrials(){
	return _trials;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
String PldaTest::getModelName(unsigned long index){
	return(_modelLine.getElement(index));
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
String PldaTest::getSegmentName(unsigned long index){
	return(_segLine.getElement(index));
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaTest::lengthNorm(){

	// Compute the norm of all models
	RealVector<double> vecNorm(0);
	vecNorm.setSize(_n_models);
	vecNorm.setAllValues(0.0);

	for(unsigned long s =0;s<_n_models;s++){
		double tmp = 0.0;
		for(unsigned long k=0;k<_vectSize;k++){
			tmp += _models(k,s)*_models(k,s);
		}
		vecNorm[s] = sqrt(tmp);
	}

	// Divide all vectors by their norm
	for(unsigned long s =0;s<_n_models;s++){
		for(unsigned long k=0;k<_vectSize;k++){
			_models(k,s) /=vecNorm[s];
		}
	}

	// Compute the norm of all segments
	vecNorm.setSize(_n_segments);
	vecNorm.setAllValues(0.0);

	for(unsigned long s =0;s<_n_segments;s++){
		double tmp = 0.0;
		for(unsigned long k=0;k<_vectSize;k++){
			tmp += _segments(k,s)*_segments(k,s);
		}
		vecNorm[s] = sqrt(tmp);
	}

	// Divide all vectors by their norm
	for(unsigned long s =0;s<_n_segments;s++){
		for(unsigned long k=0;k<_vectSize;k++){
			_segments(k,s) /=vecNorm[s];
		}
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaTest::center(RealVector<double> &mu){
	// center all models and test segments
	for(unsigned long k=0;k<_vectSize;k++){
		for(unsigned long s =0;s<_n_models;s++){
			_models(k,s) -= mu[k];
		}
		for(unsigned long s =0;s<_n_segments;s++){
			_segments(k,s) -= mu[k];
		}
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaTest::rotateLeft(Matrix<double> &M){

	Matrix<double> tmpModels(_models);
	Matrix<double> tmpSegments(_segments);

	_models.setAllValues(0.0);
	_segments.setAllValues(0.0);

	// Rotate _models
	for(unsigned long m=0;m<_n_models;m++){
		for(unsigned long i=0;i<M.rows();i++){
			for(unsigned long k=0;k<_vectSize;k++){
				_models(i,m) += M(i,k)*tmpModels(k,m);
			}
		}
	}

	// Rotate _segments
	for(unsigned long s=0;s<_n_segments;s++){
		for(unsigned long i=0;i<M.rows();i++){
			for(unsigned long k=0;k<_vectSize;k++){
				_segments(i,s) += M(i,k)*tmpSegments(k,s);
			}
		}
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void PldaTest::cosineDistance(Config &config){

	RealVector<double> normM(0),normS(0);
	normM.setSize(_n_models);
	normS.setSize(_n_segments);
	normM.setAllValues(0.0);
	normS.setAllValues(0.0);

	// Compute Norm of models and segments
	for(unsigned long k=0;k<_vectSize;k++){
		for(unsigned long m=0;m<_n_models;m++){
			normM[m] += _models(k,m)*_models(k,m);
		}
		for(unsigned long s=0;s<_n_segments;s++){
			normS[s] += _segments(k,s)*_segments(k,s);
		}
	}

	for(unsigned long m=0;m<_n_models;m++){
		normM[m] = sqrt(normM[m]);
	}
	for(unsigned long s=0;s<_n_segments;s++){
		normS[s] = sqrt(normS[s]);
	}

	for(unsigned long m=0;m<_n_models;m++){
		for(unsigned long s=0;s<_n_segments;s++){
			if(_trials(m,s)){
				for(unsigned long k=0;k<_vectSize;k++){
					_scores(m,s) += _models(k,m)*_segments(k,s);
				}
				_scores(m,s) /= (normM[m]*normS[s]);
			}
		}
	}
}















#endif

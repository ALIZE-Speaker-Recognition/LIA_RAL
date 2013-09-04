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

#if !defined(ALIZE_TVAcc_cpp)
#define ALIZE_TVAcc_cpp

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

#include<AccumulateTVStat.h>
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
#ifdef LAPACK
#include "lapacke.h"
#endif

using namespace alize;
using namespace std;

//-----------------------------------------------------------------------------------------
TVAcc::TVAcc(String & featFilename,Config & config)
	:_ms(config),_ss(config){ // constructor for a single file
		
	XList TVNdx;

	if(featFilename.endsWith(".ndx")){
		TVNdx.load(featFilename,config);
	}
	else{
		if(verboseLevel >1)cout<<"Init TVAcc with only one file"<<endl;
		XLine& tmpLine = XLine::create();
		tmpLine.addElement(featFilename);
		TVNdx.addLine()=tmpLine;
	}
	_init(TVNdx,config);
}

//-----------------------------------------------------------------------------------------
TVAcc::TVAcc(XList & ndx,Config & config)
	:_ms(config),_ss(config){ // constructor
	_init(ndx,config);
}

//-----------------------------------------------------------------------------------------
TVAcc::TVAcc(Config & config)
	:_ms(config),_ss(config){ // constructor
	_init(config);
}

//-----------------------------------------------------------------------------------------
TVAcc::~TVAcc(){
	_TETt.deleteAllObjects();
}

//-----------------------------------------------------------------------------------------
String TVAcc::getClassName() const{	return "TVAcc";}

//-----------------------------------------------------------------------------------------
void TVAcc::_init(XList &ndx, Config &config){

	///Convert the NDX file
	_fileList=ndx;
	_ndxTable=TVTranslate(ndx);

	///Load the UBM
	MixtureGD& UBM = _ms.loadMixtureGD(config.getParam("inputWorldFilename"));

	_vectSize = UBM.getVectSize();
	_n_distrib = UBM.getDistribCount();
	_svSize = _vectSize*_n_distrib;

	_rankT=1;
	if(config.existsParam("totalVariabilityNumber"))
		_rankT=config.getParam("totalVariabilityNumber").toULong();

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
	_meanW.setSize(_rankT);
	_meanW.setAllValues(0.0);

	///Create and initialise statistics acumulators
	_statN= Matrix<double>(_n_speakers, _n_distrib);
	_statF= Matrix<double>(_n_speakers, _n_distrib*_vectSize);
	_statN.setAllValues(0.0);
	_statF.setAllValues(0.0);

	_T.setDimensions(_rankT, _svSize);
	_T.setAllValues(0.0);

	_W.setDimensions(_n_speakers, _rankT);
	_W.setAllValues(0.0);

	///Create the accumulators for vEvT computation
	for(unsigned long s =0; s<_n_distrib; s++){
		_TETt.addObject(*new DoubleSquareMatrix(_rankT));
		_TETt[s].setAllValues(0.0);
	}

	///Create accumulators to compute the TotalVariability matrix
	// Initialize to dimension 1 that will be midified in real time during computation
	_A.setDimensions(1,1);
	_A.setAllValues(0.0);

	_C.setDimensions(_rankT,_svSize);
	_C.setAllValues(0.0);

	_R.setSize(_rankT);
	_r.setSize(_rankT);

	_R.setAllValues(0.0);
	_r.setAllValues(0.0);
}



//-----------------------------------------------------------------------------------------
void TVAcc::_init(Config &config){

	///Convert the NDX file
//	_fileList=ndx;
//	_ndxTable=TVTranslate(ndx);

	///Load the UBM
	MixtureGD& UBM = _ms.loadMixtureGD(config.getParam("inputWorldFilename"));

	_vectSize = UBM.getVectSize();
	_n_distrib = UBM.getDistribCount();
	_svSize = _vectSize*_n_distrib;

	_rankT=1;
	if(config.existsParam("totalVariabilityNumber"))
		_rankT=config.getParam("totalVariabilityNumber").toULong();

	///Read NDX file
//	_n_speakers=_fileList.getLineCount();
//	_n_sessions=_fileList.getAllElements().getElementCount();
	_n_speakers = 0;
	_n_sessions = 0;
	
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
	_meanW.setSize(_rankT);
	_meanW.setAllValues(0.0);

	///Create and initialise statistics acumulators
	_statN= Matrix<double>(_n_speakers, _n_distrib);
	_statF= Matrix<double>(_n_speakers, _n_distrib*_vectSize);
	_statN.setAllValues(0.0);
	_statF.setAllValues(0.0);

	_T.setDimensions(_rankT, _svSize);
	_T.setAllValues(0.0);

	_W.setDimensions(_n_speakers, _rankT);
	_W.setAllValues(0.0);

	///Create the accumulators for vEvT computation
	for(unsigned long s =0; s<_n_distrib; s++){
		_TETt.addObject(*new DoubleSquareMatrix(_rankT));
		_TETt[s].setAllValues(0.0);
	}

	///Create accumulators to compute the TotalVariability matrix
	// Initialize to dimension 1 that will be midified in real time during computation
	_A.setDimensions(1,1);
	_A.setAllValues(0.0);

	_C.setDimensions(_rankT,_svSize);
	_C.setAllValues(0.0);

	_R.setSize(_rankT);
	_r.setSize(_rankT);

	_R.setAllValues(0.0);
	_r.setAllValues(0.0);
}

//-----------------------------------------------------------------------------------------
void TVAcc::computeAndAccumulateTVStat(Config& config){
	
	#ifdef THREAD          
	if (config.existsParam("numThread") && config.getParam("numThread").toULong() >0){
		computeAndAccumulateTVStatThreaded(config.getParam("numThread").toULong(),config);				//accumulate stats
	}
	else	computeAndAccumulateTVStatUnThreaded(config); 			//unthreaded version
	#else
		computeAndAccumulateTVStatUnThreaded(config);			//accumute stats
	#endif
}

//-----------------------------------------------------------------------------------------
void TVAcc::computeAndAccumulateTVStatUnThreaded(Config& config){

	if (verboseLevel >= 1) cout << "(AccumulateTVStat) Compute Statistics UnThreaded"<<endl;

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
	double *n, *f_x;
	 n=_statN.getArray();f_x=_statF.getArray();
	
	String currentSource="";unsigned long loc=0;unsigned long session=0;
	while((seg=selectedSegments.getSeg())!=NULL){

		unsigned long begin=seg->begin()+fs.getFirstFeatureIndexOfASource(seg->sourceName()); 				/// Idx of the first frame of the current file in the feature server

		TVTranslate ndxTable=TVTranslate(_fileList);
		
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
				n[loc*_n_distrib+k]   +=aPost[k];
				for (unsigned long i=0;i<_vectSize;i++) {
					f_x[loc*_svSize+(k*_vectSize+i)]   +=aPost[k]*ff[i];
				}
			}
		}
	}	
}


#ifdef THREAD
//-----------------------------------------------------------------------------------------
//				Data strucutre of thread
//-----------------------------------------------------------------------------------------
struct StatTVthread_data{

	double *N;
	double *FX;
	unsigned long firstLine;
	unsigned long lastLine;
	unsigned long svSize;
	unsigned long vectSize;
	unsigned long threadNb;
	unsigned long n_distrib;
	RefVector<MixtureGD> *world;
	RefVector<XList> *fileList;
	Config *config;
};

//-----------------------------------------------------------------------------------------
//				Thread Routine
//-----------------------------------------------------------------------------------------
void *StatTVthread(void *threadarg) {
	struct StatTVthread_data *my_data;
	my_data = (struct StatTVthread_data *) threadarg;
	
	double *N = my_data->N;
	double *F_X = my_data->FX;
	unsigned long firstLine = my_data->firstLine;
	unsigned long lastLine = my_data->lastLine;
	unsigned long svSize = my_data->svSize;
	unsigned long vectSize = my_data->vectSize;
	Config *config = my_data->config;
	unsigned long threadNb = my_data->threadNb;
	unsigned long n_distrib = my_data->n_distrib;
	MixtureGD &world=(*(my_data->world))[threadNb];
	XList &fileList=(*(my_data->fileList))[threadNb];

	StatServer _ss(*config);
	MixtureGDStat &acc=_ss.createAndStoreMixtureStat(world);

	//Create a TVTranslate object for the complete XList
	TVTranslate ndxTable(fileList);

	//Create a temporary XLinewith the current speakers
	XLine currentList;
	for(unsigned long spk=firstLine;spk<lastLine;spk++){
		for(unsigned long i=0;i<fileList.getLine(spk).getElementCount();i++){
			currentList.addElement(fileList.getLine(spk).getElement(i));
		}
	}

	///Create and initialise the feature server
	FeatureServer fs;
	fs.init(*config, currentList);

	///Create and initialise feature clusters
	SegServer segmentsServer;
	LabelServer labelServer;
	initializeClusters(currentList,segmentsServer,labelServer,*config);

	verifyClusterFile(segmentsServer,fs,*config);
	unsigned long codeSelectedFrame=labelServer.getLabelIndexByString(config->getParam("labelSelectedFrames"));

	SegCluster& selectedSegments=segmentsServer.getCluster(codeSelectedFrame);
	acc.resetOcc();
	Seg *seg; 
	selectedSegments.rewind();
	String currentSource="";unsigned long loc=0;unsigned long session=0;
	while((seg=selectedSegments.getSeg())!=NULL){

		// Idx of the first frame of the current file in the feature server
		unsigned long begin=seg->begin()+fs.getFirstFeatureIndexOfASource(seg->sourceName());
		
		if (currentSource!=seg->sourceName()) {
			currentSource=seg->sourceName();
			loc=ndxTable.locNb(currentSource);
		}

		fs.seekFeature(begin);
		Feature f;

		for (unsigned long idxFrame=0;idxFrame<seg->length();idxFrame++){
			fs.readFeature(f);
			acc.computeAndAccumulateOcc(f);
			DoubleVector aPost=acc.getOccVect();
			double *ff=f.getDataVector();

			for(unsigned long k=0;k<n_distrib;k++) {
				N[loc*n_distrib+k]   +=aPost[k];
				for (unsigned long i=0;i<vectSize;i++) {
					F_X[loc*svSize+(k*vectSize+i)]   +=aPost[k]*ff[i];
				}
			}
		}
	}
	pthread_exit((void*) 0);
	return (void*)0 ;
}

//-----------------------------------------------------------------------------------------
void TVAcc::computeAndAccumulateTVStatThreaded(unsigned long NUM_THREADS, Config &config){

	if (verboseLevel >= 1) cout << "(AccumulateTVStat) Compute Statistics Threaded"<<endl;
	if (NUM_THREADS==0) throw Exception("Num threads can not be 0",__FILE__,__LINE__);

	int rc, status;
	if (NUM_THREADS > _fileList.getLineCount()) NUM_THREADS=_fileList.getLineCount();

	MixtureGD& UBM = _ms.getMixtureGD(0);

	RefVector<MixtureGD> mGd;
	RefVector<XList> fileList;
	for(unsigned long t=0;t<NUM_THREADS;t++){
		mGd.addObject(*new MixtureGD(UBM));
		fileList.addObject(*new XList(_fileList));
	}

	double *N=_statN.getArray();
	double *F_X=_statF.getArray();

	//Compute index of the first and last speaker to process for each thread
	RealVector<unsigned long> firstLine; firstLine.setSize(NUM_THREADS); firstLine.setAllValues(0);
	RealVector<unsigned long> lastLine; lastLine.setSize(NUM_THREADS); lastLine.setAllValues(0);
	unsigned long nbLinePerThread = (unsigned long)floor((double)_fileList.getLineCount()/(double)NUM_THREADS);
	for(unsigned long i=0;i<NUM_THREADS-1;i++){
		firstLine[i] = i*nbLinePerThread;
		lastLine[i] = (i+1)*nbLinePerThread;
	}
	firstLine[NUM_THREADS-1] = (NUM_THREADS-1)*nbLinePerThread;
	lastLine[NUM_THREADS-1] = _fileList.getLineCount();

	struct StatTVthread_data *thread_data_array = new StatTVthread_data[NUM_THREADS];
	pthread_t *threads = new pthread_t[NUM_THREADS];

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	//Create threads
	for(unsigned long t=0; t<NUM_THREADS; t++){

		thread_data_array[t].N=N;
		thread_data_array[t].FX=F_X;
		thread_data_array[t].firstLine = firstLine[t];
		thread_data_array[t].lastLine = lastLine[t];
		thread_data_array[t].svSize = _svSize;
		thread_data_array[t].vectSize = _vectSize;
		thread_data_array[t].n_distrib = _n_distrib;
		thread_data_array[t].threadNb = t;
		thread_data_array[t].world = &(mGd);
		thread_data_array[t].fileList = &(fileList);
		thread_data_array[t].config = &config;

		if (verboseLevel>1) cout<<"(AccumulateTVStat) Creating thread n["<< t<< "] for speakers["<<firstLine[t]<<"-->"<<lastLine[t]-1<<"]"<<endl;
		rc = pthread_create(&threads[t], &attr, StatTVthread, (void *)&thread_data_array[t]);
		if (rc) throw Exception("ERROR; return code from pthread_create() is ",__FILE__,rc);
	}

	pthread_attr_destroy(&attr);
	for(unsigned long t=0; t<NUM_THREADS; t++) {
		rc = pthread_join(threads[t], (void **)&status);
		if (rc)  throw Exception("ERROR; return code from pthread_join() is ",__FILE__,rc);
		if (verboseLevel >1) cout <<"(AccumulateTVStat) Completed join with thread ["<<t<<"] status["<<status<<"]"<<endl;
	}

	free(thread_data_array);
	free(threads);
}
#endif


//-----------------------------------------------------------------------------------------
void TVAcc::computeAndAccumulateTVStat(FeatureServer &fs,Config & config){
	SegServer segmentsServer;
	LabelServer labelServer;
	initializeClusters(_fileList,segmentsServer,labelServer,config);
	verifyClusterFile(segmentsServer,fs,config);
	unsigned long codeSelectedFrame=labelServer.getLabelIndexByString(config.getParam("labelSelectedFrames"));
	SegCluster& selectedSegments=segmentsServer.getCluster(codeSelectedFrame); 
	this->computeAndAccumulateTVStat(selectedSegments,fs,config);
};

//-----------------------------------------------------------------------------------------
void TVAcc::computeAndAccumulateTVStat(SegCluster &selectedSegments,FeatureServer &fs,Config & config){

	MixtureGD& UBM = _ms.getMixtureGD(0);
	MixtureGDStat &acc=_ss.createAndStoreMixtureStat(UBM);

	///Fast access to vector and matrices array
	double *n, *f_x;
	n=_statN.getArray(); f_x=_statF.getArray();

	///Compute Occupations and Statistics	
	acc.resetOcc();
	Seg *seg; 
	selectedSegments.rewind();
	
	String currentSource="";unsigned long loc=0;unsigned long session=0;
	while((seg=selectedSegments.getSeg())!=NULL){

		unsigned long begin=seg->begin()+fs.getFirstFeatureIndexOfASource(seg->sourceName()); 				/// Idx of the first frame of the current file in the feature server

		TVTranslate ndxTable=TVTranslate(_fileList);
		
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
				n[loc*_n_distrib+k]   +=aPost[k];
				for (unsigned long i=0;i<_vectSize;i++) {
					f_x[loc*_svSize+(k*_vectSize+i)]   +=aPost[k]*ff[i];
				}
			}
		}
	}	
}

//-----------------------------------------------------------------------------------------
void TVAcc::resetAcc(){
	_statF.setAllValues(0.0);
	_statN.setAllValues(0.0);	
	if (verboseLevel >= 1) cout << "# TV Accumulators reset" << endl;
}

//-----------------------------------------------------------------------------------------
void TVAcc::resetTmpAcc(){
	
	_C.setAllValues(0.0);
	_A.setAllValues(0.0);

	///Reinitialise accumulators for L matrices computation
	for(unsigned long s =0; s<_n_distrib; s++){
		_TETt[s].setAllValues(0.0);
	}
}

//-----------------------------------------------------------------------------------------
void TVAcc::loadT(const String& tvFilename, Config& config){	///load an TotalVariability Matrix
	String filename = config.getParam("matrixFilesPath") + tvFilename +  config.getParam("loadMatrixFilesExtension");
	_T.load (filename, config);

	//Transpose the matrix if the number of line is higher than the numbe of columns
	if(_T.cols()<_T.rows()){
		cout<<"Load TV : number of lines ( "<<_T.rows() <<" ) higher than the number of columns ( "<<_T.cols()<<" ("<<endl;
		_T.transpose();
	}
	_rankT=_T.rows();

	if(_rankT != config.getParam("totalVariabilityNumber").toULong()){
		throw Exception("Incorrect dimension of TotalVariability Matrix",__FILE__,__LINE__);
	}

	cout << "(AccumulateTVStat) Init TV matrix from "<< filename <<"  for TotalVariability Matrix, rank: ["<<_T.rows() << "] sv size: [" << _T.cols() <<"]"<<endl;
}

//-----------------------------------------------------------------------------------------
void TVAcc::loadT(Matrix<double> & V, Config& config){
	_rankT=V.rows();
	_T=V;

	//Transpose the matrix if the number of line is higher than the numbe of columns
	if(_T.cols()<_T.rows()){
		cout<<"Load T : number of lines higher than the number of columns"<<endl;
		_T.transpose();
	}

	unsigned long rankEV = 1;
	if(config.existsParam("totalVariabilityNumber"))
		rankEV = config.getParam("totalVariabilityNumber").toULong();

	if(_rankT != config.getParam("totalVariabilityNumber").toULong()){
		throw Exception("Incorrect dimension of TotalVariability Matrix",__FILE__,__LINE__);
	}
}

//-----------------------------------------------------------------------------------------
void TVAcc::loadMeanEstimate(DoubleVector& meanEstimate){
	_ubm_means = meanEstimate;
	if(_ubm_means.size() != _svSize){
		throw Exception("Incorrect dimension of meanEstimate vector",__FILE__,__LINE__);
	}
}

//-----------------------------------------------------------------------------------------
void TVAcc::loadN(Config& config){
	String filname=config.getParam("matrixFilesPath") + config.getParam("nullOrderStatSpeaker") + config.getParam("loadMatrixFilesExtension");
	_statN.load (filname, config);

	if((_statN.rows() != _fileList.getLineCount()) || (_statN.cols() != _n_distrib)){
		throw Exception("Incorrect dimension of N Matrix",__FILE__,__LINE__);
	}
	cout<<"(AccumulateTVStat) ---------load statistics N [ "<<_statN.rows()<<" ] [ "<<_statN.cols()<<" ]"<<endl;
}

//-----------------------------------------------------------------------------------------
void TVAcc::loadF_X(Config& config){
	String filname=config.getParam("matrixFilesPath") + config.getParam("firstOrderStatSpeaker") + config.getParam("loadMatrixFilesExtension");
	_statF.load (filname, config);

	if((_statF.rows() != _fileList.getLineCount()) || (_statF.cols() != _svSize)){
		throw Exception("Incorrect dimension of F_X Matrix",__FILE__,__LINE__);
	}
	cout<<"(AccumulateTVStat) ---------load statistics F_X [ "<<_statF.rows()<<" ] [ "<<_statF.cols()<<" ]"<<endl;
}

//-----------------------------------------------------------------------------------------
void TVAcc::initT(Config& config){		///random initialisation of the TotalVariability Matrix
		_rankT=config.getParam("totalVariabilityNumber").toULong();
		_T.setDimensions(_rankT,_svSize);

		//Two type of random initializations
		//	normal  : random with a normal law by using a Box-Muller Generator
		//	uniform : random with a uniform distribution
try{
		String randomLaw = "normal";
		if(config.existsParam("randomInitLaw"))	randomLaw = config.getParam("randomInitLaw");

		//Initialize the matrix by generating random values following a uniform law
		if(config.getParam("randomInitLaw") == "uniform"){

			srand48(_svSize*_rankT);
			_T.randomInit();

			double norm=0;
			for(unsigned long k=0; k<_svSize; k++){
				norm += _ubm_invvar[k];
			}
			norm = norm/_svSize;
			for(unsigned long i=0; i<_T.rows(); i++){
				for(unsigned long j=0; j<_T.cols(); j++){
					_T(i,j) = _T(i,j)*norm;
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
			for(unsigned long i=0; i<_T.rows(); i++){
				for(unsigned long j=0; j<_T.cols(); j++){
					double val = boxMullerGenerator(0.0, 1.0);
					while((ISNAN(val)) || (ISINF(val))){
						val = boxMullerGenerator(0.0, 1.0);
					}
					_T(i,j) = val * norm * 0.001;
				}
			}
		}
		else{
			throw Exception("Selected random initialization law does not exist",__FILE__,__LINE__);
		}
		if (verboseLevel >=1) cout << "(AccumulateTVStat) Random Init for TotalVariability Matrix with "<<randomLaw<<" law "<<", rank: ["<<_T.rows() << "] sv size: [" << _T.cols() <<"]"<<endl;
}
catch (Exception& e){ 
	cout << e.toString().c_str() << endl;
}
}

//-----------------------------------------------------------------------------------------
void TVAcc::saveT(const String& filename, Config& config){
	String vName = config.getParam("matrixFilesPath") + filename +  config.getParam("saveMatrixFilesExtension");
	_T.save(vName, config);
}

//-----------------------------------------------------------------------------------------
void TVAcc::estimateTETt(Config &config){
	
	#ifdef THREAD          
	if (config.existsParam("numThread") && config.getParam("numThread").toULong() >0)	estimateTETtThreaded(config.getParam("numThread").toULong());
	else estimateTETtUnThreaded();
	#else
	estimateTETtUnThreaded();
	#endif
}

//-----------------------------------------------------------------------------------------
void TVAcc::estimateTETtUnThreaded(){	
        if (verboseLevel >= 1) cout << "(AccumulateTVStat) compute V * Sigma-1 * V "<<endl;

	for(unsigned long d=0; d<_n_distrib; d++){

		Matrix<double>ssV= _T.crop(0,d*_vectSize, _rankT,_vectSize);

		///Compute  _TETt matrices
		double *vEvT, *v, *E;
		_TETt[d].setAllValues(0.0);
		vEvT= _TETt[d].getArray();
		v= ssV.getArray();
		E= _ubm_invvar.getArray();

		for(unsigned long i = 0; i<_rankT; i++){
			for(unsigned long j = 0; j<=i; j++){
				for(unsigned long k=0; k<_vectSize;k++){
					vEvT[i*_rankT+j] += v[i*_vectSize+k] * E[d*_vectSize+k] * v[j*_vectSize+k];
				}
			}
		}
		
		for(unsigned long i=0;i<_rankT;i++){
			for(unsigned long j=i+1;j<_rankT;j++) {
				vEvT[i*_rankT+j] = vEvT[j*_rankT+i];
			}
		}
	}
}

#ifdef THREAD
//-----------------------------------------------------------------------------------------
//				Data strucutre of thread
//-----------------------------------------------------------------------------------------
struct TETthread_data{

	double *V;
	double *ubm_invvar;
	unsigned long disBottom;
	unsigned long disUp;	
	unsigned long rankEV;
	unsigned long svSize;
	unsigned long vectSize;
	RefVector <DoubleSquareMatrix>* vevT;
};

//-----------------------------------------------------------------------------------------
//				Thread Routine
//-----------------------------------------------------------------------------------------
void *TETthread(void *threadarg) {
	struct TETthread_data *my_data;
	my_data = (struct TETthread_data *) threadarg;
	
	double *v = my_data->V;
	unsigned long disBottom = my_data->disBottom;	
	unsigned long disUp = my_data->disUp;
	double *E=my_data->ubm_invvar;
	unsigned long _rankT=my_data->rankEV;
	unsigned long _svSize=my_data->svSize;
	unsigned long _vectSize=my_data->vectSize;

	for (unsigned long d=disBottom;d <disUp;d++){
		
		DoubleSquareMatrix &vEvT=(*(my_data->vevT))[d];
		vEvT.setAllValues(0.0);
		double *vevt = vEvT.getArray();
		
		for(unsigned long i = 0; i<_rankT; i++){
			for(unsigned long j = 0; j<=i; j++){
				for(unsigned long k=0; k<_vectSize;k++){
					vevt[i*_rankT+j] += v[(i*_svSize)+(d*_vectSize)+k] * E[d*_vectSize+k] * v[(j*_svSize)+(d*_vectSize)+k];
				}
			}
		}

		for(unsigned long i=0;i<_rankT;i++){
			for(unsigned long j=i+1;j<_rankT;j++) {
				vevt[i*_rankT+j]=vevt[j*_rankT+i];
			}
		}
	}

	pthread_exit((void*) 0);
	return (void*)0 ;
}

//-----------------------------------------------------------------------------------------
void TVAcc::estimateTETtThreaded(unsigned long NUM_THREADS){
	
	if (verboseLevel >=1) cout << "(AccumulateTVStat) Compute TETt Threaded"<<endl;
	if (NUM_THREADS==0) throw Exception("Num threads can not be 0",__FILE__,__LINE__);
	
	double *V=_T.getArray();
	double *ubm_invvar=_ubm_invvar.getArray();
	
	int rc, status;
	if (NUM_THREADS > _n_distrib) NUM_THREADS=_n_distrib;
	
	struct TETthread_data *thread_data_array = new TETthread_data[NUM_THREADS];
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
		thread_data_array[t].rankEV=_rankT;
		thread_data_array[t].svSize=_svSize;
		thread_data_array[t].vectSize=_vectSize;
		thread_data_array[t].vevT=&(_TETt);

		if (verboseLevel > 1) cout<<"(AccumulateTVStat) Creating thread n["<< t<< "] for distributions["<<disBottom<<"-->"<<disUp<<"]"<<endl;
		rc = pthread_create(&threads[t], &attr, TETthread, (void *)&thread_data_array[t]);
		if (rc) throw Exception("ERROR; return code from pthread_create() is ",__FILE__,rc);

		disBottom = disUp;
	}

	pthread_attr_destroy(&attr);

	for(unsigned long t=0; t<NUM_THREADS; t++) {
		rc = pthread_join(threads[t], (void **)&status);
		if (rc)  throw Exception("ERROR; return code from pthread_join() is ",__FILE__,rc);
		if (verboseLevel >1) cout <<"(AccumulateTVStat) Completed join with thread ["<<t<<"] status["<<status<<"]"<<endl;
	}


	free(thread_data_array);
	free(threads);

	if (verboseLevel >= 1) cout << "(AccumulateTVStat) Done " << endl;
}
#endif

//-----------------------------------------------------------------------------------------
void  TVAcc::getTW(DoubleVector &vy, String &file){		///Compute VY for the speaker corresponding to the given file given file
	vy.setAllValues(0.0);	
	
	double *v, *y; y=_W.getArray(); v=_T.getArray();

	unsigned long idx=_ndxTable.locNb(file);
	for (unsigned long i=0;i<_svSize;i++){
		for (unsigned long j=0;j<_rankT;j++){
			vy[i] += v[j*_svSize + i] * y[idx*_rankT + j ];
		}
	}	
}

//-----------------------------------------------------------------------------------------
void  TVAcc::getTW(DoubleVector &vy, unsigned long spk) {		///Compute VY for speaker spk
	vy.setAllValues(0.0);
	
	double *v, *y; y=_W.getArray(); v=_T.getArray();
	for (unsigned long i=0;i<_svSize;i++){
		for (unsigned long j=0;j<_rankT;j++){
			vy[i] += v[j*_svSize + i] * y[spk*_rankT + j ];
		}
	}
}

//-----------------------------------------------------------------------------------------
// Compute supervector of client M_s_h=M+ VY_s
void TVAcc::getMplusTW(DoubleVector &Sp, String& file){
	Sp.setAllValues(0.0);
	DoubleVector vy(_svSize,_svSize);
	getTW(vy,file);		
	for (unsigned long i=0;i<_svSize;i++){
		Sp[i] = _ubm_means[i] + vy[i];
	}
}

//-----------------------------------------------------------------------------------------
// Compute supervector of client M_s=M+ TW
void TVAcc::getMplusTW(DoubleVector &Sp, unsigned long spk){
	Sp.setAllValues(0.0);
	DoubleVector vy(_svSize,_svSize);
	getTW(vy,spk);
	for (unsigned long i=0;i<_svSize;i++){
		Sp[i]=_ubm_means[i] + vy[i];
	}
}

//-----------------------------------------------------------------------------------------
void TVAcc::updateTestimate(){

	if(verboseLevel>0) cout << "(AccumulateTVStat) Update T Matrix"<<endl;

	Matrix<double> tmpC; tmpC.setDimensions(_rankT, _svSize); tmpC.setAllValues(0.0);

	DoubleSquareMatrix tmpA(_rankT);
	for(unsigned long d=0;d<_n_distrib;d++){
		DoubleSquareMatrix invA;
		invA.setSize(_rankT);
	
		for(unsigned long i=0;i<_rankT;i++)
			for(unsigned long j=0;j<_rankT;j++)
				tmpA(i,j) = _A(d, i*(_rankT) +j);
		tmpA.invert(invA);

		double *tmpc, * inva, *cev;
		tmpc=tmpC.getArray(); inva=invA.getArray(); cev=_C.getArray();
	
		for(unsigned long i=0; i<_rankT;i++){
			for(unsigned long j=0; j<_vectSize;j++){
				for(unsigned long k=0; k<_rankT;k++){
					tmpC(i,d*_vectSize+j) += invA(i,k) *  _C(k,d*_vectSize+j);
				}
			}
		}
	}
	_C=tmpC;
	_T=_C;

	if (verboseLevel>0) cout << "(AccumulateTVStat) Done " << endl;
}

//-----------------------------------------------------------------------------------------
XList& TVAcc::getXList(){
	return(_fileList);
}

//-----------------------------------------------------------------------------------------
unsigned long TVAcc::getNSpeakers(){
	return(_n_speakers);
}

//-----------------------------------------------------------------------------------------
unsigned long TVAcc::getNDistrib(){
	return(_n_distrib);
}

//-----------------------------------------------------------------------------------------
unsigned long TVAcc::getVectSize(){
	return(_vectSize);
}

//-----------------------------------------------------------------------------------------
unsigned long TVAcc::getSvSize(){
	return(_svSize);
}

//-----------------------------------------------------------------------------------------
unsigned long TVAcc::getRankT(){
	return(_rankT);
}

//-----------------------------------------------------------------------------------------
DoubleVector& TVAcc::getUbmMeans(){
	return(_ubm_means);
}

//-----------------------------------------------------------------------------------------
DoubleVector& TVAcc::getUbmInvVar(){
	return(_ubm_invvar);
}

//-----------------------------------------------------------------------------------------
Matrix<double> TVAcc::getT(){
	return(_T);
}

//-----------------------------------------------------------------------------------------
Matrix<double> TVAcc::getW(){
	return(_W);
}

//-----------------------------------------------------------------------------------------
void TVAcc::saveW(String yFile,Config &config){
	_W.save(yFile,config);
}

//-----------------------------------------------------------------------------------------
Matrix <double>& TVAcc::getN() {
	return _statN;
}	
	
//-----------------------------------------------------------------------------------------
Matrix <double>& TVAcc::getF() {
	return _statF;
}	

//-----------------------------------------------------------------------------------------
DoubleSquareMatrix& TVAcc::getTETt(unsigned long idx){
	return _TETt[idx];
}

//-----------------------------------------------------------------------------------------
void TVAcc::substractM(Config & config){
	#ifdef THREAD          
	if (config.existsParam("numThread") && config.getParam("numThread").toULong() >0)	substractMThreaded(config.getParam("numThread").toULong());
	else substractMUnThreaded();
	#else
	substractMUnThreaded();
	#endif
}

//-----------------------------------------------------------------------------------------
void TVAcc::substractMUnThreaded(){
	
	if (verboseLevel >= 1) cout <<"(AccumulateTVStat) Substract Mean from Speaker Statistics" << endl;	

	double *F_X = _statF.getArray();
	
	for(unsigned long spk=0; spk<_n_speakers; spk++){

		double *n=_statN.getArray();
		
		//Subtract M
		for(unsigned long i=0; i<_n_distrib; i++){
			for(unsigned long j = 0; j< _vectSize;j++){
				F_X[spk*_svSize+i*_vectSize+j] -= _ubm_means[i*_vectSize+j]*n[spk*_n_distrib+i];
			}
		}
	}
}

#ifdef THREAD
//-----------------------------------------------------------------------------------------
//				Data strucutre of thread
//-----------------------------------------------------------------------------------------
struct S_Mthread_data{
	double *N;
	double *F_X;
	double *ubm_means;
	unsigned long spkBottom;
	unsigned long spkUp;	
	unsigned long vectSize;
	unsigned long svSize;
	unsigned long n_distrib;
	unsigned long nt;
};

//-----------------------------------------------------------------------------------------
//				Thread Routine
//-----------------------------------------------------------------------------------------
void *S_Mthread(void *threadarg) {
	struct S_Mthread_data *my_data;
	my_data = (struct S_Mthread_data *) threadarg;

	double *n = my_data->N;
	double *f_x = my_data->F_X;

	double *ubm_means = my_data->ubm_means;
	unsigned long spkBottom = my_data->spkBottom;	
	unsigned long spkUp = my_data->spkUp;
	unsigned long _vectSize=my_data->vectSize;
	unsigned long _svSize=my_data->svSize;
	unsigned long _n_distrib=my_data->n_distrib;

	for(unsigned long spk=spkBottom; spk<spkUp; spk++){
		//Substract M
		for(unsigned long i=0; i<_n_distrib; i++){
			for(unsigned long j = 0; j< _vectSize;j++){
				f_x[spk*_svSize+i*_vectSize+j] -= ubm_means[i*_vectSize+j]*n[spk*_n_distrib+i];
			}
		}
	}
	pthread_exit((void*) 0);
	return (void*)0 ;
}

//-----------------------------------------------------------------------------------------
void TVAcc::substractMThreaded(unsigned long NUM_THREADS){
	
	if (verboseLevel >= 1) cout <<"((AccumulateTVStat) Substract Mean from Speaker Statistics Threaded" << endl;
	if (NUM_THREADS==0) throw Exception("Num threads can not be 0",__FILE__,__LINE__);

	int rc, status;
	if (NUM_THREADS > _n_speakers) NUM_THREADS=_n_speakers;

	struct S_Mthread_data *thread_data_array = new S_Mthread_data[NUM_THREADS];
	pthread_t *threads = new pthread_t[NUM_THREADS];

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	unsigned long offset=_n_speakers/NUM_THREADS;
	
	double *N=_statN.getArray(); 
	double *F_X=_statF.getArray();

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

		thread_data_array[t].ubm_means=ubm_means;
		thread_data_array[t].spkBottom=spkBottom;
		thread_data_array[t].spkUp=spkUp;
		thread_data_array[t].vectSize=_vectSize;
		thread_data_array[t].svSize=_svSize;
		thread_data_array[t].n_distrib=_n_distrib;

		if (verboseLevel>1) cout<<"(AccumulateTVStat) Creating thread n["<< t<< "] for speakers ["<<spkBottom<<"-->"<<spkUp-1<<"]"<<endl;
		rc = pthread_create(&threads[t], &attr, S_Mthread, (void *)&thread_data_array[t]);
		if (rc) throw Exception("ERROR; return code from pthread_create() is ",__FILE__,rc);
		
		spkBottom = spkUp;
	}
	
	pthread_attr_destroy(&attr);
	for(unsigned long t=0; t<NUM_THREADS; t++) {
		rc = pthread_join(threads[t], (void **)&status);
		if (rc)  throw Exception("ERROR; return code from pthread_join() is ",__FILE__,rc);
		if (verboseLevel>1) cout <<"(AccumulateTVStat) Completed join with thread ["<<t<<"] status["<<status<<"]"<<endl;
	}

	free(thread_data_array);
	free(threads);

	if (verboseLevel >= 1) cout << "(AccumulateTVStat) Done " << endl;
}

#endif

//-----------------------------------------------------------------------------------------
void TVAcc::normStatistics(Config & config){
	#ifdef THREAD          
	if (config.existsParam("numThread") && config.getParam("numThread").toULong() >0)	normStatisticsThreaded(config.getParam("numThread").toULong());
	else normStatisticsUnThreaded();
	#else
	normStatisticsUnThreaded();
	#endif
}

//-----------------------------------------------------------------------------------------
void TVAcc::normStatisticsUnThreaded(){	// Threaded version not implemented yet

	if (verboseLevel >= 1) cout <<"(AccumulateTVStat) Normalize Statistics" << endl;	

	double *F_X = _statF.getArray();
	double *n=_statN.getArray();

	//Subtract mean and normalize by using the variance of the UBM
	for(unsigned long i=0; i<_n_distrib; i++){
		for(unsigned long j = 0; j< _vectSize;j++){
			double sqrtInvCov = sqrt(_ubm_invvar[i*_vectSize+j]);
			for(unsigned long spk=0; spk<_n_speakers; spk++){
				F_X[spk*_svSize+i*_vectSize+j] -= _ubm_means[i*_vectSize+j]*n[spk*_n_distrib+i];
				F_X[spk*_svSize+i*_vectSize+j] *= sqrtInvCov;
			}
		}
	}
}

#ifdef THREAD
//-----------------------------------------------------------------------------------------
//				Data strucutre of thread
//-----------------------------------------------------------------------------------------
struct normStat_thread_data{
	double *N;
	double *F_X;
	double *ubm_means;
	double *ubm_invvar;
	unsigned long spkBottom;
	unsigned long spkUp;	
	unsigned long vectSize;
	unsigned long svSize;
	unsigned long n_distrib;
	unsigned long nt;
};

//-----------------------------------------------------------------------------------------
//				Thread Routine
//-----------------------------------------------------------------------------------------
void *normStatthread(void *threadarg) {
	
	struct normStat_thread_data *my_data;
	my_data = (struct normStat_thread_data *) threadarg;

	double *n = my_data->N;
	double *f_x = my_data->F_X;

	double *ubm_means = my_data->ubm_means;
	double *ubm_invvar = my_data->ubm_invvar;
	unsigned long spkBottom = my_data->spkBottom;	
	unsigned long spkUp = my_data->spkUp;
	unsigned long _vectSize=my_data->vectSize;
	unsigned long _svSize=my_data->svSize;
	unsigned long _n_distrib=my_data->n_distrib;


//	for(unsigned long spk=spkBottom; spk<spkUp; spk++){
//		//Subtract mean and normalize by using the variance of the UBM
//		for(unsigned long i=0; i<_n_distrib; i++){
//			for(unsigned long j = 0; j< _vectSize;j++){
////				f_x[spk*_svSize+i*_vectSize+j] -= ubm_means[i*_vectSize+j]*n[spk*_n_distrib+i];
//			}
//		}
//	}

	for(unsigned long i=0; i<_n_distrib; i++){
		for(unsigned long j = 0; j< _vectSize;j++){
			double sqrtInvCov = sqrt(ubm_invvar[i*_vectSize+j]);
			for(unsigned long spk=spkBottom; spk<spkUp; spk++){
				f_x[spk*_svSize+i*_vectSize+j] -= ubm_means[i*_vectSize+j]*n[spk*_n_distrib+i];
				f_x[spk*_svSize+i*_vectSize+j] *= sqrtInvCov;
			}
		}
	}

	pthread_exit((void*) 0);
	return (void*)0 ;
}

//-----------------------------------------------------------------------------------------
void TVAcc::normStatisticsThreaded(unsigned long NUM_THREADS){
	
	if (verboseLevel >= 1) cout <<"(AccumulateTVStat) Normalize Statistics Threaded" << endl;
	if (NUM_THREADS==0) throw Exception("Num threads can not be 0",__FILE__,__LINE__);

	int rc, status;
	if (NUM_THREADS > _n_speakers) NUM_THREADS=_n_speakers;

	struct normStat_thread_data *thread_data_array = new normStat_thread_data[NUM_THREADS];
	pthread_t *threads = new pthread_t[NUM_THREADS];

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	unsigned long offset=_n_speakers/NUM_THREADS;
	
	double *N=_statN.getArray(); 
	double *F_X=_statF.getArray();

	double *ubm_means=_ubm_means.getArray();
	double *ubm_invvar=_ubm_invvar.getArray();
	unsigned long spkBottom = 0;
	unsigned long spkUp=0;
	unsigned long re=_n_speakers - NUM_THREADS*offset;
	
	//Create threads
	for(unsigned long t=0; t<NUM_THREADS; t++){
		spkUp = spkBottom +offset;
		if(t<re) spkUp +=1;

		thread_data_array[t].N=N;
		thread_data_array[t].F_X=F_X;

		thread_data_array[t].ubm_means=ubm_means;
		thread_data_array[t].ubm_invvar=ubm_invvar;
		thread_data_array[t].spkBottom=spkBottom;
		thread_data_array[t].spkUp=spkUp;
		thread_data_array[t].vectSize=_vectSize;
		thread_data_array[t].svSize=_svSize;
		thread_data_array[t].n_distrib=_n_distrib;

		if (verboseLevel>1) cout<<"(AccumulateTVStat) Creating thread n["<< t<< "] for speakers ["<<spkBottom<<"-->"<<spkUp-1<<"]"<<endl;
		rc = pthread_create(&threads[t], &attr, normStatthread, (void *)&thread_data_array[t]);
		if (rc) throw Exception("ERROR; return code from pthread_create() is ",__FILE__,rc);
		
		spkBottom = spkUp;
	}
	
	pthread_attr_destroy(&attr);
	for(unsigned long t=0; t<NUM_THREADS; t++) {
		rc = pthread_join(threads[t], (void **)&status);
		if (rc)  throw Exception("ERROR; return code from pthread_join() is ",__FILE__,rc);
		if (verboseLevel>1) cout <<"(AccumulateTVStat) Completed join with thread ["<<t<<"] status["<<status<<"]"<<endl;
	}

	free(thread_data_array);
	free(threads);

	if (verboseLevel >= 1) cout << "(AccumulateTVStat) Done " << endl;
}

#endif

//-----------------------------------------------------------------------------------------
void TVAcc::substractMplusTW(Config &config){
	#ifdef THREAD          
	if (config.existsParam("numThread") && config.getParam("numThread").toULong() >0)	substractMplusTWThreaded(config.getParam("numThread").toULong());
	else substractMplusTWUnThreaded();
	#else
	substractMplusTWUnThreaded();
	#endif
}

//-----------------------------------------------------------------------------------------
void TVAcc::substractMplusTWUnThreaded(){
	
	if (verboseLevel >= 1) cout <<"(AccumulateTVStat) Substract Mean + VY from Speaker Statistics" << endl;	

	double *f_x = _statF.getArray();
	
	for(unsigned long spk=0; spk<_n_speakers; spk++){
		//compute M+VY for the current speaker
		DoubleVector MplusVY; MplusVY.setSize(_svSize); MplusVY.setAllValues(0.0);
		double *mplusvy=MplusVY.getArray();
		this->getMplusTW(MplusVY,spk);
		double *n=_statN.getArray();

		//substract M+VY
		for(unsigned long i=0; i<_n_distrib; i++){
			for(unsigned long j = 0; j< _vectSize;j++){
				f_x[spk*_svSize+i*_vectSize+j] -= mplusvy[i*_vectSize+j]*n[spk*_n_distrib+i];
			}
		}
	}
}

#ifdef THREAD
//-----------------------------------------------------------------------------------------
//				Data strucutre of thread
//-----------------------------------------------------------------------------------------
struct S_MplusTWthread_data{
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

//-----------------------------------------------------------------------------------------
//				Thread Routine
//-----------------------------------------------------------------------------------------
void *S_MplusTWthread(void *threadarg) {
	
	struct S_MplusTWthread_data *my_data;
	my_data = (struct S_MplusTWthread_data *) threadarg;

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
	unsigned long _rankT=my_data->rankEV;

	DoubleVector MplusVY;MplusVY.setSize(_svSize);
	double* mplusvy=MplusVY.getArray();	

	for(unsigned long spk=spkBottom; spk<spkUp;spk++){

		//compute M+VY for the current speaker
		MplusVY.setAllValues(0.0);

		//Calcul de mplusvy
		for (unsigned long i=0;i<_svSize;i++){
			for (unsigned long j=0;j<_rankT;j++){
				mplusvy[i] += v[j*_svSize + i] * y[spk*_rankT + j ];
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

//-----------------------------------------------------------------------------------------
void TVAcc::substractMplusTWThreaded(unsigned long NUM_THREADS){
	
	if (verboseLevel >= 1) cout <<"(AccumulateTVStat) Substract Mean + VY from Speaker Statistics Threaded" << endl;
	if (NUM_THREADS==0) throw Exception("Num threads can not be 0",__FILE__,__LINE__);

	int rc, status;
	if (NUM_THREADS > _n_speakers) NUM_THREADS=_n_speakers;

	struct S_MplusTWthread_data *thread_data_array = new S_MplusTWthread_data[NUM_THREADS];
	pthread_t *threads = new pthread_t[NUM_THREADS];

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	unsigned long offset=_n_speakers/NUM_THREADS;
	
	double *N=_statN.getArray(); 
	double *F_X=_statF.getArray();
	double *V=_T.getArray();
	double *Y=_W.getArray();
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
		thread_data_array[t].rankEV=_rankT;

		if (verboseLevel>1) cout<<"(AccumulateTVStat) Creating thread n["<< t<< "] for speakers ["<<spkBottom<<"-->"<<spkUp-1<<"]"<<endl;
		rc = pthread_create(&threads[t], &attr, S_MplusTWthread, (void *)&thread_data_array[t]);
		if (rc) throw Exception("ERROR; return code from pthread_create() is ",__FILE__,rc);
		
		spkBottom = spkUp;
	}
	
	pthread_attr_destroy(&attr);
	for(unsigned long t=0; t<NUM_THREADS; t++) {
		rc = pthread_join(threads[t], (void **)&status);
		if (rc)  throw Exception("ERROR; return code from pthread_join() is ",__FILE__,rc);
		if (verboseLevel>1) cout <<"(AccumulateTVStat) Completed join with thread ["<<t<<"] status["<<status<<"]"<<endl;
	}

	free(thread_data_array);
	free(threads);

	if (verboseLevel >= 1) cout << "(AccumulateTVStat) Done " << endl;
}

#endif

//-----------------------------------------------------------------------------------------
void TVAcc::getSpeakerModel(MixtureGD &mixture, String& file){
	
	DoubleVector Sp, MplusVY, UX;
	Sp.setSize(_svSize); Sp.setAllValues(0.0);
	MplusVY.setSize(_svSize);

	this->getMplusTW(MplusVY, file);
	
	for(unsigned long i=0; i<_svSize; i++){
		Sp[i] = MplusVY[i];
	}
	svToModel(Sp,mixture);
}

//-----------------------------------------------------------------------------------------
void TVAcc::orthonormalizeT(){
	// Gram-Schmidt algorithm

	Matrix<double> Q,R;
	Q.setDimensions(_rankT, _svSize);
	Q.setAllValues(0.0);
	R.setDimensions(_rankT, _svSize);
	R.setAllValues(0.0);

	for(unsigned long j=0;j<_rankT;j++){

		DoubleVector v(_svSize,_svSize);
		for(unsigned long k=0; k<_svSize;k++){
			v[k] = _T(j,k);
		}
		for(unsigned long i=0;i<j;i++){
			for(unsigned long k=0; k<_svSize;k++){
				R(i,j) += Q(i,k)*_T(j,k);
			}
			for(unsigned long k=0; k<_svSize;k++){
				v[k] -= R(i,j)* Q(i,k);
			}
		}

		double nV = 0;
		for(unsigned long k=0; k<_svSize;k++){
			nV += v[k]*v[k];
		}
		R(j,j) = sqrt(nV);

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

	//Copy the new matrix in _T
	for(unsigned long i=0;i<_rankT;i++){
		for(unsigned long j=0;j<_svSize;j++){
			_T(i,j) = Q(i,j);
		}
	}
}

//-----------------------------------------------------------------------------------------
// Normalize the Total Variability matrix by using the square root of UBM inverse co-varariance
void TVAcc::normTMatrix(){

	if(verboseLevel > 1)	cout<<"(AccumulateTVStat) Normalize matrix T"<<endl;
	for(unsigned long i=0;i<_svSize;i++){
		double sqrtInvCov = sqrt(_ubm_invvar[i]);
		for(unsigned long j=0;j<_rankT;j++){
			_T(j,i) = _T(j,i) * sqrtInvCov;
		}
	}
}


//-----------------------------------------------------------------------------------------
///Save Accumulators on disk
void TVAcc::saveAccs(Config &config) {

	String fxName, fxhName, nName, nhName;
	fxName = "F_X.mat"; nName = "N.mat";
	if(config.existsParam("nullOrderStatSpeaker"))	nName = config.getParam("matrixFilesPath") + config.getParam("nullOrderStatSpeaker")+config.getParam("saveMatrixFilesExtension");
	if(config.existsParam("firstOrderStatSpeaker")) fxName = config.getParam("matrixFilesPath") + config.getParam("firstOrderStatSpeaker")+config.getParam("saveMatrixFilesExtension");

	_statF.save(fxName,config);
	_statN.save(nName,config);
	if (verboseLevel>=1) cout << "(AccumulateTVStat) Statisitc Accumulators saved" << endl;
}

//-----------------------------------------------------------------------------------------
double TVAcc::getLLK(SegCluster &selectedSegments,MixtureGD &model,FeatureServer&fs,Config & config){

	// TO DO: SPEED UP BY USING APPROXIMATION

	if (verboseLevel >= 1) cout << "(TotalVariability) Compute Likelihood" << endl;
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

//-----------------------------------------------------------------------------------------
void TVAcc::verifyEMLK(Config& config){

	XList ndx(_fileList);
	XLine *pline; String *pFile; ndx.rewind();

	double total=0.0;
	unsigned long maxLLKcomputed=1;
	maxLLKcomputed=config.getParam("computeLLK").toULong();	

	unsigned long cnt=0;
	while((pline=ndx.getLine())!=NULL && cnt < maxLLKcomputed) { 
		while((pFile=pline->getElement())!=NULL && cnt < maxLLKcomputed) {

			/// Compute TV model
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
			if (verboseLevel>1) cout << "(TotalVariability) LLK["<<*pFile<<"]="<<llk<<endl;
			cnt++;
			total+=llk;
		}
	}
	if (verboseLevel >=1) cout << "*** (Verify LLK) Total LLK="<<total<<" ***"<<endl;
}

//-----------------------------------------------------------------------------------------
void TVAcc::estimateAandC(Config& config){
	
	#ifdef THREAD          
	if (config.existsParam("numThread") && config.getParam("numThread").toULong() >0) 	estimateAandCThreaded(config.getParam("numThread").toULong());
	else estimateAandCUnthreaded(config);
	#else
	estimateAandCUnthreaded(config);
	#endif
}

//-----------------------------------------------------------------------------------------
void TVAcc::estimateAandCUnthreaded(Config& config){

	if (verboseLevel >= 1) cout << "(AccumulateTVStat) Compute and Inverse L Matrix for TotalVariability "<<endl;
	
	DoubleSquareMatrix L(_rankT);
	_W.setAllValues(0.0);
	Matrix<double> AUX(1,_rankT);

	_R.setAllValues(0.0);
	_r.setAllValues(0.0);
	_meanW.setAllValues(0.0);
	
	double *y, *v, *f_x, *aux, *invVar, *n, *c, *R, *r, *meanY;

	y=_W.getArray(); v=_T.getArray(); f_x=_statF.getArray(); aux=AUX.getArray(); invVar=_ubm_invvar.getArray(); c=_C.getArray();n=_statN.getArray(); R=_R.getArray(); r=_r.getArray(); meanY=_meanW.getArray();

	_A.setDimensions(_n_distrib,_rankT*_rankT);
	_A.setAllValues(0.0);

	for(unsigned long spk=0; spk<_n_speakers; spk++){	
		L.setAllValues(0.0);
		AUX.setAllValues(0.0);
		for(unsigned long i=0; i<_rankT; i++){	L(i,i)=1.0;}

		double *l;
		l=L.getArray();

		for(unsigned long dis=0; dis<_n_distrib;dis++){
			double *vevt=_TETt[dis].getArray();
			for(unsigned long i=0; i<_rankT; i++){
				for(unsigned long j=0; j<=i; j++){
					l[i*_rankT+j] =+ l[i*_rankT+j] + vevt[i*_rankT+j]*n[spk*_n_distrib+dis];
				}
			}
		}
		//As L is symetric, copy the second half of coefficients
		for(unsigned long i=0;i<_rankT;i++){
			for(unsigned long j=i+1;j<_rankT;j++) {
				l[i*_rankT+j] = l[j*_rankT+i];
			}
		}
		
		//Inverse L
		DoubleSquareMatrix Linv(_rankT);
		L.invert(Linv);
		double *invl = Linv.getArray();

		for(unsigned long i=0;i<_rankT;i++){
			for(unsigned long k=0;k<_svSize;k++) {
				aux[i] += f_x[spk*_svSize+k] * invVar[k] * v[i*_svSize+k];
			}
		}
	
		//multiplication by invL
		for(unsigned long i=0; i<_rankT;i++){
			for(unsigned long k=0; k<_rankT; k++){
				y[spk*_rankT+i] += aux[k] * invl[i*_rankT+k];
			}
		}

		for(unsigned long k=0;k<_rankT;k++){
			meanY[k] += y[spk*_rankT+k];
		}

		for(unsigned long i=0;i<_rankT;i++){
			for(unsigned long j=0;j<_rankT;j++){
					invl[i*_rankT+j] += y[spk*_rankT+i]*y[spk*_rankT+j];
					//Update the Minimum Divergence Accumulator
					R[i*_rankT+j] += invl[i*_rankT+j];
			}
			r[i] += y[spk*_rankT+i];
		}

		for(unsigned long dis = 0; dis<_n_distrib;dis++){
			for(unsigned long i=0;i<_rankT;i++){
				for(unsigned long j = 0;j<_rankT;j++){
					_A(dis, i*(_rankT) +j) += invl[i*_rankT+j] * n[spk*_n_distrib+dis];
				
				}
			}
		}

		for(unsigned long i=0;i<_rankT;i++){
			for(unsigned long j=0;j<_svSize;j++){
					c[i*_svSize+j] += y[spk*_rankT+i] * f_x[spk*_svSize+j];
			}
		}
	}

	// Compute the mean of iVectors
	for(unsigned long k=0;k<_rankT;k++){
		meanY[k] /= _n_speakers;
	}
}

#ifdef THREAD
pthread_mutex_t mutexA=PTHREAD_MUTEX_INITIALIZER;				// Mutex for A
pthread_mutex_t mutexC=PTHREAD_MUTEX_INITIALIZER;				// Mutex for C
	
//-----------------------------------------------------------------------------------------
//				Data strucutre of thread
//-----------------------------------------------------------------------------------------
struct estimateAandCTthread_data{

	double *N;
	unsigned long spkBottom;
	unsigned long spkUp;	
	unsigned long rankEV;
	unsigned long n_distrib;
	unsigned long svSize;
	unsigned long nt;
	RefVector <DoubleSquareMatrix>* vevT;
	RefVector <DoubleSquareMatrix>* tmpR;
	RefVector <DoubleVector>* tmpr;
	RefVector <DoubleVector>* meanY;

	double *Y;
	double *V;
	double *F_X;
	double *ubm_invvar;

	double* C;						// Mutex C
	double* tmpA;
	unsigned long numThread;
};

//-----------------------------------------------------------------------------------------
//				Thread Routine
//-----------------------------------------------------------------------------------------
void *estimateAandCTthread(void *threadarg) {

	pthread_mutex_init(&mutexA, NULL);		// Mutex for A
	pthread_mutex_init(&mutexC, NULL);			// Mutex for C

	struct estimateAandCTthread_data *my_data;
	my_data = (struct estimateAandCTthread_data *) threadarg;

	double *n = my_data->N;
	unsigned long spkBottom = my_data->spkBottom;	
	unsigned long spkUp = my_data->spkUp;
	unsigned long _rankT=my_data->rankEV;
	unsigned long _n_distrib=my_data->n_distrib;
	unsigned long svSize =my_data->svSize;
	unsigned long nt = my_data->nt;

	unsigned long numThread = my_data->numThread;
	double *v = my_data->V;
	double *y = my_data->Y;
	double *f_x = my_data->F_X;
	double *ubm_invvar = my_data->ubm_invvar;
	double* C = my_data->C;
	double* tmpA = my_data->tmpA;

	DoubleSquareMatrix &tmpR=(*(my_data->tmpR))[nt];
	double *R = tmpR.getArray();
	DoubleVector &tmpr=(*(my_data->tmpr))[nt];
	double *r = tmpr.getArray();
	DoubleVector &meanY=(*(my_data->meanY))[nt];
	double *mY = meanY.getArray();

	Matrix<double> AUX(1,_rankT);
	double* aux = AUX.getArray();

	for(unsigned long spk=spkBottom; spk<spkUp; spk++){

		DoubleSquareMatrix L(_rankT);
		L.setAllValues(0.0);	
		double *l=L.getArray();
		for(unsigned long i=0; i<_rankT; i++){	l[i*_rankT+i]=1.0;}

		for(unsigned long d=0; d<_n_distrib;d++){
			DoubleSquareMatrix &vEvT=(*(my_data->vevT))[d];
			double *vevt = vEvT.getArray();
			for(unsigned long i=0; i<_rankT; i++){
				for(unsigned long j=0; j<=i; j++){
					l[i*_rankT+j] =+ l[i*_rankT+j] + vevt[i*_rankT+j]*n[spk*_n_distrib+d];
				}
			}
		}
		//As L is symetric, copy the second half of coefficients
		for(unsigned long i=0;i<_rankT;i++){
			for(unsigned long j=i+1;j<_rankT;j++) {
				l[i*_rankT+j] = l[j*_rankT+i];
			}
		}

		//Inverse L and stock it in _l_spk_inv[spk]
		DoubleSquareMatrix Linv(_rankT);
		L.invert(Linv);
		double *invl = Linv.getArray();

		AUX.setAllValues(0.0);
		for(unsigned long i=0;i<_rankT;i++){
			for(unsigned long k=0;k<svSize;k++) {
				aux[i] += f_x[spk*svSize+k] * ubm_invvar[k] * v[i*svSize+k];
			}
		}
	
		//multiplication by invL
		for(unsigned long i=0; i<_rankT;i++){
			for(unsigned long k=0; k<_rankT; k++){
				y[spk*_rankT+i] += aux[k] * invl[i*_rankT+k];
			}
		}

		for(unsigned long k=0;k<_rankT;k++){
			mY[k] += y[spk*_rankT+k];
		}

		for(unsigned long i=0;i<_rankT;i++){
			for(unsigned long j=0;j<_rankT;j++){
				invl[i*_rankT+j] += y[spk*_rankT+i]*y[spk*_rankT+j];
				//Update the Minimum Divergence Accumulator
				R[i*_rankT+j] += invl[i*_rankT+j];
			}
			r[i] += y[spk*_rankT+i];
		}

	pthread_mutex_lock(&mutexA);		// Lock Mutex
	for(unsigned long dis = 0; dis<_n_distrib;dis++){
		for(unsigned long i=0;i<_rankT;i++){
			for(unsigned long j = 0;j<_rankT;j++){
				tmpA[dis*(_rankT*_rankT)+ i*(_rankT) +j] += invl[i*_rankT+j] * n[spk*_n_distrib+dis];
			}
		}
	}
	pthread_mutex_unlock(&mutexA);	// Unlock Mutex


	pthread_mutex_lock(&mutexC);	// Lock Mutex
	for(unsigned long i=0;i<_rankT;i++){
		for(unsigned long j=0;j<svSize;j++){
			C[i*svSize+j] += y[spk*_rankT+i] * f_x[spk*svSize+j];
		}
	}
	pthread_mutex_unlock(&mutexC);	// Unlock Mutex
	}
	pthread_exit((void*) 0);
	return (void*)0 ;
}

//-----------------------------------------------------------------------------------------
void TVAcc::estimateAandCThreaded(unsigned long NUM_THREADS){
	
	if (verboseLevel >= 1) cout << "(AccumulateTVStat) Estimate A and C Matrices for TotalVariability Threaded"<<endl;
	if (NUM_THREADS==0) throw Exception("Num threads can not be 0",__FILE__,__LINE__);

	int rc, status;
	if (NUM_THREADS > _n_speakers) NUM_THREADS=_n_speakers;

	struct estimateAandCTthread_data *thread_data_array = new estimateAandCTthread_data[NUM_THREADS];
	pthread_t *threads = new pthread_t[NUM_THREADS];

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	unsigned long offset=_n_speakers/NUM_THREADS;

	_W.setAllValues(0.0);
	_R.setAllValues(0.0);
	_r.setAllValues(0.0);
	_meanW.setAllValues(0.0);

	double *N =_statN.getArray(); 
	double *V = _T.getArray();
	double *Y = _W.getArray();
	double *F_X = _statF.getArray();
	double *ubm_invvar = _ubm_invvar.getArray();

	double *c = _C.getArray();
	_A.setDimensions(_n_distrib,_rankT*_rankT);
	_A.setAllValues(0.0);
	double *tmpA = _A.getArray();

	//Temporary accumulatores for Minimum Divergence
	RefVector<DoubleSquareMatrix> tmpR;
	RefVector<DoubleVector> tmpr;
	RefVector<DoubleVector> meanY;
	for(unsigned long nt=0;nt<NUM_THREADS;nt++){
		tmpR.addObject(*new DoubleSquareMatrix(_rankT));
		tmpR[nt].setAllValues(0.0);
		tmpr.addObject(*new DoubleVector(_rankT,_rankT));
		tmpr[nt].setAllValues(0.0);
		meanY.addObject(*new DoubleVector(_rankT,_rankT));
		meanY[nt].setAllValues(0.0);
	}

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
		thread_data_array[t].rankEV=_rankT;
		thread_data_array[t].n_distrib=_n_distrib;
		thread_data_array[t].svSize=_svSize;
		thread_data_array[t].vevT=&(_TETt);
		thread_data_array[t].nt=t;
		thread_data_array[t].numThread=NUM_THREADS;
		thread_data_array[t].V=V;
		thread_data_array[t].Y=Y;
		thread_data_array[t].F_X=F_X;
		thread_data_array[t].ubm_invvar=ubm_invvar;
		thread_data_array[t].C = c;
		thread_data_array[t].tmpA=tmpA;								
		thread_data_array[t].tmpR=&(tmpR);
		thread_data_array[t].tmpr=&(tmpr);
		thread_data_array[t].meanY=&(meanY);

		if (verboseLevel>1) cout<<"(AccumulateTVStat) Creating thread n["<< t<< "] for speakers["<<spkBottom<<"-->"<<spkUp-1<<"]"<<endl;
		rc = pthread_create(&threads[t], &attr, estimateAandCTthread, (void *)&thread_data_array[t]);
		if (rc) throw Exception("ERROR; return code from pthread_create() is ",__FILE__,rc);
		
		spkBottom = spkUp;
	}
	
	pthread_attr_destroy(&attr);
	for(unsigned long t=0; t<NUM_THREADS; t++) {
		rc = pthread_join(threads[t], (void **)&status);
		if (rc)  throw Exception("ERROR; return code from pthread_join() is ",__FILE__,rc);
		if (verboseLevel >1) cout <<"(AccumulateTVStat) Completed join with thread ["<<t<<"] status["<<status<<"]"<<endl;
	}
	
	free(thread_data_array);
	free(threads);

	//Sum Minimum Divergence Accumulators after multithreading
	for(unsigned long mt=0; mt<NUM_THREADS;mt++){
		for(unsigned long i=0;i<_rankT;i++){
			for(unsigned long j = 0;j<_rankT;j++){
				_R(i,j) += tmpR[mt](i,j);
			}
			_r[i] +=tmpr[mt][i];
			_meanW[i] += meanY[mt][i];
		}
	}

	// Compute the mean of iVectors
	for(unsigned long k=0;k<_rankT;k++){
		_meanW[k] /= _n_speakers;
	}

	if (verboseLevel >= 1) cout << "(AccumulateTVStat) Done " << endl;
}
#endif

//-----------------------------------------------------------------------------------------
void TVAcc::minDivergence(){
	
	if (verboseLevel>0) cout << "(AccumulateTVStat) Minimum Divergence step " << endl;

	for(unsigned long i=0;i<_rankT;i++){
		_r[i] /= _n_sessions;
	}

	for(unsigned long i=0;i<_rankT;i++){
		for(unsigned long j=0;j<_rankT;j++){
			_R(i,j) = _R(i,j)/_n_sessions - (_r[i]*_r[j]);
		}
	}

	DoubleSquareMatrix Ch;
	Ch.setSize(_rankT);
	_R.upperCholesky(Ch);
	
	DoubleVector newMean(_svSize,_svSize);
	newMean = _ubm_means;
	for(unsigned long j=0;j<_svSize;j++){
		for(unsigned long k=0;k<_rankT;k++){
			newMean[j] += _meanW[k]*_T(k,j);
		}
	}
	_ubm_means = newMean;

	Matrix<double> tmpV;
	tmpV.setDimensions(_T.rows(),_T.cols());
	tmpV.setAllValues(0.0);
	for(unsigned long i=0;i<_rankT;i++){
		for(unsigned long j=0;j<_T.cols();j++){
			for(unsigned long k=0;k<_rankT;k++){
				tmpV(i,j) += Ch(i,k)*_T(k,j);
			}
		}
	}
	for(unsigned long i=0;i<_rankT;i++){
		for(unsigned long j=0;j<_T.cols();j++){
			_T(i,j) = tmpV(i,j);
		}
	}
	if (verboseLevel>0) cout << "(AccumulateTVStat) Done " << endl;
}
#endif

//-----------------------------------------------------------------------------------------
void TVAcc::estimateW(Config& config){
	
	#ifdef THREAD          
	if (config.existsParam("numThread") && config.getParam("numThread").toULong() >0)	estimateWThreaded(config.getParam("numThread").toULong());
	else estimateWUnThreaded(config);
	#else
	estimateWUnThreaded(config);
	#endif
}

//-----------------------------------------------------------------------------------------
void TVAcc::estimateWUnThreaded(Config& config){

	if (verboseLevel >= 1) cout << "(AccumulateTVStat) Estimate Y from Statistics "<<endl;
	
	DoubleSquareMatrix L(_rankT);
	_W.setAllValues(0.0);
	Matrix<double> AUX(1,_rankT);

	double *y, *t, *f_x, *aux, *invVar, *n;

	y=_W.getArray(); t=_T.getArray(); f_x=_statF.getArray(); aux=AUX.getArray(); invVar=_ubm_invvar.getArray(); n=_statN.getArray();

	for(unsigned long spk=0; spk<_n_speakers; spk++){
		L.setAllValues(0.0);
		for(unsigned long i=0; i<_rankT; i++){	L(i,i)=1.0;}

		double *l;
		l=L.getArray();

		for(unsigned long dis=0; dis<_n_distrib;dis++){
			double *tett=_TETt[dis].getArray();
			for(unsigned long i=0; i<_rankT; i++){
				for(unsigned long j=0; j<=i; j++){
					l[i*_rankT+j] =+ l[i*_rankT+j] + tett[i*_rankT+j]*n[spk*_n_distrib+dis];
				}
			}
		}

		//As L is symmetric, copy the second half of coefficients
		for(unsigned long i=0;i<_rankT;i++){
			for(unsigned long j=i+1;j<_rankT;j++) {
				l[i*_rankT+j] = l[j*_rankT+i];
			}
		}
	
		//Inverse L
		DoubleSquareMatrix Linv(_rankT);
		L.invert(Linv);
		double *invl = Linv.getArray();

		AUX.setAllValues(0.0);
		for(unsigned long i=0;i<_rankT;i++){
			for(unsigned long k=0;k<_svSize;k++) {
				aux[i] += f_x[spk*_svSize+k] * invVar[k] * t[i*_svSize+k];
			}
		}

		//multiplication by invL
		for(unsigned long i=0; i<_rankT;i++){
			for(unsigned long k=0; k<_rankT; k++){
				y[spk*_rankT+i] += aux[k] * invl[i*_rankT+k];
			}
		}
		
	}
}

#ifdef THREAD
//-----------------------------------------------------------------------------------------
//				Data structure of thread
//-----------------------------------------------------------------------------------------
struct estimateWTthread_data{

	double *N;
	unsigned long spkBottom;
	unsigned long spkUp;	
	unsigned long rankEV;
	unsigned long n_distrib;
	unsigned long svSize;
	unsigned long nt;
	RefVector <DoubleSquareMatrix>* vevT;

	double *Y;
	double *V;
	double *F_X;
	double *ubm_invvar;

	unsigned long numThread;

};

//-----------------------------------------------------------------------------------------
//				Thread Routine
//-----------------------------------------------------------------------------------------
void *estimateWTthread(void *threadarg) {
	struct estimateWTthread_data *my_data;
	my_data = (struct estimateWTthread_data *) threadarg;

	double *n = my_data->N;
	unsigned long spkBottom = my_data->spkBottom;	
	unsigned long spkUp = my_data->spkUp;
	unsigned long _rankT=my_data->rankEV;
	unsigned long _n_distrib=my_data->n_distrib;
	unsigned long svSize =my_data->svSize;
	unsigned long nt = my_data->nt;
	unsigned long numThread = my_data->numThread;

	double *v = my_data->V;
	double *y = my_data->Y;
	double *f_x = my_data->F_X;
	double *ubm_invvar = my_data->ubm_invvar;

	DoubleSquareMatrix Linv(_rankT);
	DoubleVector tmpAux(svSize,svSize);

	for(unsigned long spk=spkBottom; spk<spkUp; spk++){

		tmpAux.setAllValues(0.0);

		DoubleSquareMatrix L(_rankT);	
		L.setAllValues(0.0);	
		double *l=L.getArray();
		for(unsigned long i=0; i<_rankT; i++){	l[i*_rankT+i]=1.0;}

		for(unsigned long d=0; d<_n_distrib;d++){
			DoubleSquareMatrix &vEvT=(*(my_data->vevT))[d];
			double *vevt = vEvT.getArray();
			for(unsigned long i=0; i<_rankT; i++){
				for(unsigned long j=0; j<=i; j++){
					l[i*_rankT+j] =+ l[i*_rankT+j] + vevt[i*_rankT+j]*n[spk*_n_distrib+d];
				}
			}
		}

		//As L is symetric, copy the second half of coefficients
		for(unsigned long i=0;i<_rankT;i++){
			for(unsigned long j=i+1;j<_rankT;j++) {
				l[i*_rankT+j] = l[j*_rankT+i];
			}
		}
	
		//Inverse L and stock it
		L.invert(Linv);
		double *invl = Linv.getArray();

		tmpAux.setAllValues(0.0);
		for(unsigned long i=0;i<_rankT;i++){
			for(unsigned long k=0;k<svSize;k++) {
				tmpAux[i] += f_x[spk*svSize+k] * ubm_invvar[k] * v[i*svSize+k];
			}
		}

		//multiplication by invL
		for(unsigned long i=0; i<_rankT;i++){
			for(unsigned long k=0; k<_rankT; k++){
				y[spk*_rankT+i] += tmpAux[k] * invl[i*_rankT+k];
			}
		}
	}
	pthread_exit((void*) 0);
	return (void*)0 ;
}

//-----------------------------------------------------------------------------------------
void TVAcc::estimateWThreaded(unsigned long NUM_THREADS){
	
	if (verboseLevel >= 1) cout << "(AccumulateTVStat) Estimate Y from Statistics Threaded"<<endl;
	if (NUM_THREADS==0) throw Exception("Num threads can not be 0",__FILE__,__LINE__);

	int rc, status;
	if (NUM_THREADS > _n_speakers) NUM_THREADS=_n_speakers;

	struct estimateWTthread_data *thread_data_array = new estimateWTthread_data[NUM_THREADS];
	pthread_t *threads = new pthread_t[NUM_THREADS];

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	unsigned long offset=_n_speakers/NUM_THREADS;

	_W.setAllValues(0.0);
	
	double *N =_statN.getArray(); 
	double *V = _T.getArray();
	double *Y = _W.getArray();
	double *F_X = _statF.getArray();
	double *ubm_invvar = _ubm_invvar.getArray();

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
		thread_data_array[t].rankEV=_rankT;
		thread_data_array[t].n_distrib=_n_distrib;
		thread_data_array[t].svSize=_svSize;
		thread_data_array[t].vevT=&(_TETt);
		thread_data_array[t].nt=t;
		thread_data_array[t].numThread=NUM_THREADS;
		thread_data_array[t].V=V;
		thread_data_array[t].Y=Y;
		thread_data_array[t].F_X=F_X;
		thread_data_array[t].ubm_invvar=ubm_invvar;

		if (verboseLevel>1) cout<<"(AccumulateTVStat) Creating thread n["<< t<< "] for speakers["<<spkBottom<<"-->"<<spkUp-1<<"]"<<endl;
		rc = pthread_create(&threads[t], &attr, estimateWTthread, (void *)&thread_data_array[t]);
		if (rc) throw Exception("ERROR; return code from pthread_create() is ",__FILE__,rc);
		
		spkBottom = spkUp;
	}
	
	pthread_attr_destroy(&attr);
	for(unsigned long t=0; t<NUM_THREADS; t++) {
		rc = pthread_join(threads[t], (void **)&status);
		if (rc)  throw Exception("ERROR; return code from pthread_join() is ",__FILE__,rc);
		if (verboseLevel >1) cout <<"(AccumulateTVStat) Completed join with thread ["<<t<<"] status["<<status<<"]"<<endl;
	}
	free(thread_data_array);
	free(threads);

	if (verboseLevel >= 1) cout << "(AccumulateTVStat) Done " << endl;
}
#endif

//-----------------------------------------------------------------------------------------
void TVAcc::estimateWUbmWeight(DoubleSquareMatrix &W, Config& config){

	#ifdef THREAD          
	if (config.existsParam("numThread") && config.getParam("numThread").toULong() >0)	estimateWUbmWeightThreaded(W,config.getParam("numThread").toULong());
	else estimateWUbmWeightUnThreaded(W,config);
	#else
	estimateWUbmWeightUnThreaded(W,config);
	#endif
}

//-----------------------------------------------------------------------------------------
void TVAcc::estimateWUbmWeightUnThreaded(DoubleSquareMatrix &W, Config &config){	// Threaded version not implemented yet

	if (verboseLevel >= 1) cout << "(AccumulateTVStat) Estimate Y from Statistics "<<endl;
	
	DoubleSquareMatrix L(_rankT);
	_W.setAllValues(0.0);
	Matrix<double> AUX(1,_rankT);

	double *y, *t, *f_x, *aux, *n;

	y=_W.getArray(); t=_T.getArray(); f_x=_statF.getArray(); aux=AUX.getArray(); n=_statN.getArray();

	for(unsigned long spk=0; spk<_n_speakers; spk++){

		double *l, *w;
		l=L.getArray();
		w= W.getArray();

		double n_sum = 0.0;
		for(unsigned long c=0; c<_n_distrib; c++)
			n_sum += n[spk*_n_distrib+c];

		for(unsigned long i=0;i<_rankT;i++){
			for(unsigned long j=0;j<_rankT;j++){
				l[i*_rankT+j] = n_sum*w[i*_rankT+j];
			}
			l[i*_rankT+i] += 1.0;
		}

		//Inverse L
		DoubleSquareMatrix Linv(_rankT);
		L.invert(Linv);
		double *invl = Linv.getArray();

		AUX.setAllValues(0.0);
		for(unsigned long i=0;i<_rankT;i++){
			for(unsigned long k=0;k<_svSize;k++) {
				aux[i] += f_x[spk*_svSize+k] * t[i*_svSize+k];
			}
		}

		//multiplication by invL
		for(unsigned long i=0; i<_rankT;i++){
			for(unsigned long k=0; k<_rankT; k++){
				y[spk*_rankT+i] += aux[k] * invl[i*_rankT+k];
			}
		}
	}
}

#ifdef THREAD
//-----------------------------------------------------------------------------------------
//				Data structure of thread
//-----------------------------------------------------------------------------------------
struct estimateWUbmWeightTthread_data{

	double *N;
	unsigned long spkBottom;
	unsigned long spkUp;	
	unsigned long rankTV;
	unsigned long n_distrib;
	unsigned long svSize;
	unsigned long nt;

	double *Y;
	double *W;
	double *T;
	double *F_X;

	unsigned long numThread;

};

//-----------------------------------------------------------------------------------------
//				Thread Routine
//-----------------------------------------------------------------------------------------
void *estimateWUbmWeightTthread(void *threadarg) {
	struct estimateWUbmWeightTthread_data *my_data;
	my_data = (struct estimateWUbmWeightTthread_data *) threadarg;

	double *n = my_data->N;
	unsigned long spkBottom = my_data->spkBottom;	
	unsigned long spkUp = my_data->spkUp;
	unsigned long _rankT=my_data->rankTV;
	unsigned long _n_distrib=my_data->n_distrib;
	unsigned long svSize =my_data->svSize;
	unsigned long nt = my_data->nt;
	unsigned long numThread = my_data->numThread;

	double *t = my_data->T;
	double *y = my_data->Y;
	double *w = my_data->W;
	double *f_x = my_data->F_X;

	DoubleSquareMatrix Linv(_rankT);
	DoubleVector tmpAux(svSize,svSize);

	for(unsigned long spk=spkBottom; spk<spkUp; spk++){

		tmpAux.setAllValues(0.0);

		DoubleSquareMatrix L(_rankT);	
		L.setAllValues(0.0);	
		double *l=L.getArray();

		double n_sum = 0.0;
		for(unsigned long c=0; c<_n_distrib; c++)
			n_sum += n[spk*_n_distrib+c];

		for(unsigned long i=0;i<_rankT;i++){
			for(unsigned long j=0;j<_rankT;j++){
				l[i*_rankT+j] = n_sum*w[i*_rankT+j];
			}
			l[i*_rankT+i] += 1.0;
		}

		//Inverse L
		DoubleSquareMatrix Linv(_rankT);
		L.invert(Linv);
		double *invl = Linv.getArray();

		tmpAux.setAllValues(0.0);
		for(unsigned long i=0;i<_rankT;i++){
			for(unsigned long k=0;k<svSize;k++) {
				tmpAux[i] += f_x[spk*svSize+k] * t[i*svSize+k];
			}
		}

		//multiplication by invL
		for(unsigned long i=0; i<_rankT;i++){
			for(unsigned long k=0; k<_rankT; k++){
				y[spk*_rankT+i] += tmpAux[k] * invl[i*_rankT+k];
			}
		}
	}
	pthread_exit((void*) 0);
	return (void*)0 ;
}

//-----------------------------------------------------------------------------------------
void TVAcc::estimateWUbmWeightThreaded(DoubleSquareMatrix &W, unsigned long NUM_THREADS){
	
	if (verboseLevel >= 1) cout << "(AccumulateTVStat) Approximate W  by using UBM weights Threaded"<<endl;
	if (NUM_THREADS==0) throw Exception("Num threads can not be 0",__FILE__,__LINE__);

	int rc, status;
	if (NUM_THREADS > _n_speakers) NUM_THREADS=_n_speakers;

	struct estimateWUbmWeightTthread_data *thread_data_array = new estimateWUbmWeightTthread_data[NUM_THREADS];
	pthread_t *threads = new pthread_t[NUM_THREADS];

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	unsigned long offset=_n_speakers/NUM_THREADS;

	_W.setAllValues(0.0);
	
	double *N =_statN.getArray(); 
	double *T = _T.getArray();
	double *Y = _W.getArray();
	double *w = W.getArray();
	double *F_X = _statF.getArray();
	
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
		thread_data_array[t].rankTV=_rankT;
		thread_data_array[t].n_distrib=_n_distrib;
		thread_data_array[t].svSize=_svSize;
		thread_data_array[t].nt=t;
		thread_data_array[t].numThread=NUM_THREADS;
		thread_data_array[t].T=T;
		thread_data_array[t].W=w;
		thread_data_array[t].Y=Y;
		thread_data_array[t].F_X=F_X;

		if (verboseLevel>1) cout<<"(AccumulateTVStat) Creating thread n["<< t<< "] for speakers["<<spkBottom<<"-->"<<spkUp-1<<"]"<<endl;
		rc = pthread_create(&threads[t], &attr, estimateWUbmWeightTthread, (void *)&thread_data_array[t]);
		if (rc) throw Exception("ERROR; return code from pthread_create() is ",__FILE__,rc);
		
		spkBottom = spkUp;
	}
	
	pthread_attr_destroy(&attr);
	for(unsigned long t=0; t<NUM_THREADS; t++) {
		rc = pthread_join(threads[t], (void **)&status);
		if (rc)  throw Exception("ERROR; return code from pthread_join() is ",__FILE__,rc);
		if (verboseLevel >1) cout <<"(AccumulateTVStat) Completed join with thread ["<<t<<"] status["<<status<<"]"<<endl;
	}
	free(thread_data_array);
	free(threads);

	if (verboseLevel >= 1) cout << "(AccumulateTVStat) Done " << endl;
}
#endif

//-----------------------------------------------------------------------------------------
void TVAcc::estimateWEigenDecomposition(Matrix<double> D, Matrix<double> Q, Config& config){
	#ifdef THREAD          
	if (config.existsParam("numThread") && config.getParam("numThread").toULong() >0)	estimateWEigenDecompositionThreaded(D,Q,config.getParam("numThread").toULong());
	else estimateWEigenDecompositionUnThreaded(D,Q,config);
	#else
	estimateWEigenDecompositionUnThreaded(D,Q,config);
	#endif
}

//-----------------------------------------------------------------------------------------
void TVAcc::estimateWEigenDecompositionUnThreaded(Matrix<double> D, Matrix<double> Q, Config& config){

	if (verboseLevel >= 1) cout << "(AccumulateTVStat) Approximate I-Vector by using Eigen Decomposition "<<endl;

	for(unsigned long spk=0; spk<_n_speakers; spk++){

		// Approximate matrix L's inverse
		DoubleVector invL(_rankT,_rankT); invL.setAllValues(0.0);
		for(unsigned long i=0 ; i<_rankT ; i++){
			double tmp = 1.0;
			for(unsigned long cc=0 ; cc<_n_distrib ; cc++){
				tmp += _statN(spk,cc)*D(cc,i);
			}
			invL[i] = 1/tmp;
		}

		// Compute i-vector
		DoubleVector aux(_rankT,_rankT); aux.setAllValues(0.0);
		double *w, *t, *f_x, *invVar, *invl;
		w=_W.getArray(); t=_T.getArray(); f_x=_statF.getArray(); invVar=_ubm_invvar.getArray(); invl=invL.getArray();

		for(unsigned long i=0;i<_rankT;i++){
			for(unsigned long k=0;k<_svSize;k++) {
				aux[i] += f_x[spk*_svSize+k] * t[i*_svSize+k];
			}
		}

		// Compute Q*invL*Q'
		Matrix<double> appL(_rankT,_rankT); appL.setAllValues(0.0);
		double *q, *appl;
		q = Q.getArray(); appl = appL.getArray();
		for(unsigned long i=0; i<_rankT ; i++)
			for(unsigned long j=0; j<_rankT ; j++)
				for(unsigned long k=0; k<_rankT ; k++)
					appL(i,j) += Q(i,k)*invL[k]*Q(j,k);

		//multiplication by invL
		for(unsigned long i=0; i<_rankT;i++){
			for(unsigned long k=0; k<_rankT; k++){
				w[spk*_rankT+i] += aux[k] * appl[i*_rankT+k];
			}
		}
	}
}

#ifdef THREAD
//-----------------------------------------------------------------------------------------
//				Data structure of thread
//-----------------------------------------------------------------------------------------
struct estimateWEigenDecompositionTthread_data{

	double *N;
	double *D;
	double *Q;
	unsigned long spkBottom;
	unsigned long spkUp;	
	unsigned long rankTV;
	unsigned long n_distrib;
	unsigned long svSize;
	unsigned long nt;

	double *W;
	double *T;
	double *F_X;
	double *ubm_invvar;

	unsigned long numThread;
};

//-----------------------------------------------------------------------------------------
//				Thread Routine
//-----------------------------------------------------------------------------------------
void *estimateWEigenDecompositionTthread(void *threadarg) {
	struct estimateWEigenDecompositionTthread_data *my_data;
	my_data = (struct estimateWEigenDecompositionTthread_data *) threadarg;

	double *n = my_data->N;
	unsigned long spkBottom = my_data->spkBottom;	
	unsigned long spkUp = my_data->spkUp;
	unsigned long _rankT=my_data->rankTV;
	unsigned long _n_distrib=my_data->n_distrib;
	unsigned long svSize =my_data->svSize;
	unsigned long nt = my_data->nt;
	unsigned long numThread = my_data->numThread;

	double *t = my_data->T;
	double *d = my_data->D;
	double *q = my_data->Q;
	double *w = my_data->W;
	double *f_x = my_data->F_X;
	double *ubm_invvar = my_data->ubm_invvar;

	//DoubleSquareMatrix Linv(_rankT);
	//DoubleVector tmpAux(svSize,svSize);

	for(unsigned long spk=spkBottom; spk<spkUp; spk++){

		// Approximate matrix L's inverse
		DoubleVector invL(_rankT,_rankT); invL.setAllValues(0.0);
		for(unsigned long i=0 ; i<_rankT ; i++){
			double tmp = 1.0;
			for(unsigned long cc=0 ; cc<_n_distrib ; cc++){
				tmp += n[spk*_n_distrib+cc]*d[cc*_rankT+i];
				//tmp += _statN(spk,cc)*D(cc,i);
			}
			invL[i] = 1/tmp;
		}

		// Compute i-vector
		DoubleVector aux(_rankT,_rankT); aux.setAllValues(0.0);
//		double *w, *t, *f_x, *invVar, *invl;
		double *invl = invL.getArray();
//		w=_W.getArray(); t=_T.getArray(); f_x=_statF.getArray(); invVar=_ubm_invvar.getArray(); invl=invL.getArray();

		for(unsigned long i=0;i<_rankT;i++){
			for(unsigned long k=0;k<svSize;k++) {
				aux[i] += f_x[spk*svSize+k] * t[i*svSize+k];
			}
		}

		// Compute Q*invL*Q'
		Matrix<double> appL(_rankT,_rankT); appL.setAllValues(0.0);
//		double *q, *appl;
//		q = Q.getArray(); 
		double *appl = appL.getArray();
		for(unsigned long i=0; i<_rankT ; i++)
			for(unsigned long j=0; j<_rankT ; j++)
				for(unsigned long k=0; k<_rankT ; k++)
					appl[i*_rankT+j] += q[i*_rankT+k]*invL[k]*q[j*_rankT+k];
					//appL(i,j) += Q(i,k)*invL[k]*Q(j,k);

		//multiplication by invL
		for(unsigned long i=0; i<_rankT;i++){
			for(unsigned long k=0; k<_rankT; k++){
				w[spk*_rankT+i] += aux[k] * appl[i*_rankT+k];
			}
		}
	}
	pthread_exit((void*) 0);
	return (void*)0 ;
}

//-----------------------------------------------------------------------------------------
void TVAcc::estimateWEigenDecompositionThreaded(Matrix<double> &D, Matrix<double> &Q, unsigned long NUM_THREADS){
	
	if (verboseLevel >= 1) cout << "(AccumulateTVStat) Approximate W  by using Eigen Decomposition Threaded"<<endl;
	if (NUM_THREADS==0) throw Exception("Num threads can not be 0",__FILE__,__LINE__);

	int rc, status;
	if (NUM_THREADS > _n_speakers) NUM_THREADS=_n_speakers;

	struct estimateWEigenDecompositionTthread_data *thread_data_array = new estimateWEigenDecompositionTthread_data[NUM_THREADS];
	pthread_t *threads = new pthread_t[NUM_THREADS];

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	unsigned long offset=_n_speakers/NUM_THREADS;

	_W.setAllValues(0.0);
	
	double *N =_statN.getArray(); 
	double *T = _T.getArray();
	double *d = D.getArray();
	double *q = Q.getArray();
	double *w = _W.getArray();
	double *F_X = _statF.getArray();
	double *ubm_invvar = _ubm_invvar.getArray();

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
		thread_data_array[t].rankTV=_rankT;
		thread_data_array[t].n_distrib=_n_distrib;
		thread_data_array[t].svSize=_svSize;
		thread_data_array[t].nt=t;
		thread_data_array[t].numThread=NUM_THREADS;
		thread_data_array[t].T=T;
		thread_data_array[t].W=w;
		thread_data_array[t].D=d;
		thread_data_array[t].Q=q;
		thread_data_array[t].F_X=F_X;
		thread_data_array[t].ubm_invvar=ubm_invvar;

		if (verboseLevel>1) cout<<"(AccumulateTVStat) Creating thread n["<< t<< "] for speakers["<<spkBottom<<"-->"<<spkUp-1<<"]"<<endl;
		rc = pthread_create(&threads[t], &attr, estimateWEigenDecompositionTthread, (void *)&thread_data_array[t]);
		if (rc) throw Exception("ERROR; return code from pthread_create() is ",__FILE__,rc);
		
		spkBottom = spkUp;
	}
	
	pthread_attr_destroy(&attr);
	for(unsigned long t=0; t<NUM_THREADS; t++) {
		rc = pthread_join(threads[t], (void **)&status);
		if (rc)  throw Exception("ERROR; return code from pthread_join() is ",__FILE__,rc);
		if (verboseLevel >1) cout <<"(AccumulateTVStat) Completed join with thread ["<<t<<"] status["<<status<<"]"<<endl;
	}
	free(thread_data_array);
	free(threads);

	if (verboseLevel >= 1) cout << "(AccumulateTVStat) Done " << endl;
}
#endif



















//-----------------------------------------------------------------------------------------
void TVAcc::saveWbyFile(Config &config){
	
	String svPath=config.getParam("saveVectorFilesPath");
	String yExtension = ".y";
	if(config.existsParam("vectorFilesExtension"))	yExtension = config.getParam("vectorFilesExtension");

	String inputClientListFileName = config.getParam("targetIdList");
	XList inputClientList(inputClientListFileName,config);

	XLine * linep;
	unsigned long session = 0;
	while ((linep=inputClientList.getLine()) != NULL){             	// linep gives the XLine with the Id of a given client and the list of files
		String id=linep->getElement(0); 
		String yFile=svPath+id+yExtension;
		
		Matrix<double> sessionY;
		sessionY.setDimensions(1,_rankT);
		for(unsigned long i=0;i<_rankT;i++){
			sessionY(0,i) = _W(session,i);
		}
		sessionY.save(yFile,config);
		session++;
	}
}


//-----------------------------------------------------------------------------------------
void TVAcc::getWeightedCov(DoubleSquareMatrix &W, DoubleVector& weight, Config& config){

	#ifdef THREAD          
	if (config.existsParam("numThread") && config.getParam("numThread").toULong() >0)	getWeightedCovThreaded(W,weight,config.getParam("numThread").toULong());
	else getWeightedCovUnThreaded(W,weight);
	#else
	getWeightedCovUnThreaded(W,weight);
	#endif
}

//-----------------------------------------------------------------------------------------
void TVAcc::getWeightedCovUnThreaded(DoubleSquareMatrix &W, DoubleVector& weight){

	if(verboseLevel > 1)	cout<<"(AccumulateTVStat) Compute weighted covariance matrix"<<endl;
	W.setAllValues(0.0);

	for (unsigned long cc=0 ; cc<_n_distrib ; cc++){
		Matrix<double> Tc = _T.crop(0,cc*_vectSize,_rankT,_vectSize);
		// Compute half of the coefficient because the matrix is symetric
		for(unsigned long i=0 ; i<_rankT ; i++)
			for(unsigned long j=0 ; j<i+1 ; j++)
				for(unsigned long k=0 ; k<_vectSize ; k++)
					W(i,j) += weight[cc]*Tc(i,k)*Tc(j,k);
	}

	// Copy the other half of the coefficients
	for(unsigned long i=0 ; i<_rankT ; i++)
			for(unsigned long j=i ; j<_rankT ; j++)
				W(i,j) = W(j,i);
}

#ifdef THREAD
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Data structure of thread
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
struct WeightedCovthread_data{

	double *T;
	double *weight;
	unsigned long disBottom;
	unsigned long disUp;	
	unsigned long rankTV;
	unsigned long svSize;
	unsigned long vectSize;
	RefVector <DoubleSquareMatrix>* W;
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Thread Routine
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void *WeightedCovthread(void *threadarg){

	struct WeightedCovthread_data *my_data;
	my_data = (struct WeightedCovthread_data *) threadarg;
	
	double *t = my_data->T;
	unsigned long disBottom = my_data->disBottom;	
	unsigned long disUp = my_data->disUp;
	double *weight=my_data->weight;
	unsigned long _rankT=my_data->rankTV;
	unsigned long _svSize=my_data->svSize;
	unsigned long _vectSize=my_data->vectSize;

	for (unsigned long cc=disBottom; cc<disUp; cc++){
	
		DoubleSquareMatrix &W=(*(my_data->W))[cc];
		W.setAllValues(0.0);
		double *w = W.getArray();

		// Compute half of the coefficient because the matrix is symetric
		for(unsigned long i=0 ; i<_rankT ; i++)
			for(unsigned long j=0 ; j<i+1 ; j++)
				for(unsigned long k=0 ; k<_vectSize ; k++)
					w[i*_rankT+j] += weight[cc] * t[(i*_svSize)+(cc*_vectSize)+k] * t[(j*_svSize)+(cc*_vectSize)+k];
	}

	pthread_exit((void*) 0);
	return (void*)0 ;
}

//-----------------------------------------------------------------------------------------
void TVAcc::getWeightedCovThreaded(DoubleSquareMatrix &W, DoubleVector& weight, unsigned long NUM_THREADS){
	
	if (verboseLevel >=1) cout << "(AccumulateTVStat) Compute weighted covariance matrix Threaded"<<endl;
	if (NUM_THREADS==0) throw Exception("Num threads can not be 0",__FILE__,__LINE__);
	
	double *T=_T.getArray();
	double *ubm_weight=weight.getArray();
	
	int rc, status;
	if (NUM_THREADS > _n_distrib) NUM_THREADS=_n_distrib;
	
	struct WeightedCovthread_data *thread_data_array = new WeightedCovthread_data[NUM_THREADS];
	pthread_t *threads = new pthread_t[NUM_THREADS];

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	unsigned long offset=_n_distrib/NUM_THREADS;
	
	unsigned long disBottom = 0;
	unsigned long disUp = 0;
	unsigned long re=_n_distrib - NUM_THREADS*offset;

	// Create temporary RefVector
	RefVector<DoubleSquareMatrix> tmpW;
	for(unsigned long s =0; s<_n_distrib; s++){
		tmpW.addObject(*new DoubleSquareMatrix(_rankT));
		tmpW[s].setAllValues(0.0);
	}

	//Create threads : one per distribution as a maximum
	for(unsigned long t=0; t<NUM_THREADS; t++){
	
		disUp = disBottom + offset;
		if(t<re) disUp +=1;

		thread_data_array[t].T = T;
		thread_data_array[t].weight=ubm_weight;
		thread_data_array[t].disBottom=disBottom;
		thread_data_array[t].disUp=disUp;
		thread_data_array[t].rankTV=_rankT;
		thread_data_array[t].svSize=_svSize;
		thread_data_array[t].vectSize=_vectSize;
		thread_data_array[t].W=&(tmpW);

		if (verboseLevel > 1) cout<<"(AccumulateTVStat) Creating thread n["<< t<< "] for distributions["<<disBottom<<"-->"<<disUp<<"]"<<endl;
		rc = pthread_create(&threads[t], &attr, WeightedCovthread, (void *)&thread_data_array[t]);
		if (rc) throw Exception("ERROR; return code from pthread_create() is ",__FILE__,rc);

		disBottom = disUp;
	}

	pthread_attr_destroy(&attr);

	for(unsigned long t=0; t<NUM_THREADS; t++) {
		rc = pthread_join(threads[t], (void **)&status);
		if (rc)  throw Exception("ERROR; return code from pthread_join() is ",__FILE__,rc);
		if (verboseLevel >1) cout <<"(AccumulateTVStat) Completed join with thread ["<<t<<"] status["<<status<<"]"<<endl;
	}

	// Accumulate all temporary matrices and fill the second half of the symetric matrix
	for(unsigned long i=0 ; i<_rankT ; i++){
		for(unsigned long j=i ; j<_rankT ; j++){
			for(unsigned long s =0; s<_n_distrib; s++){
				//vevt[i*_rankT+j]=vevt[j*_rankT+i];
				W(i,j) += tmpW[s](i,j);
			}
			W(j,i) = W(i,j);
		}
	}
	tmpW.deleteAllObjects();

	free(thread_data_array);
	free(threads);

	if (verboseLevel >= 1) cout << "(AccumulateTVStat) Done " << endl;
}
#endif


//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void TVAcc::computeEigenProblem(Matrix<double> &EP,Matrix<double> &eigenVect, long rank, Config& config)
{
	// EP has to be a square matrix
	Matrix <double> eigenVal;
	eigenVal.setDimensions (EP.rows(),EP.rows());
	eigenVal.setAllValues(0.0);
	computeEigenProblem(EP,eigenVect,eigenVal,rank,config);
}

#ifdef LAPACK
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void TVAcc::computeEigenProblem(Matrix<double> &EP,Matrix<double> &eigenVect,Matrix<double> &eigenVal, long rank, Config& config)
{
	unsigned long matSize = EP.rows();
	lapack_int n = matSize;

	double vl[matSize*matSize]; // left eigen vector - not used
	double vr[matSize*matSize]; // right eigen vector - the ones we want
	double wr[matSize];           // right eigen value
	double wi[matSize];           // left eigen valures - not needed

	lapack_int info;

	if(verboseLevel>2)      cout<<"         (PldaDev) Compute Eigen Problem using Lapack"<<endl;
	
	double * EPdata=EP.getArray();
	// call to lapackr 
	// 'N' : not interested in the left eigen vectors
	// 'V' : we want the right  eigen vectors
	info = LAPACKE_dgeev( LAPACK_ROW_MAJOR, 'N', 'V',matSize, EPdata, matSize, wr, wi,vl, matSize, vr, matSize );

	if (verboseLevel >2) cout << 		"--- Eigen Problem solved" <<endl;

	// get and check eigen values
	LKVector EV(0,0);
	for(int i=0;i<matSize;i++){ 
		if(wi[i]!=0)	cout << "WARNING eigenvalue [ " << i << "] has an imaginary part" << endl;
		LKVector::type s;
		s.idx = i;
		s.lk = wr[i];
		EV.addValue(s);
	}
	eigenVal.setAllValues(0.0);

	// Order the EigenValues
	EV.descendingSort();
	for(unsigned long k=0; k<matSize;k++){
		for(int j=0;j<rank;j++){
			eigenVect(k,j)=vr[k*matSize+EV[j].idx];
		}
	}

	for(int j=0;j<rank;j++){ 
		eigenVal(j,j) = EV[j].lk;
	}
	
	if(verboseLevel>3){
		cerr<<"EigenValues"<<endl;
		for(int i=0;i<rank;i++){ 
			cerr<<eigenVal(i,0)<<endl;
		}
	}
}

#else
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void TVAcc::computeEigenProblem(Matrix<double> &EP,Matrix<double> &eigenVect,Matrix<double> &eigenVal, long rank, Config& config)
{
	if(verboseLevel>2)	cout<<"		(PldaDev) Compute Eigen Problem"<<endl;

	// convert ALIZE matrix into Eigen::MatrixXd
	Eigen::MatrixXd A(EP.rows(),EP.cols());
	for(unsigned long i=0;i<EP.rows();i++){
		for(unsigned long j=0;j<EP.cols();j++)
		A(i,j) = EP(i,j);
	}	

	// Compute Eigen Decomposition
	Eigen::EigenSolver<Eigen::MatrixXd> es(A);
	complex<double> lambda = es.eigenvalues()[0];
	Eigen::MatrixXcd V = es.eigenvectors();

	// get and check eigen values
	LKVector EV(0,0);
	unsigned long imagPart = 0;
	for(unsigned long i=0;i<EP.rows();i++){ 
		if(imag(es.eigenvalues()[i])!=0)	imagPart++;
		LKVector::type s;
		s.idx = i;
		s.lk = real(es.eigenvalues()[i]);
		EV.addValue(s);
	}
	if(imagPart>0) cout << "WARNING "<<imagPart<<" eigenvalues have an imaginary part" << endl;
	eigenVal.setAllValues(0.0);

	// Order the EigenValues
	EV.descendingSort();
	for(unsigned long k=0; k<EP.rows();k++){
		for(int j=0;j<rank;j++){
			eigenVect(k,j)= real(V(k,j));
		}
	}
	for(int j=0;j<rank;j++){ 
		eigenVal(j,j) = EV[j].lk;
	}
	
	if(verboseLevel>3){
		cerr<<"EigenValues"<<endl;
		for(int i=0;i<rank;i++){ 
			cerr<<eigenVal(i,i)<<endl;
		}
	}
}
#endif

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void TVAcc::approximateTcTc(Matrix<double> &D, Matrix<double> &Q, Config &config){
	#ifdef THREAD          
	if (config.existsParam("numThread") && config.getParam("numThread").toULong() >0)	approximateTcTcThreaded(D,Q,config.getParam("numThread").toULong());
	else approximateTcTcUnThreaded(D,Q);
	#else
	approximateTcTcUnThreaded(D,Q);
	#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void TVAcc::approximateTcTcUnThreaded(Matrix<double> &D, Matrix<double> &Q){
	
	Matrix<double> A(_vectSize,_rankT);
		
	for (unsigned long cc=0 ; cc<_n_distrib ; cc++){

		A.setAllValues(0.0);
		double *a, *t, *q, *d;
		a = A.getArray(); t = _T.getArray(); q = Q.getArray(); d = D.getArray();

		for(unsigned long i=0 ; i<_vectSize ; i++)
			for(unsigned long j=0 ; j<_rankT ; j++)
				for(unsigned long k=0 ; k<_rankT; k++)
					a[i*_rankT+j] += t[(k*_svSize)+(cc*_vectSize)+i] * q[k*_rankT+j];

		// Compute Diagonal terms of D
		for(unsigned long i=0 ; i<_rankT ; i++)
			for(unsigned long k=0 ; k<_vectSize; k++)
				d[cc*_rankT+i] += a[k*_rankT+i]*a[k*_rankT+i];
	}
}

#ifdef THREAD
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Data structure of thread
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
struct approximateTcTcthread_data{

	double *T;
	double *Q;
	double *D;
	unsigned long disBottom;
	unsigned long disUp;	
	unsigned long rankTV;
	unsigned long svSize;
	unsigned long vectSize;
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Thread Routine
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void *approximateTcTcthread(void *threadarg){

	struct approximateTcTcthread_data *my_data;
	my_data = (struct approximateTcTcthread_data *) threadarg;
	
	double *t = my_data->T;
	double *d = my_data->D;
	double *q = my_data->Q;
	unsigned long disBottom = my_data->disBottom;	
	unsigned long disUp = my_data->disUp;
	
	unsigned long _rankT=my_data->rankTV;
	unsigned long _svSize=my_data->svSize;
	unsigned long _vectSize=my_data->vectSize;

	Matrix<double> A(_vectSize,_rankT);
	double *a = A.getArray();

	for (unsigned long cc=disBottom; cc<disUp; cc++){

		A.setAllValues(0.0);

		for(unsigned long i=0 ; i<_vectSize ; i++)
			for(unsigned long j=0 ; j<_rankT ; j++)
				for(unsigned long k=0 ; k<_rankT; k++)
					a[i*_rankT+j] += t[(k*_svSize)+(cc*_vectSize)+i] * q[k*_rankT+j];
		
		// Compute diagonal terms of D
		for(unsigned long i=0 ; i<_rankT ; i++)
			for(unsigned long k=0 ; k<_vectSize; k++)
				d[cc*_rankT+i] += a[k*_rankT+i]*a[k*_rankT+i];
	}

	pthread_exit((void*) 0);
	return (void*)0 ;
}

//-----------------------------------------------------------------------------------------
void TVAcc::approximateTcTcThreaded(Matrix<double> &D, Matrix<double> &Q, unsigned long NUM_THREADS){
	
	if (verboseLevel >=1) cout << "(AccumulateTVStat) Compute matrix D, approximation of Tc'*Tc Threaded"<<endl;
	if (NUM_THREADS==0) throw Exception("Num threads can not be 0",__FILE__,__LINE__);
	
	double *T=_T.getArray();
	double *d=D.getArray();
	double *q=Q.getArray();

	int rc, status;
	if (NUM_THREADS > _n_distrib) NUM_THREADS=_n_distrib;
	
	struct approximateTcTcthread_data *thread_data_array = new approximateTcTcthread_data[NUM_THREADS];
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

		thread_data_array[t].T = T;
		thread_data_array[t].D = d;
		thread_data_array[t].Q = q;
		thread_data_array[t].disBottom=disBottom;
		thread_data_array[t].disUp=disUp;
		thread_data_array[t].rankTV=_rankT;
		thread_data_array[t].svSize=_svSize;
		thread_data_array[t].vectSize=_vectSize;
	
		if (verboseLevel > 1) cout<<"(AccumulateTVStat) Creating thread n["<< t<< "] for distributions["<<disBottom<<"-->"<<disUp<<"]"<<endl;
		rc = pthread_create(&threads[t], &attr, approximateTcTcthread, (void *)&thread_data_array[t]);
		if (rc) throw Exception("ERROR; return code from pthread_create() is ",__FILE__,rc);

		disBottom = disUp;
	}

	pthread_attr_destroy(&attr);

	for(unsigned long t=0; t<NUM_THREADS; t++) {
		rc = pthread_join(threads[t], (void **)&status);
		if (rc)  throw Exception("ERROR; return code from pthread_join() is ",__FILE__,rc);
		if (verboseLevel >1) cout <<"(AccumulateTVStat) Completed join with thread ["<<t<<"] status["<<status<<"]"<<endl;
	}

	free(thread_data_array);
	free(threads);

	if (verboseLevel >= 1) cout << "(AccumulateTVStat) Done " << endl;
}
#endif
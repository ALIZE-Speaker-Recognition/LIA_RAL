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


using namespace alize;
using namespace std;

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
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

			//XList TVNdx;
			TVNdx.addLine()=tmpLine;
		}
		_init(TVNdx,config);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
TVAcc::TVAcc(XList & ndx,Config & config)
	:_ms(config),_ss(config){ // constructor
	_init(ndx,config);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
TVAcc::~TVAcc(){
	_vEvT.deleteAllObjects();
	_Aev.deleteAllObjects();
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
String TVAcc::getClassName() const{	return "TVAcc";}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void TVAcc::_init(XList &ndx, Config &config){

	///Convert the NDX file
	_fileList=ndx;
	_ndxTable=TVTranslate(ndx);

	///Load the UBM
	MixtureGD& UBM = _ms.loadMixtureGD(config.getParam("inputWorldFilename"));

	_vectSize = UBM.getVectSize();
	_n_distrib = UBM.getDistribCount();
	_svSize = _vectSize*_n_distrib;

	_rankEV=1;
	if(config.existsParam("totalVariabilityNumber"))
		_rankEV=config.getParam("totalVariabilityNumber").toULong();

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

	_meanY.setSize(_rankEV);
	_meanY.setAllValues(0.0);


	///Create and initialise statistics acumulators
	_matN= Matrix<double>(_n_speakers, _n_distrib);
	_F_X= Matrix<double>(_n_speakers, _n_distrib*_vectSize);
	_matN.setAllValues(0.0);
	_F_X.setAllValues(0.0);

	_V.setDimensions(_rankEV, _svSize);
	_V.setAllValues(0.0);

	_Y.setDimensions(_n_speakers, _rankEV);
	_Y.setAllValues(0.0);

	///Create the accumulators for L matrices computation
	for(unsigned long s =0; s<_n_distrib; s++){
		_vEvT.addObject(*new DoubleSquareMatrix(_rankEV));
		_vEvT[s].setAllValues(0.0);
	}

	///Create accumulators to compute the TotalVariability matrix
	// Initialize to dimension 1 that will be midified in real time during computation
	for(unsigned long d =0; d<_n_distrib; d++){
		_Aev.addObject(*new DoubleSquareMatrix(1));
		_Aev[d].setAllValues(0.0);
	}
	_Cev.setDimensions(_rankEV,_svSize);
	_Cev.setAllValues(0.0);

	_R.setSize(_rankEV);
	_r.setSize(_rankEV);

	_R.setAllValues(0.0);
	_r.setAllValues(0.0);

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
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

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void TVAcc::computeAndAccumulateTVStatUnThreaded(Config& config){

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
	 n=_matN.getArray();f_x=_F_X.getArray();
	
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
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Data strucutre of thread
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
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

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Thread Routine
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
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

		unsigned long begin=seg->begin()+fs.getFirstFeatureIndexOfASource(seg->sourceName()); 				/// Idx of the first frame of the current file in the feature server
		
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

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
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

	double *N=_matN.getArray();
	double *F_X=_F_X.getArray();

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


//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void TVAcc::computeAndAccumulateTVStat(FeatureServer &fs,Config & config){
	SegServer segmentsServer;
	LabelServer labelServer;
	initializeClusters(_fileList,segmentsServer,labelServer,config);
	verifyClusterFile(segmentsServer,fs,config);
	unsigned long codeSelectedFrame=labelServer.getLabelIndexByString(config.getParam("labelSelectedFrames"));
	SegCluster& selectedSegments=segmentsServer.getCluster(codeSelectedFrame); 
	this->computeAndAccumulateTVStat(selectedSegments,fs,config);
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void TVAcc::computeAndAccumulateTVStat(SegCluster &selectedSegments,FeatureServer &fs,Config & config){

	MixtureGD& UBM = _ms.getMixtureGD(0);
	MixtureGDStat &acc=_ss.createAndStoreMixtureStat(UBM);

	///Fast access to vector and matrices array
	double *n, *f_x;
	n=_matN.getArray(); f_x=_F_X.getArray();

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

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void TVAcc::resetAcc(){
	_F_X.setAllValues(0.0);
	_matN.setAllValues(0.0);	
	if (verboseLevel >= 1) cout << "# TV Accumulators reset" << endl;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void TVAcc::resetTmpAcc(){
	
	_Cev.setAllValues(0.0);

	///Reinitialise accumulators for L matrices computation
	for(unsigned long s =0; s<_n_distrib; s++){
		_vEvT[s].setAllValues(0.0);
		_Aev[s].setAllValues(0.0);

	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void TVAcc::loadEV(const String& evFilename, Config& config){	///load an TotalVariability Matrix
	String filename = config.getParam("matrixFilesPath") + evFilename +  config.getParam("loadMatrixFilesExtension");
	_V.load (filename, config);

	//Transpose the matrix if the number of line is higher than the numbe of columns
	if(_V.cols()<_V.rows()){
		cout<<"Load EV : number of lines ( "<<_V.rows() <<" ) higher than the number of columns ( "<<_V.cols()<<" ("<<endl;
		_V.transpose();
	}
	_rankEV=_V.rows();

	if(_rankEV != config.getParam("totalVariabilityNumber").toULong()){
		throw Exception("Incorrect dimension of TotalVariability Matrix",__FILE__,__LINE__);
	}

	cout << "(AccumulateTVStat) Init TV matrix from "<< filename <<"  for TotalVariability Matrix: "<<", rank: ["<<_V.rows() << "] sv size: [" << _V.cols() <<"]"<<endl;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void TVAcc::loadEV(Matrix<double> & V, Config& config){
	_rankEV=V.rows();
	_V=V;

	//Transpose the matrix if the number of line is higher than the numbe of columns
	if(_V.cols()<_V.rows()){
		cout<<"Load EV : number of lines higher than the number of columns"<<endl;
		_V.transpose();
	}

	unsigned long rankEV = 1;
	if(config.existsParam("totalVariabilityNumber"))
		rankEV = config.getParam("totalVariabilityNumber").toULong();

	if(_rankEV != config.getParam("totalVariabilityNumber").toULong()){
		throw Exception("Incorrect dimension of TotalVariability Matrix",__FILE__,__LINE__);
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void TVAcc::loadMeanEstimate(DoubleVector& meanEstimate){
	_ubm_means = meanEstimate;
	if(_ubm_means.size() != _svSize){
		throw Exception("Incorrect dimension of meanEstimate vector",__FILE__,__LINE__);
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void TVAcc::loadN(Config& config){
	String filname=config.getParam("matrixFilesPath") + config.getParam("nullOrderStatSpeaker") + config.getParam("loadMatrixFilesExtension");
	_matN.load (filname, config);

	if((_matN.rows() != _fileList.getLineCount()) || (_matN.cols() != _n_distrib)){
		throw Exception("Incorrect dimension of N Matrix",__FILE__,__LINE__);
	}
	cout<<"(AccumulateTVStat) ---------load statistics N [ "<<_matN.rows()<<" ] [ "<<_matN.cols()<<" ]"<<endl;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void TVAcc::loadF_X(Config& config){
	String filname=config.getParam("matrixFilesPath") + config.getParam("firstOrderStatSpeaker") + config.getParam("loadMatrixFilesExtension");
	_F_X.load (filname, config);

	if((_F_X.rows() != _fileList.getLineCount()) || (_F_X.cols() != _svSize)){
		throw Exception("Incorrect dimension of F_X Matrix",__FILE__,__LINE__);
	}
	cout<<"(AccumulateTVStat) ---------load statistics F_X [ "<<_F_X.rows()<<" ] [ "<<_F_X.cols()<<" ]"<<endl;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void TVAcc::initEV(Config& config){		///random initialisation of the TotalVariability Matrix
		_rankEV=config.getParam("totalVariabilityNumber").toULong();
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
		if (verboseLevel >=1) cout << "(AccumulateTVStat) Random Init for TotalVariability Matrix with "<<randomLaw<<" law "<<", rank: ["<<_V.rows() << "] sv size: [" << _V.cols() <<"]"<<endl;
}
catch (Exception& e){ 
	cout << e.toString().c_str() << endl;
}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void TVAcc::saveV(const String& filename, Config& config){
	String vName = config.getParam("matrixFilesPath") + filename +  config.getParam("saveMatrixFilesExtension");
	_V.save(vName, config);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void TVAcc::estimateVEVT(Config &config){
	
	#ifdef THREAD          
	if (config.existsParam("numThread") && config.getParam("numThread").toULong() >0)	estimateVEVTThreaded(config.getParam("numThread").toULong());
	else estimateVEVTUnThreaded();
	#else
	estimateVEVTUnThreaded();
	#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void TVAcc::estimateVEVTUnThreaded(){	
        if (verboseLevel >= 1) cout << "(AccumulateTVStat) compute V * Sigma-1 * V "<<endl;

	for(unsigned long d=0; d<_n_distrib; d++){

		Matrix<double>ssV= _V.crop(0,d*_vectSize, _rankEV,_vectSize);

		///Compute  _vEvT matrices
		double *vEvT, *v, *E;
		_vEvT[d].setAllValues(0.0);
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

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Thread Routine
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void *TETthread(void *threadarg) {
	struct TETthread_data *my_data;
	my_data = (struct TETthread_data *) threadarg;
	
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
void TVAcc::estimateVEVTThreaded(unsigned long NUM_THREADS){
	
	if (verboseLevel >=1) cout << "(AccumulateTVStat) Compute vEvT Threaded"<<endl;
	if (NUM_THREADS==0) throw Exception("Num threads can not be 0",__FILE__,__LINE__);
	
	double *V=_V.getArray();
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
		thread_data_array[t].rankEV=_rankEV;
		thread_data_array[t].svSize=_svSize;
		thread_data_array[t].vectSize=_vectSize;
		thread_data_array[t].vevT=&(_vEvT);

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

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void  TVAcc::getVY(DoubleVector &vy, String &file){		///Compute VY for the speaker corresponding to the given file given file
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
void  TVAcc::getVY(DoubleVector &vy, unsigned long spk) {		///Compute VY for speaker spk
	vy.setAllValues(0.0);
	
	double *v, *y; y=_Y.getArray(); v=_V.getArray();
	for (unsigned long i=0;i<_svSize;i++){
		for (unsigned long j=0;j<_rankEV;j++){
			vy[i] += v[j*_svSize + i] * y[spk*_rankEV + j ];
		}
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Compute supervector of client M_s_h=M+ VY_s
void TVAcc::getMplusVY(DoubleVector &Sp, String& file){
	Sp.setAllValues(0.0);
	DoubleVector vy(_svSize,_svSize);
	getVY(vy,file);		
	for (unsigned long i=0;i<_svSize;i++){
		Sp[i] = _ubm_means[i] + vy[i];
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Compute supervector of client M_s_h=M+ VY_s
void TVAcc::getMplusVY(DoubleVector &Sp, unsigned long spk){
	Sp.setAllValues(0.0);
	DoubleVector vy(_svSize,_svSize);
	getVY(vy,spk);
	for (unsigned long i=0;i<_svSize;i++){
		Sp[i]=_ubm_means[i] + vy[i];
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void TVAcc::estimateAndInverseL_EV(Config& config){
	
	#ifdef THREAD          
	if (config.existsParam("numThread") && config.getParam("numThread").toULong() >0)	estimateAndInverseLThreaded_EV(config.getParam("numThread").toULong(),config);
	else estimateAndInverseLUnThreaded_EV(config);
	#else
	estimateAndInverseLUnThreaded_EV(config);
	#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void TVAcc::estimateAndInverseLUnThreaded_EV(Config& config){
	
	if (verboseLevel >= 1) cout << "(AccumulateTVStat) Compute and Inverse L Matrix for TotalVariability "<<endl;
	
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
		
		//Inverse L
		DoubleSquareMatrix Linv(_rankEV);
		L.invert(Linv);

		// Save the matrix invL in Binary format
		String s,LinvFile;
		LinvFile = config.getParam("matrixFilesPath")+"Linv"+s.valueOf(spk)+config.getParam("saveMatrixFilesExtension");
		Matrix<double> ll(Linv);
		ll.save(LinvFile,config);
	}
}

#ifdef THREAD
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Data strucutre of thread
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
struct LTthread_data{

	double *N;
	unsigned long spkBottom;
	unsigned long spkUp;	
	unsigned long rankEV;
	unsigned long n_distrib;
	unsigned long nt;
	RefVector <DoubleSquareMatrix>* vevT;
	Config	config;
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Thread Routine
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void *LTthread(void *threadarg) {
	struct LTthread_data *my_data;
	my_data = (struct LTthread_data *) threadarg;

	double *n = my_data->N;
	unsigned long spkBottom = my_data->spkBottom;	
	unsigned long spkUp = my_data->spkUp;
	unsigned long _rankEV=my_data->rankEV;
	unsigned long _n_distrib=my_data->n_distrib;
	Config config = my_data->config;

	DoubleSquareMatrix L(_rankEV);
	double *l=L.getArray();
	DoubleSquareMatrix Linv(_rankEV);

	for(unsigned long spk=spkBottom; spk<spkUp; spk++){
		
		L.setAllValues(0.0);	
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
		
		//Inverse L and stock it in _l_spk_inv[spk]
		L.invert(Linv);

		// Save the matrix invL in Binary format
		String s,LinvFile;
		LinvFile = config.getParam("matrixFilesPath")+"Linv"+s.valueOf(spk)+config.getParam("saveMatrixFilesExtension");
		Matrix<double> ll(Linv);
		ll.save(LinvFile,config);
	}
	pthread_exit((void*) 0);
	return (void*)0 ;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void TVAcc::estimateAndInverseLThreaded_EV(unsigned long NUM_THREADS,Config& config){
	
	if (verboseLevel >= 1) cout << "(AccumulateTVStat) Compute and Inverse L Matrix for TotalVariability Threaded"<<endl;
	if (NUM_THREADS==0) throw Exception("Num threads can not be 0",__FILE__,__LINE__);

	int rc, status;
	if (NUM_THREADS > _n_speakers) NUM_THREADS=_n_speakers;

	struct LTthread_data *thread_data_array = new LTthread_data[NUM_THREADS];
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
		thread_data_array[t].config=config;
		thread_data_array[t].nt=t;

		if (verboseLevel>1) cout<<"(AccumulateTVStat) Creating thread n["<< t<< "] for speakers["<<spkBottom<<"-->"<<spkUp-1<<"]"<<endl;
		rc = pthread_create(&threads[t], &attr, LTthread, (void *)&thread_data_array[t]);
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


//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void TVAcc::estimateYandV(Config& config){
	#ifdef THREAD          
	if (config.existsParam("numThread") && config.getParam("numThread").toULong() >0)	estimateYandVThreaded(config.getParam("numThread").toULong(),config);
	else estimateYandVUnThreaded(config);
	#else
	estimateYandVUnThreaded(config);
	#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void TVAcc::estimateYandVUnThreaded(Config& config){	//estimate Y for all speakers and V
	
	if (verboseLevel >=1) cout << "(AccumulateTVStat) Compute Y  and V Estimate "<<endl;
	_Y.setAllValues(0.0);
	Matrix<double> AUX(_n_speakers,_rankEV);
	AUX.setAllValues(0.0);
	
	_R.setAllValues(0.0);
	_r.setAllValues(0.0);
	_meanY.setAllValues(0.0);

	double *y, *v, *f_x, *aux, *invVar, *n, *c, *R, *r, *meanY;
	y=_Y.getArray(); v=_V.getArray(); f_x=_F_X.getArray(); aux=AUX.getArray(); invVar=_ubm_invvar.getArray(); n=_matN.getArray(); c=_Cev.getArray();; R=_R.getArray(); r=_r.getArray(); meanY=_meanY.getArray();

	for(unsigned long d =0; d<_n_distrib; d++){
		_Aev[d].setSize(_rankEV);
		_Aev[d].setAllValues(0.0);
	}

	//For each speaker
	for(unsigned long spk=0; spk<_n_speakers; spk++){

		//Load L matrix of the current speaker
		String s, LinvFile;
		LinvFile = config.getParam("matrixFilesPath")+"Linv"+s.valueOf(spk)+config.getParam("saveMatrixFilesExtension");
		Matrix<double> Linv;
		Linv.load (LinvFile, config);
		double *invl = Linv.getArray();

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

		for(unsigned long k=0;k<_rankEV;k++){
			meanY[k] += y[spk*_rankEV+k];
		}

		for(unsigned long i=0;i<_rankEV;i++){
			for(unsigned long j=0;j<_rankEV;j++){
					invl[i*_rankEV+j] += y[spk*_rankEV+i]*y[spk*_rankEV+j];
					//Update the Minimum Divergence Accumulator
					R[i*_rankEV+j] += invl[i*_rankEV+j];
			}
			r[i] += y[spk*_rankEV+i];
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
	// Compute the mean of iVectors
	for(unsigned long k=0;k<_rankEV;k++){
		meanY[k] /= _n_speakers;
	}
}

#ifdef THREAD
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Data strucutre of thread
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
struct estimateWandTthread_data{
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
	RefVector <DoubleVector>* meanY;
	Config config;
};


struct TVcomputeC_data{
	double* tmpC;
	double *Y;
	double *F_X;
	unsigned long numThread;
	unsigned long spkBottom;
	unsigned long spkUp;	
	unsigned long rankEV;
	unsigned long svSize;
};


struct TVcomputeINVLandA_data{
	double *N;
	double *Y;
	double* R;
	double* r;
	double* L;
	unsigned long spk;
	unsigned long _n_distrib;
	unsigned long numThread;
	unsigned long rankBottom;
	unsigned long rankUp;	
	unsigned long rankEV;
	unsigned long svSize;
	unsigned long nt;
	RefVector <DoubleSquareMatrix>* A;
//	Config config;
};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Thread Routine
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void *estimateWandTthread(void *threadarg) {
	
	struct estimateWandTthread_data *my_data;
	my_data = (struct estimateWandTthread_data *) threadarg;

	double *v = my_data->V;
	double *y = my_data->Y;
	double *f_x = my_data->F_X;
	double *aux = my_data->AUX;
	double *ubm_invvar = my_data->ubm_invvar;
	unsigned long spkBottom = my_data->spkBottom;	
	unsigned long spkUp = my_data->spkUp;
	unsigned long _rankEV=my_data->rankEV;
	unsigned long _svSize=my_data->svSize;
	unsigned long nt = my_data->nt;
	DoubleVector &meanY=(*(my_data->meanY))[nt];
	Config config = my_data->config;
	
	double *mY = meanY.getArray();
	
	//For each session
	for(unsigned long spk= spkBottom; spk<spkUp; spk++){

		//Load current session L Matrix
		String s, LinvFile;
		LinvFile = config.getParam("matrixFilesPath")+"Linv"+s.valueOf(spk)+config.getParam("saveMatrixFilesExtension");
		Matrix<double> Linv;
		Linv.load (LinvFile, config);
		double *invl = Linv.getArray();


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
		for(unsigned long k=0;k<_rankEV;k++){
			mY[k] += y[spk*_rankEV+k];
		}
	}
	pthread_exit((void*) 0);
	return (void*)0 ;
}


void *TVcomputeC(void *threadarg){
	
	struct TVcomputeC_data *my_data;
	my_data = (struct TVcomputeC_data *) threadarg;

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

void *TVcomputeINVLandA(void *threadarg){

	struct TVcomputeINVLandA_data *my_data;
	my_data = (struct TVcomputeINVLandA_data *) threadarg;

	double *n = my_data->N;
	double *y = my_data->Y;
	double *R = my_data->R;
	double *r = my_data->r;
	double *invl = my_data->L;
	unsigned long spk = my_data->spk;
	unsigned long _n_distrib = my_data->_n_distrib;
	unsigned long rankBottom = my_data->rankBottom;	
	unsigned long rankUp = my_data->rankUp;
	unsigned long _rankEV=my_data->rankEV;
	unsigned long nt = my_data->nt;
	
	for(unsigned long i= rankBottom; i<rankUp; i++){
		for(unsigned long j=0;j<_rankEV;j++){
			invl[i*_rankEV+j] += y[spk*_rankEV+i]*y[spk*_rankEV+j];
			//Update the Minimum Divergence Accumulator
			R[i*_rankEV+j] += invl[i*_rankEV+j];
			for(unsigned long dis = 0; dis<_n_distrib;dis++){
				DoubleSquareMatrix &A=(*(my_data->A))[dis];
				double *a = A.getArray();
				a[i*_rankEV+j] += invl[i*_rankEV+j] * n[spk*_n_distrib+dis];
			}
			r[i] += y[spk*_rankEV+i];
		}
	}
	pthread_exit((void*) 0);
	return (void*)0 ;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void TVAcc::estimateYandVThreaded(unsigned long NUM_THREADS,Config& config){

	if (verboseLevel >= 1) cout << "(AccumulateTVStat) Estimate Y and V for each session Threaded"<<endl;
	if (NUM_THREADS==0) throw Exception("Num threads can not be 0",__FILE__,__LINE__);

	int rc, status;
	if (NUM_THREADS > _n_speakers) NUM_THREADS=_n_speakers;

	struct estimateWandTthread_data *thread_data_array = new estimateWandTthread_data[NUM_THREADS];
	pthread_t *threads = new pthread_t[NUM_THREADS];

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	unsigned long offset=_n_speakers/NUM_THREADS;

	_Y.setAllValues(0.0);
	Matrix<double> AUX(_n_speakers,_rankEV);
	AUX.setAllValues(0.0);

	_R.setAllValues(0.0);
	_r.setAllValues(0.0);
	_meanY.setAllValues(0.0);
	
	double *y, *f_x, *n, *c;
	y=_Y.getArray(); f_x=_F_X.getArray(); n=_matN.getArray(); c=_Cev.getArray();

	//Temporary accumulatores for Minimum Divergence
	RefVector<DoubleVector> meanY;
	for(unsigned long nt=0;nt<NUM_THREADS;nt++){
		meanY.addObject(*new DoubleVector(_rankEV,_rankEV));
		meanY[nt].setAllValues(0.0);
	}

	double *V =_V.getArray();
	double *Y =_Y.getArray();
	double *F_X =_F_X.getArray();
	double *aux = AUX.getArray(); 
	double *ubm_invvar=_ubm_invvar.getArray();
	double *R = _R.getArray();
	double *r = _r.getArray();

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
		thread_data_array[t].nt=t;
		thread_data_array[t].meanY=&(meanY);
		thread_data_array[t].config = config;

		if (verboseLevel >1) cout<<"(AccumulateTVStat) Creating thread n [ "<< t<< " ] for speakers[ "<<spkBottom<<" --> "<<spkUp-1<<" ]"<<endl;
		rc = pthread_create(&threads[t], &attr, estimateWandTthread, (void *)&thread_data_array[t]);
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


	//Compute the i-vector mean
	for(unsigned long mt=0; mt<NUM_THREADS;mt++){
		for(unsigned long i=0;i<_rankEV;i++){
			_meanY[i] += meanY[mt][i];
		}
	}
	meanY.deleteAllObjects();

	//Create a temporary matrix in by concatenating C matrices
	Matrix<double> _tmpC;
	_tmpC.setDimensions(_rankEV*NUM_THREADS,_svSize);
	_tmpC.setAllValues(0.0);

	struct TVcomputeC_data *thread_data_array2 = new TVcomputeC_data[NUM_THREADS];
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


		if (verboseLevel >1) cout<<"(AccumulateTVStat) Compute C creating thread n [ "<< t<< " ] for speakers[ "<<spkBottom<<" --> "<<spkUp-1<<" ]"<<endl;
			rc = pthread_create(&threads2[t], &attr2, TVcomputeC, (void *)&thread_data_array2[t]);
		if (rc) throw Exception("ERROR; return code from pthread_create() is ",__FILE__,rc);

		spkBottom = spkUp;
	}
	
	pthread_attr_destroy(&attr2);
	for(unsigned long t=0; t<NUM_THREADS; t++) {
		rc = pthread_join(threads2[t], (void **)&status);
		if (rc)  throw Exception("ERROR; return code from pthread_join() is ",__FILE__,rc);
		if (verboseLevel>1) cout <<"(AccumulateTVStat) Completed join with thread ["<<t<<"] status["<<status<<"]"<<endl;
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
	_tmpC.setDimensions(1,1);

	// Modify the dimension of A accumulators
	for(unsigned long d =0; d<_n_distrib; d++){
		_Aev[d].setSize(_rankEV);
		_Aev[d].setAllValues(0.0);
	}


	for(unsigned long spk=0; spk<_n_speakers; spk++){
		struct TVcomputeINVLandA_data *thread_data_array3 = new TVcomputeINVLandA_data[NUM_THREADS];
		pthread_t *threads3 = new pthread_t[NUM_THREADS];

		pthread_attr_t attr3;
		pthread_attr_init(&attr3);
		pthread_attr_setdetachstate(&attr3, PTHREAD_CREATE_JOINABLE);
		unsigned long offset3=_rankEV/NUM_THREADS;

		unsigned long rankBottom = 0;
		unsigned long rankUp=0;
		re=_rankEV - NUM_THREADS*offset3;

//*******************************
// Load L matrix before thread

	String s, LinvFile;
	LinvFile = config.getParam("matrixFilesPath")+"Linv"+s.valueOf(spk)+config.getParam("saveMatrixFilesExtension");
	Matrix<double> Linv;
	Linv.load (LinvFile, config);
	double *L = Linv.getArray();
//*******************************

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
			thread_data_array[t].nt=t;
			thread_data_array3[t].A = &(_Aev);
			thread_data_array3[t].R= R;
			thread_data_array3[t].r= r;

//************************** modif to load  matrix L before thread
			thread_data_array3[t].L = L;
//			thread_data_array3[t].config = config;

			if (verboseLevel >2) cout<<"(AccumulateTVStat) ComputeLinvandA creating thread n [ "<< t<< " ] for speakers[ "<<rankBottom<<" --> "<<rankUp-1<<" ]"<<endl;
				rc = pthread_create(&threads3[t], &attr3, TVcomputeINVLandA, (void *)&thread_data_array3[t]);
			if (rc) throw Exception("ERROR; return code from pthread_create() is ",__FILE__,rc);
			rankBottom = rankUp;
		}

		pthread_attr_destroy(&attr3);
		for(unsigned long t=0; t<NUM_THREADS; t++) {
			rc = pthread_join(threads3[t], (void **)&status);
			if (rc)  throw Exception("ERROR; return code from pthread_join() is ",__FILE__,rc);
			if (verboseLevel>2) cout <<"(AccumulateTVStat) Completed join with thread ["<<t<<"] status["<<status<<"]"<<endl;
		}
		free(thread_data_array3);
		free(threads3);
	}

	if (verboseLevel >= 1) cout << "(AccumulateTVStat) Done " << endl;
}
#endif

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void TVAcc::updateVestimate(){

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

	for(unsigned long d =0; d<_n_distrib; d++){
		_Aev[d].setSize(1);
		_Aev[d].setAllValues(0.0);
	}

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
XList& TVAcc::getXList(){
	return(_fileList);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned long TVAcc::getNSpeakers(){
	return(_n_speakers);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned long TVAcc::getNSessions(){
	return(_n_sessions);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned long TVAcc::getNDistrib(){
	return(_n_distrib);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned long TVAcc::getVectSize(){
	return(_vectSize);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned long TVAcc::getSvSize(){
	return(_svSize);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
unsigned long TVAcc::getRankEV(){
	return(_rankEV);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
DoubleVector& TVAcc::getUbmMeans(){
	return(_ubm_means);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
DoubleVector& TVAcc::getUbmInvVar(){
	return(_ubm_invvar);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
Matrix<double> TVAcc::getV(){
	return(_V);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
Matrix<double> TVAcc::getY(){
	return(_Y);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void TVAcc::saveY(String yFile,Config &config){
	_Y.save(yFile,config);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
Matrix <double>& TVAcc::getN() {
	return _matN;
}	
	
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
Matrix <double>& TVAcc::getF_X() {
	return _F_X;
}	
	
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void TVAcc::substractM(Config & config){
	#ifdef THREAD          
	if (config.existsParam("numThread") && config.getParam("numThread").toULong() >0)	substractMThreaded(config.getParam("numThread").toULong());
	else substractMUnThreaded();
	#else
	substractMUnThreaded();
	#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void TVAcc::substractMUnThreaded(){
	
	if (verboseLevel >= 1) cout <<"(AccumulateTVStat) Substract Speaker FA Stats M " << endl;	

	double *F_X = _F_X.getArray();
	
	for(unsigned long spk=0; spk<_n_speakers; spk++){

		double *n=_matN.getArray();
		
		//Substract M+DZ
		for(unsigned long i=0; i<_n_distrib; i++){
			for(unsigned long j = 0; j< _vectSize;j++){
				F_X[spk*_svSize+i*_vectSize+j] -= _ubm_means[i*_vectSize+j]*n[spk*_n_distrib+i];
			}
		}
	}
}

#ifdef THREAD
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Data strucutre of thread
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
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

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Thread Routine
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
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
		//Substract M+DZ
		for(unsigned long i=0; i<_n_distrib; i++){
			for(unsigned long j = 0; j< _vectSize;j++){
				f_x[spk*_svSize+i*_vectSize+j] -= ubm_means[i*_vectSize+j]*n[spk*_n_distrib+i];
			}
		}
	}
	pthread_exit((void*) 0);
	return (void*)0 ;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void TVAcc::substractMThreaded(unsigned long NUM_THREADS){
	
	if (verboseLevel >= 1) cout <<"(AccumulateTVStat) Substract Speaker FA Stats M Threaded" << endl;
	if (NUM_THREADS==0) throw Exception("Num threads can not be 0",__FILE__,__LINE__);

	int rc, status;
	if (NUM_THREADS > _n_speakers) NUM_THREADS=_n_speakers;

	struct S_Mthread_data *thread_data_array = new S_Mthread_data[NUM_THREADS];
	pthread_t *threads = new pthread_t[NUM_THREADS];

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	unsigned long offset=_n_speakers/NUM_THREADS;
	
	double *N=_matN.getArray(); 
	double *F_X=_F_X.getArray();

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

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void TVAcc::substractMplusVY(Config &config){
	#ifdef THREAD          
	if (config.existsParam("numThread") && config.getParam("numThread").toULong() >0)	substractMplusVYThreaded(config.getParam("numThread").toULong());
	else substractMplusVYUnThreaded();
	#else
	substractMplusVYUnThreaded();
	#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void TVAcc::substractMplusVYUnThreaded(){
	
	if (verboseLevel >= 1) cout <<"(AccumulateTVStat) Compute and Substract Speaker FA Stats M + VY" << endl;	

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

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Thread Routine
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
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
void TVAcc::substractMplusVYThreaded(unsigned long NUM_THREADS){
	
	if (verboseLevel >= 1) cout <<"(AccumulateTVStat) Compute and Substract Speaker FA Stats M + VY Threaded" << endl;
	if (NUM_THREADS==0) throw Exception("Num threads can not be 0",__FILE__,__LINE__);

	int rc, status;
	if (NUM_THREADS > _n_speakers) NUM_THREADS=_n_speakers;

	struct S_MplusTWthread_data *thread_data_array = new S_MplusTWthread_data[NUM_THREADS];
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

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void TVAcc::getSpeakerModel(MixtureGD &mixture, String& file){
	
	DoubleVector Sp, MplusVY, UX;
	Sp.setSize(_svSize); Sp.setAllValues(0.0);
	MplusVY.setSize(_svSize);

	this->getMplusVY(MplusVY, file);
	
	for(unsigned long i=0; i<_svSize; i++){
		Sp[i] = MplusVY[i];
	}
	svToModel(Sp,mixture);
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void TVAcc::orthonormalizeV(){
	// Gram-Schmidt algorithm

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
			for(unsigned long k=0; k<_svSize;k++){
				R(i,j) += Q(i,k)*_V(j,k);
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

	//Copy the new matrix in _V
	for(unsigned long i=0;i<_rankEV;i++){
		for(unsigned long j=0;j<_svSize;j++){
			_V(i,j) = Q(i,j);
		}
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
///Save Accumulators on disk
void TVAcc::saveAccs(Config &config) {

	String fxName, fxhName, nName, nhName;
	fxName = "F_X.mat"; nName = "N.mat";
	if(config.existsParam("nullOrderStatSpeaker"))	nName = config.getParam("matrixFilesPath") + config.getParam("nullOrderStatSpeaker")+config.getParam("saveMatrixFilesExtension");
	if(config.existsParam("firstOrderStatSpeaker")) fxName = config.getParam("matrixFilesPath") + config.getParam("firstOrderStatSpeaker")+config.getParam("saveMatrixFilesExtension");

	_F_X.save(fxName,config);
	_matN.save(nName,config);
	if (verboseLevel>=1) cout << "(AccumulateTVStat) FA Accs states saved" << endl;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
double TVAcc::getLLK(SegCluster &selectedSegments,MixtureGD &model,FeatureServer&fs,Config & config){

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

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
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

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void TVAcc::estimateAandC(Config& config){
	
	#ifdef THREAD          
	if (config.existsParam("numThread") && config.getParam("numThread").toULong() >0)	estimateAandCThreaded(config.getParam("numThread").toULong());
	else estimateAandCUnthreaded(config);
	#else
	estimateAandCUnthreaded(config);
	#endif
}

//-------------------------------------------------------------------------------------------------------------------------
void TVAcc::estimateAandCUnthreaded(Config& config){

	if (verboseLevel >= 1) cout << "(AccumulateTVStat) Compute and Inverse L Matrix for TotalVariability "<<endl;
	
	DoubleSquareMatrix L(_rankEV);
	_Y.setAllValues(0.0);
	Matrix<double> AUX(_n_speakers,_rankEV);
	AUX.setAllValues(0.0);

	_R.setAllValues(0.0);
	_r.setAllValues(0.0);
	_meanY.setAllValues(0.0);
	
	double *y, *v, *f_x, *aux, *invVar, *n, *c, *R, *r, *meanY;

	y=_Y.getArray(); v=_V.getArray(); f_x=_F_X.getArray(); aux=AUX.getArray(); invVar=_ubm_invvar.getArray(); c=_Cev.getArray();n=_matN.getArray(); R=_R.getArray(); r=_r.getArray(); meanY=_meanY.getArray();

	for(unsigned long d =0; d<_n_distrib; d++){
		_Aev[d].setSize(_rankEV);
		_Aev[d].setAllValues(0.0);
	}

	for(unsigned long spk=0; spk<_n_speakers; spk++){	
		L.setAllValues(0.0);
		for(unsigned long i=0; i<_rankEV; i++){	L(i,i)=1.0;}

		double *l;
		l=L.getArray();

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
		
		//Inverse L
		DoubleSquareMatrix Linv(_rankEV);
		L.invert(Linv);
		double *invl = Linv.getArray();

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

		for(unsigned long k=0;k<_rankEV;k++){
			meanY[k] += y[spk*_rankEV+k];
		}

		for(unsigned long i=0;i<_rankEV;i++){
			for(unsigned long j=0;j<_rankEV;j++){
					invl[i*_rankEV+j] += y[spk*_rankEV+i]*y[spk*_rankEV+j];
					//Update the Minimum Divergence Accumulator
					R[i*_rankEV+j] += invl[i*_rankEV+j];
			}
			r[i] += y[spk*_rankEV+i];
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

	// Compute the mean of iVectors
	for(unsigned long k=0;k<_rankEV;k++){
		meanY[k] /= _n_speakers;
	}

}

#ifdef THREAD
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Data strucutre of thread
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
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

	double *AUX;
	double *Y;
	double *V;
	double *F_X;
	double *ubm_invvar;

	double* tmpC;
	double* tmpA;
	unsigned long numThread;

};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Thread Routine
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void *estimateAandCTthread(void *threadarg) {
	struct estimateAandCTthread_data *my_data;
	my_data = (struct estimateAandCTthread_data *) threadarg;

	double *n = my_data->N;
	unsigned long spkBottom = my_data->spkBottom;	
	unsigned long spkUp = my_data->spkUp;
	unsigned long _rankEV=my_data->rankEV;
	unsigned long _n_distrib=my_data->n_distrib;
	unsigned long svSize =my_data->svSize;
	unsigned long nt = my_data->nt;

	unsigned long numThread = my_data->numThread;
	double *aux =my_data->AUX;
	double *v = my_data->V;
	double *y = my_data->Y;
	double *f_x = my_data->F_X;
	double *ubm_invvar = my_data->ubm_invvar;
	double* tmpC = my_data->tmpC;
	double* tmpA = my_data->tmpA;

	DoubleSquareMatrix &tmpR=(*(my_data->tmpR))[nt];
	double *R = tmpR.getArray();
	DoubleVector &tmpr=(*(my_data->tmpr))[nt];
	double *r = tmpr.getArray();
	DoubleVector &meanY=(*(my_data->meanY))[nt];
	double *mY = meanY.getArray();

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

		//DoubleSquareMatrix &linv=(*(my_data->linv))[spk];		
		DoubleSquareMatrix Linv(_rankEV);
		//Inverse L and stock it in _l_spk_inv[spk]
		L.invert(Linv);
		double *invl = Linv.getArray();

		for(unsigned long i=0;i<_rankEV;i++){
			for(unsigned long k=0;k<svSize;k++) {
				aux[spk*_rankEV+i] += f_x[spk*svSize+k] * ubm_invvar[k] * v[i*svSize+k];
			}
		}
	
		//multiplication by invL
		for(unsigned long i=0; i<_rankEV;i++){
			for(unsigned long k=0; k<_rankEV; k++){
				y[spk*_rankEV+i] += aux[spk*_rankEV+k] * invl[i*_rankEV+k];
			}
		}

		for(unsigned long k=0;k<_rankEV;k++){
			mY[k] += y[spk*_rankEV+k];
		}

		for(unsigned long i=0;i<_rankEV;i++){
			for(unsigned long j=0;j<_rankEV;j++){
				invl[i*_rankEV+j] += y[spk*_rankEV+i]*y[spk*_rankEV+j];
				//Update the Minimum Divergence Accumulator
				R[i*_rankEV+j] += invl[i*_rankEV+j];
			}
			r[i] += y[spk*_rankEV+i];
		}

		for(unsigned long dis = 0; dis<_n_distrib;dis++){
			for(unsigned long i=0;i<_rankEV;i++){
				for(unsigned long j = 0;j<_rankEV;j++){
					tmpA[nt * (_n_distrib*_rankEV*_rankEV) + dis*(_rankEV*_rankEV)+ i*(_rankEV) +j] += invl[i*_rankEV+j] * n[spk*_n_distrib+dis]; 
				}
			}
		}
		for(unsigned long i=0;i<_rankEV;i++){
			for(unsigned long j=0;j<svSize;j++){
				tmpC[(nt*_rankEV+i)*svSize+j] += y[spk*_rankEV+i] * f_x[spk*svSize+j];	
			}
		}
	}
	pthread_exit((void*) 0);
	return (void*)0 ;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void TVAcc::estimateAandCThreaded(unsigned long NUM_THREADS){
	
	if (verboseLevel >= 1) cout << "(AccumulateTVStat) Compute and Inverse L Matrix for TotalVariability Threaded"<<endl;
	if (NUM_THREADS==0) throw Exception("Num threads can not be 0",__FILE__,__LINE__);

	int rc, status;
	if (NUM_THREADS > _n_speakers) NUM_THREADS=_n_speakers;

	struct estimateAandCTthread_data *thread_data_array = new estimateAandCTthread_data[NUM_THREADS];
	pthread_t *threads = new pthread_t[NUM_THREADS];

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	unsigned long offset=_n_speakers/NUM_THREADS;

	_Y.setAllValues(0.0);
	Matrix<double> AUX(_n_speakers,_rankEV);
	AUX.setAllValues(0.0);

	_R.setAllValues(0.0);
	_r.setAllValues(0.0);
	_meanY.setAllValues(0.0);
	
	double *N =_matN.getArray(); 
	double *V = _V.getArray();
	double *aux = AUX.getArray();
	double *Y = _Y.getArray();
	double *F_X = _F_X.getArray();
	double *ubm_invvar = _ubm_invvar.getArray();
	double *c = _Cev.getArray();

	//Temporary tmpC matrix
	Matrix<double> _tmpC;
	_tmpC.setDimensions(_rankEV*NUM_THREADS,_svSize);
	_tmpC.setAllValues(0.0);
	double *tmpC = _tmpC.getArray();
	//Temporary tmpA matrix
	Matrix<double> _tmpA;
	_tmpA.setDimensions(NUM_THREADS,_n_distrib*_rankEV*_rankEV);
	_tmpA.setAllValues(0.0);
	double *tmpA = _tmpA.getArray();
	//Temporary accumulatores for Minimum Divergence
	RefVector<DoubleSquareMatrix> tmpR;
	RefVector<DoubleVector> tmpr;
	RefVector<DoubleVector> meanY;
	for(unsigned long nt=0;nt<NUM_THREADS;nt++){
		tmpR.addObject(*new DoubleSquareMatrix(_rankEV));
		tmpR[nt].setAllValues(0.0);
		tmpr.addObject(*new DoubleVector(_rankEV,_rankEV));
		tmpr[nt].setAllValues(0.0);
		meanY.addObject(*new DoubleVector(_rankEV,_rankEV));
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
		thread_data_array[t].rankEV=_rankEV;
		thread_data_array[t].n_distrib=_n_distrib;
		thread_data_array[t].svSize=_svSize;
		thread_data_array[t].vevT=&(_vEvT);
		thread_data_array[t].nt=t;
		thread_data_array[t].numThread=NUM_THREADS;
		thread_data_array[t].AUX=aux;
		thread_data_array[t].V=V;
		thread_data_array[t].Y=Y;
		thread_data_array[t].F_X=F_X;
		thread_data_array[t].ubm_invvar=ubm_invvar;
		thread_data_array[t].tmpC=tmpC;
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

	//Sum matrix C after multithreading
	for(unsigned long mt=0; mt<NUM_THREADS;mt++){
		for(unsigned long i=0; i<_rankEV;i++){
			for(unsigned long j=0; j<_svSize;j++){
				c[i*_svSize+j] += _tmpC(mt*_rankEV+i,j);
			}
		}
	}

	for(unsigned long d =0; d<_n_distrib; d++){
		_Aev[d].setSize(_rankEV);
		_Aev[d].setAllValues(0.0);
	}

	//Sum matrices A after multithreading
	for(unsigned long dis = 0; dis<_n_distrib;dis++){
		for(unsigned long i=0;i<_rankEV;i++){
			for(unsigned long j = 0;j<_rankEV;j++){
				for(unsigned long mt=0; mt<NUM_THREADS;mt++){
					_Aev[dis](i,j) += _tmpA(mt,dis*(_rankEV*_rankEV)+ i*(_rankEV) +j);
				}
			}
		}
	}

	//Sum Minimum Divergence Accumulators after multithreading
	for(unsigned long mt=0; mt<NUM_THREADS;mt++){
		for(unsigned long i=0;i<_rankEV;i++){
			for(unsigned long j = 0;j<_rankEV;j++){
				_R(i,j) += tmpR[mt](i,j);
			}
			_r[i] +=tmpr[mt][i];
			_meanY[i] += meanY[mt][i];
		}
	}

	// Compute the mean of iVectors
	for(unsigned long k=0;k<_rankEV;k++){
		_meanY[k] /= _n_speakers;
	}

	if (verboseLevel >= 1) cout << "(AccumulateTVStat) Done " << endl;
}
#endif


void TVAcc::minDivergence(){
	
	for(unsigned long i=0;i<_rankEV;i++){
		_r[i] /= _n_sessions;
	}

	for(unsigned long i=0;i<_rankEV;i++){
		for(unsigned long j=0;j<_rankEV;j++){
			_R(i,j) = _R(i,j)/_n_sessions - (_r[i]*_r[j]);
		}
	}

	DoubleSquareMatrix Ch;
	Ch.setSize(_rankEV);
	_R.upperCholesky(Ch);
	
	DoubleVector newMean(_svSize,_svSize);
	newMean = _ubm_means;
	for(unsigned long j=0;j<_svSize;j++){
		for(unsigned long k=0;k<_rankEV;k++){
			newMean[j] += _meanY[k]*_V(k,j);
		}
	}
	_ubm_means = newMean;

	Matrix<double> tmpV;
	tmpV.setDimensions(_V.rows(),_V.cols());
	tmpV.setAllValues(0.0);
	for(unsigned long i=0;i<_rankEV;i++){
		for(unsigned long j=0;j<_V.cols();j++){
			for(unsigned long k=0;k<_rankEV;k++){
				tmpV(i,j) += Ch(i,k)*_V(k,j);
			}
		}
	}
	for(unsigned long i=0;i<_rankEV;i++){
		for(unsigned long j=0;j<_V.cols();j++){
			_V(i,j) = tmpV(i,j);
		}
	}
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void TVAcc::storeAccs(){ // save JFA state in temporary variables
		_cF_X=_F_X;
		_cN=_matN;
		if (verboseLevel>=1) cout << "(AccumulateJFAStat) JFA Accs states stored" << endl;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void TVAcc::restoreAccs() {	// Restore Accumulators in temporary variables
	_matN=_cN;
	_F_X=_cF_X;
	if (verboseLevel>=1) cout << "(AccumulateJFAStat) JFA Accs states restored" << endl;			
}

#endif

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void TVAcc::estimateY(Config& config){
	
	#ifdef THREAD          
	if (config.existsParam("numThread") && config.getParam("numThread").toULong() >0)	estimateYThreaded(config.getParam("numThread").toULong());
	else estimateYUnThreaded(config);
	#else
	estimateYUnThreaded(config);
	#endif
}

//-------------------------------------------------------------------------------------------------------------------------
void TVAcc::estimateYUnThreaded(Config& config){

	if (verboseLevel >= 1) cout << "(AccumulateTVStat) Estimate Y from Statistics "<<endl;
	
	DoubleSquareMatrix L(_rankEV);
	_Y.setAllValues(0.0);
	Matrix<double> AUX(_n_speakers,_rankEV);
	AUX.setAllValues(0.0);

	double *y, *v, *f_x, *aux, *invVar, *n;

	y=_Y.getArray(); v=_V.getArray(); f_x=_F_X.getArray(); aux=AUX.getArray(); invVar=_ubm_invvar.getArray(); n=_matN.getArray();

	for(unsigned long spk=0; spk<_n_speakers; spk++){
		L.setAllValues(0.0);
		for(unsigned long i=0; i<_rankEV; i++){	L(i,i)=1.0;}

		double *l;
		l=L.getArray();

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
		
		//Inverse L
		DoubleSquareMatrix Linv(_rankEV);
		L.invert(Linv);
		double *invl = Linv.getArray();

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
struct estimateYTthread_data{

	double *N;
	unsigned long spkBottom;
	unsigned long spkUp;	
	unsigned long rankEV;
	unsigned long n_distrib;
	unsigned long svSize;
	unsigned long nt;
	RefVector <DoubleSquareMatrix>* vevT;

//	double *AUX;
	double *Y;
	double *V;
	double *F_X;
	double *ubm_invvar;

	unsigned long numThread;

};

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
//				Thread Routine
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void *estimateYTthread(void *threadarg) {
	struct estimateYTthread_data *my_data;
	my_data = (struct estimateYTthread_data *) threadarg;

	double *n = my_data->N;
	unsigned long spkBottom = my_data->spkBottom;	
	unsigned long spkUp = my_data->spkUp;
	unsigned long _rankEV=my_data->rankEV;
	unsigned long _n_distrib=my_data->n_distrib;
	unsigned long svSize =my_data->svSize;
	unsigned long nt = my_data->nt;
	unsigned long numThread = my_data->numThread;

	double *v = my_data->V;
	double *y = my_data->Y;
	double *f_x = my_data->F_X;
	double *ubm_invvar = my_data->ubm_invvar;

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

		//DoubleSquareMatrix &linv=(*(my_data->linv))[spk];		
		DoubleSquareMatrix Linv(_rankEV);
		//Inverse L and stock it
		L.invert(Linv);
		double *invl = Linv.getArray();

		DoubleVector tmpAux(svSize,svSize);
		tmpAux.setAllValues(0.0);
		for(unsigned long i=0;i<_rankEV;i++){
			for(unsigned long k=0;k<svSize;k++) {
//				aux[spk*_rankEV+i] += f_x[spk*svSize+k] * ubm_invvar[k] * v[i*svSize+k];
				tmpAux[i] += f_x[spk*svSize+k] * ubm_invvar[k] * v[i*svSize+k];
			}
		}

		//multiplication by invL
		for(unsigned long i=0; i<_rankEV;i++){
			for(unsigned long k=0; k<_rankEV; k++){
//				y[spk*_rankEV+i] += aux[spk*_rankEV+k] * invl[i*_rankEV+k];
				y[spk*_rankEV+i] += tmpAux[k] * invl[i*_rankEV+k];
			}
		}
	}
	pthread_exit((void*) 0);
	return (void*)0 ;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
void TVAcc::estimateYThreaded(unsigned long NUM_THREADS){
	
	if (verboseLevel >= 1) cout << "(AccumulateTVStat) Estimate Y from Statistics Threaded"<<endl;
	if (NUM_THREADS==0) throw Exception("Num threads can not be 0",__FILE__,__LINE__);

	int rc, status;
	if (NUM_THREADS > _n_speakers) NUM_THREADS=_n_speakers;

	struct estimateYTthread_data *thread_data_array = new estimateYTthread_data[NUM_THREADS];
	pthread_t *threads = new pthread_t[NUM_THREADS];

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	unsigned long offset=_n_speakers/NUM_THREADS;

	_Y.setAllValues(0.0);
	
	double *N =_matN.getArray(); 
	double *V = _V.getArray();
	double *Y = _Y.getArray();
	double *F_X = _F_X.getArray();
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
		thread_data_array[t].rankEV=_rankEV;
		thread_data_array[t].n_distrib=_n_distrib;
		thread_data_array[t].svSize=_svSize;
		thread_data_array[t].vevT=&(_vEvT);
		thread_data_array[t].nt=t;
		thread_data_array[t].numThread=NUM_THREADS;
		thread_data_array[t].V=V;
		thread_data_array[t].Y=Y;
		thread_data_array[t].F_X=F_X;
		thread_data_array[t].ubm_invvar=ubm_invvar;

		if (verboseLevel>1) cout<<"(AccumulateTVStat) Creating thread n["<< t<< "] for speakers["<<spkBottom<<"-->"<<spkUp-1<<"]"<<endl;
		rc = pthread_create(&threads[t], &attr, estimateYTthread, (void *)&thread_data_array[t]);
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


void TVAcc::saveYbyFile(Config &config){
	
	String svPath=config.getParam("vectorFilesPath");
	String yExtension = ".y";
	if(config.existsParam("yExtension"))	yExtension = config.getParam("yExtension");

	XLine * linep;
	unsigned long session = 0;
	while ((linep=_fileList.getLine()) != NULL){             	// linep gives the XLine with the Id of a given client and the list of files

		String *id=linep->getElement(); 
		String yFile=svPath+*id+yExtension;
		
		Matrix<double> sessionY;
		sessionY.setDimensions(1,_rankEV);
		for(unsigned long i=0;i<_rankEV;i++){
			sessionY(0,i) = _Y(session,i);
		}
		sessionY.save(yFile,config);
		session++;
	}
}














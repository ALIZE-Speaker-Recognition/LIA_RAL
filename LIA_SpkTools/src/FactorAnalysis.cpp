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

#if !defined(ALIZE_FactorAnalysis_cpp)
#define ALIZE_FactorAnalysis_cpp

#include<FactorAnalysis.h>
#include<iostream>
#include<fstream>
#include<cstdio>
#include<cassert>
#include<cmath>
#include "SuperVectors.h"
#include "RealVector.h"
#if defined(THREAD)
#include <pthread.h>
#endif

using namespace alize;
using namespace std;

FactorAnalysisStat::FactorAnalysisStat(String & featFilename,FeatureServer & fs,Config & config):_ms(config),_ss(config){ // constructor for a single file
	XList faNdx;	
	XLine featLine;
	featLine.addElement(featFilename);
        faNdx.addLine()=featLine;
        _init(faNdx,fs,config);
}

FactorAnalysisStat::FactorAnalysisStat(XList & ndx,FeatureServer & fs,Config & config):_ms(config),_ss(config){ // constructeur
	_init(ndx,fs,config);
}

void FactorAnalysisStat::_init(XList & ndx,FeatureServer & fs,Config & config){
	fileList=ndx;
	_nb_speakers=fileList.getLineCount();
	_nb_sent=fileList.getAllElements().getElementCount();
	
	_ndxTable=Translate(ndx);
	_topGauss=config.existsParam("topGauss");
	if (verbose) cout << "(FactorAnalysisStat) Estimate FA stats on "<<_nb_sent<< " example(s) with "<<_nb_speakers<<" speaker(s)"<<endl;
		
	_ms.loadMixtureGD(config.getParam("inputWorldFilename")); // garde
	_ms.loadMixtureGD(config.getParam("inputWorldFilename")); // modifie
		
	MixtureGD & UBM=_ms.getMixtureGD((unsigned long) 0);
		
	if (verbose) cout <<"(FactorAnalysisStat) Model parameter Dim:" <<UBM.getVectSize()<<"___nbC:" <<UBM.getDistribCount()<<endl;

	_vsize=UBM.getVectSize();
	_mixsize=UBM.getDistribCount();
	_supervsize=_vsize*_mixsize;
	_D.setSize(_supervsize);

	// Random Init for U or init with an existing one (used in traintarget for instance)
	_matU=Matrix<double>();
	if (config.existsParam("initChannelMatrix")) {
		_matU.load(config.getParam("initChannelMatrix"),config);
		if (_matU.cols()>_matU.rows()) _matU.transpose();
		if (verbose) cout << "(FactorAnalysisStat) Init for Channel Matrix: "<< config.getParam("initChannelMatrix")<< ", rank: ["<<_matU.cols() << "] sv size: [" << _matU.rows() <<"]"<<endl;
		_rang=_matU.cols();
		if (UBM.getDistribCount()*UBM.getVectSize()!=_matU.rows()) throw Exception("Supervector size does not match with UBM model",__FILE__,__LINE__);
	}
	else {
		_rang=config.getParam("channelMatrixRank").toLong();		
		_matU.setDimensions(_supervsize,_rang);
		srand48(_supervsize*_rang);
		_matU.randomInit();			
		if (verbose) cout << "(FactorAnalysisStat) Random Init for Channel Matrix: "<<", rank: ["<<_matU.cols() << "] sv size: [" << _matU.rows() <<"]"<<endl;			
	}		
		
	// FA stats
	_matY=Matrix<double>(_nb_speakers,_supervsize);
	_matX=Matrix<double>(_nb_sent,_rang);
	_matS_X=Matrix<double>(_nb_speakers,_supervsize);
	_matS_X_h=Matrix<double>(_nb_sent,_supervsize);	
	_matN_h=Matrix<double>(_nb_sent,_mixsize);	
	_matN=Matrix<double>(_nb_speakers,_mixsize);		
	
	_matY.setAllValues(0.0);
	_matX.setAllValues(0.0);	
	_matS_X.setAllValues(0.0);	
	_matS_X_h.setAllValues(0.0);	
	_matN_h.setAllValues(0.0);
	_matN.setAllValues(0.0);
			
	_super_mean.setSize(_supervsize);
	_super_invvar.setSize(_supervsize);
	_tau=config.getParam("regulationFactor").toDouble();
	if (verbose) cout << "(FactorAnalysisStat) Regulation Factor ["<<_tau<<"]"<<endl;
	// D construction
	for(unsigned long i=0;i<_mixsize;i++){
		DistribGD & aux=UBM.getDistrib(i);
		RealVector <double> & v=aux.getCovInvVect();
		RealVector <double> & u=aux.getMeanVect();
		for(unsigned long j=0;j<_vsize;j++){
			 _D[i*_vsize+j]=sqrt(1.0/(v[j]*_tau));
			 _super_mean[i*_vsize+j]=u[j];
			 _super_invvar[i*_vsize+j]=v[j];
		}
	}
		

	if (verboseLevel > 1) cout <<"(FactorAnalysisStat) Number of examples " << _nb_sent << endl;
	for(unsigned long i=0;i<_nb_sent;i++) {
		_l_h_inv.addObject(*(new DoubleSquareMatrix()),i);
		_l_h_inv[i].setSize(_rang);
	}		
	if (verboseLevel >1) cout << "(FactorAnalysisStat) Y dimensions are: " << _matY.rows() <<","<<_matY.cols()<<"__ X dimensions are: "<< _matX.rows() <<","<< _matX.cols() <<endl << "(FactorAnalysisStat) FA Accumator Built"<<endl;     
};  

void FactorAnalysisStat::computeAndAccumulateGeneralFAStats(FeatureServer &fs,Config & config){
	SegServer segmentsServer;
	LabelServer labelServer;
	initializeClusters(fileList,segmentsServer,labelServer,config);
	verifyClusterFile(segmentsServer,fs,config);
	unsigned long codeSelectedFrame=labelServer.getLabelIndexByString(config.getParam("labelSelectedFrames"));
	SegCluster& selectedSegments=segmentsServer.getCluster(codeSelectedFrame); 
	this->computeAndAccumulateGeneralFAStats(selectedSegments,fs,config);
};

void FactorAnalysisStat::computeAndAccumulateGeneralFAStats(SegCluster &selectedSegments,FeatureServer &fs,Config & config){
	if (verbose) cout <<"(FactorAnalysisStat) Compute General FA Stats (Complete)" << endl;
	double *N_h, *N, *S_X_h, *S_X,*ff;	
	_matN_h.setAllValues(0.0);
	_matN.setAllValues(0.0);
	_matS_X_h.setAllValues(0.0);
	_matS_X.setAllValues(0.0);
	N_h=_matN_h.getArray(); N=_matN.getArray(); S_X_h=_matS_X_h.getArray();S_X=_matS_X.getArray();
	
	MixtureGD & UBM=_ms.getMixtureGD((unsigned long) 1);
	MixtureGDStat &acc=_ss.createAndStoreMixtureStat(UBM);

	// Compute Occupations and Statistics
	acc.resetOcc();
	Seg *seg; 
	selectedSegments.rewind();
	String currentSource="";unsigned long loc=0;unsigned long sent=0;
	while((seg=selectedSegments.getSeg())!=NULL){	
		unsigned long begin=seg->begin()+fs.getFirstFeatureIndexOfASource(seg->sourceName()); 				// Idx of the first frame of the current file in the feature server
		if (currentSource!=seg->sourceName()) {
		currentSource=seg->sourceName();
		loc=_ndxTable.locNb(currentSource);
		sent=_ndxTable.sessionNb(currentSource);	
		if (verbose)cout << "Processing speaker["<<currentSource<<"]"<< endl;	
		}

		fs.seekFeature(begin);
		Feature f;
		if (!_topGauss) {
			for (unsigned long idxFrame=0;idxFrame<seg->length();idxFrame++){
				fs.readFeature(f);
				acc.computeAndAccumulateOcc(f);
				RealVector <double> aPost=acc.getOccVect();
				ff=f.getDataVector();
				for(unsigned long k=0;k<_mixsize;k++) {
					N_h[sent*_mixsize+k]+=aPost[k];
					N[loc*_mixsize+k]   +=aPost[k];
					for (unsigned long i=0;i<_vsize;i++) {
						S_X_h[sent*_supervsize+(k*_vsize+i)]+=aPost[k]*ff[i];
						S_X[loc*_supervsize+(k*_vsize+i)]   +=aPost[k]*ff[i];
						}
				}	
			}
		} 
		else throw Exception("ComputeGeneralStats TopGauss not done at this level",__FILE__,__LINE__);
	}					
};
/*
void FactorAnalysisStat::computeAndAccumulateGeneralFAStatsTopGauss(FeatureServer &fs,Config & config){ // this does not work great
	if (verbose) cout <<"(FactorAnalysisStat) Compute General FA Stats (TopGauss)" << endl;
	unsigned long i,k,begin, nbg, k1;
	RealVector <double> AUX2; AUX2.setSize(_supervsize); AUX2.setAllValues(0.0);double* aux2=AUX2.getArray();	
	RealVector <double> _Prob; _Prob.setSize(_mixsize); double *Prob=_Prob.getArray();
	RealVector <double> m_xh; m_xh.setSize(_supervsize);
	double *super_mean, *N_h, *N, *S_X_h, *S_X,*ff;	
	super_mean=_super_mean.getArray(); N_h=_matN_h.getArray(); N=_matN.getArray(); S_X_h=_matS_X_h.getArray(); S_X=_matS_X.getArray();
	static unsigned long it=0;
	MixtureGD & UBM=_ms.getMixtureGD((unsigned long) 1);
	
	// Compute Occupations and Statistics
	unsigned long loc=0; 
	unsigned long sent=0; 	
	XLine *pline; String *pFile; fileList.rewind();double sum=0.0;

	while((pline=fileList.getLine())!=NULL) { 
		if (verbose) cout << "(FactorAnalysisStat) Process speaker [" << loc << "] ----> " <<endl;
		for(i=0;i<_mixsize;i++) 
			N[loc*_mixsize+i]=0.0;
		for(i=0;i<_supervsize;i++) 
			S_X[loc*_supervsize+i]=0.0;
		while((pFile=pline->getElement())!=NULL) {
			if (1) this->getSpeakerModel(UBM,loc,sent);

			if (verboseLevel > 1) cout << "Session : [" << sent <<"], ";
			for(i=0;i<_mixsize;i++) 
				N_h[sent*_mixsize+i]=0;//.000000001;
			for(i=0;i<_mixsize;i++) 
				N[loc*_mixsize+i]+=0;//.000000001;
			
			for(i=0;i<_supervsize;i++) 
				S_X_h[sent*_supervsize+i]=0.0;
			
			this->getUX(AUX2,sent);
				
			for(k=0;k<_supervsize;k++) 
				aux2[k]+=super_mean[k];

			String featureFileName=(*pFile); 			// Current file basename	
			begin=fs.getFirstFeatureIndexOfASource(featureFileName);
			fs.seekFeature(begin);
			SegServer segmentsServer;
			LabelServer labelServer;
			initializeClusters(featureFileName,segmentsServer,labelServer,config);
			verifyClusterFile(segmentsServer,fs,config);
			unsigned long codeSelectedFrame=labelServer.getLabelIndexByString(config.getParam("labelSelectedFrames"));
			SegCluster& selectedSegments=segmentsServer.getCluster(codeSelectedFrame);  
			Seg *seg; 
			selectedSegments.rewind();
			TopGauss tg;
			tg.read(featureFileName,config);
			unsigned long*_idx=tg.idx().getArray();
			unsigned long*_nbg=tg.nbg().getArray();

			unsigned long nt=0; //frame counter
			unsigned long nan=0;			unsigned long floor=0;
			while((seg=selectedSegments.getSeg())!=NULL){
				unsigned long begin=seg->begin()+fs.getFirstFeatureIndexOfASource(seg->sourceName()); 				// Idx of the first frame of the current file in the feature server      
				fs.seekFeature(begin);
				Feature f;
				unsigned long idxBegin=tg.frameToIdx(nt); // get the begin index of gaussians in idx vectors			
				for (unsigned long idxFrame=0;idxFrame<seg->length();idxFrame++){
					fs.readFeature(f);
					ff=f.getDataVector();
					nbg=_nbg[nt]; // get nb gauss for this frame
					sum=0.0;
					for(unsigned long k=0;k<nbg;k++) {
						Prob[k]=UBM.weight(_idx[idxBegin+k])*UBM.getDistrib(_idx[idxBegin+k]).computeLK(f);
						if (isnan(Prob[k])) {Prob[k]=1e-20;nan++;}
						sum+=Prob[k];												
					}
					if (log(sum)<-200) {floor++;_Prob.setAllValues(0.0);}
					else {
						for(unsigned long k=0;k<nbg;k++) 
							Prob[k]/=sum; 
					}
					for(unsigned long k=0;k<nbg;k++) {			
						k1=_idx[k+idxBegin];
						N_h[sent*_mixsize+k1]+=Prob[k];
						N[loc*_mixsize+k1]   +=Prob[k];
						for (unsigned long i=0;i<_vsize;i++) {
							S_X_h[sent*_supervsize+(k1*_vsize+i)]+=Prob[k]*ff[i];
							S_X[loc*_supervsize+(k1*_vsize+i)]   +=Prob[k]*ff[i];
							}
					}
					idxBegin+=nbg;
					nt++;	
				}	
			}
			if (verboseLevel > 2) cout << "nan:" <<nan<<" floor:" <<floor<<endl;
			sent++;
		}
		loc++;
	}
	it++;
	if (verboseLevel > 1) cout << "(FactorAnalysisStat) Done " << endl;
};*/



void FactorAnalysisStat::substractSpeakerStats(){
	if (verbose) cout <<"(FactorAnalysisStat) Compute and Substract Speaker FA Stats... " << endl;	
	RealVector <double> AUX1;AUX1.setSize(_supervsize); AUX1.setAllValues(0.0); double *aux1=AUX1.getArray();  
	
	// Compute Occupations and Statistics
	unsigned long loc=0; 
	unsigned long sent=0; 		
	XLine *pline; String *pFile; fileList.rewind();

	double *N_h=_matN_h.getArray(); 
	double *S_X_h=_matS_X_h.getArray();

	while((pline=fileList.getLine())!=NULL) { 
		while((pFile=pline->getElement())!=NULL) {
		this->getMplusDY(AUX1,*pFile);			
			for(unsigned long k=0;k<_mixsize;k++) 
				for (unsigned long i=0;i<_vsize;i++) 
					S_X_h[sent*_supervsize+(k*_vsize+i)]-= N_h[sent*_mixsize+k]*aux1[i+k*_vsize];
			sent++;
		}
		loc++;
	}
	if (verboseLevel >1) cout << "(FactorAnalysisStat) Done "<<endl;
};

void FactorAnalysisStat::substractChannelStats(){
	if (verbose) cout <<"(FactorAnalysisStat) Compute and Substract Channel FA Stats... "<<endl;	
	RealVector <double> UX;UX.setSize(_supervsize); UX.setAllValues(0.0);double* ux=UX.getArray();		
		
	// Compute Occupations and Statistics
	unsigned long loc=0;
	unsigned long sent=0; 
	XLine *pline; String *pFile; 
	fileList.rewind();
	
	double *super_mean=_super_mean.getArray();
	double *N_h=_matN_h.getArray(); 
	double *S_X=_matS_X.getArray();

	while((pline=fileList.getLine())!=NULL) { 		
		while((pFile=pline->getElement())!=NULL) {
			this->getUX(UX,*pFile);
			for(unsigned long k=0;k<_supervsize;k++) 
				ux[k]+=super_mean[k];
			
			for(unsigned long k=0;k<_mixsize;k++)
				for (unsigned long i=0;i<_vsize;i++)
					S_X[loc*_supervsize+(k*_vsize+i)]-= N_h[sent*_mixsize+k]*ux[i+k*_vsize];
			sent++;
		}
		loc++;
	}
};


void FactorAnalysisStat::getXEstimate(){
	if (verbose) cout << "(FactorAnalysisStat) Compute X Estimate "<<endl;
	RealVector <double> AUX;
	AUX.setSize(_rang);
	XLine *pline;
	String *pFile;

	_matX.setAllValues(0.0);

	double *X=_matX.getArray();	
	double *U=_matU.getArray();
	double *S_X_h=_matS_X_h.getArray();
	double *aux=AUX.getArray();
	double *super_invvar=_super_invvar.getArray();
	
	fileList.rewind();
	while((pline=fileList.getLine())!=NULL) {
		while((pFile=pline->getElement())!=NULL) {
		    unsigned long sent=_ndxTable.sessionNb(*pFile);
			AUX.setAllValues(0.0);
			for(unsigned long i=0;i<_rang;i++)
				for(unsigned long k=0;k<_supervsize;k++) 
					aux[i]+= U[k*_rang+i]*super_invvar[k]*S_X_h[sent*_supervsize+k];	
			double *l_h_inv=_l_h_inv[sent].getArray();	
			for(unsigned long i=0;i<_rang;i++)
				for(unsigned long k=0;k<_rang;k++) 
					X[sent*_rang+i]+=l_h_inv[i*_rang+k]*aux[k];
			sent++;
		}
	}
};

void FactorAnalysisStat::estimateAndInverseL(Config & config){
    #if defined(THREAD)          
    if (config.existsParam("numThread") && config.getParam("numThread").toLong() >0) estimateAndInverseLThreaded(config.getParam("numThread").toLong());
    else estimateAndInverseLUnThreaded();
    #else
    estimateAndInverseLUnThreaded();
    #endif
};

void FactorAnalysisStat::estimateAndInverseLUnThreaded(){
	if (verbose) cout << "(FactorAnalysisStat) Inverse L Matrix ... "<<endl;	
	unsigned long mk;
	DoubleSquareMatrix L(_rang);	
	L.setAllValues(0.0);
	RealVector <double> AUX;
	AUX.setSize(_rang);
	unsigned long sent=0;
	XLine *pline;
	fileList.rewind();

	double *N_h=_matN_h.getArray(); 
	double *U=_matU.getArray();
	double *LV=L.getArray();
	double *super_invvar=_super_invvar.getArray();
	
	while((pline=fileList.getLine())!=NULL) {
		while(pline->getElement()!=NULL) {
			L.setAllValues(0.0);
			AUX.setAllValues(0.0);
			for(unsigned long i=0;i<_rang;i++){
				for(unsigned long j=0;j<=i;j++){
					for(unsigned long k=0;k<_supervsize;k++){
						mk=k/_vsize;
						LV[i*_rang+j]+=N_h[sent*_mixsize+mk]*super_invvar[k]*U[k*_rang+i]*U[k*_rang+j];
						}
					}
				}
			for(unsigned long i=0;i<_rang;i++)
				for(unsigned long j=i+1;j<_rang;j++) 
					LV[i*_rang+j]=LV[j*_rang+i];
				
			for(unsigned long i=0;i<_rang;i++) 
				LV[i*_rang+i]+=1.0;
			L.invert(_l_h_inv[sent]);	
			sent++;
		}
	}
};

/*****************************************************************************************************
********************* Threaed Version of L Matrices Inversion (Nico Scheffer)
*****************************************************************************************************/
#if defined(THREAD)
struct Lthread_data{
	double *N_h;
	double *U;
	double *super_invvar;
	unsigned long sentBottom;
	unsigned long sentUp;	
	unsigned long rang;
	unsigned long sv;
	unsigned long vsize;
	unsigned long mixsize;
	RefVector <DoubleSquareMatrix>* linv;
};

void *Lthread(void *threadarg) {
	struct Lthread_data *my_data;
	my_data = (struct Lthread_data *) threadarg;
	double *U = my_data->U;
	unsigned long sentBottom = my_data->sentBottom;	
	unsigned long sentUp = my_data->sentUp;
	double *N_h=my_data->N_h;
	double *super_invvar=my_data->super_invvar;
	unsigned long _rang=my_data->rang;
	unsigned long _supervsize=my_data->sv;
	unsigned long _vsize=my_data->vsize;
	unsigned long _mixsize=my_data->mixsize;
	for (unsigned long sent=sentBottom;sent <sentUp;sent++) {
		DoubleSquareMatrix L(_rang);	
		L.setAllValues(0.0);	
		double *LV=L.getArray();
		unsigned long mk;
		RealVector <double> AUX;
		AUX.setSize(_rang);
		AUX.setAllValues(0.0);
		for(unsigned long i=0;i<_rang;i++){
			for(unsigned long j=0;j<=i;j++){
				for(unsigned long k=0;k<_supervsize;k++){
					mk=k/_vsize;
					LV[i*_rang+j]+=N_h[sent*_mixsize+mk]*super_invvar[k]*U[k*_rang+i]*U[k*_rang+j];
				}
			}
		}
		for(unsigned long i=0;i<_rang;i++)
			for(unsigned long j=i+1;j<_rang;j++) 
				LV[i*_rang+j]=LV[j*_rang+i];
					
		for(unsigned long i=0;i<_rang;i++) 
			LV[i*_rang+i]+=1.0;	
		DoubleSquareMatrix &linv=(*(my_data->linv))[sent];		
		L.invert(linv);
	}
	pthread_exit((void*) 0);
	return (void*)0;
}

void FactorAnalysisStat::estimateAndInverseLThreaded(unsigned long NUM_THREADS){
	if (verbose) cout << "(FactorAnalysisStat) Inverse L Matrix Threads ... "<<endl;	
	if (NUM_THREADS==0) throw Exception("Num threads can be 0",__FILE__,__LINE__);
	double *N_h=_matN_h.getArray(); 
	double *U=_matU.getArray();
	double *super_invvar=_super_invvar.getArray();

	int rc, status;
	if (NUM_THREADS > _nb_sent) NUM_THREADS=_nb_sent;
//	struct Lthread_data thread_data_array[NUM_THREADS];
	struct Lthread_data *thread_data_array = new Lthread_data[NUM_THREADS];
//	pthread_t threads[NUM_THREADS];	
	pthread_t *threads = new pthread_t[NUM_THREADS];

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	unsigned long offset=_nb_sent/NUM_THREADS;
	for(unsigned long t=0, sentBottom=0; t<NUM_THREADS; t++,sentBottom+=offset){
		thread_data_array[t].N_h=N_h;
		thread_data_array[t].U=U;
		thread_data_array[t].super_invvar=super_invvar;
		thread_data_array[t].sentBottom=sentBottom;
		unsigned long sentUp=sentBottom+offset;
		if (_nb_sent%NUM_THREADS!=0) sentUp++;
		thread_data_array[t].sentUp=sentUp;
		thread_data_array[t].rang=_rang;
		thread_data_array[t].sv=_supervsize;
		thread_data_array[t].vsize=_vsize;
		thread_data_array[t].mixsize=_mixsize;
		thread_data_array[t].linv=&(_l_h_inv);

		if (verbose) cout<<"(FactorAnalysisStat) Creating thread n["<< t<< "] for sessions["<<sentBottom<<"-->"<<sentUp<<"]"<<endl;		
		rc = pthread_create(&threads[t], &attr, Lthread, (void *)&thread_data_array[t]);		
		if (rc) throw Exception("ERROR; return code from pthread_create() is ",__FILE__,rc);
	}
	pthread_attr_destroy(&attr);
	for(unsigned long t=0; t<NUM_THREADS; t++) {
		rc = pthread_join(threads[t], (void **)&status);
		if (rc)  throw Exception("ERROR; return code from pthread_join() is ",__FILE__,rc);
		if (verbose) cout <<"(FactorAnalysisStat) Completed join with thread ["<<t<<"] status["<<status<<"]"<<endl;
	}

	free(thread_data_array);
	free(threads);

	if (verboseLevel > 1) cout << "(FactorAnalysisStat) Done " << endl;
};
#endif
/****************************************** End *****************************************************/

void FactorAnalysisStat::getYEstimate(){
	if (verbose) cout << "(FactorAnalysisStat) Compute Y Estimate "<<endl;
	_matY.setAllValues(0.0);

	double *N=_matN.getArray(); 
	double *Y=_matY.getArray();
	double *S_X=_matS_X.getArray();
	double *D=_D.getArray();
	double *super_invvar=_super_invvar.getArray();

	for(unsigned long loc=0;loc<_nb_speakers;loc++) {
		for(unsigned long i=0;i<_supervsize;i++){
			unsigned long mi=i/_vsize;
			Y[loc*_supervsize+i]= (_tau/(_tau+N[loc*_mixsize+mi]))*D[i]*super_invvar[i]*S_X[loc*_supervsize+i];
			}
		}
};

void FactorAnalysisStat::getUEstimate(Config & config){
    #if defined(THREAD)          
    if (config.existsParam("numThread") && config.getParam("numThread").toLong() >0) getUEstimateThreaded(config.getParam("numThread").toLong());
    else getUEstimateUnThreaded();
    #else
    getUEstimateUnThreaded();
    #endif
};

void FactorAnalysisStat::getUEstimateUnThreaded(){
	if (verbose) cout << "(FactorAnalysisStat) Compute U Estimate "<<endl;	
	DoubleSquareMatrix L(_rang);
	DoubleSquareMatrix L_Inv(_rang);
	RealVector <double> R;
	R.setSize(_rang);

	double *N_h=_matN_h.getArray(); 
	double *X=_matX.getArray();
	double *S_X_h=_matS_X_h.getArray();
	double *U=_matU.getArray();
	double *RV=R.getArray();
	double *LV=L.getArray();

	/*for(i=0;i<_supervsize;i++){
		unsigned long im=i/_vsize;
		L.setAllValues(0.0);
		R.setAllValues(0.0);
		for(j=0;j<_nb_sent;j++){
			double *l_h_inv=_l_h_inv[j].getArray();
			for(k=0;k<_rang;k++){
				for(l=0;l<_rang;l++) 
					LV[k*_rang+l]+=(l_h_inv[k*_rang+l]+X[j*_rang+k]*X[j*_rang+l])*N_h[j*_mixsize+im];
				RV[k]+=S_X_h[j*_supervsize+i]*X[j*_rang+k];
				}
			}		
		L.invert(L_Inv);
		double *L_InvV=L_Inv.getArray();
		for(j=0;j<_rang;j++)
			U[i*_rang+j]=0.0;
		for(j=0;j<_rang;j++)
			for(k=0;k<_rang;k++) 
				U[i*_rang+j]+=L_InvV[j*_rang+k]*RV[k];
	}*/
	_matU.setAllValues(0.0);
	for(unsigned long g=0;g<_mixsize;g++){		
		L.setAllValues(0.0);
		// Estimate LU for the gaussian g
		for(unsigned long j=0;j<_nb_sent;j++){
			double *l_h_inv=_l_h_inv[j].getArray();
			for(unsigned long k=0;k<_rang;k++)
				for(unsigned long l=0;l<_rang;l++) 
					LV[k*_rang+l]+=(l_h_inv[k*_rang+l]+X[j*_rang+k]*X[j*_rang+l])*N_h[j*_mixsize+g];
		}
		L.invert(L_Inv); 
		double *L_InvV=L_Inv.getArray();	
		
		// Estimate LU for the line lineNum
		for (unsigned long i=0;i<_vsize;i++) {
			unsigned long lineNum=_vsize*g+i;			
			// estime RU i
			R.setAllValues(0.0);
			for(unsigned long j=0;j<_nb_sent;j++)
				for(unsigned long k=0;k<_rang;k++)
					RV[k]+=S_X_h[j*_supervsize+lineNum]*X[j*_rang+k];
			// estime U i
			for(unsigned long j=0;j<_rang;j++)
				for(unsigned long k=0;k<_rang;k++) 
					U[lineNum*_rang+j]+=L_InvV[j*_rang+k]*RV[k];
		}
	}	
};

/*****************************************************************************************************
********************* Threaded Version of U Matrix Computation (Nico Scheffer)
*****************************************************************************************************/
#if defined(THREAD)
struct Uthread_data{
	double *N_h;
	double *U;
	double *S_X_h;
	double *X;
	unsigned long rang;
	unsigned long sv;
	unsigned long vsize;
	unsigned long mixsize;		
	unsigned long nbsent;	
	unsigned long idxBottom;	
	unsigned long idxUp;
	RefVector <DoubleSquareMatrix>* L;
};

void *Uthread(void *threadarg) {
	struct Uthread_data *my_data;
	my_data = (struct Uthread_data *) threadarg;
	// Get data
	double *U = my_data->U;
	double *N_h=my_data->N_h;	
	double *S_X_h=my_data->S_X_h;	
	double *X=my_data->X;
	unsigned long _rang=my_data->rang;
	unsigned long _supervsize=my_data->sv;	
	unsigned long _nb_sent=my_data->nbsent;
	unsigned long _vsize=my_data->vsize;
	unsigned long _mixsize=my_data->mixsize;
	unsigned long idxBottom=my_data->idxBottom;
	unsigned long idxUp=my_data->idxUp;
	RefVector <DoubleSquareMatrix>& _l_h_inv=*(my_data->L);

	// Routine
	DoubleSquareMatrix L(_rang);
	DoubleSquareMatrix L_Inv(_rang);
	RealVector <double> R;
	R.setSize(_rang);
	double *RV=R.getArray();
	double *LV=L.getArray();	
	for (unsigned long g=idxBottom;g<idxUp;g++) {
		L.setAllValues(0.0);
		// Estimate LU for the gaussian g
		for(unsigned long j=0;j<_nb_sent;j++){
			double *l_h_inv=_l_h_inv[j].getArray();
			for(unsigned long k=0;k<_rang;k++)
				for(unsigned long l=0;l<_rang;l++) 
					LV[k*_rang+l]+=(l_h_inv[k*_rang+l]+X[j*_rang+k]*X[j*_rang+l])*N_h[j*_mixsize+g];
		}
		L.invert(L_Inv); 
		double *L_InvV=L_Inv.getArray();	
		
		// Estimate LU for the line lineNum
		for (unsigned long i=0;i<_vsize;i++) {
			unsigned long lineNum=_vsize*g+i;			
			// estime RU i
			R.setAllValues(0.0);
			for(unsigned long j=0;j<_nb_sent;j++)
				for(unsigned long k=0;k<_rang;k++)
					RV[k]+=S_X_h[j*_supervsize+lineNum]*X[j*_rang+k];
			// estime U i
			for(unsigned long j=0;j<_rang;j++)
				for(unsigned long k=0;k<_rang;k++) 
					U[lineNum*_rang+j]+=L_InvV[j*_rang+k]*RV[k];
		}
	}
	//
	pthread_exit((void*) 0);
	return(void*)0;
}


void FactorAnalysisStat::getUEstimateThreaded(unsigned long NUM_THREADS){
	if (verbose) cout << "(FactorAnalysisStat) Compute U Estimate Threaded "<<endl;	
	double *N_h=_matN_h.getArray(); 
	double *X=_matX.getArray();
	double *S_X_h=_matS_X_h.getArray();
	double *U=_matU.getArray();
	_matU.setAllValues(0.0);
	if (NUM_THREADS==0) throw Exception("Num threads can not be 0",__FILE__,__LINE__);
	int rc, status;
	if (NUM_THREADS > _mixsize) NUM_THREADS=_mixsize;
//	struct Uthread_data thread_data_array[NUM_THREADS];
	struct Uthread_data *thread_data_array = new Uthread_data[NUM_THREADS];

	unsigned long offset=_mixsize/NUM_THREADS;
//	pthread_t threads[NUM_THREADS];	
	pthread_t *threads = new pthread_t[NUM_THREADS];

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);	
	for(unsigned long t=0, idxBottom=0; t<NUM_THREADS; t++,idxBottom+=offset){
		thread_data_array[t].N_h=N_h;
		thread_data_array[t].U=U;
		thread_data_array[t].S_X_h=S_X_h;
		thread_data_array[t].X=X;			
		thread_data_array[t].rang=_rang;
		thread_data_array[t].sv=_supervsize;
		thread_data_array[t].vsize=_vsize;
		thread_data_array[t].mixsize=_mixsize;		
		thread_data_array[t].nbsent=_nb_sent;
		thread_data_array[t].idxBottom=idxBottom;
		unsigned long idxUp=idxBottom+offset;
		if (_supervsize%NUM_THREADS!=0) idxUp++;		
		thread_data_array[t].idxUp=idxUp;		
		thread_data_array[t].L=&(_l_h_inv);
		if (verbose) cout<<"(FactorAnalysisStat) Creating thread n["<< t<< "] for U, Gaussian ["<<idxBottom<<"-->"<<idxUp<<"]"<<endl;			
		rc = pthread_create(&threads[t], &attr, Uthread, (void *)&thread_data_array[t]);
		if (rc) throw Exception("ERROR; return code from pthread_create() is ",__FILE__,rc);
	}
	pthread_attr_destroy(&attr);
	for(unsigned long t=0; t<NUM_THREADS; t++) {
		rc = pthread_join(threads[t], (void **)&status);
		if (rc)  throw Exception("ERROR; return code from pthread_join() is ",__FILE__,rc);
		if (verboseLevel >1) cout <<"(FactorAnalysisStat)Completed join with thread ["<<t<<"] status["<<status<<"]"<<endl;
	}

	free(thread_data_array);
	free(threads);

	if (verboseLevel >1) cout << "(FactorAnalysisStat) Done " << endl;	
};
#endif
/****************************************** End *********************************************/

/*void FactorAnalysisStat::getLLK(Config & config){
	if (verbose) cout << "(FactorAnalysisStat) # Compute Likelihoods" << endl;		
	unsigned long loc,sent;
	unsigned long cnt=0;
	unsigned long maxLLKcomputed=1;
	if (config.existsParam("computeLLK")) maxLLKcomputed=config.getParam("computeLLK").toLong();	
	MixtureGD & UBM=_ms.getMixtureGD((unsigned long) 1);
	MixtureGDStat &acc=_ss.createAndStoreMixtureStat(UBM);		
	DoubleVector m_xh;
	m_xh.setSize(_supervsize);	
	
	// Compute Occupations and Statistics
	loc=0; sent=0; XLine *pline; String *pFile; fileList.rewind();
	double glob_llk=0.0, llk_sent;

	while((pline=fileList.getLine())!=NULL && cnt < maxLLKcomputed) { 
		while((pFile=pline->getElement())!=NULL && cnt < maxLLKcomputed) {
			cnt++;
			String featureFileName=(*pFile); 			// Current file basename				
			FeatureServer fs(config,*pFile);
			this->getSpeakerModel(UBM,*pFile);
			
			if (_topGauss) {
				TopGauss tg;
				tg.read(*pFile,config);
				llk_sent=tg.get(UBM,*pFile,config);//This does not work on top gauss
			}
			else {
				unsigned long begin=fs.getFirstFeatureIndexOfASource(featureFileName);
				fs.seekFeature(begin);
				
				// Usual SegServer stuff
				SegServer segmentsServer;
				LabelServer labelServer;
				initializeClusters(featureFileName,segmentsServer,labelServer,config);
				verifyClusterFile(segmentsServer,fs,config);
				unsigned long codeSelectedFrame=labelServer.getLabelIndexByString(config.getParam("labelSelectedFrames"));
				SegCluster& selectedSegments=segmentsServer.getCluster(codeSelectedFrame);  
				acc.resetOcc();
				Seg *seg;          // current selected segment
				selectedSegments.rewind(); 			// reset the reader at the begin of the input stream

				acc.resetLLK();
				while((seg=selectedSegments.getSeg())!=NULL){                           	
					unsigned long begin=seg->begin()+fs.getFirstFeatureIndexOfASource(seg->sourceName()); 
					fs.seekFeature(begin);
					Feature f;
					for (unsigned long idxFrame=0;idxFrame<seg->length();idxFrame++){
						fs.readFeature(f); 
						acc.computeAndAccumulateLLK(f,1.0,DETERMINE_TOP_DISTRIBS);
					}		
				}
				llk_sent= acc.getMeanLLK();
			}
			if (verbose) cout << "(FactorAnalysisStat) # Mean LLK on [" << " " << featureFileName << "](" << sent << ")=" << llk_sent << endl ;
			glob_llk +=llk_sent;
			sent++;
		}
		loc++;
	}
	_ss.deleteMixtureStat(acc);	
	if (verbose) cout << "(FactorAnalysisStat) # Total Likelihood: " << glob_llk << endl;
};*/

/// Compute Factor Analysis model to perform likelihood computation in the classical way
void FactorAnalysisStat::getFactorAnalysisModel(MixtureGD& FA,String& file) {
	if (verbose) cout << "(FactorAnalysisStat) Compute Variance adapted Speaker Model"<<endl;		
	this->getTrueSpeakerModel(FA,file);
	unsigned long loc=_ndxTable.locNb(file);	

	/// Compute sigma_s
	RealVector <double> sigma_s;
	sigma_s.setSize(_supervsize);
	for (unsigned long i=0;i<_mixsize;i++) 	
		for (unsigned long j=0;j<_vsize;j++) 
			sigma_s[i*_vsize+j]=FA.getDistrib(i).getCov(j)/(_tau+_matN(loc,i));
		
	/// Compute sigma_s+*sigma_w
	RealVector <double> sum,prod;
	sum.setSize(_supervsize);prod.setSize(_supervsize);
	for (unsigned long i=0;i<_mixsize;i++) {
		for (unsigned long j=0;j<_vsize;j++) {
			sum[i*_vsize+j]=sigma_s[i*_vsize+j]+FA.getDistrib(i).getCov(j);
			//cout << "i,j"<<i<<","<<j<<" sigma_s "<<sigma_s[i*_vsize+j]<<" cov:"<<FA.getDistrib(i).getCov(j)<<" N"<<_N(loc,i)<<" "<<_tau<<endl;
		}
	}

	for (unsigned long i=0;i<_mixsize;i++)
		for (unsigned long j=0;j<_vsize;j++)
			FA.getDistrib(i).setCov(sum[i*_vsize+j],j);
	FA.computeAll();
}

/// Compute Log Likelihood of the Factor Analysis model
double FactorAnalysisStat::getLLK(SegCluster &selectedSegments,MixtureGD &model,FeatureServer&fs,Config & config){
	if (verbose) cout << "(FactorAnalysisStat) Compute Likelihood" << endl;		
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

// Set Y given a mixture model (substract world and divide by D)
void FactorAnalysisStat::setYFromModel(MixtureGD& M,String & file){
	if (verbose) cout << "(FactorAnalysisStat) Set Y from Model" << endl;
	unsigned long loc=_ndxTable.locNb(file);	
	RealVector<double> v(_supervsize,_supervsize);
	modelToSv(M,v);
	RealVector<double> w(_supervsize,_supervsize);
	modelToSv(_ms.getMixtureGD(0),w); // get UBM
	for (unsigned long i=0;i<v.size();i++) 
		_matY(loc,i)=(v[i]-w[i])/_D[i];
};

void FactorAnalysisStat::estimateXForKnownSpeaker(MixtureGD & M,String & file,Config &config) {
	RealVector <double> v(_supervsize,_supervsize);
	modelToSv(M,v);	
	this->resetXY();
	this->setYFromModel(M,file);      
        for(int i=0;i<config.getParam("nbTrainIt").toLong();i++){
		if (verbose) cout << "(FactorAnalysisStat) ------ Iteration ["<<i<<"] ------"<<endl;
		this->substractSpeakerStats();
		this->getXEstimate();
	}
	double * X=_matX.getArray();
	for (unsigned long i=0;i<_supervsize;i++) {
		double ux=0.0;
		for (unsigned long j=0;j<_rang;j++) {
			ux+=_matU(i,j)*X[j];
		}
		v[i]+=ux;
	}
	svToModel(v,M);
};

void FactorAnalysisStat::getUX(RealVector <double> &ux,String& file) {
	ux.setAllValues(0.0);	
	unsigned long idx=_ndxTable.sessionNb(file);
	for (unsigned long i=0;i<_supervsize;i++)
		for (unsigned long j=0;j<_rang;j++) 
			ux[i]+=_matU(i,j)*_matX(idx,j);
};

// Compute supervector of client M_s_h=M+Dy_s and get a model
void FactorAnalysisStat::getTrueSpeakerModel(MixtureGD & M,String &file) {
	if (verbose) cout << "(FactorAnalysisStat) Compute True Speaker Model"<<endl;		
	RealVector <double> Sp(_supervsize,_supervsize);
	getMplusDY(Sp,file);
	svToModel(Sp,M);
};

// Compute supervector of client M_s_h=M+Dy_s and get a model
void FactorAnalysisStat::getSessionModel(MixtureGD & M,String &file) {
	if (verbose) cout << "(FactorAnalysisStat) Compute Channel Model"<<endl;		
	RealVector <double> Sp(_supervsize,_supervsize);
	getUX(Sp,file);
	for (unsigned long i=0;i<_supervsize;i++)
		Sp[i]+=_super_mean[i];
	svToModel(Sp,M);
};

// Compute supervector of client M_s_h=M+Dy_s
void FactorAnalysisStat::getMplusDY(RealVector <double> &Sp, String& file) {
	Sp.setAllValues(0.0);
	unsigned long loc=_ndxTable.locNb(file);		
	for (unsigned long i=0;i<_supervsize;i++)
		Sp[i]=_super_mean[i]+_D[i]*_matY(loc,i);
};

// Compute supervector of client M_s_h=M+Dy_s+Ux and get a model
void FactorAnalysisStat::getSpeakerModel(MixtureGD & M,String& file) {	
	RealVector <double> Sp;Sp.setSize(_supervsize);	
	RealVector <double> ux;ux.setSize(_supervsize);				
	this->getUX(ux,file);
	this->getMplusDY(Sp,file);
	for (unsigned long i=0;i<_supervsize;i++)
		Sp[i]+=ux[i];
	svToModel(Sp,M);
};

/// Normalize features with a smooth mixture transformation o't=ot-sum(P(c|ot)Uc.x)
void FactorAnalysisStat::normalizeFeatures(SegCluster &selectedSegments,FeatureServer &fs,Config & config){
	if (verbose) cout << "(FactorAnalysisStat) Normalize Features" << endl;	
	MixtureGD & clientMixture=_ms.getMixtureGD(1); // copy the UBM mixture		
	unsigned long nt=0;	
	RealVector <double> m_xh_1; m_xh_1.setSize(_supervsize); 	
	double *_m_xh_1=m_xh_1.getArray();
	Seg *seg;          // current selectd segment
	selectedSegments.rewind();
	String currentSource="";
	while((seg=selectedSegments.getSeg())!=NULL){                	
		unsigned long begin=seg->begin()+fs.getFirstFeatureIndexOfASource(seg->sourceName()); 
		if (currentSource!=seg->sourceName()) {
			currentSource=seg->sourceName();
			this->getUX(m_xh_1,currentSource);
			this->getSpeakerModel(clientMixture,currentSource);			
			if (verbose)cout << "Processing speaker["<<currentSource<<"]"<< endl;	
		}		
		fs.seekFeature(begin);
		Feature f;
		if (!_topGauss) {
			for (unsigned long idxFrame=0;idxFrame<seg->length();idxFrame++){
				fs.readFeature(f,0);
				double *ff=f.getDataVector();				
				double sum=0.0;
				RealVector <double> P;
				P.setSize(_mixsize);
				double *Prob=P.getArray();
				for(unsigned long k=0;k<_mixsize;k++) {
					Prob[k]=clientMixture.weight(k)*clientMixture.getDistrib(k).computeLK(f);
					sum+=Prob[k];
					}
				for(unsigned long k=0;k<_mixsize;k++) 
					Prob[k]/=sum; 
				for(unsigned long k=0;k<_mixsize;k++) {
					for (unsigned long i=0;i<_vsize;i++) 
						ff[i]-= Prob[k]*_m_xh_1[k*_vsize+i];
					}
				fs.writeFeature(f);
				nt++;		
			}	
		}
		else {
			throw Exception("no topgauss yet",__FILE__,__LINE__);
		}
	}
};	

/*void FactorAnalysisStat::normalizeFeatures(TopGauss &tg,SegCluster &selectedSegments,FeatureServer &fs,Config & config){
	if (verbose) cout << "(FactorAnalysisStat) Normalize Features" << endl;			
	unsigned long nbg, k1;
	unsigned long nt=0;
	RealVector <double> m_xh_1; m_xh_1.setSize(_supervsize); 	
	double *_m_xh_1=m_xh_1.getArray();
	Seg *seg;          // current selectd segment
	selectedSegments.rewind();
	
	unsigned long*_idx=tg.idx().getArray();
	unsigned long*_nbg=tg.nbg().getArray();		

	while((seg=selectedSegments.getSeg())!=NULL){                	
		unsigned long begin=seg->begin()+fs.getFirstFeatureIndexOfASource(seg->sourceName()); 
		fs.seekFeature(begin);
		Feature f;
		unsigned long idxBegin=tg.frameToIdx(nt); // get the begin index of gaussians in idx vectors			
		for (unsigned long idxFrame=0;idxFrame<seg->length();idxFrame++){
			fs.readFeature(f,0);
			double* ff=f.getDataVector();
			nbg=_nbg[nt]; // get nb gauss for this frame				
			double sum=0.0;
			RealVector <double> P;
			P.setSize(nbg);
			double *Prob=P.getArray();
			for(unsigned long k=0;k<nbg;k++) {
				Prob[k]=clientMixture.weight(_idx[k+idxBegin])*clientMixture.getDistrib(_idx[k+idxBegin]).computeLK(f);
				sum+=Prob[k];
				}
			for(unsigned long k=0;k<nbg;k++) 
				Prob[k]/=sum; 
			for(unsigned long k=0;k<nbg;k++) {
				k1=_idx[k+idxBegin];
				for (unsigned long i=0;i<_vsize;i++) 
					ff[i]-= Prob[k]*_m_xh_1[k1*_vsize+i];
				}
			fs.writeFeature(f);
			nt++;
			idxBegin+=nbg;
		}
	}		
};*/

void FactorAnalysisStat::normalizeFeatures(FeatureServer &fs,Config & config){			
	SegServer segmentsServer;
	LabelServer labelServer;
	initializeClusters(fileList,segmentsServer,labelServer,config);
	verifyClusterFile(segmentsServer,fs,config);
	unsigned long codeSelectedFrame=labelServer.getLabelIndexByString(config.getParam("labelSelectedFrames"));
	SegCluster& selectedSegments=segmentsServer.getCluster(codeSelectedFrame);  
	normalizeFeatures(selectedSegments,fs,config);
};

void FactorAnalysisStat::estimateXYAndNorm(FeatureServer &fs,Config &config) {
	this->resetXY();

	for(int i=0;i<config.getParam("nbTrainIt").toLong();i++){
		if (verbose) cout << "------ Iteration ["<<i<<"] ------"<<endl;
		this->computeAndAccumulateGeneralFAStats(fs,config);
		this->substractSpeakerStats();
		this->substractChannelStats(); 
		this->estimateAndInverseL(config);
		this->getXEstimate();
		this->getYEstimate();
	  }
	this->normalizeFeatures(fs,config);
};

void FactorAnalysisStat::estimateXYAndNorm(SegCluster &selectedSegments,FeatureServer &fs,Config &config) {
	this->resetXY();
	for(int i=0;i<config.getParam("nbTrainIt").toLong();i++){
		if (verbose) cout << "------ Iteration ["<<i<<"] ------"<<endl;
		this->computeAndAccumulateGeneralFAStats(selectedSegments,fs,config);
		this->substractSpeakerStats();
		this->substractChannelStats(); 
		this->estimateAndInverseL(config);
		this->getXEstimate();
		this->getYEstimate();
	  }
	this->normalizeFeatures(selectedSegments,fs,config);
};
#endif

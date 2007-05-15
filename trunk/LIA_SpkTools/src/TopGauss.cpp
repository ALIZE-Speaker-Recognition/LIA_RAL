// TopGauss.cpp
// Main Author: Driss Matrouf
// Integration in LIA_RAL: Nico
// Used to compute posterior and write top indexes to files, this needs to be changed to a proper way to do that in ALIZE

#if !defined(ALIZE_TopGauss_cpp)
#define ALIZE_TopGauss_cpp

#define myTINY 0.0000001

#include <iostream>
#include <fstream>  // pour outFile
#include <cstdio>   // pour printf()
#include <cassert> // pour le debug pratique
#include <cmath>
#include "liatools.h"
#include "TopGauss.h"
#define EPS_LK 1e-200
unsigned long TopGauss::frameToIdx(unsigned long & f) {
	unsigned long cnt=0;
	unsigned long*_pnbg=_nbg.getArray();
	for (unsigned long t=0;t<f;t++)
		cnt+=_pnbg[t];
return cnt;
}

void TopGauss::read(String & featureFilename,Config &config) {
	if (verbose) cout << "(TopGauss::read) File ["<<featureFilename<<"] ... ";
	String nbGaussPerFrameFilename=config.getParam("nbGaussianFilesDir")+featureFilename;
	if (verboseLevel > 2) cout << "(TopGauss) Reading indexes in file ["<<nbGaussPerFrameFilename<<"]"<<endl;	
	ifstream nbGaussFile(nbGaussPerFrameFilename.c_str(),ios::in);
	if(!nbGaussFile){throw Exception("Cannot find nbGaussian file", __FILE__, __LINE__);}	
	// Read number of frames
	nbGaussFile.read((char *)&_nt,sizeof(unsigned long));
	nbGaussFile.read((char *)&_nbgcnt,sizeof(unsigned long));	
	_idx.setAllValues(0);_idx.setSize(_nbgcnt);
	_nbg.setAllValues(0);_nbg.setSize(_nt);
	_snsl.setAllValues(0);_snsl.setSize(_nt);
	_snsw.setAllValues(0);_snsw.setSize(_nt);	
	nbGaussFile.read((char *)_nbg.getArray(),sizeof(unsigned long)*_nt);
	nbGaussFile.read((char *)_idx.getArray(),sizeof(unsigned long)*_nbgcnt);
	nbGaussFile.read((char *)_snsw.getArray(),sizeof(double)*_nt);
	nbGaussFile.read((char *)_snsl.getArray(),sizeof(double)*_nt);	
	nbGaussFile.close();
}

void TopGauss::write(String & featureFilename,Config & config) {
	if (verbose) cout << "(TopGauss::write) File ["<<featureFilename<<"] "<<endl;	
	String nbGaussFilename=config.getParam("nbGaussianFilesDir")+featureFilename;
	ofstream nbGaussFile (nbGaussFilename.c_str(),ios::out);
	if(!nbGaussFile){throw Exception("Cannot find nbGaussian file", __FILE__, __LINE__);}	
	nbGaussFile.write((char *)&_nt,sizeof(unsigned long));		// write number of frames in cluster
	nbGaussFile.write((char *)&_nbgcnt,sizeof(unsigned long));		// write number of frames in cluster	
	// Write number of Gaussians + Gaussian Indexes;
	nbGaussFile.write((char *)_nbg.getArray(),sizeof(unsigned long)*_nbg.size());
	nbGaussFile.write((char *)_idx.getArray(),sizeof(unsigned long)*_idx.size());
	nbGaussFile.write((char *)_snsw.getArray(),sizeof(double)*_snsw.size());
	nbGaussFile.write((char *)_snsl.getArray(),sizeof(double)*_snsl.size());	
	nbGaussFile.close();
}

// Init with a Filename 
double TopGauss::compute(MixtureGD & UBM,String & featureFilename, Config & config) {
	if (verbose) cout << "(TopGauss::compute on a single file) File ["<<featureFilename<<"] ... ";
	FeatureServer fs(config,featureFilename);  
	double llk=compute(UBM,fs,featureFilename,config);	
	if (verbose) cout <<_nbgcnt/_nt <<" per frames LLK="<<llk<<endl;	
return llk;
}

//Init with a FeatureServer
RealVector <double> TopGauss::compute(MixtureGD & UBM,FeatureServer &fs,Config & config) {
	RealVector <double> llkv;
	for (unsigned long i=0;i<fs.getSourceCount();i++) {
		String filename=fs.getNameOfASource(i);
		if (verbose) cout << "(TopGauss::compute on FeatureServer) File ["<<filename<<"] ... ";			
		llkv.addValue(compute(UBM,fs,filename,config));
		cout << "_nt"<<_nt<<" idx"<<_idx.size()<<" _snsw"<<_snsw.size()<<endl;
		this->write(filename,config); // this should be removed!! but when passing a feature server things get worse		
		if (verbose) cout <<_nbgcnt/_nt <<" per frames LLK="<<llkv[i]<<endl;	
	}
return llkv;
}


// Main init function
double TopGauss::compute(MixtureGD & UBM,FeatureServer &fs,String & featureFilename,Config & config){
	StatServer ss(config);
	MixtureGDStat &acc=ss.createAndStoreMixtureStat(UBM);	
	unsigned long _mixsize=UBM.getDistribCount();
	String labelSelectedFrames =config.getParam("labelSelectedFrames");
	unsigned long begin=fs.getFirstFeatureIndexOfASource(featureFilename);
	fs.seekFeature(begin);
	SegServer segmentsServer;
	LabelServer labelServer;
	initializeClusters(featureFilename,segmentsServer,labelServer,config);
	verifyClusterFile(segmentsServer,fs,config);
	unsigned long codeSelectedFrame=labelServer.getLabelIndexByString(labelSelectedFrames);	
	SegCluster& selectedSegments=segmentsServer.getCluster(codeSelectedFrame);  
	acc.resetLLK();
	double topD=config.getParam("topGauss").toDouble();
	if (verbose) {if(topD<1.0) cout << "LLK %="<< topD << "% ";else cout << "Top-"<<topD<<" ";}
	
	// Class values
	_nt=totalFrame(selectedSegments);	
	_nbg.setSize(_nt); _idx.setSize(0);_snsw.setSize(0); _snsl.setSize(0);
	_nbg.setAllValues(0); _idx.setAllValues(0);_snsw.setAllValues(0.0);_snsl.setAllValues(0.0);
	_nbgcnt=0;
	Seg *seg;          // current selected segment
	selectedSegments.rewind();		
	unsigned long t=0; //cnt frames
	while((seg=selectedSegments.getSeg())!=NULL){                       	
		unsigned long begin=seg->begin()+fs.getFirstFeatureIndexOfASource(seg->sourceName()); 
		fs.seekFeature(begin);
		Feature f;
		for (unsigned long idxFrame=0;idxFrame<seg->length();idxFrame++){
			fs.readFeature(f); 
			double llk=acc.computeAndAccumulateLLK(f,1.0,DETERMINE_TOP_DISTRIBS);
			const LKVector &topV=ss.getTopDistribIndexVector();
			double lk_tot=exp(llk);
			
			double val=0.0;
			if (topD<1.0) {
				for(unsigned long j=0;j<_mixsize;j++){
					if (val > topD*lk_tot) break;
					val+=(topV[j].lk);
					_nbg[t]++;
				}
			} else _nbg[t]=(unsigned long)topD;
			_nbgcnt+=_nbg[t];
			 
			double snsw=1.0;
			double snsl=lk_tot;					
			for(unsigned long j=0;j<_nbg[t];j++) {
				_idx.addValue(topV[j].idx);    		
				snsw -=UBM.weight(topV[j].idx);
				snsl -=topV[j].lk;
			}

			_snsw.addValue(snsw);
			if (snsl < EPS_LK)
				_snsl.addValue(EPS_LK);
			else _snsl.addValue(snsl);
			t++;
		}		
	}
	if (t!=_nt) cout << "W: t("<<t<<") != _nt(" <<_nt<<")"<<endl;
return acc.getMeanLLK();
}



RealVector <double> TopGauss::get(MixtureGD & UBM,FeatureServer &fs,Config & config){
	RealVector <double> llkv;
	for (unsigned long i=0;i<fs.getSourceCount();i++) {
		String filename=fs.getNameOfASource(i);
		if (verbose) cout << "(TopGauss::get on FeatureServer) File ["<<filename<<"] ... ";	
		this->read(filename,config); // this should be removed!! but when passing a feature server things get worse			
		llkv.addValue(get(UBM,fs,filename,config));
		if (verbose) cout <<_nbgcnt/_nt <<" per frames LLK="<<llkv[i]<<endl;	
	}
return llkv;
}

// Can use this function to get likelihood with a topgauss
double TopGauss::get(MixtureGD & UBM,String & featureFilename,Config & config){
	FeatureServer fs(config,featureFilename);  
	if (verbose) cout << "(TopGauss::get on a single file) File ["<<featureFilename<<"] ... ";	
	double llk=get(UBM,fs,featureFilename,config);	
	if (verbose) cout <<_nbgcnt/_nt <<" per frames LLK="<<llk<<endl;	
return llk;	
}

// Can use this function to get likelihood with a topgauss
/*double TopGauss::get(MixtureGD & UBM,FeatureServer &fs,String & featureFilename,Config & config){
	double llkcluster=0.0;
	StatServer ss(config);
	String labelSelectedFrames =config.getParam("labelSelectedFrames");
	unsigned long begin=fs.getFirstFeatureIndexOfASource(featureFilename);
	fs.seekFeature(begin);
	SegServer segmentsServer;
	LabelServer labelServer;
	initializeClusters(featureFilename,segmentsServer,labelServer,config);
	verifyClusterFile(segmentsServer,fs,config);
	unsigned long codeSelectedFrame=labelServer.getLabelIndexByString(labelSelectedFrames);	
	SegCluster& selectedSegments=segmentsServer.getCluster(codeSelectedFrame);  
	
	Seg *seg;          // current selected segment
	selectedSegments.rewind();		
	unsigned long t=0; //cnt frames
	double minLLK=config.getParam("minLLK").toDouble();
	double maxLLK=config.getParam("maxLLK").toDouble();	
	while((seg=selectedSegments.getSeg())!=NULL){  
		double llkseg=0.0;
		unsigned long begin=seg->begin()+fs.getFirstFeatureIndexOfASource(seg->sourceName()); 
		fs.seekFeature(begin);
		Feature f;
		for (unsigned long idxFrame=0;idxFrame<seg->length();idxFrame++){
			double lkframe=0.0;
			fs.readFeature(f); 
			unsigned long idx=this->frameToIdx(t);
			unsigned long nbg=_nbg[t];
			for (unsigned long i=0;i<nbg;i++) {
				double lkdistrib=0.0;
				lkdistrib=(UBM.weight(_idx[idx+i])*UBM.getDistrib(_idx[idx+i]).computeLK(f));
				if (std::isnan(lkdistrib)) lkdistrib=exp(minLLK);
				lkframe+=lkdistrib;
			}
			if (log(lkframe) <= minLLK) llkseg+=minLLK;
			else if(log(lkframe) >= maxLLK) llkseg+=maxLLK;
			else llkseg+=log(lkframe);			
			t++;
		}	
		llkcluster+=llkseg;
	}
	if (t!=_nt) cout << "W: t("<<t<<") != _nt(" <<_nt<<")"<<endl;
return llkcluster/totalFrame(selectedSegments);
}*/

// Can use this function to get likelihood with a topgauss
double TopGauss::get(MixtureGD & UBM,FeatureServer &fs,String & featureFilename,Config & config){
	StatServer ss(config);
	String labelSelectedFrames =config.getParam("labelSelectedFrames");
	unsigned long begin=fs.getFirstFeatureIndexOfASource(featureFilename);
	fs.seekFeature(begin);
	SegServer segmentsServer;
	LabelServer labelServer;
	initializeClusters(featureFilename,segmentsServer,labelServer,config);
	verifyClusterFile(segmentsServer,fs,config);
	unsigned long codeSelectedFrame=labelServer.getLabelIndexByString(labelSelectedFrames);	
	SegCluster& selectedSegments=segmentsServer.getCluster(codeSelectedFrame);  
	MixtureGDStat &acc=ss.createAndStoreMixtureStat(UBM);
	
	Seg *seg;          // current selected segment
	selectedSegments.rewind();		
	unsigned long t=0; //cnt frames
	acc.resetLLK();
	unsigned long idxBegin=0;
	while((seg=selectedSegments.getSeg())!=NULL){  
		unsigned long begin=seg->begin()+fs.getFirstFeatureIndexOfASource(seg->sourceName()); 
		fs.seekFeature(begin);
		Feature f;
		idxBegin=this->frameToIdx(t);
		for (unsigned long idxFrame=0;idxFrame<seg->length();idxFrame++){
			fs.readFeature(f); 
			//unsigned long idx=this->frameToIdx(t);
			unsigned long nbg=_nbg[t];	
			ULongVector index;
			double sumNonSelectedWeights=_snsw[t];
			double sumNonSelectedLLK=_snsl[t];
			for (unsigned long i=0;i<nbg;i++) {
				index.addValue(_idx[idxBegin+i]);
			}		
			char c[100];
			sprintf(c,"%d",(int)index.size());
			config.setParam("topDistribsCount",c); // this should be high enough	
			if (t==0) {acc.computeAndAccumulateLLK(f,1.0,DETERMINE_TOP_DISTRIBS);acc.resetLLK();} // to remove in ALIZE, this is to init the LKvector
			ss.setTopDistribIndexVector(index, sumNonSelectedWeights, sumNonSelectedLLK);
			acc.computeAndAccumulateLLK(f,1.0,USE_TOP_DISTRIBS);
			idxBegin+=nbg;
			t++;
		}	
	}	
	//ss.deleteMixtureStat(acc);
	if (t!=_nt || idxBegin !=_nbgcnt) cout << "W: t("<<t<<") != _nt(" <<_nt<<")"<<"W: idxBegin("<<idxBegin<<") != _nbgcnt(" <<_nbgcnt<<")"<<endl;
return acc.getMeanLLK();
}

#endif

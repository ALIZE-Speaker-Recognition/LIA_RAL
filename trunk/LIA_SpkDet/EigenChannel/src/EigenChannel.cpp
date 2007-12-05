#if !defined(ALIZE_EigenChannel_cpp)
#define ALIZE_EigenChannel_cpp

#include <fstream>
#include <cstdio>		
#include <cassert>
#include <cmath>
#include "liatools.h"
#include <EigenChannel.h>

using namespace alize;
using namespace std;

void verifyEMLK(FactorAnalysisStat & FA,XList &ndx,FeatureServer &fs,Config &config) {
	XLine *pline; String *pFile; ndx.rewind();	
	double total=0.0;
	unsigned long maxLLKcomputed=1;
	maxLLKcomputed=config.getParam("computeLLK").toLong();	
	bool FALLK=false;
	if (config.existsParam("FALLK")) {FALLK=true;if(verbose) cout<<"(EigenChannel) Computing Factor Analysis Likelihoods"<<endl;}
	unsigned long cnt=0;
	while((pline=ndx.getLine())!=NULL && cnt < maxLLKcomputed) { 
		while((pFile=pline->getElement())!=NULL && cnt < maxLLKcomputed) {
			/// Compute FA model
			MixtureServer ms(config);
			MixtureGD &model=ms.loadMixtureGD(config.getParam("inputWorldFilename"));
			if (FALLK) FA.getFactorAnalysisModel(model,*pFile);
			else FA.getSpeakerModel(model,*pFile);
			
			/// Get LLK
			FeatureServer fs(config,*pFile);
			SegServer segmentsServer;
			LabelServer labelServer;
			initializeClusters(*pFile,segmentsServer,labelServer,config);
			verifyClusterFile(segmentsServer,fs,config);
			unsigned long codeSelectedFrame=labelServer.getLabelIndexByString(config.getParam("labelSelectedFrames"));
			SegCluster& selectedSegments=segmentsServer.getCluster(codeSelectedFrame);  
			double llk=FA.getLLK(selectedSegments,model,fs,config); 
			if (verbose) cout << "(EigenChannel) LLK["<<*pFile<<"]="<<llk<<endl;
			cnt++;
			total+=llk;
		}
	}
	if (verbose) cout << "(EigenChannel) Total LLK="<<total<<endl;
}

int EigenChannel(Config & config){
	unsigned long nbIt=config.getParam("nbIt").toLong();
	bool computeLLK=false;
	if (config.existsParam("computeLLK")) computeLLK=true;	
	bool init=false;
	if (config.existsParam("loadAccs")) init=true;

	XList ndx(config.getParam("ndxFilename"));
	XLine allFiles=ndx.getAllElements();			
	FeatureServer fs;
	
	if (verbose) cout << "** Create Factor Analysis Statistics Accumulator" << endl;
	FactorAnalysisStat FA((XList&)ndx,fs,config);	
	
	// Init, compute Stats on UBM and save for further experiements
	if (!init) {
		fs.init(config,allFiles);		
		FA.computeAndAccumulateGeneralFAStats(fs,config);
		FA.saveAccs(config);		
	}		
	else FA.loadAccs(config);
		
	
	FA.storeAccs(); // save FA state
	for(unsigned long i=0;i<nbIt;i++){
		if (computeLLK) verifyEMLK(FA,ndx,fs,config);
		if (verbose) cout << "(EigenChannel) --------- Iteration "<<i<<"--------- "<<endl;
		FA.estimateAndInverseL(config);
		FA.substractSpeakerStats();
		FA.getXEstimate();
		
		FA.substractChannelStats(); 
		FA.getYEstimate();

		FA.getUEstimate(config);
		if (1) { //should be debug
			String mat=config.getParam("channelMatrix")+".lastIt";
			FA.getU().save(mat,config); // save last it matrix in debug mode
		}
		FA.restoreAccs(); // restore state
	}
	FA.getU().save(config.getParam("channelMatrix"),config);
return 0;
}
#endif 

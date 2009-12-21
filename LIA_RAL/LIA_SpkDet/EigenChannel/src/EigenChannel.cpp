// EigenChannel.cpp
// 
// This file is a part of Mistral Package and LIA Software 
// LIA_SpkDet, based on Mistral_Ral toolkit 
// LIA_SpkDet is a free, open tool for speaker recognition
// LIA_SpkDet is a development project initiated and funded by the LIA lab.
//
// See mistral.univ-avignon.fr 
// 
// ALIZE is needed for LIA_SpkDet
// for more information about ALIZE, see http://alize.univ-avignon.fr
//
// Copyright (C) 2004 - 2005 - 2006 - 2007 -2008
//  Laboratoire d'informatique d'Avignon [www.lia.univ-avignon.fr]
//  Jean-Francois Bonastre [jean-francois.bonastre@univ-avignon.fr]
//      
// Mistral is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
// You should have received a copy of the GNU General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// The LIA team as well as the ALIZE project want to highlight the limits of voice authentication
// in a forensic context. 
// The following paper proposes a good overview of this point:
// [Bonastre J.F., Bimbot F., Boe L.J., Campbell J.P., Douglas D.A., Magrin-chagnolleau I.,
//  Person  Authentification by Voice: A Need of Caution,
//  Eurospeech 2003, Genova]
// The conclusion of the paper of the paper is proposed bellow:
// [Currently, it is not possible to completely determine whether the
//  similarity between two recordings is due to the speaker or to other
//  factors, especially when: (a) the speaker does not cooperate, (b) there
//  is no control over recording equipment, (c) recording conditions are not 
//  known, (d) one does not know whether the voice was disguised and, to a
//  lesser extent, (e) the linguistic content of the message is not
//  controlled. Caution and judgment must be exercised when applying speaker
//  recognition techniques, whether human or automatic, to account for these
//  uncontrolled factors. Under more constrained or calibrated situations,
//  or as an aid for investigative purposes, judicious application of these
//  techniques may be suitable, provided they are not considered as infallible.
//  At the present time, there is no scientific process that enables one to
//  uniquely characterize a person=92s voice or to identify with absolute
//  certainty an individual from his or her voice.]
//
// Contact Jean-Francois Bonastre (jean-francois.bonastre@lia.univ-avignon.fr) for
// more information about the licence or the use of LIA_SpkDet
// First version 15/07/2004
// New version 23/02/2005
// 
// Last review 4 nov 2008
//

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
	bool _computeLLK=false;
	if (config.existsParam("computeLLK")) _computeLLK=true;	
	bool init=false;
	if (config.existsParam("loadAccs")) init=config.getParam("loadAccs").toBool();;

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
		if (_computeLLK) verifyEMLK(FA,ndx,fs,config);
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

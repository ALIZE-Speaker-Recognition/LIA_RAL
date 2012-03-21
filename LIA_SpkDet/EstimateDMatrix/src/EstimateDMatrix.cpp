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

#if !defined(ALIZE_EstimateDMatrix_cpp)
#define ALIZE_EstimateDMatrix_cpp

#include <fstream>
#include <cstdio>		
#include <cassert>
#include <cmath>
#include "liatools.h"
#include <EstimateDMatrix.h>

using namespace alize;
using namespace std;


//-----------------------------------------------------------------------------------------------------------------------------------------------------------
/*void verifyEMLK(JFAAcc & JFA,XList &ndx,Config &config) {

	XLine *pline; String *pFile; ndx.rewind();	
	double total=0.0;
	unsigned long maxLLKcomputed=1;
	maxLLKcomputed=config.getParam("computeLLK").toLong();	
	bool JFALLK=false;
	if (config.existsParam("JFALLK")) {JFALLK=true;if(verbose) cout<<"(EigenVoice) Computing Joint Factor Analysis Likelihoods"<<endl;}
	unsigned long cnt=0;
	while((pline=ndx.getLine())!=NULL && cnt < maxLLKcomputed) { 
		while((pFile=pline->getElement())!=NULL && cnt < maxLLKcomputed) {
			/// Compute JFA model
			MixtureServer ms(config);
			MixtureGD &model=ms.loadMixtureGD(config.getParam("inputWorldFilename"));
			if (JFALLK) JFA.getJFAModel(model,*pFile);
			else JFA.getSpeakerModel(model,*pFile);
			
			/// Get LLK
			FeatureServer fs(config,*pFile);
			SegServer segmentsServer;
			LabelServer labelServer;
			initializeClusters(*pFile,segmentsServer,labelServer,config);
			verifyClusterFile(segmentsServer,fs,config);
			unsigned long codeSelectedFrame=labelServer.getLabelIndexByString(config.getParam("labelSelectedFrames"));
			SegCluster& selectedSegments=segmentsServer.getCluster(codeSelectedFrame);  
			double llk=FA.getLLK(selectedSegments,model,fs,config); 
			if (verbose) cout << "(EigenVoice) LLK["<<*pFile<<"]="<<llk<<endl;
			cnt++;
			total+=llk;
		}
	}
	if (verbose) cout << "(EigenVoice) Total LLK="<<total<<endl;
}*/

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
int EstimateDMatrix(Config & config){
	
	//Read the NDX file
	String ndxFilename = config.getParam("ndxFilename");
	
	//Create and initialise the accumulator
	JFAAcc jfaAcc(ndxFilename, config, "EstimateD");

	//Statistics
	if((config.existsParam("loadAccs")) && config.getParam("loadAccs").toBool()){	//load pre-computed statistics
		jfaAcc.loadN(config);
		jfaAcc.loadN_h(config);
		jfaAcc.loadF_X(config);
		jfaAcc.loadF_X_h(config);
	}
	else{															//Compute statistics if they don't exists
		jfaAcc.computeAndAccumulateJFAStat(config);
		jfaAcc.saveAccs(config);
	}

	//Initialise the D Matrix
	if((config.existsParam("loadInitDMatrix"))&&(config.getParam("loadInitDMatrix").toBool())){	//Load the D matrix when existing
		jfaAcc.loadD(config.getParam("initDMatrix"), config);
	}
	else{	//Initialise the D matrix randomly if does not exists
		jfaAcc.initD(config);
	}

	//Save the initial D matrix to be able restart the process with the same initialisation
	if(config.existsParam("saveInitD") && config.getParam("saveInitD").toBool()){
		String saveFilename = config.getParam("DMatrix") + "_init";
		jfaAcc.saveD(saveFilename, config);
		cout<<"	(EstimateDMatrix) Save the initial D Matrix in "<<saveFilename<<endl;
	}
	
	//Initialise the EV Matrix
	if(config.existsParam("eigenVoiceMatrix")){	//Load the EV matrix when existing
		jfaAcc.loadEV(config.getParam("eigenVoiceMatrix"), config);
	}
	else{	//Initialise the EV matrix to zero if does not exist
		jfaAcc.getV().setAllValues(0.0);
		cout<<"(EstimateDMatrix) Initialise NULL EigenVoice Matrix"<<endl;
	}
	
	//Estimate Y factors for each speaker
	jfaAcc.storeAccs();
	jfaAcc.estimateVEVT(config);
	jfaAcc.estimateAndInverseL_EV(config);
	jfaAcc.substractMplusDZ(config);
	jfaAcc.substractUX(config);
	jfaAcc.estimateY(config);
	//Reinitialise the accumulators
	jfaAcc.restoreAccs();

	//Initialise the EC Matrix
	if(config.existsParam("eigenChannelMatrix")){	//Load the EC matrix when existing
		jfaAcc.loadEC(config.getParam("eigenChannelMatrix"), config);
	}
	else{	//Initialise the EC matrix to zero if does not exist
		jfaAcc.getU().setAllValues(0.0);
		cout<<"	(EstimateDMatrix) Initialise NULL EigenChannel Matrix"<<endl;
	}

	//Estimate X factors for each session
	jfaAcc.storeAccs();
	jfaAcc.estimateUEUT(config);
	jfaAcc.estimateAndInverseL_EC(config);
	jfaAcc.substractMplusVYplusDZ(config);
	jfaAcc.estimateX(config);
	//Reinitialise the accumulators
	jfaAcc.restoreAccs();	

	//iteratively retrain the D matrix
	unsigned long nbIt = config.getParam("nbIt").toULong();

	jfaAcc.storeAccs();
	for(unsigned long it=0; it<nbIt; it++){

		cout<<"	(EstimateDMatrix) --------- start iteration "<<it<<" --------"<<endl;

		//Substract speaker statistics M + VY
		jfaAcc.substractMplusVY(config);
		
		//Substract channel statistics UX
		jfaAcc.substractUX(config);
	
		//Estimate Z for each speaker
		jfaAcc.estimateZandD();

		//Reinitialise the accumulators
		jfaAcc.resetTmpAcc();
		jfaAcc.restoreAccs();
		
		//Save the D matrix at the end of the iteration
		if(config.getParam("saveAllDMatrices").toBool()){
			String s;
			String output = config.getParam("DMatrix") + s.valueOf(it);
			jfaAcc.saveD(output, config);
		}
	}

	cout<<"	(EstimateDMatrix) --------- save D Matrix --------"<<endl;
	jfaAcc.saveD(config.getParam("DMatrix"), config);
	cout<<"	(EstimateDMatrix) --------- end of process --------"<<endl;
	
return 0;
}
#endif 

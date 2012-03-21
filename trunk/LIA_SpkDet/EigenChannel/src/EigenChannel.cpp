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

#if !defined(ALIZE_EigenChannel_cpp)
#define ALIZE_EigenChannel_cpp

#include <iostream>
#include <fstream>  
#include <cstdio>		
#include <cassert>
#include <cmath>
#include "liatools.h"
#include "EigenChannel.h"

using namespace alize;
using namespace std;

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
int EigenChannelJFA(Config & config){
	
	//Read the NDX file
	String ndxFilename = config.getParam("ndxFilename");
	
	//Create and initialise the accumulator
	JFAAcc jfaAcc(ndxFilename, config, "EigenChannelJFA");

	//Option used to check the Likelihood at each iteration
	bool _checkLLK = false;
	if (config.existsParam("checkLLK")) _checkLLK= config.getParam("checkLLK").toBool();
	else if (verboseLevel >=1) _checkLLK= true;

	
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
	
	//Initialise the EC Matrix
	if(config.existsParam("loadInitChannelMatrix") && config.getParam("loadInitChannelMatrix").toBool()){	//Load the EC matrix when existing
		jfaAcc.loadEC(config.getParam("initEigenChannelMatrix"), config);
	}
	else{	//Initialise the EC matrix randomly if does not exists
		jfaAcc.initEC(config);
	}
	
	//Save the initial U matrix to be able restart the process with the same initialisation
	bool saveInitMatrix = false;
	if(config.existsParam("saveInitChannelMatrix"))	saveInitMatrix = config.getParam("saveInitChannelMatrix").toBool();
	if(saveInitMatrix){
		jfaAcc.saveU(config.getParam("initEigenChannelMatrix"), config);
		cout<<"	(EigenChannel) Save the initial EigenChannel Matrix in "<<config.getParam("initEigenChannelMatrix")<<endl;
	}

	//Initialise the EV Matrix
	if(config.existsParam("eigenVoiceMatrix")){	//Load the EV matrix when existing
		jfaAcc.loadEV(config.getParam("eigenVoiceMatrix"), config);
	}
	else{	//Initialise the EV matrix to zero
		jfaAcc.getV().setAllValues(0.0);
		cout<<"	(EigenChannel) Initialise NULL EigenVoice Matrix"<<endl;
	}

	//iteratively retrain the EC matrix
	unsigned long nbIt = config.getParam("nbIt").toULong();
	
	//Estimate Y factors for each speaker
	jfaAcc.storeAccs();
	jfaAcc.estimateVEVT(config);
	jfaAcc.estimateAndInverseL_EV(config);
	jfaAcc.substractMplusDZ(config);
	jfaAcc.substractUX(config);
	jfaAcc.estimateY(config);
	
	//Reinitialise the accumulators
	jfaAcc.restoreAccs();

	jfaAcc.storeAccs();
	for(unsigned long it=0; it<nbIt; it++){
		
		cout<<"	(EigenChannel) --------- start iteration "<<it<<" --------"<<endl;
		if (_checkLLK) jfaAcc.verifyEMLK(config);

		//Compute the vEvT matrices
		jfaAcc.estimateUEUT(config);

		//Estimate and inverse L matrices
		jfaAcc.estimateAndInverseL_EC(config);

		//Substract the speaker statistics
		jfaAcc.substractMplusVYplusDZ(config);

		//On update X pour toutes les sessions
		jfaAcc.estimateXandU(config);
			
		//Update _U
		jfaAcc.updateUestimate();

		//Reinitialise the accumulators
		jfaAcc.resetTmpAcc();
		jfaAcc.restoreAccs();

		//Save the U matrix at the end of the iteration
		bool saveAllMatrices = false;
		if(config.existsParam("saveAllECMatrices")) saveAllMatrices=config.getParam("saveAllECMatrices").toBool();
		if(saveAllMatrices)
		{
			String s;
			String output = config.getParam("eigenChannelMatrix") + s.valueOf(it);
			jfaAcc.saveU(output, config);
		}
	}
	
	cout<<"	(EigenChannel) --------- save EigenChannel Matrix --------"<<endl;
	jfaAcc.saveU(config.getParam("eigenChannelMatrix") ,config);
	cout<<"	(EigenChannel) --------- end of process --------"<<endl;

return 0;
}


//-----------------------------------------------------------------------------------------------------------------------------------------------------------
int EigenChannelLFA(Config & config){
	
	if(verbose) cout<<"	Enter EigenChannel mode LFA"<<endl;

	//Read the NDX file
	String ndxFilename = config.getParam("ndxFilename");
	
	//Create and initialise the accumulator
	JFAAcc jfaAcc(ndxFilename, config);

	//Option used to check the Likelihood at each iteration
	bool _checkLLK = false;
	if (config.existsParam("checkLLK")) _checkLLK= config.getParam("checkLLK").toBool();
	
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
	
	//Initialise the EC Matrix
	if(config.existsParam("loadInitChannelMatrix") && config.getParam("loadInitChannelMatrix").toBool()){	//Load the EC matrix when existing
		jfaAcc.loadEC(config.getParam("initEigenChannelMatrix"), config);
	}
	else{	//Initialise the EC matrix randomly if does not exists
		jfaAcc.initEC(config);
	}
	
	//Save the initial U matrix to be able restart the process with the same initialisation
	if(config.existsParam("saveInitChannelMatrix") && config.getParam("saveInitChannelMatrix").toBool()){
		jfaAcc.saveU(config.getParam("initEigenChannelMatrix"), config);
		cout<<"	(EigenChannel) Save the initial EigenChannel Matrix in "<<config.getParam("initEigenChannelMatrix")<<endl;
	}

	//Initialise the EV Matrix
	if(config.existsParam("eigenVoiceMatrix")){	//Load the EV matrix when existing
		jfaAcc.loadEV(config.getParam("eigenVoiceMatrix"), config);
	}
	else{	//Initialise the EV matrix to zero
		jfaAcc.getV().setAllValues(0.0);
		cout<<"	(EigenChannel) Initialise NULL EigenVoice Matrix"<<endl;
	}


	//Initialise the D Matrix with MAP paradigm
	jfaAcc.initD(config);

	//iteratively retrain the EC matrix
	unsigned long nbIt = config.getParam("nbIt").toULong();

	//Estimate Y factors for each speaker
	jfaAcc.storeAccs();
	jfaAcc.estimateVEVT(config);
	jfaAcc.estimateAndInverseL_EV(config);
	jfaAcc.substractMplusDZ(config);
	jfaAcc.substractUX(config);
	jfaAcc.estimateY(config);
	
	//Reinitialise the accumulators
	jfaAcc.restoreAccs();

	jfaAcc.storeAccs();
	for(unsigned long it=0; it<nbIt; it++){
		
		cout<<"	(EigenChannel) --------- start iteration "<<it<<" --------"<<endl;
		if (_checkLLK) jfaAcc.verifyEMLK(config);

		//Compute the vEvT matrices
		jfaAcc.estimateUEUT(config);

		//Estimate and inverse L matrices
		jfaAcc.estimateAndInverseL_EC(config);

		//Substract the speaker statistics
		jfaAcc.substractMplusVYplusDZ(config);

		//Estimate X
		jfaAcc.estimateX(config);

		//Substract channel statistics
		jfaAcc.substractMplusUX();

		//Estimate Z (Y for JFA ??? or Y and Z dans le cas ou on utilise la matrice D du MAP)
		jfaAcc.estimateZMAP(config.getParam("regulationFactor").toLong());

		//estimate U
		jfaAcc.estimateU();

		//Update _U
		jfaAcc.updateUestimate();

		//Reinitialise the accumulators
		jfaAcc.resetTmpAcc();
		jfaAcc.restoreAccs();

		//Save the U matrix at the end of the iteration
		if(config.getParam("saveAllECMatrices").toBool()){
			String s;
			String output = config.getParam("eigenChannelMatrix") + s.valueOf(it);
			jfaAcc.saveU(output, config);
		}
	}
	
	cout<<"	(EigenChannel) --------- save EigenChannel Matrix --------"<<endl;
	jfaAcc.saveU(config.getParam("eigenChannelMatrix") ,config);
	cout<<"	(EigenChannel) --------- end of process --------"<<endl;

return 0;
}












#endif 

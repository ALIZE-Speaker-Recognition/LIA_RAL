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


#if !defined(ALIZE_TotalVariability_cpp)
#define ALIZE_TotalVariability_cpp

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <liatools.h>
#include "TotalVariability.h"

using namespace std;
using namespace alize;

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
int TotalVariability(Config & config){

	//Read the NDX file
	String ndxFilename = config.getParam("ndxFilename");

	//Create and initialise the accumulator
	TVAcc tvAcc(ndxFilename, config);

	//Option used to check the Likelihood at each iteration
	bool _checkLLK = false;
	if (config.existsParam("checkLLK")) _checkLLK= config.getParam("checkLLK").toBool();
	else if (verboseLevel >=1) _checkLLK= true;
	//Statistics
	if((config.existsParam("loadAccs")) && config.getParam("loadAccs").toBool()){	//load pre-computed statistics
		cout<<"	()Load Accumulators"<<endl;
		tvAcc.loadN(config);
		tvAcc.loadF_X(config);
	}
	else{															//Compute statistics if they don't exists
		tvAcc.computeAndAccumulateTVStat(config);
		tvAcc.saveAccs(config);
	}

	//Initialize the Total Variability matrix
	bool loadInitEigenVoiceMatrix = false;
	if(config.existsParam("loadInitTotalVariabilityMatrix")) loadInitEigenVoiceMatrix = config.getParam("loadInitTotalVariabilityMatrix").toBool();
	if(loadInitEigenVoiceMatrix){	//Load the Total Variability matrix when existing
		tvAcc.loadT(config.getParam("initTotalVariabilityMatrix"), config);
	}
	else{	//Initialize the EV matrix randomly if does not exists
		tvAcc.initT(config);
	}
	
	//Save the initial Total Variability matrix to be able restart the process with the same initialisation
	if(config.existsParam("saveInitTotalVariabilityMatrix") && config.getParam("saveInitTotalVariabilityMatrix").toBool()){
		String initV = config.getParam("totalVariabilityMatrix")+"_init";
		tvAcc.saveT(initV, config);
		cout<<"	(TotalVariability) Save the initial TotalVariability Matrix in "<<initV<<endl;
	}
	
	//Select minimum divergence criteria
	bool minDiv = false;
	if(config.existsParam("minDivergence")) minDiv = config.getParam("minDivergence").toBool();

	//Iteratively retrain the EV matrix
	unsigned long nbIt = config.getParam("nbIt").toULong();

	for(unsigned long it=0; it<nbIt; it++){
		
		cout<<"	(TotalVariability) --------- start iteration "<<it<<" --------"<<endl;

		//Subtract mean from the statistics
		tvAcc.substractM(config);

		//Compute TETt for each distribution
		tvAcc.estimateTETt(config);

		// Compute inverse(L) and estimate TotalVariability matrix
		tvAcc.estimateAandC(config);

		//Compute Likelihood over the data if required
		if (_checkLLK) tvAcc.verifyEMLK(config);

		//Update _T
		tvAcc.updateTestimate();

		//Minimum Divergence step
		if(minDiv){
			tvAcc.minDivergence();
		}

		//If the option is on, orthonormalize the matrix V
		if(config.existsParam("orthonormalizeT") && (config.getParam("orthonormalizeT").toBool())){
			if(verboseLevel > 0) cout<<"Orthonormalize TV matrix"<<endl;
			tvAcc.orthonormalizeT();
		}

		//Reinitialize the accumulators
		tvAcc.resetTmpAcc();

		//Reload statistics
		tvAcc.loadN(config);
		tvAcc.loadF_X(config);

		//Save the T matrix at the end of the iteration
		bool saveAllEVMatrices = false;
		if(config.existsParam("saveAllTVMatrices")) saveAllEVMatrices=config.getParam("saveAllTVMatrices").toBool();
		if(saveAllEVMatrices){
			String s;
			String output = config.getParam("totalVariabilityMatrix") + s.valueOf(it);
			tvAcc.saveT(output, config);
			
			// Save the new mean computed by Minimum Divergence if required
			if(minDiv){
				String sm;
				String outputm = config.getParam("matrixFilesPath") + config.getParam("meanEstimate") + sm.valueOf(it) + config.getParam("saveMatrixFilesExtension");
				((Matrix<double>)tvAcc.getUbmMeans()).save(outputm, config);
			}
		}
	}

	cout<<"	(TotalVariability) --------- save Matrix --------"<<endl;
	tvAcc.saveT(config.getParam("totalVariabilityMatrix"), config);
	if(minDiv){
		String outputm = config.getParam("matrixFilesPath") + config.getParam("meanEstimate") + config.getParam("saveMatrixFilesExtension");
		((Matrix<double>)tvAcc.getUbmMeans()).save(outputm, config);
	}
	cout<<"	(TotalVariability) --------- end of process --------"<<endl;

return 0;
}


#endif 

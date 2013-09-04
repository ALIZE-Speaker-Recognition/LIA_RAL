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


#if !defined(ALIZE_IvNorm_cpp)
#define ALIZE_IvNorm_cpp

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <liatools.h>
#include "IvNorm.h"

using namespace std;
using namespace alize;


//-----------------------------------------------------------------------------------------------------------------------------------------------------------
int IvNorm(Config & config){

	if(verboseLevel>0) cout<<"(IvNorm) Normalize i-vectors"<<endl;

	// Estimate normalization parameters if required
	if(!config.getParam("ivNormLoadParam").toBool()){

		//Initialize development data
		String backgroundNdxFilename = config.getParam("backgroundNdxFilename");
		PldaDev dev(backgroundNdxFilename,config);

		if(config.getParam("ivNormIterationNb").toULong() > 0)	//Estimate normalization parameters
			dev.sphericalNuisanceNormalization(config);

		if(config.existsParam("LDA") && config.getParam("LDA").toBool()){	//if LDA is required
			Matrix<double> ldaMat;
			unsigned long ldaRank = config.getParam("ldaRank").toULong();
			String ldaFilename = config.getParam("matrixFilesPath")+config.getParam("ldaMatrix")+config.getParam("loadMatrixFilesExtension");

			//Compute LDA matrix
			dev.computeLDA(ldaMat,ldaRank,config);
			ldaMat.save(ldaFilename,config);
		}
	}

	// Apply normalization to input vectors if required
	if(config.existsParam("inputVectorFilename") || config.existsParam("ndxFilename")){

		PldaTest test(config);
	
		if(config.getParam("ivNormIterationNb").toULong() > 0)	//Apply normalization
			test.sphericalNuisanceNormalization(config);
	
		if(config.existsParam("LDA") && config.getParam("LDA").toBool()){	// Apply only LDA reduction
			String ldaFilename = config.getParam("matrixFilesPath")+config.getParam("ldaMatrix")+config.getParam("loadMatrixFilesExtension");
			Matrix<double> LDAmat(ldaFilename,config);
			test.rotateLeft(LDAmat);
		}

		//Save test data
		String outputDir = config.getParam("saveVectorFilesPath");

		//if input is a list of files, then save only test segments and not models (they are the same)
		if(config.existsParam("inputVectorFilename") && !config.existsParam("targetIdList"))
			test.saveSegments(outputDir, config);

		//if there are models and test, then save everything
		else{
			test.saveVectors(outputDir, config);
		}
	}

return 0;
}


#endif 

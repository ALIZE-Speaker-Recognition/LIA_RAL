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
 
#if !defined(ALIZE_IvTest_cpp)
#define ALIZE_IvTest_cpp
 
#include <iostream>
#include <fstream> 
#include <cassert> 
#include <cmath> 
#include <liatools.h>
#include "IvTest.h"
#include <sys/stat.h>

using namespace alize;
using namespace std;


//-------------------------------------------------------------------------------------------------------
//	Compute Test for I-Vector front end
//-------------------------------------------------------------------------------------------------------
int IvTest(Config& config){

 try{

	String scoring = config.getParam("scoring");
	String outputNISTFileName = config.getParam("outputFilename");			// Result file in NIST style (.nist) result file format
	String gender=config.getParam("gender");								// gives the gender for compatibility reasons with NIST standard output file

	real_t decisionThreshold;
	if (config.existsParam("decisionThreshold"))							// Define the threshold if needed
		decisionThreshold=config.getParam("decisionthreshold").toDouble();
	else decisionThreshold=0;
	unsigned long maxClientLine=CST_MAX_CLIENT_LINE;						// Max of target Id for a ndx line
	if (config.existsParam("maxTargetLine")) maxClientLine=config.getParam("maxTargetLine").toLong();

	if (config.getParam_debug())debug=true;  else debug=false;

	bool WCCN = false;
	bool computeWCCN = false;
	if(config.existsParam("wccn") && config.getParam("wccn").toBool()){
		WCCN = config.getParam("wccn").toBool();
		if(config.existsParam("loadWccnMatrix") && !config.getParam("loadWccnMatrix").toBool())
			computeWCCN = !config.getParam("loadWccnMatrix").toBool();
			
	}

	bool loadMahalanobisMatrix = false;
	if((config.getParam("scoring")=="mahalanobis") && (config.getParam("loadMahalanobisMatrix").toBool()))
		loadMahalanobisMatrix = config.getParam("loadMahalanobisMatrix").toBool();

	bool loadWccnMatrix = false;
	if((config.getParam("scoring")=="cosine") && config.getParam("wccn").toBool())
		loadWccnMatrix = config.getParam("loadWccnMatrix").toBool();

	bool load2covMatrix = false;
	if((config.getParam("scoring")=="2cov") &&  config.getParam("load2covMatrix").toBool())
		load2covMatrix = config.getParam("load2covMatrix").toBool();

	// Load test data
	PldaTest pldaTest(config);

	// Load backgorund data in case one of the following condition is fullfilled but not in the case of PLDA scoring
	//	config.existsParam("ivNorm")&&config.getParam("ivNorm").toBool()&& !config.getParam("ivNormLoadParam").toBool()
	// (scoring=="mahalanobis") && (!loadMahalanobisMatrix)
	// (scoring=="cosine") && (!loadWccnMatrix)
	// (scoring=="2cov") && (!load2covMatrix)
	//
	if(( (config.existsParam("ivNorm")&&config.getParam("ivNorm").toBool()&& !config.getParam("ivNormLoadParam").toBool()) ||	
		((scoring=="mahalanobis") && (!loadMahalanobisMatrix)) ||
		(WCCN && (!loadWccnMatrix)) ||
		((scoring=="2cov") && (!load2covMatrix))
		) &&
		scoring != "plda"
		){

		//Initialize development data
		if(verboseLevel>0) cout<<"(IvTest)	Load development data"<<endl;
		String backgroundNdxFilename = config.getParam("backgroundNdxFilename");

		PldaDev dev(backgroundNdxFilename,config);

		//Compute normalization parameters including LDA
		if(config.existsParam("ivNorm")&&config.getParam("ivNorm").toBool()){

			if(verboseLevel>0) cout<<"(IvTest) Normalize i-vectors"<<endl;

			// Estimate normalization parameters if required
			if(!config.getParam("ivNormLoadParam").toBool()){

				if(config.getParam("ivNormIterationNb").toULong() > 0){	//Estimate normalization parameters
					if(verboseLevel>0) cout<<"(IvTest)	Estimate EFR parameters"<<endl;
					dev.sphericalNuisanceNormalization(config);
				}

				if(config.existsParam("LDA") && config.getParam("LDA").toBool()){	//if LDA is required

					if(verboseLevel>0) cout<<"(IvTest)	Estimate LDA parameters"<<endl;

					Matrix<double> ldaMat;
					unsigned long ldaRank = config.getParam("ldaRank").toULong();
					String ldaFilename = config.getParam("matrixFilesPath")+config.getParam("ldaMatrix")+config.getParam("loadMatrixFilesExtension");

					//Compute LDA matrix
					dev.computeLDA(ldaMat,ldaRank,config);
					
					//Apply LDA on development data
					if(verboseLevel>0)	cout<<"(IvTest)		Apply LDA on dev data"<<endl;
					dev.rotateLeft(ldaMat);

					//Save LDA matrix
					ldaMat.save(ldaFilename,config);
				}
			}

//If normalization parameters already exist and should be apply on dev data
else{

   	if(config.getParam("ivNormIterationNb").toULong() > 0){ //Estimate normalization parameters
                if(verboseLevel>0) cout<<"(IvTest)      Estimate EFR parameters"<<endl;
        	dev.applySphericalNuisanceNormalization(config);
     	}

	if(config.existsParam("LDA") && config.getParam("LDA").toBool()){       //if LDA is required
		//Load LDA matrix
        	Matrix<double> ldaMat;
        	String ldaFilename = config.getParam("matrixFilesPath")+config.getParam("ldaMatrix")+config.getParam("loadMatrixFilesExtension");
        	cout<<"         Load LDA matrix from: "<<ldaFilename<<endl;
        	ldaMat.load (ldaFilename,config);	

        	//Apply LDA on development data
        	if(verboseLevel>0)      cout<<"(IvTest)         Apply LDA on dev data"<<endl;
		dev.rotateLeft(ldaMat);
	}

}

		}

		//Compute Within Class Covariance Matrix for Cosine scoring
		if(computeWCCN){
			if(verboseLevel>0) cout<<"(IvTest)	Estimate WCCN parameters"<<endl;
			
			String wccnFilenameTmp = "WCCN";
			if(config.existsParam("wccnMatrix")){
				wccnFilenameTmp = config.getParam("wccnMatrix");
			}

			String wccnFilename = config.getParam("matrixFilesPath")+ wccnFilenameTmp +config.getParam("loadMatrixFilesExtension");

			Matrix<double> W;
			DoubleSquareMatrix WCCN;
			WCCN.setSize(1);
			WCCN.setAllValues(0.0);
			dev.computeWccnChol(WCCN,config);

			Matrix<double> tmpW(WCCN);
			W = tmpW;

			// Save WCCN matrix
			W.save(wccnFilename,config);
		}

		//Compute Mahalanobis matrix
		if(scoring=="mahalanobis"){
			if(verboseLevel>0) cout<<"	Compute Mahalanobis matrix"<<endl;

			//Get Mahalanobis matrix filename
			String mahalanobisFilename = "Mahalanobis";
			if(config.existsParam("mahalanobisMatrix")){
				mahalanobisFilename = config.getParam("matrixFilesPath")+config.getParam("mahalanobisMatrix")+config.getParam("saveMatrixFilesExtension");
			}
		
			//Compute Mahalanobis matrix
			Matrix<double> Mah;

			//Compute Mahalanobis matrix
			DoubleSquareMatrix M;
			dev.computeMahalanobis(M,config);
			Matrix<double> tmpM(M);
			Mah = tmpM;
			Mah.save(mahalanobisFilename,config);
		}

		//Compute covariance matrices for 2cov scoring
		if(scoring=="2cov"){

			DoubleSquareMatrix W, B, Sigma;
			dev.computeCovMat(Sigma,W,B,config);

			//Save covariance matrices
			String W2covFilename = "2Cov_W";
			String B2covFilename = "2Cov_B";
			if(config.existsParam("TwoCovFilename")){
				W2covFilename = config.getParam("TwoCovFilename")+"_W";
				B2covFilename = config.getParam("TwoCovFilename")+"_B";
			}
			W2covFilename = config.getParam("matrixFilesPath")+W2covFilename+config.getParam("saveMatrixFilesExtension");
			B2covFilename = config.getParam("matrixFilesPath")+B2covFilename+config.getParam("saveMatrixFilesExtension");

			Matrix<double> tmpW(W);
			Matrix<double> tmpB(B);
			tmpW.save(W2covFilename,config);
			tmpB.save(B2covFilename,config);
		}
	}

	//For PLDA scoring, normalization is done separately using PldaDev from the PldaModel and PLDA training data
	if((scoring=="plda") && (!config.getParam("pldaLoadModel").toBool())){

		//Create PLDA model
		PldaModel plda("train",config);

		//If normalization is required, estimate parameters here
		if((config.existsParam("ivNorm")&&config.getParam("ivNorm").toBool()&& !config.getParam("ivNormLoadParam").toBool())){

			if(verboseLevel>0) cout<<"(IvTest) Normalize i-vectors"<<endl;

			if(config.getParam("ivNormIterationNb").toULong() > 0){	//Estimate normalization parameters
				if(verboseLevel>0) cout<<"(IvTest)	Estimate EFR parameters"<<endl;
				plda.getDev().sphericalNuisanceNormalization(config);
			}

			if(config.existsParam("LDA") && config.getParam("LDA").toBool()){	//if LDA is required

				if(verboseLevel>0) cout<<"(IvTest)	Estimate LDA parameters"<<endl;

				Matrix<double> ldaMat;
				unsigned long ldaRank = config.getParam("ldaRank").toULong();
				String ldaFilename = config.getParam("matrixFilesPath")+config.getParam("ldaMatrix")+config.getParam("loadMatrixFilesExtension");

				//Compute LDA matrix
				plda.getDev().computeLDA(ldaMat,ldaRank,config);
					
				//Apply LDA on development data
				if(verboseLevel>0)	cout<<"(IvTest)		Apply LDA on dev data"<<endl;
				plda.getDev().rotateLeft(ldaMat);

				//Save LDA matrix
				ldaMat.save(ldaFilename,config);
			}
		}

		//Estimate PLDA model parameters
		// Center data
		plda.updateModel(config);
		plda.centerData();

		// EM iterations
		unsigned long nbIt = config.getParam("pldaNbIt").toULong();
		for(unsigned long it=0;it<nbIt;it++){
			if(verboseLevel>0) cout<<"	(PldaModel)	EM iteration [ "<<it<<" ]"<<endl;
			plda.em_iteration(config,it);
		}

		// save PLDA model
		plda.saveModel(config);
	}


	//Normalize test data before scoring
	if(config.existsParam("ivNorm")&&config.getParam("ivNorm").toBool()){

		if(config.getParam("ivNormIterationNb").toULong() > 0){
			pldaTest.sphericalNuisanceNormalization(config);
		}
		if(config.existsParam("LDA") && config.getParam("LDA").toBool()){

			//Load LDA matrix
			Matrix<double> ldaMat;
			String ldaFilename = config.getParam("matrixFilesPath")+config.getParam("ldaMatrix")+config.getParam("loadMatrixFilesExtension");
			cout<<"		Load LDA matrix from: "<<ldaFilename<<endl;
			ldaMat.load (ldaFilename,config);

			//Apply LDA to test data
			pldaTest.rotateLeft(ldaMat);
		}
	}

	//SCORING
	if(scoring=="cosine"){

		if(WCCN){
			//rotate test data
				String wccnFilenameTmp = "WCCN";
				if(config.existsParam("wccnMatrix")){
					wccnFilenameTmp = config.getParam("wccnMatrix");
				}

				String wccnFilename = config.getParam("matrixFilesPath")+ wccnFilenameTmp +config.getParam("loadMatrixFilesExtension");
				Matrix<double> Wccn(wccnFilename,config);
				pldaTest.rotateLeft(Wccn);
		}

		// Scoring
		pldaTest.cosineDistance(config);
	}

	else if(scoring=="mahalanobis"){

		//Get Mahalanobis matrix
		String mahalanobisFilename = "Mahalanobis";
		if(config.existsParam("mahalanobisMatrix")){
			mahalanobisFilename = config.getParam("matrixFilesPath")+config.getParam("mahalanobisMatrix")+config.getParam("loadMatrixFilesExtension");
		}
		
		Matrix<double> Mah;
		//Load Mahalanobis
		Mah.load (mahalanobisFilename,config);

		//Scoring
		pldaTest.mahalanobisDistance(Mah,config);
	}

	else if(scoring=="2cov"){

		//Load covariance matrices
		String W2covFilename = "2Cov_W";
		String B2covFilename = "2Cov_B";
		if(config.existsParam("TwoCovFilename")){
			W2covFilename = config.getParam("TwoCovFilename") + "_W";
			B2covFilename = config.getParam("TwoCovFilename") + "_B";
		}
		W2covFilename = config.getParam("matrixFilesPath")+ W2covFilename + config.getParam("loadMatrixFilesExtension");
		B2covFilename = config.getParam("matrixFilesPath")+ B2covFilename + config.getParam("loadMatrixFilesExtension");
		

		Matrix<double> tmpW(W2covFilename,config);
		Matrix<double> tmpB(B2covFilename,config);
		DoubleSquareMatrix W,B;
		W.setSize(tmpW.rows()); B.setSize(tmpB.rows());
		for(unsigned long ii=0;ii<tmpW.rows();ii++)
			for(unsigned long jj=0;jj<tmpW.rows();jj++){
				W(ii,jj) = tmpW(ii,jj);
				B(ii,jj) = tmpB(ii,jj);
			}

		//Scoring
		pldaTest.twoCovScoring(W, B, config);
	}

	else if(scoring=="plda"){

		PldaModel plda("test",config);

		// Select plda scoring mode: Multiple-Segment LLR, Mean, ScoreFusion
		String pldaScoring = "native";
		if(config.existsParam("pldaScoring")) pldaScoring = config.getParam("pldaScoring");
		
		if(pldaScoring == "enrollMean"){
			pldaTest.pldaMeanScoring(plda,config);
		}
		else if((pldaScoring == "native") || (pldaTest.getMaxEnrollmentSession() == 1)){
			pldaTest.pldaNativeScoring(plda,config);
		}
	}
	else{
		cout<<"Scoring option is invalid, must be: cosine OR mahalanobis OR 2cov OR plda"<<endl;
	}

	String outputScoreFormat = "ascii";
	if(config.existsParam("outputScoreFormat")) outputScoreFormat = config.getParam("outputScoreFormat");

	if(outputScoreFormat == "ascii"){
		ofstream outNist(outputNISTFileName.c_str(),ios::out | ios::trunc);    				// Initialise the output file
		if(verboseLevel >0) cout<<"Score computation done, writing the file in ASCII file"<<endl;
	
		double *_s;
		bool *_t;
		Matrix<double> _S(pldaTest.getScores());
		BoolMatrix Trials(pldaTest.getTrials());
		_s = _S.getArray();
		_t = Trials.getArray();

		unsigned long segNb = pldaTest.getSegmentsNumber();
		unsigned long modNb = pldaTest.getModelsNumber();

		// Write scores in ASCII outputFile
		for(unsigned long s=0;s<segNb;s++){
			for(unsigned long m=0;m<modNb;m++){
				if(_t[m*segNb+s]){
					char decision=setDecision(_s[m*segNb+s],decisionThreshold);                       // take a decision
					outputResultLine(_s[m*segNb+s], pldaTest.getModelName(m),pldaTest.getSegmentName(s) ,gender ,decision,outNist);
				}
			}
		}
		outNist.close(); 
	}

	else if(outputScoreFormat == "binary"){

		unsigned long segNb = pldaTest.getSegmentsNumber();
		unsigned long modNb = pldaTest.getModelsNumber();

		Matrix<double> _S(pldaTest.getScores());

		//save model list in one file
		String modelFilename = config.getParam("outputFilename") + "_model.txt";
		ofstream outModel(modelFilename.c_str(),ios::out | ios::trunc);    				// Initialise the output file
		for(unsigned long m=0;m<modNb;m++)
			outModel<<pldaTest.getModelName(m)<<endl;
		outModel.close();

		//save segment list in one file
		String segFilename = config.getParam("outputFilename") + "_testSeg.txt";
		ofstream outSegment(segFilename.c_str(),ios::out | ios::trunc);    				// Initialise the output file
		for(unsigned long ts=0;ts<segNb;ts++)
			outSegment<<pldaTest.getSegmentName(ts)<<endl;
		outSegment.close();

		//save matrix of score in Binary format (DB)
		String scoreFilename = config.getParam("outputFilename") + config.getParam("saveMatrixFilesExtension");
		_S.save(scoreFilename,config);
	}
}// fin try
catch (Exception& e){ 
	cout << e.toString().c_str() << endl;
}
return 0;
}



























#endif //!defined(ALIZE_ComputeTest_cpp)

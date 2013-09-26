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

#if !defined(ALIZE_IvExtractor_cpp)
#define ALIZE_IvExtractor_cpp

#include <iostream>
#include <fstream>  // pour outFile
#include <cstdio>   // pour printf()
#include <cassert> // pour le debug pratique
#include <cmath>
#include <liatools.h>
#include "IvExtractor.h"

using namespace alize;
using namespace std;

// Information on the quantity of data available by client
// could output a list with the selected files for a defined quantity of data
// Uses the same input file than IvExtractor and output a new list
int InfoTarget(Config& config)
{
  String inputClientListFileName = config.getParam("targetIdList");
  bool fixedLabelSelectedFrame;
  String labelSelectedFrames;
  if (config.existsParam("useIdForSelectedFrame"))      // the ID of each speaker is used as labelSelectedFrame
    fixedLabelSelectedFrame=false;
  else{                                                // the label is decided by the command line and is unique for the run
    labelSelectedFrames=config.getParam("labelSelectedFrames");
    if (verbose) cout << "Computing on" << labelSelectedFrames << " label" << endl;
    fixedLabelSelectedFrame=true;
  }
  unsigned long maxFrame=config.getParam("maxFrame").toLong();
  String outputFilename=config.getParam("outputFilename");
  
  
  ofstream outputFile(outputFilename.c_str(),ios::out| ios::trunc);
  try{
    XList inputClientList(inputClientListFileName,config);          // read the Id + filenames for each client
    XLine * linep;
    if (verbose) cout << "InfoTarget" << endl;
    // *********** Target loop ***************** 
    while ((linep=inputClientList.getLine()) != NULL){             // linep gives the XLine with the Id of a given client and the list of files
      String *id=linep->getElement();                              // Get the Client ID (id)
      outputFile<<*id;
      String currentFile="";
      XLine featureFileListp=linep->getElements();	           // Get the list of feature file for the client (end of the line)
      if (verbose) cout << "Info model ["<<*id<<"]"<<endl;   
      if (!fixedLabelSelectedFrame){                                // the ID is used as label for selecting the frame
	labelSelectedFrames=*id;
	if (debug) cout <<*id<<" is used for label selected frames"<<endl;
      }
      // label files reading - It creates, for each file and each label, a cluster of segments - will be integrated witth the featre s - asap
      SegServer segmentsServer;                                    // Reading the segmentation files for each feature input file
      LabelServer labelServer;                                     
      initializeClusters(featureFileListp,segmentsServer,labelServer,config);           // Reading the segmentation files for each feature input file
      unsigned long codeSelectedFrame=labelServer.getLabelIndexByString(labelSelectedFrames);            // Get the index of the cluster with in interest audio segments
      SegCluster& selectedSegments=segmentsServer.getCluster(codeSelectedFrame); // Gives the cluster of the selected/used segments   
      Seg *seg;                                                                  // Will give the current segment
      unsigned long frameCount=0;
      selectedSegments.rewind();                                                 // at the begin of the selected segments list
      while(((seg=selectedSegments.getSeg())!=NULL) && (frameCount<maxFrame)){   // For each of the selected segments until the amount of data is get
	frameCount+=seg->length();
	cout << seg->sourceName()<<" "<<seg->begin()<<" "<<seg->length()<<" Total time="<<frameCount<<endl;
	if (seg->sourceName()!=currentFile){                                     
	  outputFile<<" "<<seg->sourceName();
	  currentFile=seg->sourceName();
	}
      }                                                                          // end of the initial Train Iteration loop
      outputFile<<endl;
      if (verbose) cout << "Save info client ["<<*id<<"]" << endl;
    }                                                                            // end of the the target loop 
  } // fin try

  catch (Exception& e)
    { 
      cout << e.toString().c_str() << endl;
    }
  return 0;
}
 
//-----------------------------------------------------------------------------------------------------------------------------------------------------------
int IvExtractor(Config& config)
{
	String inputWorldFilename = config.getParam("inputWorldFilename");

	// label for selected frames - Only the frames associated with this label, in the label files, will be used
	bool fixedLabelSelectedFrame=true;
	String labelSelectedFrames;
	if (config.existsParam("useIdForSelectedFrame"))    // the ID of each speaker is used as labelSelectedFrame ?
		fixedLabelSelectedFrame=(config.getParam("useIdForSelectedFrame").toBool()==false);  
	if (fixedLabelSelectedFrame)                        // the label is decided by the command line and is unique for the run
		labelSelectedFrames=config.getParam("labelSelectedFrames");
 
try{
	MixtureServer ms(config);
	if (verbose) cout << "(IvExtractor) TotalVariability - Load world model [" << inputWorldFilename<<"]"<<endl;
	MixtureGD& world = ms.loadMixtureGD(inputWorldFilename);      
	
	//Load the statistics from files or compute statistics for all segments at once
	//Read the NDX file
	String ndxFilename = config.getParam("targetIdList");

	//Remove the first element of each line which is the model name
	XList tmpFileList(ndxFilename);
	XList fileList;
	for(unsigned long ll=0;ll<tmpFileList.getLineCount();ll++){
		fileList.addLine();
		for(unsigned long i=1;i<tmpFileList.getLine()->getElementCount();i++){
			fileList.getLine(fileList.getLineCount()-1).addElement(tmpFileList.getLine(ll).getElement(i));
		}
	}

	//Create and initialise the accumulator
	TVAcc tvAcc(fileList, config);
	
	// Load TotalVariability Matrix
	tvAcc.loadT(config.getParam("totalVariabilityMatrix"),config);


	//Statistics
	if((config.existsParam("loadAccs")) && config.getParam("loadAccs").toBool()){	//load pre-computed statistics
		cout<<"	(IvExtractor) Load Accumulators"<<endl;
		tvAcc.loadN(config);
		tvAcc.loadF_X(config);
	}
	else{															//Compute statistics if they don't exists
		tvAcc.computeAndAccumulateTVStat(config);
		tvAcc.saveAccs(config);
	}

	// Then load the meanEstimate computed by minDiv if required
	DoubleVector meanEstimate = tvAcc.getUbmMeans();
	if(config.existsParam("minDivergence")&& config.getParam("minDivergence").toBool()){
		String minDivName = config.getParam("matrixFilesPath") + config.getParam("meanEstimate") + config.getParam("loadMatrixFilesExtension");
		Matrix<double> tmpMean(minDivName,config);
		for(unsigned long i=0;i<meanEstimate.size();i++){
			meanEstimate[i] = tmpMean(0,i);
		}
	}
	//Update the mean Estimate
	cout<<"	(IvExtractor) Load Mean Estimate"<<endl;
	tvAcc.loadMeanEstimate(meanEstimate);

	//Substract mean from the statistics
	tvAcc.substractM(config);

	//Compute vEvT for each session
	tvAcc.estimateTETt(config);

	// Estimate I-Vectors
	tvAcc.estimateW(config);

	cout<<"--------- save IV by File --------"<<endl;
	tvAcc.saveWbyFile(config);
	cout<<"--------- end of process --------"<<endl;

} // fin try
catch (Exception& e) {cout << e.toString().c_str() << endl;}
return 0;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
int IvExtractorUbmWeigth(Config& config)
{
	String inputWorldFilename = config.getParam("inputWorldFilename");

	// label for selected frames - Only the frames associated with this label, in the label files, will be used
	bool fixedLabelSelectedFrame=true;
	String labelSelectedFrames;
	if (config.existsParam("useIdForSelectedFrame"))    // the ID of each speaker is used as labelSelectedFrame ?
		fixedLabelSelectedFrame=(config.getParam("useIdForSelectedFrame").toBool()==false);  
	if (fixedLabelSelectedFrame)                        // the label is decided by the command line and is unique for the run
		labelSelectedFrames=config.getParam("labelSelectedFrames");
 
try{
	MixtureServer ms(config);
	if (verbose) cout << "(IvExtractor) Approximate i-vector by using UBM's weight parameters"<<endl;
	if (verbose) cout << "(IvExtractor) TotalVariability - Load world model [" << inputWorldFilename<<"]"<<endl;
	MixtureGD& world = ms.loadMixtureGD(inputWorldFilename);      

	//Read the NDX file
	String ndxFilename = config.getParam("targetIdList");

	//Remove the first element of each line which is the model name
	XList tmpFileList(ndxFilename);
	XList fileList;
	for(unsigned long ll=0;ll<tmpFileList.getLineCount();ll++){
		fileList.addLine();
		for(unsigned long i=1;i<tmpFileList.getLine()->getElementCount();i++){
			fileList.getLine(fileList.getLineCount()-1).addElement(tmpFileList.getLine(ll).getElement(i));
		}
	}
	
	//Create and initialise the accumulator
	TVAcc tvAcc(fileList, config);

	DoubleSquareMatrix W(tvAcc.getRankT());

	if(config.existsParam("loadUbmWeightParam") && config.getParam("loadUbmWeightParam").toBool()){	// Load normalized T matrix and weighted Covariance matrix if pre-computed
		//Load TotalVariability matrix
		String normTFilename = config.getParam("totalVariabilityMatrix") + "_norm";
		tvAcc.loadT(normTFilename, config);

		//Load weighted Covariance matrix
		String wFilename = config.getParam("matrixFilesPath") + config.getParam("totalVariabilityMatrix") + "_weightedCov" + config.getParam("loadMatrixFilesExtension");
		Matrix<double> tmpW(wFilename,config);
		for(unsigned long i=0;i<tvAcc.getRankT();i++)
			for(unsigned long j=0;j<tvAcc.getRankT();j++)
				W(i,j)=tmpW(i,j);

	}
	else{
		//Load TotalVariability matrix
		tvAcc.loadT(config.getParam("totalVariabilityMatrix"), config);

		// Normalize matrix T
		tvAcc.normTMatrix();

		// Compute weighted co-variance matrix by using UBM weight coefficients
		W.setAllValues(0.0);
		tvAcc.getWeightedCov(W,world.getTabWeight(),config);
	}

	//Load the statistics from files or compute statistics for all segments at once
	if((config.existsParam("loadAccs")) && config.getParam("loadAccs").toBool()){	//load pre-computed statistics
		cout<<"	(IvExtractor) Load Accumulators"<<endl;
		tvAcc.loadN(config);
		tvAcc.loadF_X(config);
	}
	else{															//Compute statistics if they don't exists
		tvAcc.computeAndAccumulateTVStat(config);
		tvAcc.saveAccs(config);
	}

	// Then load the meanEstimate computed by minDiv if required
	DoubleVector meanEstimate = tvAcc.getUbmMeans();
	if(config.existsParam("minDivergence")&& config.getParam("minDivergence").toBool()){
		String minDivName = config.getParam("matrixFilesPath") + config.getParam("meanEstimate") + config.getParam("loadMatrixFilesExtension");
		Matrix<double> tmpMean(minDivName,config);
		for(unsigned long i=0;i<meanEstimate.size();i++){
			meanEstimate[i] = tmpMean(0,i);
		}
	}
	//Update the mean Estimate
	cout<<"	(IvExtractor) Load Mean Estimate"<<endl;
	tvAcc.loadMeanEstimate(meanEstimate);

	//Substract mean from the statistics and normalize co-variance
	tvAcc.normStatistics(config);

	// Estimate I-Vectors
	tvAcc.estimateWUbmWeight(W,config);

	cout<<"--------- save IV by File --------"<<endl;
	tvAcc.saveWbyFile(config);
	cout<<"--------- end of process --------"<<endl;


} // fin try
catch (Exception& e) {cout << e.toString().c_str() << endl;}
return 0;
}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
int IvExtractorEigenDecomposition(Config& config)
{
	String inputWorldFilename = config.getParam("inputWorldFilename");

	// label for selected frames - Only the frames associated with this label, in the label files, will be used
	bool fixedLabelSelectedFrame=true;
	String labelSelectedFrames;
	if (config.existsParam("useIdForSelectedFrame"))    // the ID of each speaker is used as labelSelectedFrame ?
		fixedLabelSelectedFrame=(config.getParam("useIdForSelectedFrame").toBool()==false);  
	if (fixedLabelSelectedFrame)                        // the label is decided by the command line and is unique for the run
		labelSelectedFrames=config.getParam("labelSelectedFrames");
 
try{
	MixtureServer ms(config);
	if (verbose) cout << "(IvExtractor) Approximate i-vector by using Eigen Decomposition"<<endl;
	if (verbose) cout << "(IvExtractor) TotalVariability - Load world model [" << inputWorldFilename<<"]"<<endl;
	MixtureGD& world = ms.loadMixtureGD(inputWorldFilename);      

	//Read the NDX file
	String ndxFilename = config.getParam("targetIdList");

	//Remove the first element of each line which is the model name
	XList tmpFileList(ndxFilename);
	XList fileList;
	for(unsigned long ll=0;ll<tmpFileList.getLineCount();ll++){
		fileList.addLine();
		for(unsigned long i=1;i<tmpFileList.getLine()->getElementCount();i++){
			fileList.getLine(fileList.getLineCount()-1).addElement(tmpFileList.getLine(ll).getElement(i));
		}
	}
	
	//Create and initialise the accumulator
	TVAcc tvAcc(fileList, config);
	Matrix<double> Q(tvAcc.getRankT(),tvAcc.getRankT());
	Matrix<double> D(tvAcc.getNDistrib(),tvAcc.getRankT());

	if(config.existsParam("loadEigenDecompositionParam") && config.getParam("loadEigenDecompositionParam").toBool()){	// Load normalized T matrix and weighted Covariance matrix if pre-computed

		//Load TotalVariability matrix
		String normTFilename = config.getParam("totalVariabilityMatrix") + "_norm";
		tvAcc.loadT(normTFilename, config);

		//Load D and Q matrices
		String dFilename = config.getParam("matrixFilesPath") + config.getParam("totalVariabilityMatrix") + "_EigDec_D" + config.getParam("loadMatrixFilesExtension");
		D.load(dFilename,config);

		String qFilename = config.getParam("matrixFilesPath") + config.getParam("totalVariabilityMatrix") + "_EigDec_Q" + config.getParam("loadMatrixFilesExtension");
		Q.load(qFilename,config);

	}
	else{
		//Load TotalVariability matrix
		tvAcc.loadT(config.getParam("totalVariabilityMatrix"), config);

		// Normalize matrix T
		tvAcc.normTMatrix();

		// Compute weighted co-variance matrix by using UBM weight coefficients
		DoubleSquareMatrix W(tvAcc.getRankT());
		W.setAllValues(0.0);
		tvAcc.getWeightedCov(W,world.getTabWeight(),config);

		// Eigen Decomposition of W to get Q
		Matrix<double> tmpW(W);
		Q.setAllValues(0.0);
		tvAcc.computeEigenProblem(tmpW,Q,tvAcc.getRankT(),config);

		// Compute D matrices (approximation of Tc'Tc matrices)
		D.setAllValues(0.0);
		tvAcc.approximateTcTc(D,Q,config);
	}

	//Load the statistics from files or compute statistics for all segments at once
	if((config.existsParam("loadAccs")) && config.getParam("loadAccs").toBool()){	//load pre-computed statistics
		cout<<"	(IvExtractor) Load Accumulators"<<endl;
		tvAcc.loadN(config);
		tvAcc.loadF_X(config);
	}
	else{															//Compute statistics if they don't exists
		tvAcc.computeAndAccumulateTVStat(config);
		tvAcc.saveAccs(config);
	}

	// Then load the meanEstimate computed by minDiv if required
	DoubleVector meanEstimate = tvAcc.getUbmMeans();
	if(config.existsParam("minDivergence")&& config.getParam("minDivergence").toBool()){
		String minDivName = config.getParam("matrixFilesPath") + config.getParam("meanEstimate") + config.getParam("loadMatrixFilesExtension");
		Matrix<double> tmpMean(minDivName,config);
		for(unsigned long i=0;i<meanEstimate.size();i++){
			meanEstimate[i] = tmpMean(0,i);
		}
	}
	//Update the mean Estimate
	cout<<"	(IvExtractor) Load Mean Estimate"<<endl;
	tvAcc.loadMeanEstimate(meanEstimate);

	//Substract mean from the statistics and normalize co-variance
	tvAcc.normStatistics(config);

	// Estimate I-Vectors
	tvAcc.estimateWEigenDecomposition(D,Q,config);

	cout<<"--------- save IV by File --------"<<endl;
	tvAcc.saveWbyFile(config);
	cout<<"--------- end of process --------"<<endl;


} // fin try
catch (Exception& e) {cout << e.toString().c_str() << endl;}
return 0;
}

#endif //!defined(ALIZE_IvExtractor_cpp)

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

#if !defined(ALIZE_Scoring_cpp)
#define ALIZE_Scoring_cpp

#include <iostream>
#include <fstream>  // pour outFile
#include <cstdio>   // pour printf()
#include <cassert> // pour le debug pratique
#include <cmath>
#include "Scoring.h"
#include <liatools.h>


using namespace alize;
using namespace std;

//-------------------------------------------------------------------------
// Take a decision, print t if above threshold f otherwise
String getDecision(double score, double threshold)
{
	String decision;
	if (score > threshold)
	{decision="t";}	
	else 
	{decision="f";}
	return decision;
}

//-------------------------------------------------------------------------
// Take a decision, print true if above threshold false otherwise
String getLongDecision(double score, double threshold)
{
	String decision;
	if (score > threshold)
	{decision="true";}	
	else 
	{decision="false";}
	return decision;
}

//-------------------------------------------------------------------------
// Retrieve info in a nist file providing the good fields
void retrieveNISTSegmentInfo(XList & inputList, String & gender, String & clientName, String & seg, unsigned long genderField, unsigned long nameField, unsigned long segField, unsigned long & i) {
	if (inputList.getLine(i).getElement(genderField)=="F") {gender="f";} else {gender="m";}
	clientName=inputList.getLine(i).getElement(nameField);
	seg=inputList.getLine(i).getElement(segField);
}

//-------------------------------------------------------------------------
// 
unsigned long getIndexOfMaxScore(XList & inputList,unsigned long scoreField, unsigned long segField, unsigned long & i, unsigned long nbListLines)
{
	String seg=inputList.getLine(i).getElement(segField);	//this is a new segment
	long max_score=-200;
	long score;
	unsigned long maxIndex=0; 

	while (inputList.getLine(i).getElement(segField)==seg)	// while same segment test
	{
		score= inputList.getLine(i).getElement(scoreField).toLong();
		if (score >= max_score) {max_score=score; maxIndex=i;}	// store max score and its index
		i++;
		if (i >= nbListLines) break; // break the loop if end of file
	}
return maxIndex;
}

//-------------------------------------------------------------------------
// Produce a tab with mean and cov by segment (mdtm and etf only)
void getSegmentalMeanCov(XList & inputList, unsigned long nbLoc, double* meanTab, double * covTab)
{
	XLine *linep;
	unsigned long cpt=0;
	unsigned long cpttab=0;
	double meanAcc=0.0;
	double covAcc=0.0;
     
	while((linep=inputList.getLine())!=NULL)
	{
		double score=linep->getElement(6).toDouble();
		meanAcc+=score;			// Accumulate mean and cov
		covAcc+=(score*score);
		if (cpt%nbLoc==(nbLoc-1))	// when changing segment
		{
			meanTab[cpttab]=meanAcc;
			covTab[cpttab]=covAcc;
			cpttab++;	
			meanAcc=0.0;
			covAcc=0.0;
		}
		cpt++;
	}
}

//-------------------------------------------------------------------------
// Get a tab with indexes of speakers with maximum likelihood (mdtm and etf only)
void getTarClientIdx(Config & config, XList & inputList, unsigned long nbLoc, unsigned long * tarTab)
{
	XLine *linep;
	unsigned long cpt=0;
	unsigned long cpttab=0;
	double minLLK=config.getParam("minLLK").toDouble();
	double maxScore=minLLK;
	unsigned long idxTar=0;
    bool verbose=config.existsParam("verbose");
    
	while((linep=inputList.getLine())!=NULL)
	{
		double score=linep->getElement(6).toDouble();
		if (score>=maxScore)
		{
			maxScore=score;		// this is the maximum score
			idxTar=cpt;			// index is just the line
			if (verbose) {cout << "giving highest score to " << linep->getElement(1) << " "<<maxScore << endl;} 
		}                         
		if (cpt%nbLoc==(nbLoc-1))	// when changing segment
		{
			tarTab[cpttab]=idxTar;   // idx of the target goes in the tab
			if (verbose) {cout << linep->getElement(1) << " max score: "<<maxScore <<"idx: "<<idxTar<<"cpt: "<<cpt<< endl;}
			cpttab++;	
			maxScore=minLLK; 	//reset maxScore
		}
		cpt++;
	}
}

//-------------------------------------------------------------------------
// Produce a tab with mean and cov by segment without the maximum score(mdtm and etf only)
void getSegmentalMeanCovWithoutMax(Config & config, XList & inputList, unsigned long nbLoc, unsigned long * tarTab, double* meanTab, double * covTab)
{
	XLine *linep;
	unsigned long cpt=0;
	unsigned long cpttab=0;
	double minLLK=config.getParam("minLLK").toDouble();
	double maxScore=minLLK;
	double meanAcc=0.0;
	double covAcc=0.0;
	unsigned long idxTar=0;
	bool verbose=config.existsParam("verbose");

	while((linep=inputList.getLine())!=NULL)
	{
		double score=linep->getElement(6).toDouble();
		if (score>=maxScore)
		{
			maxScore=score;		// this is the maximum score
			idxTar=cpt;			// index is just the line
			if (verbose) {cout << "giving highest score to " << linep->getElement(1) << " "<<maxScore << endl;} 
		}
		meanAcc+=score;			// Accumulate mean and cov
		covAcc+=(score*score);
          
		if (cpt%nbLoc==(nbLoc-1))	// when changing segment
		{	
			tarTab[cpttab]=idxTar;
			meanAcc-=maxScore;	//remove max from Stats
			covAcc-=(maxScore*maxScore);
			meanTab[cpttab]=meanAcc;
			covTab[cpttab]=covAcc;
			if (verbose) {cout << linep->getElement(1) << " max score: "<<maxScore <<"idx: "<<idxTar<<"cpt: "<<cpt<< " meanA: "<<meanAcc<<" covA: "<<covAcc<<endl;}
			cpttab++;	
			maxScore=minLLK;
			meanAcc=0.0;
			covAcc=0.0;
		}
		cpt++;
	}
}  

int Scoring(Config& config)
{

	using namespace alize;
	using namespace std;

	try{

	if (config.existsParam("debug"))debug=true; else debug=false;  
	if (config.existsParam("verbose"))verbose=true; else verbose=false;
	String Mode = config.getParam("mode");
	String inputFileName = config.getParam("inputFile");
	String outputFileName = config.getParam("outputFile");
	bool hard=config.existsParam("hardDecision"); //Force a true decision to be taken among the test (identification)
	double threshold=0.0;
	if (!hard) threshold = config.getParam("threshold").toDouble(); // Ask for a threshold in verification mode
	else if (verbose) cout << "W: Force a decision to be taken" << endl;	// else display a warning
	XList inputList(inputFileName,config);
	ofstream outFile(outputFileName.c_str(),ios::out | ios::trunc);
		
		if ( Mode == "NIST")
		{
			String segTypeTest= config.getParam("segTypeTest");
			String trainTypeTest= config.getParam("trainTypeTest");
			String adaptationMode= config.getParam("adaptationMode");
			unsigned long genderField, decisionField, segField, scoreField, nameField;
			String gender, clientName, seg, decision;
			setLIAInfoFields(genderField, nameField, decisionField, segField, scoreField);	// Set fields position for a LIA file format F name - seg score
			unsigned long i=0;
			unsigned long maxIndex;
			double LLR=0.0;
			
			while (i < inputList.getLineCount())
			{
				if (!hard) // an ACCPET decision is not compulsory
					{
					decision=getDecision(inputList.getLine(i).getElement(scoreField).toDouble(),threshold);
					retrieveNISTSegmentInfo(inputList, gender, clientName, seg, genderField, nameField, segField, i);
					LLR=inputList.getLine(i).getElement(scoreField).toDouble();
					outputResultNIST04Line(trainTypeTest,adaptationMode,segTypeTest,gender,clientName,seg,decision,LLR,outFile);
					i++;
				}
				else {
					maxIndex=getIndexOfMaxScore(inputList,scoreField,segField,i,inputList.getLineCount());
					decision="true";						// decision is true
					retrieveNISTSegmentInfo(inputList, gender, clientName, seg, genderField, nameField, segField, maxIndex);
					LLR=inputList.getLine(maxIndex).getElement(scoreField).toDouble();
					outputResultNIST04Line(trainTypeTest,adaptationMode,segTypeTest,gender,clientName,seg,decision,LLR,outFile); // just output the interesting lines
					}
				
			}
		}
		else if ( Mode == "leaveMaxOutTnorm")
		{
			unsigned long nbLoc=config.getParam("nbLoc").toLong(); 
			unsigned long dimTabs=(unsigned long)(inputList.getLineCount()/nbLoc)+1;	//number of segments is number of lines divided by nb of segments
			unsigned long * tarTab=new unsigned long [dimTabs];   		//tab for target speaker idx
			std::fill_n( tarTab, dimTabs, static_cast<unsigned long>( 0 ) );
			double* meanTab=new double [dimTabs];					 //tab for t-norm mean
			std::fill_n( meanTab, dimTabs, static_cast<double>( 0 ) );
			double* covTab=new double [dimTabs];						//tab for t-norm cov
			std::fill_n( covTab, dimTabs, static_cast<double>( 0 ) );    
			double LLR;
			unsigned long tabcpt=0;
			unsigned long cpt=0;
			unsigned long scoreField=6;

			getSegmentalMeanCovWithoutMax(config, inputList, nbLoc, tarTab, meanTab, covTab); //return target idx in tarTab		 
			inputList.rewind();

			for (unsigned long i=0; i < inputList.getLineCount(); i++) 	// Loop in the XList
			{
				String subType=inputList.getLine(i).getElement(0); 
				String event=inputList.getLine(i).getElement(1);
				String channel="1";
				String source=inputList.getLine(i).getElement(3);
				String start=inputList.getLine(i).getElement(4);
				String duration=inputList.getLine(i).getElement(5);	

				if (cpt%nbLoc==0 && cpt!=0) 
					{tabcpt++;}				// Increments when changing segment
				double tmpScore=inputList.getLine(i).getElement(scoreField).toDouble();
				double mean=0.0;
				double cov=0.0;
				if (cpt==tarTab[tabcpt]) // if this score is the highest
					{
					mean=meanTab[tabcpt]/(nbLoc-1); // don't remove score from mean & cov accumulators
					cov=covTab[tabcpt]/(nbLoc-1)-(mean*mean); 
					}
				else
					{
					mean=(meanTab[tabcpt]-tmpScore)/(nbLoc-2); // remove score from mean & cov accumulators
					cov=(covTab[tabcpt]-(tmpScore*tmpScore))/(nbLoc-2)-(mean*mean); 
					}	
				LLR=(tmpScore-mean)/(sqrt(cov));
				cpt++;
				outputResultLIARALLine(subType,event,channel,source,start,duration, LLR,outFile);
			}
			delete []meanTab;
			delete []covTab;
			delete []tarTab;
		}

		else if ( Mode == "ETF")
		{

			unsigned long nbLoc=config.getParam("nbLoc").toLong(); 
			unsigned long dimTabs=(unsigned long)(inputList.getLineCount()/nbLoc);	//number of segments is number of lines divided by nb of segmentsi
			unsigned long * tarTab=new unsigned long [dimTabs];   //tab for target speaker idx
			std::fill_n( tarTab, dimTabs, static_cast<unsigned long>( 0 ) );
			double LLR;
			String decision="unknown";
			unsigned long scoreField=6;
			unsigned long tabcpt=0;
			unsigned long cpt=0;

			ofstream TARFile("MAX",ios::out | ios::trunc);
			ofstream NONFile("MIN",ios::out | ios::trunc);

			if	(config.existsParam("hard"))				//For Hard Decision, assume tests are not cross-gendered
			{
				cerr << "W: Applying a hard decision, 1 decision by segment" << endl;
				getTarClientIdx(config, inputList, nbLoc, tarTab); //return target idx in tarTab		 
				inputList.rewind();
			}

			for (unsigned long i=0; i < inputList.getLineCount(); i++) 			// Loop in the XList
			{
				String type="spk";
				String subType;
				if (inputList.getLine(i).getElement(0)=="F") 
					{subType="female";}
				else if (inputList.getLine(i).getElement(0)=="M")
					{subType="male";}
				else {subType="unknown";}
				String event=inputList.getLine(i).getElement(1);
				String channel="1";
				String source=inputList.getLine(i).getElement(3);
				String start=inputList.getLine(i).getElement(4);
				double duration=inputList.getLine(i).getElement(5).toDouble()-inputList.getLine(i).getElement(4).toDouble()-0.01;	

				LLR=inputList.getLine(i).getElement(scoreField).toDouble(); // get Score according to the Scorefield
				if (config.existsParam("hard"))
				{
					if (cpt==tarTab[tabcpt])
					{
						decision=getLongDecision(inputList.getLine(i).getElement(scoreField).toDouble(),threshold);
						outputResultETFLine(source, channel, start, duration, type, subType, event, LLR, decision, TARFile);
						tabcpt++;
					}
					else 
					{	
						decision="false";
						outputResultETFLine(source, channel, start, duration, type, subType, event, LLR, decision, NONFile);
					}
					cpt++;
				}
				else
				{
					decision=getDecision(inputList.getLine(i).getElement(scoreField).toDouble(),threshold);
				}
				outputResultETFLine(source, channel, start, duration, type, subType, event, LLR, decision, outFile);
			}
			delete []tarTab;
		}	
		else if ( Mode == "MDTM")
		{
			unsigned long scoreField=6;
			for (unsigned long i=0; i < inputList.getLineCount(); i++)
			{
				String type="speaker";
				String subType;
				if (inputList.getLine(i).getElement(0)=="F") 
				{subType="adult_female";} 
				else if (inputList.getLine(i).getElement(0)=="M")
				{subType="adult_male";}
				else if (inputList.getLine(i).getElement(0)=="C")
				{subType="child";}
				else {subType="unknown";}
				String channel=inputList.getLine(i).getElement(1);
				String source=inputList.getLine(i).getElement(3);
				String start=inputList.getLine(i).getElement(4);
				double duration=(inputList.getLine(i).getElement(5).toDouble()-inputList.getLine(i).getElement(4).toDouble());	
				double LLR=inputList.getLine(i).getElement(scoreField).toDouble();
				outputResultMDTMLine(source, channel, start, duration, type, LLR, subType, outFile);
			}
		}
		else { cerr << "ERROR: Unknown Mode" << endl; exit(1);}
		outFile.close();

	} // fin try

catch (Exception& e)
	{ 
		cout << e.toString().c_str() << endl;
	}
return 0;
}

int WarpScores(Config & config) {
	using namespace alize;
	using namespace std;
	
if (debug) cout <<endl<< "Warp scores"<<endl;
try {
	// Defines NIST04 fields  //Nico / remove ASAP, patch in preparation with retrieve info
	//unsigned long scoreField=7;
	

	if (config.existsParam("debug"))debug=true; else debug=false;  
	if (config.existsParam("verbose"))verbose=true; else verbose=false;
	String inputFileName = config.getParam("inputNISTFile");
	String outputFileName = config.getParam("outputNISTFile");
	String inputHistoFilename = config.getParam("inputHisto");
	String destHistoFilename = config.getParam("destHisto");
	double nonObserved=0.0;
	if (config.existsParam("nonObserved")) nonObserved=config.getParam("nonObserved").toDouble();
	double area=0.0;
	if (config.existsParam("refArea")) area=config.getParam("refArea").toDouble();	 
	String  trainTypeTest, adaptationMode, segTypeTest, gender, clientName, decision, seg;
	double newLLR, LLR=0.0;
	unsigned long i=0;
	XList inputList(inputFileName,config);
if (debug) cout <<"file["<<inputFileName<<"] loaded"<<endl;
	Histo destH;
	destH.load(destHistoFilename);
	Histo inputH;
	inputH.load(inputHistoFilename);
	if (debug) cout <<"histo loaded"<<endl;
	ofstream outFile(outputFileName.c_str(),ios::out | ios::trunc);
	/*	
	if (verbose) cerr<< "Compute Input Histo" << endl;
	while (i < inputList.getLineCount()) { //building histo of input file
			LLR=inputList.getLine(i).getElement(scoreField).toDouble();
			input.accumulateValue(LLR);
			i++;
			}
	input.computeHisto();
	*/
	/*
	i=0;
	while (i < inputRef.getLineCount()) { //building histo of reference file
	  LLR=inputRef.getLine(i).getElement(scoreField).toDouble();
	  ref.accumulateValue(LLR);
	  i++;
	}
	ref.computeHisto();
	}
	*/
	if (verbose) cerr << "Warping in process...." << endl;
	i=0;
	while (i < inputList.getLineCount()) {
	  String trainTypeTest=inputList.getLine(i).getElement(0);
	  String adaptationMode=inputList.getLine(i).getElement(1);
	  String segTypeTest=inputList.getLine(i).getElement(2);
	  unsigned long genderField=3; // remove all this ASAP patch in preparation
	  unsigned long nameField=4;
	  unsigned long segField=5;
	  String decision="-";
	  retrieveNISTSegmentInfo(inputList, gender, clientName, seg, genderField,nameField, segField, i);
	  LLR=inputList.getLine(i).getElement(7).toDouble();
	  newLLR=scoreWarping(LLR, inputH,destH,nonObserved,area);
	  outputResultNIST04Line(trainTypeTest,adaptationMode,segTypeTest,gender,clientName,seg,decision,newLLR,outFile);
	  i++; 
	}
} // fin try

catch (Exception& e)
	{cout << e.toString().c_str() << endl;}
return 0;	
}

#endif //!defined(ALIZE_Scoring_cpp)

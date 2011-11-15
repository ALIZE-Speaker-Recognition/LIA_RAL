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

#if !defined(ALIZE_ComputeNorm_cpp)
#define ALIZE_ComputeNorm_cpp

#include <iostream>
#include <fstream>  
#include <cstdio> 
#include <cassert> 
#include <cmath>
#include <cstdlib>

#include "liatools.h"
#include "ComputeNorm.h"


// JF Bonastre 4/11/2009 - cleaning and add dichotomic search in Norm
// JF Bonastre 27/4/2010 -cleaning and add 
//                          * impostor (target independent) selection
//                          * High score discard
//							* Low score discard
//                          * Median computation instead of aryth mean
                          

using namespace alize;
using namespace std;


class DistribNorm{
	LKVector tabScore; // table of norm scores
	/* Use of LKVector ALIZE object, it provides sort functions */
	String *tabIdImp; // table of norm scores
	Config conf;
	public:
	void init();
	void addScore(double, const String&);
	bool computeMeanStd(double &, double&, char,double,double);
	unsigned long size();
	void print();
	DistribNorm(Config &);  
	~DistribNorm();
};

void DistribNorm::init(){
 	tabScore.clear(); 
}

unsigned long DistribNorm::size(){
	return tabScore.size();
}

void DistribNorm::addScore(double d, const String& idImp){
	LKVector::type t;
 	unsigned long ind;
	/* idx value allows to keep in memory the impostor name related to each score, even after applying the sort function */
        /* After sorting, the name of an impostor is retrieved thanks to the tabIdImp[idx] */
        
	if(tabScore.size() == 0)
		ind = (unsigned long)0;
	else
		ind = (tabScore.size());
	t.lk = d;
	t.idx = ind;
	tabIdImp[ind] = idImp;
	tabScore.addValue(t);
}


bool DistribNorm::computeMeanStd(double &mean,double &std,char mode, double percentH,double percentL){
	if(tabScore.size() == 0) return false;
	else{
		unsigned long size=tabScore.size();
		unsigned long begin=0;
		unsigned long end=size;
		if ((percentH!=0) || (percentL!=0)){
			tabScore.descendingSort();
			unsigned long discardH=(unsigned long) (((double) size) * percentH);
			unsigned long discardL=(unsigned long) (((double) size) * percentL);
			size-=(discardH+discardL);
			begin=discardH;
			end=tabScore.size()-discardL;
			if (debug) cout<<"discardH["<<discardH<<"] discardL["<<discardL<<"] begin["<<begin<<"] end["<<end<<"]"<<endl;   
		   }
		double sum=0;
		double sum2=0;
		switch (mode){
			case 0: //classical mean computation	
					for(unsigned int i=begin; i<end; i++){
						sum += tabScore[i].lk;  
						sum2 += tabScore[i].lk * tabScore[i].lk;
						}
					mean=sum/(double)(size);
					std=sqrt((sum2/(double)(size)-(mean*mean)));  
					break;
			case 1: // Median instead of mean
			 		mean=tabScore[begin+(size/2)].lk;
			 		for(unsigned int i=begin; i<end; i++)
						sum += abs(tabScore[i].lk-mean); 
			 		std=sum/(double)(size);
			 		break;
			default: cout <<"mean compute mode["<<mode<<"] unknown"<<endl;
					 exit(0);
		}  				
	
	    return true;
	}
}

/*
double DistribNorm::std(){
	double m;	

	if(tabScore.size() == 0){
		cout << "No normalization score available for std computation: std=1.0" << endl;
		return 1.0;	 
	}
	m = mean();
	return 
}*/


void DistribNorm::print(){
	cout << "Distrib Score Print (nbScore=" << tabScore.size() << ")" << endl;
 	for(unsigned int i=0; i<tabScore.size(); i++){
		cout << tabIdImp[tabScore[i].idx] << " ";
	 	printf("%.8lE ", tabScore[i].lk);
	}
	cout << endl;
}

DistribNorm::DistribNorm(Config &config){
	conf = config;
	unsigned long nbScoreMax = (unsigned long)(conf.getParam("maxScoreDistribNb").toLong());
	tabIdImp = new String[nbScoreMax];	
//	meanVal = 0.0;
//	stdVal = 1.0;
}

DistribNorm::~DistribNorm() {
 	delete []tabIdImp;
}

class NormSeg{
	String name;
 	double mean;
	double std;
	bool computed;
//	double sum;
//	double sum2;
	unsigned long nb;
	DistribNorm *distribNorm; // A pointer on an imp score distrib used once for computing the norm parameters
	public:
	double getMean(char,double,double);
	double getStd(char,double,double);
	unsigned long getNb();
	String getName();
	void setNb(unsigned long);
	void set(String,double,double); 
	void setName(String); 
	void addScore(double,const String &imp);
        void newDistribNorm(Config &);
	void deleteDistribNorm();
	NormSeg();
	~NormSeg();
};
void NormSeg::addScore(double score,const String & imp){
	distribNorm->addScore(score,imp);
}
double NormSeg::getMean(char computeMode,double discardH,double discardL){
	if (!computed){
			if (!distribNorm->computeMeanStd(mean,std,computeMode,discardH,discardL)){
			cout << "Problem: empty impostor cohort"<<endl;
			exit(0);
		    }
		computed=true;
		}
	return mean; 
}

double NormSeg::getStd(char computeMode,double discardH,double discardL){
	if (!computed){
		if (!distribNorm->computeMeanStd(mean,std,computeMode,discardH,discardL)){
			cout << "Problem: empty impostor cohort"<<endl;
			exit(0);
		    }
		computed=true;
		}
	return std;	
}

unsigned long NormSeg::getNb(){
	return nb; 
}

String NormSeg::getName(){
	return name; 
}


void NormSeg::setNb(unsigned long m){
	nb = m; 
}

void NormSeg::set(String n,double m,double s){
	name = n; 
	mean = m; 
	std = s;
	computed=true; 
}
void NormSeg::setName(String n){
	name = n; 
}
void NormSeg::newDistribNorm(Config &config){
	distribNorm=new DistribNorm(config);
}
void NormSeg::deleteDistribNorm(){
	delete distribNorm;
	distribNorm=NULL;
}
NormSeg::NormSeg(){
	 distribNorm=NULL; // initialize at NULL the distribNorm, initialized independently of NormSeg but freezed automatically
}
NormSeg::~NormSeg(){ // Clean, if needed, the DistribNorm
	if (distribNorm) deleteDistribNorm();
}

class Norm {
 	NormSeg *normTab;	 // table of normSeg => seg (or Id for znorm) name associated with mean and std values
	unsigned long nbNorm;    // number of normSeg in tabScore
	unsigned long nbNormMax; // Max number of normSeg in tabScore
	unsigned long findEntityIdxInNorm(String);             // Find where to add the seg info
        void moveEntity(unsigned long,unsigned long); // Make the room for the new Entity
	unsigned long idxLastFind;             // Save the idx of the last segment search
	char _computeMode;             // The mode for mean/std computation
								  // 0=classical
								  // 1=MedianBased							  
    double _percentH;        // % of Highest scores discarded
    double _percentL;        // % of Lowest scores discarded						  
	public:
	bool findEntityInNorm(String,unsigned long &);   
	double getMean(unsigned long);
	double getStd(unsigned long);
	//double getSum(unsigned long);
	//double getSum2(unsigned long);
	unsigned long getNb(unsigned long);
	void print();
	unsigned long addSeg(String, double, double);
//	unsigned long addSeg(String, double, double, unsigned long);
	unsigned long findAddSeg(String,Config&); // TAKE CARE, DO NOT INITIALIZE MEAN AND COV
	void addScore(unsigned long,double,const String &);
	void deleteAllDistribNorm(); // Clean the memory (done automatically - in ~NormSeg() - but could be asked manually)
	void setNormAndFreeDistrib();// Compute Mean and Cov and free the DistribNorm (per seg)
//	void setNorm();              // Same, but not clean the DistribNorm
	Norm(unsigned long);         // Create a Norm object, with x NormSeg and Classical mean/std computation
	Norm(unsigned long, char,double,double); // Id but with a defined type of computation for mean/std
	~Norm();
};
void Norm::setNormAndFreeDistrib(){
	for(unsigned long idx=0; idx<nbNorm; idx++){
		normTab[idx].getMean(_computeMode,_percentH,_percentL);
  		normTab[idx].deleteDistribNorm();	
		}
}

void Norm::addScore(unsigned long idx,double score,const String & imp)
{		
	normTab[idx].addScore(score,imp);
}
void Norm::deleteAllDistribNorm(){		// clean the distribNorm to save space
	for(unsigned long i=0; i<nbNorm; i++)
		normTab[i].deleteDistribNorm();	
}
// Make the room for the new element
// Could be time consumming. good to sort the score file (.res) before
void Norm::moveEntity(unsigned long beg,unsigned long end){ 
	
        if (debug || (verboseLevel>2)) 
		cout << "moveEntity beg["<<beg<<"] end["<<end<<"]"<<endl;
	for (unsigned idx=end;idx>beg;idx--)
		normTab[idx]=normTab[idx-1];
	//memmove((void*)&(normTab[beg+1]),(void*)&(normTab[beg]),sizeof(NormSeg)*((beg-end)+1));
}
unsigned long Norm::findEntityIdxInNorm(String name){
	unsigned long beg=0;	
	unsigned long end=nbNorm-1;
	unsigned long idx=0;	
  	if (nbNorm==0) return 0; // empty...
	while ((end-beg)>1){
		idx=(beg+end)/2;
		if (name > normTab[idx].getName())
			beg=idx;
			else end=idx;
		}
	if (name <=normTab[beg].getName()) return beg;
		else if (name>normTab[end].getName()) return end+1;
			else return end; 
}
bool Norm::findEntityInNorm(String name,unsigned long & idx){
	if(nbNorm == 0) {idx=0;return false;}
	// last one searched?
	if (normTab[idxLastFind].getName() == name){
		idx=idxLastFind;
	 	return true;		
	}
        idx=findEntityIdxInNorm(name);
	if(normTab[idx].getName() == name){
		idxLastFind=idx; 
		return true;
		}
	else return false; 
}

void Norm::print(){
	cout << "nbnorm: " << nbNorm << endl; 
 	for(unsigned long i=0; i<nbNorm; i++){
		cout << normTab[i].getName() << " " << getMean(i) << " " << getStd(i) << endl; 
	} 
}

double Norm::getMean(unsigned long ind){
 	return normTab[ind].getMean(_computeMode,_percentH,_percentL);
}

double Norm::getStd(unsigned long ind){
 	return normTab[ind].getStd(_computeMode,_percentH,_percentL);
}

unsigned long Norm::getNb(unsigned long ind){
 	return normTab[ind].getNb();
}

unsigned long Norm::findAddSeg(String name,Config & config){ // TAKE CARE, JUST ADD THE SEG, DO NOT INITIALIZE MEAN/COV
        if(nbNorm == (nbNormMax)){
		cout << "Table norm Capacity out of bound: " << nbNorm << endl; 
		exit(-1);
	}
	unsigned long idx;
	if (!findEntityInNorm(name,idx)){       // the seg is missing
		if (debug) cout <<"addSeg:seg["<<name<<"] idx["<<idx<<"]"<<endl;
		moveEntity(idx,nbNorm);                      // Make the room for the new Entity
		normTab[idx].setName(name);
		nbNorm++;
		normTab[idx].newDistribNorm(config);
	}
    return idx;
}

unsigned long Norm::addSeg(String name, double mean, double std){
        if(nbNorm == (nbNormMax)){
		cout << "Table norm Capacity out of bound: " << nbNorm << endl; 
		exit(-1);
	}
	unsigned long idx=findEntityIdxInNorm(name);   // Find where to add the seg info
	if (debug) cout <<"addSeg:seg["<<name<<"] idx["<<idx<<"]"<<endl;
	moveEntity(idx,nbNorm);                      // Make the room for the new Entity
	normTab[idx].set(name,mean,std); 
	nbNorm++;	
    return (idx);
}

Norm::Norm(unsigned long nbMax){
	 _computeMode=0;
	 _percentH=0;_percentL=0;
	 nbNorm = 0;
	 nbNormMax = nbMax;
	 normTab = new NormSeg[nbNormMax];
	 idxLastFind=0;
}
Norm::Norm(unsigned long nbMax, char compute,double percentH,double percentL){
	 _computeMode=compute;
	 _percentH=percentH;
	 _percentL=percentL;
	 nbNorm = 0;
	 nbNormMax = nbMax;
	 normTab = new NormSeg[nbNormMax];
	 idxLastFind=0;
}
Norm::~Norm(){
  	delete [] normTab; 
}

// Get the imp scores, for one seg or once for all
//-------------------------------------------------------------------------------------------------------
bool selectImp(String &fieldOne,String& fieldTwo,XList &impList,char selectMode)
{
	bool select;
	switch (selectMode){
		case 0:select=true;break; // No selection
		case 1:select=(impList.findLine(fieldTwo) !=NULL);break; // target independent selection 
		default:select=true;cout <<"selectMode unknown. No impostor scores selection is applied"<<endl;
		}     
	return select;		       	     
}																
void getAllScores(XList &list, unsigned long indFieldOne, unsigned long indFieldTwo, unsigned long indScore, Norm &norm,XList &impList,char selectMode,Config &config){
	XLine *linep;
	// ind* are used to locate information in NIST File in the case of format change
	list.getLine(0);
	while ((linep=list.getLine()) != NULL){
		String fieldOne=linep->getElement(indFieldOne);
		String fieldTwo=linep->getElement(indFieldTwo);
		double score=linep->getElement(indScore).toDouble();		      	 
		if (selectImp(fieldOne,fieldTwo,impList,selectMode)){ 
		   // If no selection of impostor scores is needed or if there is a selection and the impostor is in the list 
		   // Find the index in norm or add a seg/id in norm if needed, return the idx of the seg
		   // And initialize (memory ans others) the corresponding DistribNorm if needed
           unsigned long idx=norm.findAddSeg(fieldOne,config);
		   norm.addScore(idx,score, fieldTwo);
		   }
		}	
}	
		
// getAllScoresFirstNormed computes the score distributions for t or znorm AFTER applying firstNorm	normalization
// used for ztnorm and tznorm												
void getAllScoresFirstNormed(XList &list, unsigned long indFieldOne, unsigned long indFieldTwo,unsigned long indScore, Norm &firstNorm, 
 			Norm &outNorm,XList &impList,char selectMode,Config &config){
	XLine *linep;
	// TODO Replace 0 and 3 by the indElt, indIdImp
	list.getLine(0);
	while ((linep=list.getLine()) != NULL){
		String fieldOne=linep->getElement(indFieldOne);
		String fieldTwo=linep->getElement(indFieldTwo);
		double score=linep->getElement(indScore).toDouble();	 	 
		unsigned long ind; 		
		if (selectImp(fieldOne,fieldTwo,impList,selectMode)){ 
		   	// If no selection of impostor scores is needed or if there is a selection and the impostor is in the list 
			if(firstNorm.findEntityInNorm(fieldTwo,ind)){
			    unsigned long idx=outNorm.findAddSeg(fieldOne,config);
				outNorm.addScore(idx,(score - firstNorm.getMean(ind)) / firstNorm.getStd(ind), fieldTwo);
		        }
			else{
				cout << "distribution for["<<fieldTwo<<"] not found!" << endl;
				cout <<"line:"<<linep<<endl;
				exit(-1); 
				}
		}
	}
}		
//-------------------------------------------------------------------------------------------------------
int ComputeNorm(Config& config)
{

	using namespace alize;
	using namespace std;

	int decision=0; // TODO change ASAP
    String outputNISTFileName = config.getParam("outputFileBaseName");                        // Result file BAsenamein NIST style (.nist) result file format
    String znormFilesExtension = config.getParam("znormFilesExtension");                      // Result file znorm extension
    String tnormFilesExtension = config.getParam("tnormFilesExtension");                      // Result file tnorm extension
    String tznormFilesExtension = config.getParam("tznormFilesExtension");                    // Result file tztnorm extension
    String ztnormFilesExtension = config.getParam("ztnormFilesExtension");                    // Result file tztnorm extension
  
	// Norm Type = znorm, tnorm ou ztnorm
	unsigned long maxIdNb = (unsigned long)(config.getParam("maxIdNb").toLong());
	unsigned long maxSegNb = (unsigned long)(config.getParam("maxSegNb").toLong());
	String normType = config.getParam("normType");
	String testNistFile = config.getParam("testNistFile"); 
    char selectMode=0;       // 0 if no selection, 1 if target independent selection
    XList impList;
    if (config.existsParam("impostorIDList")){
    	selectMode=1; // Target independent Impostor score selection mode  	
    	impList.load(config.getParam("impostorIDList"),config);
    }  
	char computeMode = (char)(config.getParam("meanMode").toLong()); // 0 classical, 1 Median
	double percentH=(double)(config.getParam("percentH").toDouble()); // % of hihest scores discarded
	double percentL=(double)(config.getParam("percentL").toDouble());// % of lowest scores discarded
	
	// Define the field of the score files
	char fieldGender=(char) config.getParam("fieldGender").toLong();
	char fieldName=(char)config.getParam("fieldName").toLong();
	//char fieldDecision=(char)config.getParam("fieldDecision").toLong(); // not used
	char fieldSeg=(char)config.getParam("fieldSeg").toLong();
	char fieldLLR=(char)config.getParam("fieldLLR").toLong();
	
	XList testList(testNistFile, config);
	testList.getLine(0); // goto first line still to debug
	
	try{		
	if(normType == "tnorm"){	
		String outputFilename=outputNISTFileName+tnormFilesExtension;	   // Create the complete output file filename
		ofstream outFile(outputFilename.c_str(),ios::out | ios::trunc);    // Initialise the output file
	 	Norm tnorm(maxSegNb,computeMode,percentH,percentL);// storage of mean and std for all segments
		String tnormNistFile = config.getParam("tnormNistFile"); 
		if (verbose) cout << "Tnorm, reading tnormList"<<endl;
		XList tnormList(tnormNistFile, config);
		getAllScores(tnormList, fieldSeg, fieldName, fieldLLR, tnorm,impList,selectMode,config);
		tnorm.setNormAndFreeDistrib(); // finalize the mean/cov computation and free the score (the per seg DistribNorm)

		if (verbose) cout << "Tnorm, begin the test list"<<endl;
	 	XLine * linep;
		while ((linep=testList.getLine()) != NULL){
			// The test segment filename
			String seg=linep->getElement(fieldSeg);
			unsigned long ind;
			if(debug) cout << endl << "findEntityInNorm SEG[" << seg << "]"; 				
			if (!tnorm.findEntityInNorm(seg,ind)){// Problem, the seg parameters are not present...
				cout << "tnorm impostor score not found for seg["<<seg<<endl;
				exit(-1); 
				}
			if(debug) cout << " idx["<< ind<< "]"<<endl; 										
			double score = ((linep->getElement(fieldLLR)).toDouble() - tnorm.getMean(ind)) / tnorm.getStd(ind);
	        outputResultLine(score, linep->getElement(fieldName), seg, linep->getElement(fieldGender), decision, outFile);	
		}
		outFile.close();    			
	    if(debug||verbose){
			cout << "tnorm distrib" << endl;
	       	tnorm.print();       	 
		}
		      
	}
    else if(normType == "znorm"){	
    	String outputFilename=outputNISTFileName+znormFilesExtension;	   // Create the complete output file filename
		ofstream outFile(outputFilename.c_str(),ios::out | ios::trunc);    // Initialise the output file
		if(verbose) cout << endl << "Compute Znorm" << endl; 
		// storage of mean and std for all segments
	 	Norm znorm(maxIdNb,computeMode,percentH,percentL);
		String znormNistFile = config.getParam("znormNistFile"); 	
		XList znormList(znormNistFile, config);
	 	XLine * linep;		
	 	if (verbose) cout << "Computing once the Znorm distributions"<<endl<<
							 "Take Care, really faster if the score files are sorted!"<<endl;
		getAllScores(znormList, fieldName, fieldSeg, fieldLLR, znorm,impList,selectMode,config);
		znorm.setNormAndFreeDistrib(); // finalize the mean/cov computation and free the score (the per seg DistribNorm)
		
		while ((linep=testList.getLine()) != NULL){
			// The test segment filename
			String id=linep->getElement(fieldName);
			unsigned long ind; 	
			if(verboseLevel>2) cout << endl << "Compute Znorm [" << id << "]"<< endl; 							
			if(debug) cout << endl << "findEntityInNorm Id[" << id << "]"; 				
			if (!znorm.findEntityInNorm(id,ind)){// Problem, the seg parameters are not present...
				cout << "znorm impostor score not found for id["<<id<<endl;
				exit(-1); 
			    }
			if(debug) cout << " idx["<< ind<< "]"<<endl; 						
			double score = ((linep->getElement(fieldLLR)).toDouble() - znorm.getMean(ind)) / znorm.getStd(ind);
	        outputResultLine(score, id, linep->getElement(fieldSeg), linep->getElement(fieldGender), decision, outFile);
		    }	 
		outFile.close();		    
	    if(debug||verbose){
			cout << "znorm distrib" << endl;
	       	znorm.print();       	 
		    }	
	    }
     else if(normType == "ztnorm"){ // Ztnorm outputs ztnorm AND tnorm results
     	String outputTnormFilename=outputNISTFileName+tnormFilesExtension;	    // Create the complete output file filename
     	String outputZtnormFilename=outputNISTFileName+ztnormFilesExtension;	// Create the complete output file filename
     	
		ofstream outFileTnorm(outputTnormFilename.c_str(),ios::out | ios::trunc);      // Initialise the output file
		ofstream outFileZtnorm(outputZtnormFilename.c_str(),ios::out | ios::trunc);    // Initialise the output file
		if(verbose) cout << endl << "Compute ZTnorm (and tnorm)" << endl; 
		// storage of mean and std for all segments
	 	Norm tnorm(maxSegNb,computeMode,percentH,percentL);
	 	Norm znorm(maxIdNb,computeMode,percentH,percentL);
		Norm ztnorm(maxSegNb,computeMode,percentH,percentL);
		
		String tnormNistFile = config.getParam("tnormNistFile"); 
		String znormNistFile = config.getParam("znormNistFile"); 
		String ztnormNistFile = config.getParam("ztnormNistFile"); 
		XList tnormList(tnormNistFile, config);
		XList znormList(znormNistFile, config);
		XList ztnormList(ztnormNistFile, config);
				
		// tnorm distribution for the impostor segments (from imp_imp) in ztnorm		
		if (verbose) cout << "Computing once the Tnorm distribution for impostor segments"<<endl<<
							 "Take Care, really faster if the score files are sorted!"<<endl;
		getAllScores(ztnormList, fieldSeg, fieldName, fieldLLR, ztnorm,impList,selectMode,config);
		ztnorm.setNormAndFreeDistrib(); // finalize the mean/cov computation and free the score (the per seg DistribNorm)
		// compute the tnorm distributions for the (test) seg in tnorm
	    if (verbose) cout << "Computing once all the tnorm distributions for the test segments "<<endl
					      << "Take Care, really faster if the score files are sorted!"<<endl;
		getAllScores(tnormList, fieldSeg, fieldName, fieldLLR, tnorm,impList,selectMode,config);
		tnorm.setNormAndFreeDistrib(); // finalize the mean/cov computation and free the score (the per seg DistribNorm)
		
		// compute the znorm distrib for the tnormed scores in znorm
		if (verbose) cout << "Computing once all the znorm (tnormed) distributions "<<endl
					      << "Take Care, really faster if the score files are sorted!"<<endl;
		getAllScoresFirstNormed(znormList, fieldName, fieldSeg, fieldLLR, ztnorm, znorm,impList,selectMode,config);
		znorm.setNormAndFreeDistrib();// finalize the mean/cov computation and free the score (the per seg DistribNorm)
	
		// for all the test lines (one test seg against several target speakers)
		XLine * linep;
		while ((linep=testList.getLine()) != NULL){
			// The test segment filename
			String seg=linep->getElement(fieldSeg);    // the test seg 
			String id= linep->getElement(fieldName);   // the Id (target spk)
			if (verboseLevel>2) cout << "ztnorm seg["<<seg<<"] Id["<<id<<"]"<<endl;
			unsigned long ind; 		
			// segment tnorm distribution already computed ?
			if(debug) cout << endl << "findEntityInNorm SEG[" << seg << "]"; 				
			if (!tnorm.findEntityInNorm(seg,ind)){// Problem, the seg parameters are not present...
				cout << "tnorm impostor distribution not found for seg["<<seg<<"]"<<endl;
				exit(-1); 
				}
			if(debug) cout << " idx["<< ind<< "]"<<endl; 	
			double scoreTNorm = ((linep->getElement(fieldLLR)).toDouble() - tnorm.getMean(ind)) / tnorm.getStd(ind);

			// ID mean and std already computed ?
			if(!znorm.findEntityInNorm(id,ind)){// Problem, the seg parameters are not present...
				cout << "znorm(tnormed) distribution not found for id["<<id<<"]"<<endl;
				exit(-1); 
				}
			double score = (scoreTNorm - znorm.getMean(ind)) / znorm.getStd(ind);
			
		    outputResultLine(score, linep->getElement(fieldName), seg, linep->getElement(fieldGender), decision, outFileZtnorm);
		    outputResultLine(scoreTNorm, linep->getElement(fieldName), seg, linep->getElement(fieldGender), decision, outFileTnorm);
		}	      
		outFileZtnorm.close();
		outFileTnorm.close();
	    if(debug||verbose){
		    cout << "tnorm distrib" << endl;
		   	tnorm.print();       	 
			cout << "znorm distrib" << endl;
		   	znorm.print();       	 
		}
	}     
	   else if(normType == "tznorm"){ // tznorm outputs tznorm and znorm scores
     	String outputTznormFilename=outputNISTFileName+tznormFilesExtension;	// Create the complete output file filename  	
     	String outputZnormFilename=outputNISTFileName+znormFilesExtension;	    // Create the complete output file filename  	
		ofstream outFileTznorm(outputTznormFilename.c_str(),ios::out | ios::trunc);    // Initialise the output file
		ofstream outFileZnorm(outputZnormFilename.c_str(),ios::out | ios::trunc);    // Initialise the output file
		if(verbose) cout << endl << "Compute TZnorm" << endl; 
		// storage of mean and std for all segments
	 	Norm tnorm(maxSegNb,computeMode,percentH,percentL);
	 	Norm znorm(maxIdNb,computeMode,percentH,percentL);
		Norm tznorm(maxSegNb,computeMode,percentH,percentL);
		
		String tnormNistFile = config.getParam("tnormNistFile"); 
		String znormNistFile = config.getParam("znormNistFile"); 
		String tznormNistFile = config.getParam("ztnormNistFile"); 
		XList tnormList(tnormNistFile, config);
		XList znormList(znormNistFile, config);
		XList tznormList(tznormNistFile, config);
		
		// compute znorm distribs for the target speakers ID in znorm 
		if (verbose) cout << "Compute Znorm distribs for the target speakers (using znorm scores)"<<endl
						  <<"Computing once all the impostor distribution "<<endl
					      << "Take Care, really faster if the score files are sorted!"<<endl;
		getAllScores(znormList, fieldName, fieldSeg, fieldLLR, znorm,impList,selectMode,config);
		znorm.setNormAndFreeDistrib(); // finalize the mean/cov computation and free the score (the per seg DistribNorm)		
		
		// compute znorm distribs for the impostor ID in tznorm 
		if (verbose) cout << "Compute znorm distribs for the impostor speakers (using ztnorm scores)"
						  << "Computing once all the impostor distribution "<<endl
					      << "Take Care, really faster if the score files are sorted!"<<endl;
		getAllScores(tznormList, fieldName, fieldSeg, fieldLLR, tznorm,impList,selectMode,config);
		tznorm.setNormAndFreeDistrib(); // finalize the mean/cov computation and free the score (the per seg DistribNorm)		
	
		// compute tnorm distribs for the znormed scores
		if (verbose) cout << "Compute tnorm distribs of the znormed scores "<<endl
						  << "Computing once all the impostor distribution "<<endl
					      << "Take Care, really faster if the score files are sorted!"<<endl;
		getAllScoresFirstNormed(tnormList, fieldSeg, fieldName, fieldLLR, tznorm, tnorm,impList,selectMode,config);
		tnorm.setNormAndFreeDistrib();// finalize the mean/cov computation and free the score (the per seg DistribNorm)
		if(debug|| verbose){
		    cout << "tnorm distrib" << endl;
		   	tnorm.print();       	 
			cout << "znorm distrib" << endl;
		   	znorm.print();    
		    cout << "tznorm distrib" << endl;
		   	tznorm.print();       	   	 
		}		
		// for all the test lines (one test seg against several target speakers)
		if (verbose) cout << "begin the test list"<<endl;
		XLine * linep;
		while ((linep=testList.getLine()) != NULL){
			// The test segment filename
			String seg=linep->getElement(fieldSeg);    // the test seg 
			String id= linep->getElement(fieldName);   // the Id (target spk)
			if (verboseLevel>2) cout << "tznorm seg["<<seg<<"] Id["<<id<<"]"<<endl;
			unsigned long ind; 		
			// id mean and std already computed ?
			if(debug) cout << endl << "findEntityInNorm SEG[" << seg << "]"; 				
			if (!znorm.findEntityInNorm(id,ind)){// Problem, the id parameters are not present...
				cout << "znorm distribution not found for id["<<id<<"]"<<endl;
				exit(-1); 
				}
			if(debug) cout << " idx["<< ind<< "]"<<endl; 	
			double scoreZNorm = ((linep->getElement(fieldLLR)).toDouble() - znorm.getMean(ind)) / znorm.getStd(ind);
			// segment mean and std already computed ?
			if(!tnorm.findEntityInNorm(seg,ind)){// Problem, the seg parameters are not present...
				cout << "tnorm(znormed) distribution not found for seg["<<seg<<"]"<<endl;
				exit(-1); 
				}
			double score = (scoreZNorm - tnorm.getMean(ind)) / tnorm.getStd(ind);
			
		    outputResultLine(score, linep->getElement(fieldName), seg, linep->getElement(fieldGender), decision, outFileTznorm);
		    outputResultLine(scoreZNorm, linep->getElement(fieldName), seg, linep->getElement(fieldGender), decision, outFileZnorm);
		    if (verboseLevel>2) outputResultLine(score, linep->getElement(fieldName), seg, linep->getElement(fieldGender), decision, cout);

		}	      
		outFileTznorm.close();
		outFileZnorm.close();
	    if(debug || verbose){
		    cout << "tnorm distrib" << endl;
		   	tnorm.print();       	 
			cout << "znorm distrib" << endl;
		   	znorm.print();       	 
		}
	}   
	else {
		  cout << "unknown normalization mode:"<<normType<<endl;
		  exit(-1); 
		  }   
	} //try
	catch (Exception& e){ 
		cout << e.toString().c_str() << endl;
	}

              
    return 0;
}

#endif //!defined(ALIZE_ComputeNorm_cpp)

#if !defined(ALIZE_ComputeNorm_cpp)
#define ALIZE_ComputeNorm_cpp

#include <iostream>
#include <fstream>  
#include <cstdio> 
#include <cassert> 
#include <cmath>
#include <cstdlib>
#include <liatools.h>
#include "ComputeNorm.h"

using namespace alize;
using namespace std;

class NormSeg{
	String name;
 	double mean;
	double std;
	double sum;
	double sum2;
	unsigned long nb;
	Histo histo;
	
	public:
	double getMean();
	double getStd();
	Histo getHisto();
	double getSum();
	double getSum2();
	unsigned long getNb();
	String getName();
	void setMean(double);
	void setStd(double);
	void setHisto(Histo);
	void setSum(double);
	void setSum2(double);
	void setNb(unsigned long);
	void setName(String); 
};

double NormSeg::getMean(){
	return mean; 
}

double NormSeg::getStd(){
	return std; 
}

double NormSeg::getSum(){
	return sum; 
}

double NormSeg::getSum2(){
	return sum2; 
}
unsigned long NormSeg::getNb(){
	return nb; 
}

Histo NormSeg::getHisto(){
	return histo; 
}

String NormSeg::getName(){
	return name; 
}

void NormSeg::setMean(double m){
	mean = m; 
}

void NormSeg::setStd(double s){
	std = s; 
}

void NormSeg::setSum(double m){
	sum = m; 
}

void NormSeg::setSum2(double m){
	sum2 = m; 
}

void NormSeg::setNb(unsigned long m){
	nb = m; 
}

void NormSeg::setHisto(Histo h){
	histo = h; 
}

void NormSeg::setName(String n){
	name = n; 
}


class Norm {

 	NormSeg *normTab;	 // table of normSeg => seg name associated with mean and std values
	unsigned long nbNorm;    // number of normSeg in tabScore
	unsigned long nbNormMax; // Max number of normSeg in tabScore
	
	public:
	long findEntityInNorm(String);   
	double getMean(unsigned long);
	double getStd(unsigned long);
	double getSum(unsigned long);
	double getSum2(unsigned long);
	unsigned long getNb(unsigned long);
	Histo getHisto(unsigned long);	
	void print();
	unsigned long addSeg(String, double, double);
	unsigned long addSeg(String, double, double, unsigned long);
	unsigned long addSeg(String, Histo);
	Norm(unsigned long);
	~Norm();
};

long Norm::findEntityInNorm(String name){
	unsigned long i=0;

	// last one searched?
	if((nbNorm != 0) && (normTab[nbNorm-1].getName() == name))
	 	return nbNorm-1;
	
	while ((i<nbNorm) && (normTab[i].getName() != name)) i++;
	if(i == nbNorm) 
	 	return -1;
	else
	 	return i; 
}

void Norm::print(){
	cout << "nbnorm: " << nbNorm << endl; 
 	for(unsigned long i=0; i<nbNorm; i++){
		cout << normTab[i].getName() << " " << normTab[i].getMean() << " " << normTab[i].getStd() << endl; 
	} 
}

double Norm::getMean(unsigned long ind){
 	return normTab[ind].getMean();
}

double Norm::getStd(unsigned long ind){
 	return normTab[ind].getStd();
}

double Norm::getSum(unsigned long ind){
 	return normTab[ind].getSum();
}

double Norm::getSum2(unsigned long ind){
 	return normTab[ind].getSum2();
}

unsigned long Norm::getNb(unsigned long ind){
 	return normTab[ind].getNb();
}


Histo Norm::getHisto(unsigned long ind){
 	return normTab[ind].getHisto();
}

unsigned long Norm::addSeg(String name, double mean, double std){
	normTab[nbNorm].setName(name);
	normTab[nbNorm].setMean(mean);
	normTab[nbNorm].setStd(std); 
	nbNorm++;
	if(nbNorm == nbNormMax){
		cout << "Table norm Capacity out of bound: " << nbNorm << endl; 
		exit(-1);
	}
       return (nbNorm-1);
}

unsigned long Norm::addSeg(String name, double sum, double sum2, unsigned long nb){
	normTab[nbNorm].setName(name);
	normTab[nbNorm].setSum(sum);
	normTab[nbNorm].setSum2(sum2); 
	normTab[nbNorm].setNb(nb); 
	
	nbNorm++;
	if(nbNorm == nbNormMax){
		cout << "Table norm Capacity out of bound: " << nbNorm << endl; 
		exit(-1);
	}
       return (nbNorm-1);
}


unsigned long Norm::addSeg(String name, Histo h){
	normTab[nbNorm].setName(name);
	normTab[nbNorm].setHisto(h);
	nbNorm++;
	if(nbNorm == nbNormMax){
		cout << "Table norm Capacity out of bound: " << nbNorm << endl; 
		exit(-1);
	}
       return (nbNorm-1);
}


Norm::Norm(unsigned long nbMax){
	 nbNorm = 0;
	 nbNormMax = nbMax;
	 normTab = new NormSeg[nbNormMax];
}

Norm::~Norm(){
  	delete [] normTab; 
}

class DistribNorm{
	LKVector tabScore; // table of norm scores
	/* Use of LKVector ALIZE object, provided sort functions */
	String *tabIdImp; // table of norm scores
	double meanVal;   // storage of mean value to avoid recomputation
	double stdVal;    // storage of std value to avoid recomputation
	Config conf;
	bool debug;
	bool verbose;
		
	public:
	void init();
	void addScore(double, String&);
	void selectScore(XLine *, String);
	long findImp(String);
	double mean();
	double meanVrais();
	double std();
	double sum();
	double sum2(); 
	unsigned long size();
	Histo histo();	
	void print();
	static int compare(const void*, const void*);
	void sortScore();
	DistribNorm(Config &);  
	~DistribNorm();
};

void DistribNorm::init(){
 	tabScore.clear();
	meanVal = 0.0;
	stdVal = 1.0;	 
}

unsigned long DistribNorm::size(){
	return tabScore.size();
}

void DistribNorm::addScore(double d, String& idImp){
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


void DistribNorm::selectScore(XLine *l, String selectType){
	
	if(selectType == "noSelect"){		
	 	return;	 
	}	 
	
        if(selectType == "selectNBestByTestSegment"){
		// Sort of impostor scores	
	 	if(debug){
			cout << "Sort the score" << endl;		 
		}
	 	sortScore();		 
	 	
		// Size of impostor cohort
		String sMax = conf.getParam("cohortNb");
		unsigned long maxCohortNb = (unsigned long)sMax.toLong();	  	
	 	
		if(debug){
			cout << "Size of the cohort " << maxCohortNb <<endl;		 
		}
		tabScore.setSize(maxCohortNb);
		
		if(debug){
			cout << "Score distribution after selection" << endl;
		 	print();
		}
		return;
	}
       if(selectType == "selectTargetDependentCohortInFile"){
		/* Scores are selected according to an impostor name list, which is dependent on target model */
	 	if(debug){
			cout << "Retrieve client dependent cohort score" << endl;		 
		}
		String impFile = conf.getParam("cohortFilePath")+l->getElement(1)+conf.getParam("cohortFileExt");
		cout << "ouverture file : " << impFile << endl;
		XList impFileList(impFile, conf);
	 	XLine * linep;
		unsigned long impNb=0;
		while ((linep=impFileList.getLine()) != NULL){
			String impName = linep->getElement(0);
			/* necessary to find where is stored the informationof this impostor */
			long indImp = findImp(impName);			
			
		        if(indImp == -1){
				cout << "dependent Impostor name does not correspond to any impostor" << endl; 
				exit(0);
			 }
			else{
			 	/* permutation of structures in tabScore */
				LKVector::type t = tabScore[impNb];
				tabScore[impNb] = tabScore[indImp];
				tabScore[indImp] = t;
				impNb++;				
			}
		}
	        tabScore.setSize(impNb);
		if(debug){
			cout << "Score distribution after Cohort selection" << endl;
		 	print();
		}
	       	
	}	
}

long DistribNorm::findImp(String impName){
	unsigned int i=0;
	
	while((i < tabScore.size()) && (tabIdImp[tabScore[i].idx] != impName)) i++;
	
	if(i == tabScore.size()) return -1;
	return i;	
}

double DistribNorm::sum(){
 	double sum=0.0;
	
	if(tabScore.size() == 0){
		cout << "No normalization score available for mean computation: mean=0.0" << endl;
		return 0.0;
	}

	for(unsigned int i=0; i<tabScore.size(); i++)
		sum += tabScore[i].lk; 
       	
        return sum;
}


double DistribNorm::sum2(){
 	double sum2=0.0;
	
	if(tabScore.size() == 0){
		cout << "No normalization score available for mean computation: mean=0.0" << endl;
		return 0.0;
	}

	for(unsigned int i=0; i<tabScore.size(); i++)
		sum2 += tabScore[i].lk * tabScore[i].lk; 
       	
        return sum2;
}

double DistribNorm::mean(){
 	double sum=0.0;
	
	if(meanVal != 0.0)
		return meanVal; 


	if(tabScore.size() == 0){
		cout << "No normalization score available for mean computation: mean=0.0" << endl;
		return 0.0;
	}

	for(unsigned long i=0; i<tabScore.size(); i++)
		sum += tabScore[i].lk; 
       	
        meanVal = (sum/(double)tabScore.size());
        return meanVal;
}

double DistribNorm::meanVrais(){
 	double sum=0.0;
	
	if(meanVal != 0.0)
		return meanVal; 


	if(tabScore.size() == 0){
		cout << "No normalization score available for mean computation: mean=0.0" << endl;
		return 0.0;
	}
       
	for(unsigned int i=0; i<tabScore.size(); i++)
		sum += exp(tabScore[i].lk); 
       	
	meanVal = log((sum/(double)tabScore.size()));
        return meanVal;
}


double DistribNorm::std(){
	double sum2 = 0.0, m;
		
	if(stdVal != 1.0)
	 	return stdVal;
	
	if(tabScore.size() == 0){
		cout << "No normalization score available for std computation: std=1.0" << endl;
		return 1.0;	 
	}
 
       	m = mean();
	for(unsigned int i=0; i<tabScore.size(); i++){
		sum2 += tabScore[i].lk*tabScore[i].lk; 
	}
        stdVal = sqrt((sum2/(double)(tabScore.size())-(m*m)));
	
        return stdVal;
}

Histo DistribNorm::histo(){
	Histo h(10);
	
	if(tabScore.size() == 0){
		cout << "No normalization score available for histo computation" << endl;
	}
       
	for(unsigned int i=0; i<tabScore.size(); i++){
		h.accumulateValue(tabScore[i].lk);
	}
	h.computeHisto();
	
	return h;
}

void DistribNorm::sortScore(){
  	tabScore.descendingSort();
}

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
	
	String do_verbose = conf.getParam("verbose");
	static bool verbose=false;
	if (do_verbose == "true")
		verbose=true;
	String do_debug = conf.getParam("debug");
	static bool debug=false;
	if (do_debug == "true")
		debug=true;

	meanVal = 0.0;
	stdVal = 1.0;
}

DistribNorm::~DistribNorm() {
 	delete []tabIdImp;
}


//-------------------------------------------------------------------------------------------------------

void getScore(XList &list, unsigned long indIdClient, unsigned long indIdImp, String name, unsigned long indScore, DistribNorm &norm){
	XLine *linep;
	
	// ind* are used to locate information in NIST File in the case of format change
	list.getLine(0);
	while ((linep=list.getLine()) != NULL){
		if(linep->getElement(indIdClient) == name){
		 	 //cout << "score: " << linep->getElement(indRet).toDouble() << endl;
			 norm.addScore(linep->getElement(indScore).toDouble(), linep->getElement(indIdImp));
		}
	}
}																		

//-------------------------------------------------------------------------------------------------------
// Retrieve znorm scores after applying tnorm normalization !!!!
 
void getScoreTnormed(XList &list, unsigned long indElt, unsigned long indIdImp, String name, unsigned long indRet, Norm &tnorm, 
 			DistribNorm &norm){
	XLine *linep;
	
	list.getLine(0);
	while ((linep=list.getLine()) != NULL){
		if(linep->getElement(indElt) == name){		 	 
			 int ind = tnorm.findEntityInNorm(linep->getElement(3));		
			 if(ind != -1)
				norm.addScore((linep->getElement(indRet).toDouble() - tnorm.getMean(ind)) / tnorm.getStd(ind), linep->getElement(indIdImp));
			 else{
			 	cout << "tnormed segment not found!" << endl;
				exit(-1); 
			 }
		}
	}
}																		

//-------------------------------------------------------------------------------------------------------

int decision=0; //change ASAP
 
int ComputeNorm(Config& config)
{

	using namespace alize;
	using namespace std;

	
	String do_verbose = config.getParam("verbose");
	static bool verbose=false;
	if (do_verbose == "true")
		verbose=true;
	String do_debug = config.getParam("debug");
	static bool debug=false;
	if (do_debug == "true")
		debug=true;
	
        String outputNISTFileName = config.getParam("outputFile");                       // Result file in NIST style (.nist) result file format
	// Norm Type = znorm, tnorm ou ztnorm
	unsigned long maxIdNb = (unsigned long)(config.getParam("maxIdNb").toLong());
	unsigned long maxSegNb = (unsigned long)(config.getParam("maxSegNb").toLong());
	String normType = config.getParam("normType");
	String testNistFile = config.getParam("testNistFile"); 

	String outputFilename=config.getParam("outputFile");
	ofstream outFile(outputFilename.c_str(),ios::out | ios::trunc);    // Initialise the output file

	String selectType = "noSelect";
	if(config.existsParam("selectType")){
		selectType = config.getParam("selectType");
	}
	if((normType != "tnorm") && (selectType == "selectTargetDependentCohortInFile")){
		cout << "Caution: The Target dependent cohort selection works only with tnorm normalization technique" << endl;
		cout << "selectType => noSelect" << endl;		
		selectType = "noSelect";
	}

	static bool warpingNorm = false;
	double nonObserved=0.0;
	double area=0.0;
	Histo normalLawHisto;
	if(config.existsParam("warpingNorm")){
		String do_warpingNorm = config.getParam("warpingNorm");
		if(do_warpingNorm == "true"){	
			String normalLawHistoFilePath = "./";
			if(config.existsParam("normalLawHistoFilePath"))
				normalLawHistoFilePath = config.getParam("normalLawHistoFilePath");

			String normalLawHistoFileExt = ".txt";
			if(config.existsParam("normalLawHistoFileExt"))
				normalLawHistoFileExt = config.getParam("normalLawHistoFileExt");
			
			String normalLawHistoFileName = normalLawHistoFilePath+config.getParam("normalLawHistoFileName")+normalLawHistoFileExt;
			if(verbose)
				cout << "Tnorm Warping applied. normal law file: " << normalLawHistoFileName << endl;		
			normalLawHisto.load(normalLawHistoFileName);
			nonObserved=0.1;
			if (config.existsParam("nonObserved")) nonObserved=config.getParam("nonObserved").toDouble();
			area=0.5;
			if (config.existsParam("refArea")) area=config.getParam("refArea").toDouble();
			warpingNorm = true;
		}
	}
	
	// Define the field of the score files
	char fieldGender=0;
	char fieldName=1;
	//char fieldDecision=2;
	char fieldSeg=3;
	char fieldLLR=4;

	XList testList(testNistFile, config);
	testList.getLine(0); // goto first line still to debug
	
	try{
        ofstream outNist(outputNISTFileName.c_str(),ios::out | ios::trunc);    // Initialise the output file
		
	if(normType == "tnorm"){
				
		// storage of mean and std for all segments
	 	Norm tnorm(maxSegNb);
		String tnormNistFile = config.getParam("tnormNistFile"); 
		String sMax = config.getParam("tnormCohortNb");
		config.setParam("cohortNb", sMax);
		DistribNorm tnormSeg(config);
		
		// Mean and std for each seg
		XList tnormList(tnormNistFile, config);
	 	XLine * linep;
		while ((linep=testList.getLine()) != NULL){
			// The test segment filename
			String seg=linep->getElement(fieldSeg);
			int ind=0;
			/* In the specific case of "selectTargetDependentCohortInFile"
			the impostor distribution may change gor the same segment
			depending on the target cohort => ind=-1 for reloading and selecting new score 
			*/
			if(selectType == "selectTargetDependentCohortInFile"){
				ind = -1;
			}
			else{
				ind = tnorm.findEntityInNorm(seg);	
			}
			if(verbose) cout << endl << "Compute Tnorm segment: " << seg << endl; 			
			// if ind != -1 => the score distribution for that segment has been already stored
			if(ind == -1){
				// segment mean and std to compute
				// storage of score distribution for a given segment
			 	tnormSeg.init();
				// Lit nist tnorm file and copy score distribution in tnormseg
				getScore(tnormList, fieldSeg, fieldName, seg, fieldLLR, tnormSeg);
				// Apply selection of scores if necessary
				tnormSeg.selectScore(linep, selectType);
				
				if(warpingNorm){
					ind = tnorm.addSeg(seg, tnormSeg.histo());
				}
				else
					ind = tnorm.addSeg(seg, tnormSeg.mean(), tnormSeg.std());
			}
			double score;
			if(warpingNorm){
				score = scoreWarping((linep->getElement(fieldLLR)).toDouble(), tnorm.getHisto(ind), normalLawHisto, nonObserved, area); 
				if(debug)
					cout << endl << "score before warping: " << (linep->getElement(fieldLLR)).toDouble() << " / score after warping: " << score << endl;
			}
			else
				score = ((linep->getElement(fieldLLR)).toDouble() - tnorm.getMean(ind)) / tnorm.getStd(ind);
				
	        	outputResultLine(score, linep->getElement(fieldName), seg, linep->getElement(fieldGender), decision, outFile);				 			 
		}
	        if((debug) && (warpingNorm == false)){
			cout << "tnorm distrib" << endl;
	       		tnorm.print();       	 
		}
		      
	}
        else if(normType == "znorm"){	
		if(verbose) cout << endl << "Compute Znorm" << endl; 
		// storage of mean and std for all segments
	 	Norm znorm(maxIdNb);
		String znormNistFile = config.getParam("znormNistFile"); 
		String sMax = config.getParam("znormCohortNb");
		config.setParam("cohortNb", sMax);
		DistribNorm znormSeg(config);
		
		XList znormList(znormNistFile, config);
	 	XLine * linep;		
		
		while ((linep=testList.getLine()) != NULL){
			// The test segment filename
			String id=linep->getElement(fieldName);
			int ind = znorm.findEntityInNorm(id);		
			// segment mean and std already computed
			if(ind == -1){
			        // segment mean and std to compute
				// storage of score distribution for a given identity
			 	znormSeg.init();
				// Lit nist znorm file and copy score distribution in znormseg
				getScore(znormList, fieldName, fieldSeg, id, fieldLLR, znormSeg);
				// Apply selection of scores if necessary
				znormSeg.selectScore(linep, selectType);
				ind = znorm.addSeg(id, znormSeg.mean(), znormSeg.std());								 
			}
			double score = ((linep->getElement(fieldLLR)).toDouble() - znorm.getMean(ind)) / znorm.getStd(ind);
	        	outputResultLine(score, id, linep->getElement(fieldSeg), linep->getElement(fieldGender), decision, outFile);				 			 
		}	 
	        if(debug){
			cout << "znorm distrib" << endl;
	       		znorm.print();       	 
		}		
	}
        else if(normType == "ztnorm"){
		if(verbose) cout << endl << "Compute ZTnorm" << endl; 
		// storage of mean and std for all segments
	 	Norm tnorm(maxSegNb);
	 	Norm znorm(maxIdNb);
		Norm ztnorm(maxSegNb);
		
		String tnormNistFile = config.getParam("tnormNistFile"); 
		String znormNistFile = config.getParam("znormNistFile"); 
		String ztnormNistFile = config.getParam("ztnormNistFile"); 
		
		String sMax = config.getParam("tnormCohortNb");
		config.setParam("cohortNb", sMax);
		DistribNorm tnormSeg(config);
		  	
		sMax = config.getParam("znormCohortNb");
		config.setParam("cohortNb", sMax);
		DistribNorm znormSeg(config);
					
		XList tnormList(tnormNistFile, config);
		XList znormList(znormNistFile, config);
		XList ztnormList(ztnormNistFile, config);		
		
		XLine *linep;
		// Mean and std for each znorm seg
		while ((linep=znormList.getLine()) != NULL){
			// The ztnorm segment filename
			String seg=linep->getElement(fieldSeg);
			int ind = ztnorm.findEntityInNorm(seg);		
			if(ind == -1){
				// storage of score distribution for a given segment
				tnormSeg.init();
			 	// Lit nist ztnorm file and copy score distribution in tnormseg
				getScore(ztnormList, fieldSeg, fieldName, seg, fieldLLR, tnormSeg);
				// Apply selection of scores if necessary
				tnormSeg.selectScore(linep, selectType);
				ind = ztnorm.addSeg(seg, tnormSeg.mean(), tnormSeg.std());								 
			}			 			 
		}
		if(debug){
		        cout << "ztnorm distrib" << endl;
			ztnorm.print();       	 
		}

		while ((linep=testList.getLine()) != NULL){
			// The test segment filename
			String seg=linep->getElement(fieldSeg);
			int ind = tnorm.findEntityInNorm(seg);		
			// segment mean and std already computed
			if(ind == -1){
			        // segment mean and std to compute
				// storage of score distribution for a given segment
				tnormSeg.init();
				// Lit nist tnorm file and copy score distribution in tnormseg
				getScore(tnormList, fieldSeg, fieldName, seg, fieldLLR, tnormSeg);
				// Apply selection of scores if necessary
				tnormSeg.selectScore(linep, selectType);
				ind = tnorm.addSeg(seg, tnormSeg.mean(), tnormSeg.std());
				 				
			}
			double scoreTNorm = ((linep->getElement(fieldLLR)).toDouble() - tnorm.getMean(ind)) / tnorm.getStd(ind);
			
		        String id= linep->getElement(fieldName);
			ind = znorm.findEntityInNorm(id);
			// segment mean and std already computed
			if(ind == -1){
			        // segment mean and std to compute
				// storage of score distribution for a given identity
			 	znormSeg.init();
				// Lit nist znorm file, tnorm the scores and copy score distribution in znormseg
				getScoreTnormed(znormList, fieldName, fieldSeg, id, fieldLLR, ztnorm, znormSeg);
				// Apply selection of scores if necessary
				znormSeg.selectScore(linep, selectType);
				ind = znorm.addSeg(id, znormSeg.mean(), znormSeg.std());								 
			}
			double score = (scoreTNorm - znorm.getMean(ind)) / znorm.getStd(ind);
			
		       	outputResultLine(score, linep->getElement(fieldName), seg, linep->getElement(fieldGender), decision, outFile);				 
		}	        	 
	        if(debug){
		        cout << "tnorm distrib" << endl;
		       	tnorm.print();       	 
			cout << "znorm distrib" << endl;
		       	znorm.print();       	 
		}
	}

	if(normType == "poolnorm"){
				
		// storage of mean and std for all segments
	 	Norm poolnorm(maxSegNb);
		String poolnormNistFile = config.getParam("poolnormNistFile"); 
		String sMax = config.getParam("poolnormCohortNb");
		config.setParam("cohortNb", sMax);
		DistribNorm poolnormSeg(config);
		DistribNorm poolnormId(config);
		
		// Mean and std for each seg
		XList poolnormList(poolnormNistFile, config);
	 	XLine * linep;
		
		while ((linep=testList.getLine()) != NULL){
			// The test segment filename
		
			String seg=linep->getElement(fieldSeg);
			if(verbose) cout << endl << "Retrieve poolnorm score for segment: " << seg << endl; 			
			int indSeg = poolnorm.findEntityInNorm(seg);	
			// if ind != -1 => the score distribution for that segment has been already stored
			// segment mean and std already computed
			if(indSeg == -1){	
				// segment mean and std to compute
				// storage of score distribution for a given segment
			 	poolnormSeg.init();
				getScore(poolnormList, fieldSeg, fieldName, seg, fieldLLR, poolnormSeg);
				poolnormSeg.selectScore(linep, selectType);
				indSeg = poolnorm.addSeg(seg, poolnormSeg.sum(), poolnormSeg.sum2(), poolnormSeg.size());
			}
			else{
				cout << "Found index: " << indSeg << " for segment: " << seg << endl;
			}	
						
			String id=linep->getElement(fieldName);
			if(verbose) cout << endl << "Retrieve poolnorm score for id: " << id << endl; 			
			int indId = poolnorm.findEntityInNorm(id);		
			// if ind != -1 => the score distribution for that segment has been already stored
			// segment mean and std already computed
			if(indId == -1){
				// segment mean and std to compute
				// storage of score distribution for a given segment
			 	poolnormId.init();
				getScore(poolnormList, fieldName, fieldSeg, id, fieldLLR, poolnormId);
				poolnormId.selectScore(linep, selectType);
				indId = poolnorm.addSeg(id, poolnormId.sum(), poolnormId.sum2(), poolnormId.size());
			}
			else{
				cout << "Found index: " << indId << " for id: " << id << endl;
			}	
			// Fusion of pooled scored
			// getMean and getStd return in this case the sum and squared sum !
			double sum = poolnorm.getSum(indId) + poolnorm.getSum(indSeg);
			double sum2 = poolnorm.getSum2(indId) + poolnorm.getSum2(indSeg);
			
			double mean = sum / (double)(poolnorm.getNb(indId) + poolnorm.getNb(indSeg));
			double std = sqrt(sum2/(double)(poolnorm.getNb(indId) + poolnorm.getNb(indSeg)) - (mean*mean));
			
			if(verbose){
				double mId = poolnorm.getSum(indId)/(double)(poolnorm.getNb(indId));
				double sId =  sqrt(poolnorm.getSum2(indId)/(double)(poolnorm.getNb(indId)) - mId*mId);
				double mSeg = poolnorm.getSum(indSeg)/(double)(poolnorm.getNb(indSeg));
				double sSeg =  sqrt(poolnorm.getSum2(indSeg)/(double)(poolnorm.getNb(indSeg)) - mSeg*mSeg);
			
				cout << linep->getElement(fieldName) << " " << seg << " => " << mean << " " << std << endl;
			
				cout << linep->getElement(fieldName) << " => " << mId << " " << sId << " " << poolnorm.getNb(indId) << " " << poolnorm.getSum(indId) << " " << poolnorm.getSum2(indId) << endl;
				cout << seg << " " << " => " << mSeg << " " << sSeg << " " << poolnorm.getNb(indSeg) << " " << poolnorm.getSum(indSeg) << " " << poolnorm.getSum2(indSeg) <<endl;
			}
			double score = ((linep->getElement(fieldLLR)).toDouble() - mean) / std;
			
				
	        	outputResultLine(score, linep->getElement(fieldName), seg, linep->getElement(fieldGender), decision, outFile);				 			 
		}
	        if(debug){
			cout << "poolnorm distrib" << endl;
	       		poolnorm.print();       	 
		}
		      
	}


	outNist.close();
	      
	} //try
	catch (Exception& e){ 
		cout << e.toString().c_str() << endl;
	}

              
    return 0;
}

#endif //!defined(ALIZE_ComputeNorm_cpp)

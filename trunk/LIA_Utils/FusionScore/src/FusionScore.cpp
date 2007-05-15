#if !defined(ALIZE_FusionScore_cpp)
#define ALIZE_FusionScore_cpp

#include <iostream>
#include <fstream>  // pour outFile
#include <cstdio>   // pour printf()
#include <cassert> // pour le debug pratique
#include <cmath>
#include "FusionScore.h"
#include <liatools.h>

using namespace alize;
using namespace std;

//-------------------------------------------------------------------------
void computeFusion(XLine & lineTab, XLine & finalScore, XList & weights,Histo &destH, Config & config)
{
	String fusionMethod=config.getParam("fusionMethod");

	if (fusionMethod=="ArithMean")
	{
		double Score=0;
		double sumWeight=0;
		for (unsigned long i=0; i<lineTab.getElementCount(); i++)
		{
			cout << weights.getLine(0).getElement(i) << " * " << lineTab.getElement(i) << endl;
			sumWeight += weights.getLine(0).getElement(i).toDouble();
			Score+=weights.getLine(0).getElement(i).toDouble()*lineTab.getElement(i).toDouble();
		}
		Score/=sumWeight; 
		finalScore.addElement(String::valueOf(Score));
	}
	else if (fusionMethod=="GeoMean")
	{
		double Score=1;
		double sumWeight=0;
		double nbElm=lineTab.getElementCount();
		for (int i=0; i<nbElm; i++){
			sumWeight += weights.getLine(0).getElement(i).toDouble();			
			Score*=pow(weights.getLine(0).getElement(i).toDouble(),lineTab.getElement(i).toDouble());
		}
		if (Score < 0)	{Score=pow(-1*Score,1.0/sumWeight); Score=-1*Score;}
		else {Score=pow(Score,1.0/sumWeight);}
		finalScore.addElement(String::valueOf(Score));
	}
	else if (fusionMethod=="multipleAryth"){
		double Score=0;
		long bin= destH(lineTab.getElement(0).toDouble(),1); 
		if (debug) cout <<"score["<<lineTab.getElement(0).toDouble()<<"], fusion with line["<<bin<<"]"<<endl;
		for (unsigned long i=0; i<lineTab.getElementCount(); i++)
		{
			Score+=weights.getLine(bin).getElement(i).toDouble()*lineTab.getElement(i).toDouble();
		}
		//Score/=lineTab.getElementCount(); give percentage in weight, it's better
		finalScore.addElement(String::valueOf(Score));
	}

	else {cerr <<"fusionMethod"<<fusionMethod<<"] unknown"<<endl; exit(1);}
}

//-------------------------------------------------------------------------
void findNbest(XLine & tabScoreLine, int Nbest, int comp_double(const void *a, const void *b))
{
	int nbElm=tabScoreLine.getElementCount();
	double tab[nbElm];

	for (int i=0; i < nbElm; i++)
		tab[i]=tabScoreLine.getElement(i).toDouble();

	qsort(tab,nbElm,sizeof(double),comp_double);

	tabScoreLine.reset();

	for (int j=nbElm-1; j >= nbElm-Nbest; j--)
		tabScoreLine.addElement(String::valueOf(tab[j]));
}

//-------------------------------------------------------------------------
int comp_double(const void *a, const void *b)
{
	if ((*(double*) a - *(double*) b)<=0)	{return -1;}
	else	{return 1;}
}

String extractModel(String &classModel){
	String tmp;
	String m = classModel;
	long ind =  m.find("_", 0);
	for(int i=0; i<ind; i++) {
		tmp += m[i];
	}
	return tmp;
}

//-------------------------------------------------------------------------
int FusionScore(Config& config)
{
	using namespace alize;
	using namespace std;

	String inputNISTListFileName = config.getParam("inputFileList");
	String outputNISTFileName = config.getParam("outputFile");
	String weightsFile = config.getParam("weights");
	String format=config.getParam("format");
	bool debug=false;
	if (config.existsParam("debug")) debug=true;
	String fusionMethod=config.getParam("fusionMethod");
	Histo destH;
	if (fusionMethod=="multipleAryth")
	  destH.load(config.getParam("destH"));
	unsigned long nbFiles;
	unsigned long genderField, nameField, decisionField, segField, scoreField; 		//Necessary field in files
	if (format=="lia") {setLIAInfoFields(genderField, nameField, decisionField, segField,scoreField);}
	else if(format=="nist") {setNIST04InfoFields(genderField, nameField, decisionField, segField,scoreField);}
	else if (format=="etf") {scoreField=7;} // To remove, but sorry we're in eval!! retrieveETFInfoFields(genderField, nameField, decisionField, segField,scoreField);}
	else {cerr << "Error: Format unknown" << endl; exit(1);}
	
	XList tabScore;
	XLine finalScore;

	try{
		XList inputNISTList(inputNISTListFileName,config);	// List of files to fuse
		XList weights(weightsFile,config);	// coeffs
		nbFiles=inputNISTList.getLineCount();

		for (unsigned long i=0; i< nbFiles; i++) {
			XList inputNISTFile(inputNISTList.getLine(i).getElement(0),config);
			for (unsigned long j=0; j < inputNISTFile.getLineCount(); j++) {
				if (i==0) {tabScore.addLine();}
				tabScore.getLine(j).addElement(inputNISTFile.getLine(j).getElement(scoreField));
			}
		}

		for (unsigned long i=0; i < tabScore.getLineCount(); i++)	{//Fusion Process
			if (config.existsParam("Nbest"))	{ // In case of a nbest fusion
				findNbest(tabScore.getLine(i), config.getParam("Nbest").toLong(), comp_double);
			}

		if (debug) {
			for (unsigned long j=0; j < tabScore.getLine(i).getElementCount(); j++)
				cout << tabScore.getLine(i).getElement(j) << " ";
				cout << endl;
		}
			computeFusion(tabScore.getLine(i),finalScore,weights,destH,config); // See above method
		}
		XList inputNISTFile(inputNISTList.getLine(0).getElement(0),config);
		ofstream outFile(outputNISTFileName.c_str(),ios::out | ios::trunc);	// Preparing output
		String subType, event, channel, source;
		unsigned long duration, start;
		double LLR;
		int decision=0; // no decision in the fusion process
 
		if (format=="lia") {
			for (unsigned long i=0; i < tabScore.getLineCount(); i++) {
				setLIAInfoFields(genderField, nameField, decisionField, segField, scoreField); // retrieve LIA fields
				subType=inputNISTFile.getLine(i).getElement(genderField);
				event=inputNISTFile.getLine(i).getElement(nameField);
				decision=decision; // there is a pb here, outputResult takes a int as decision
				source=inputNISTFile.getLine(i).getElement(segField);
				LLR=finalScore.getElement(i).toDouble();
				outputResultLine(LLR,event,source, subType,decision,outFile);	
			}
		}
		if (format=="nist") {
			for (unsigned long i=0; i < tabScore.getLineCount(); i++) {
				String trainTypeTest="1side"; // clean this ASAP
				String adaptationMode="n";
				String segTypeTest="1side";
				setNIST04InfoFields(genderField, nameField, decisionField, segField, scoreField); // retrieve NIST fields
				subType=inputNISTFile.getLine(i).getElement(genderField);
				event=inputNISTFile.getLine(i).getElement(nameField);
				String decision="-";// there is a pb here, outputResult takes a int as decision
				source=inputNISTFile.getLine(i).getElement(segField);
				LLR=finalScore.getElement(i).toDouble();
				outputResultNIST04Line(trainTypeTest,adaptationMode,segTypeTest,subType,event,source,decision,LLR,outFile);				
			}
		}
		else if (format=="etf") {
			for (unsigned long i=0; i < tabScore.getLineCount(); i++) { // retrieve files from etf format has to be done
				subType=inputNISTFile.getLine(i).getElement(0);
				event=inputNISTFile.getLine(i).getElement(1);
				decision=inputNISTFile.getLine(i).getElement(2).toLong();
				source=inputNISTFile.getLine(i).getElement(3);
				start=inputNISTFile.getLine(i).getElement(4).toLong();
				duration=inputNISTFile.getLine(i).getElement(5).toLong();
				LLR=finalScore.getElement(i).toDouble();
				outputResultLine(LLR,event,source, start, start+duration, subType,decision,outFile);	
			}	
		}	
} // fin try

catch (Exception& e)
	{ 
	  cout <<" FusionScore" << e.toString().c_str() << endl;
	}
	return 0;
}


int FusionClassScore(Config& config)
{

	using namespace alize;
	using namespace std;

	String inputNISTFileName = config.getParam("inputNistFileName");
	String outputNISTFileName = config.getParam("outputFile");
	ofstream outFile(outputNISTFileName.c_str(),ios::out | ios::trunc);	// Preparing output
	String format=config.getParam("format");
	bool debug=false;
	if (config.existsParam("debug")) debug=true;

	String fusionMethod=config.getParam("fusionMethod");
	Histo destH;
	if (fusionMethod=="multipleAryth")
	  destH.load(config.getParam("destH"));
	
	unsigned long genderField, nameField, decisionField, segField, scoreField; 		//Necessary field in files
	if (format=="lia") {setLIAInfoFields(genderField, nameField, decisionField, segField,scoreField);}
	else if(format=="nist") {setNIST04InfoFields(genderField, nameField, decisionField, segField,scoreField);}
	else {cerr << "Error: Format unknown" << endl; exit(1);}

	/* Weight file managing */
	String weightFileExt = ".txt";
	if(config.existsParam("weightFileExt")){
		weightFileExt = config.getParam("weightFileExt");	
	}
	String weightFilePath = "./";
	if(config.existsParam("weightFilePath")){
		weightFilePath = config.getParam("weightFilePath");	
	}
	

	try{
		String subType, event, channel, source;
		double LLR;
		//int decision=0; // no decision in the fusion process
		String trainTypeTest="1side"; // clean this ASAP
		String adaptationMode="n";
		String segTypeTest="1side";

		XList tabScore;
		tabScore.addLine();
		XLine finalScore;
		XList weights;
		weights.addLine();
		XList tmpW;
		
		XList nistList(inputNISTFileName, config);
	 	XLine * currentlinep = nistList.getLine();
		XLine * nextlinep;
		
		String currentModel = extractModel(currentlinep->getElement(nameField));
		if(debug)
			cout << "Model: " << currentModel << " Accumule:" << currentlinep->getElement(scoreField) << endl;
		tabScore.getLine(0).addElement(currentlinep->getElement(scoreField));
		String weightFileName = weightFilePath+currentlinep->getElement(nameField)+weightFileExt;
		tmpW.load(weightFileName, config);
		weights.getLine(0).addElement(tmpW.getLine(0).getElement(0));

		while ((nextlinep=nistList.getLine()) != NULL){
			String nextModel = extractModel(nextlinep->getElement(nameField));		
			if(nextModel != currentModel){
				computeFusion(tabScore.getLine(0),finalScore,weights,destH, config); // See above method
				if(debug)
					cout << "Score Fusion: " << finalScore.getElement(0)<< endl;
				/* Output scores in files */
				if (format=="lia") {
					LLR=finalScore.getElement(0).toDouble();
					outputResultLine(LLR,currentModel,currentlinep->getElement(segField),currentlinep->getElement(genderField),int(currentlinep->getElement(decisionField).toLong()),outFile);	
				}
				if (format=="nist") {
					LLR=finalScore.getElement(0).toDouble();
					outputResultNIST04Line(trainTypeTest,adaptationMode,segTypeTest,currentlinep->getElement(genderField),currentlinep->getElement(nameField),currentlinep->getElement(segField),currentlinep->getElement(decisionField),LLR,outFile);				
				}
				
				currentModel = nextModel;
				tabScore.getLine(0).reset();
				weights.getLine(0).reset();
				finalScore.reset();
				currentlinep = nextlinep;	
				if(debug)
					cout << "Model: " << currentModel << " Accumule:" << currentlinep->getElement(scoreField) << endl;
				tabScore.getLine(0).addElement(currentlinep->getElement(scoreField));
				String weightFileName = weightFilePath+currentlinep->getElement(nameField)+weightFileExt;
				tmpW.load(weightFileName, config);
				weights.getLine(0).addElement(tmpW.getLine(0).getElement(0));
			}
			else{
				currentlinep = nextlinep;	
				if(debug)
					cout << "Model: " << currentModel << " Accumule:" << currentlinep->getElement(scoreField) << endl;
				tabScore.getLine(0).addElement(currentlinep->getElement(scoreField));
				String weightFileName = weightFilePath+currentlinep->getElement(nameField)+weightFileExt;
				tmpW.load(weightFileName, config);
				weights.getLine(0).addElement(tmpW.getLine(0).getElement(0));
			}
			
		}
		
		/* last line ???!!!! */
		computeFusion(tabScore.getLine(0),finalScore,weights,destH,config); // See above method
		if(debug)
			cout << "Score Fusion: " << finalScore.getElement(0)<< endl;
		/* Output scores in files */
		if (format=="lia") {
			LLR=finalScore.getElement(0).toDouble();
			outputResultLine(LLR,currentModel,currentlinep->getElement(segField),currentlinep->getElement(genderField),int(currentlinep->getElement(decisionField).toLong()),outFile);	
		}
		if (format=="nist") {
			LLR=finalScore.getElement(0).toDouble();
			outputResultNIST04Line(trainTypeTest,adaptationMode,segTypeTest,currentlinep->getElement(genderField),currentlinep->getElement(nameField),currentlinep->getElement(segField),currentlinep->getElement(decisionField),LLR,outFile);				
		}
		
} // fin try

catch (Exception& e)
	{ 
	  cout << "FusionClasse :"<<e.toString().c_str() << endl;
	}
	return 0;
}

#endif //!defined(ALIZE_FusionScore_cpp)

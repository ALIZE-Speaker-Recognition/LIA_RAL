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

#if !defined(ALIZE_Svm_cpp)
#define ALIZE_Svm_cpp

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <stdlib.h>
#include "Svm.h"
//#include "svm.h"
#include "liatools.h"

#define INF HUGE_VAL
#define TAU 1e-12	

using namespace alize;
using namespace std;

/*******************************************
 This is C default in SVM Light (inverse of average of train examples vectors norm)
*********************************************/
double getC(Matrix <double> train) {
	double acc=0.0;
	for (unsigned long i=0;i<train.rows();i++) 
		for (unsigned long j=0;j<train.cols();j++) 
			acc+=train(i,j)*train(i,j);
	acc/=train.rows();
return 1/acc;	
}

/*******************************************
Fill the struct to define parameters
*********************************************/
svm_parameter definesParameter (Config & config,Matrix <double> & mat) {
	struct svm_parameter param;
	param.svm_type = C_SVC;
	if (config.existsParam("kernelType")) {
		param.kernel_type = config.getParam("kernelType").toLong();
	} else {param.kernel_type=LINEAR;}
	param.degree = 1;
	param.gamma = 0;	// 1/k
	param.coef0 = 0;
	param.nu = 0.5;
	param.cache_size = 10000;
	if (config.existsParam("C")) {
		param.C = config.getParam("C").toDouble();
	} else {param.C = getC(mat);} // c value of svm light
	if (verbose) cout << "(Svm) Regularization Factor C=" << param.C <<endl;
	param.eps = 1e-3;
	param.p = 0.1;
	param.shrinking = 1;
	param.probability = 0;
	param.nr_weight = 0;
	param.weight_label = NULL;
	param.weight = NULL;
	if (config.existsParam("targetPenalty")) { //unblanced data (only two classes target 1 and impostors -1 (penalty on target only)
		param.nr_weight = 1;
		param.weight_label = (int *)realloc(param.weight_label,sizeof(int)*param.nr_weight);
		param.weight = (double *)realloc(param.weight,sizeof(double)*param.nr_weight);
		param.weight_label[param.nr_weight-1] = 1; // penalty on false reject
		param.weight[param.nr_weight-1] = config.getParam("targetPenalty").toDouble();
	}
return param;
}

/*******************************************
Fill the struct to define problem (svm_problme,svm_node)
*********************************************/
void readProblem(Matrix <double> & inputProb,RealVector <double>& labels,struct svm_problem & prob,struct svm_parameter& param,struct svm_node *x_space) {
	if (verbose) cout << "(SVM) Reading SVM Classification pb" << endl;
	int max_index=0;
	int j=0;
	int i=0;
	for(i=0;i<prob.l;i++)
	{
		prob.x[i] = &x_space[j];
		prob.y[i] = labels[i];
		unsigned long k=0;
		while(k<inputProb.cols())
		{
			(x_space[j].index)=k+1; //idx non sparse
			(x_space[j].value)=inputProb(i,k); //val
			++j;
			k++;
		}		
		if(j>=1 && x_space[j-1].index > max_index)
			max_index = x_space[j-1].index;
		x_space[j++].index = -1;
	}
	if (verbose) cout << "(Svm) Problem in memory" << endl;	
}
/*********************************************************
Matrix functions
********************************************************/
void concat(Matrix<double> & RES,Matrix<double> & A) {
	unsigned long stop=RES.rows();
	RES.setDimensions(RES.rows()+A.rows(),RES.cols());	
	for (unsigned long i=0;i<A.rows();i++)
		for (unsigned long j=0;j<A.cols();j++)
			RES(stop+i,j)=A(i,j);
}

void copyLine(Matrix<double>&RES,Matrix<double>&X,unsigned long idx) {
	for (unsigned long i=0;i<X.cols();i++)
		RES(idx,i)=X(0,i);
}

void copy(Matrix<double>&X,Matrix<double>&D) {
	D.setAllValues(0.0);
	D.setDimensions(X.rows(),X.cols());
	double*x=X.getArray();
	double*d=D.getArray();
	for (unsigned long i=0;i<X.rows();i++)	
		for (unsigned long j=0;j<X.cols();j++)
			d[i*X.cols()+j]=x[i*X.cols()+j];
}

/*****************************************************
 Hyperplane functions
******************************************************/
struct Svm_Hyperplane {
	RealVector <double> w;
	double offset;
};
// Load hyperplane from file [val1,val2,...,valn,offset]
void loadHyperplane(String &filename,struct Svm_Hyperplane &h,Config & config) {
	Matrix <double> model;
	model.load(filename,config);
	h.w.setSize(model.cols());
	for (unsigned long i=0;i<model.cols()-1;i++) 
		h.w[i]=model(0,i);
	h.offset=model(0,model.cols()-1);	
}
// return hyperplane in a struct
void getHyperplane(struct Svm_Hyperplane &h,const svm_model *model,RealVector <double> &labels,Config & config) {
	int libsvm_bug=1;
	if(labels[0]!=1) libsvm_bug=-1;
	int l = model->l;
	const double * const *sv_coef = model->sv_coef;
	const svm_node * const *SV = model->SV;
	unsigned long vsize=config.getParam("vsize").toLong();
	h.w.setSize(vsize);
	h.w.setAllValues(0.0);
	for(int i=0;i<l;i++) {
		const svm_node *p = SV[i];
		for (unsigned long j=0;j<h.w.size();j++,p++) 
			h.w[j]+=libsvm_bug*sv_coef[0][i]*p->value;
	}
	h.offset=libsvm_bug*model->rho[0];
}

// Save hyperplane + offset on one line (offset is last coeff)
void saveHyperplane(struct Svm_Hyperplane &h,String &outputFilename,Config & config) {
	Matrix <double> w;
	w.setDimensions(1,h.w.size()+1);
	for (unsigned long i=0;i<h.w.size();i++)
		w(0,i)=h.w[i];
	w(0,h.w.size())=h.offset;
	w.save(outputFilename,config);
}
/*******************************************************************************
From a matrix with examples and a mode return a vector of scores
********************************************************************************/
void predict(Matrix <double> & testInstances, String & model, RealVector <double>& scores,Config & config) {
	if (verbose) {cout << "(Svm) Processing model [" << model<<"] "<<endl;}
	struct Svm_Hyperplane h;
	loadHyperplane(model,h,config);
	/* norm of w */
	double norm_w=0.0;	
	for (unsigned long j=0;j<h.w.size();j++) // we don't use it ....
		norm_w+=h.w[j]*h.w[j];
	norm_w=sqrt(norm_w);
	/**/
	if (verbose) cout << "(Svm) Norm of hyperplane: ["<<norm_w<<"] Offset ["<<h.offset<<"]"<<endl;
	// Compute Distance wx+b / ||w|| (i think coefficents are already normalized by norm_w!!!!!!)
	unsigned long nbScores=testInstances.rows();
	unsigned long vsize=testInstances.cols();
	scores.setSize(nbScores);
	scores.setAllValues(0.0);
	double *_t=testInstances.getArray();
	for (unsigned long j=0;j<scores.size();j++)
		for (unsigned long k=0;k<vsize;k++)
			scores[j]+=h.w[k]*_t[j*vsize+k];
		//	scores[j]+=w[k]*testInstances(j,k);
	for (unsigned long j=0;j<scores.size();j++)
		scores[j]-=h.offset;
}

/*******************************************************************************
Load all Impostors examples once for all
********************************************************************************/
XLine loadBck(Matrix<double>&Bck,Config &config) {
	XList bckList(config.getParam("inputBckList"));	
	XLine impostorId;
	String vectorFilesPath = config.getParam("bckFilesPath");
	String vectorFilesExt = config.getParam("vectorFilesExtension");
	Bck.setDimensions(bckList.getLineCount(),config.getParam("vsize").toLong());
	Bck.setAllValues(0.0);
	if (verbose) cout << "(Svm) Bck Dimensions: ("<<Bck.rows()<<","<<Bck.cols()<<")"<<endl;
	XLine *pLine;
	unsigned long cnt=0;
	while((pLine=bckList.getLine())!=NULL) {
		Matrix <double> tmp;
		String file=vectorFilesPath+pLine->getElement(0)+vectorFilesExt;
		if (verboseLevel >1) cout<<"	(Svm) Adding: ["<<file<<"]"<<endl;
		impostorId.addElement(pLine->getElement(0));
		tmp.load(file,config);
		if (tmp.cols()!=(unsigned long)config.getParam("vsize").toLong()) throw Exception("Param vectSize and vector size does not match",__FILE__,__LINE__);
		copyLine(Bck,tmp,cnt);
		cnt++;
	}
	bckList.rewind();
return impostorId;
}

/*******************************************************************************
Train function
********************************************************************************/
int svmTrain(Config &config)
{
	String outputFilename = config.getParam("outputFilename");
	String vectorFilesPath = config.getParam("vectorFilesPath");
	String modelFilesPath = config.getParam("modelFilesPath");	
	String vectorFilesExt = config.getParam("vectorFilesExtension");
	String modelFilesExt = config.getParam("modelFilesExtension");	
	String inputFilename = config.getParam("inputFilename");
	try{
		Matrix <double> Bck;
		if (verbose) cout << "(Svm) Using ["<<config.getParam("inputBckList")<<"] as additionnal exemples for training" << endl;
		XLine impostorId=loadBck(Bck,config); // loading bck exemples	
		XList inputList;
		bool tnorm=false;
		if (config.existsParam("Tnorm")) {tnorm=true;if(verbose)cout <<"(Svm) Tnorm mode: Leave one out in inputBckList"<<endl;}
		inputList.load(inputFilename,config);
		XLine *pLine;

		while((pLine=inputList.getLine())!=NULL) {
			String clientModelname=pLine->getElement(0);
			String clientModelTrain=pLine->getElement(0); // same id			
			inputFilename=vectorFilesPath+clientModelTrain+vectorFilesExt;
			outputFilename=modelFilesPath+clientModelname+modelFilesExt;
			if (verbose) cout << "(Svm) Using ["<<inputFilename<<"] as training for: [" << outputFilename<<"]"<<endl;
			Matrix <double>clientVector;
			RealVector <double> labels;
			if (!tnorm) { 
				clientVector.load(inputFilename,config);
				unsigned long nbSessions=clientVector.rows();
				if (verbose) cout << "(Svm) Number of sessions: " << nbSessions << endl; // todo
				concat(clientVector,Bck);
				labels.setSize(clientVector.rows());
				for (unsigned long i=0;i<clientVector.rows();i++) {
					if (i<nbSessions) labels[i]=1.0;
					else labels[i]=-1.0;
				}				
			}
			/** Tnorm is a special case, just change the labels and do not touch to the matrix **/
			else {
				copy(Bck,clientVector);
				labels.setSize(clientVector.rows());
				for (unsigned long i=0;i<clientVector.rows();i++) { // only one session
					if (clientModelTrain==impostorId.getElement(i)) {labels[i]=1.0;}
					else labels[i]=-1.0;
				}	
			}				

			if (verbose) cout << "(Svm) Pb Dimensions: ("<<clientVector.rows()<<","<<clientVector.cols()<<")"<<endl;
			//***************** Defines structures, and read Problem ****************************************//
			struct svm_parameter param=definesParameter(config,clientVector);
			struct svm_problem prob;	
			prob.l=clientVector.rows(); //number of exemples 				
			prob.y=new double[prob.l];
			prob.x=new svm_node*[prob.l];
			struct svm_node *x_space=NULL;
			x_space=new svm_node[(prob.l+1)*clientVector.cols()];	// the +1 was tricky to find!!!!!!!!
			readProblem(clientVector,labels,prob,param,x_space); //alize way or reading problem
			//************************************************************//	
				
			const char* error_msg = svm_check_parameter(&prob,&param);
			if(error_msg) {
				fprintf(stderr,"Error: %s\n",error_msg); exit(1);
			}
			else {				
				if (verbose) cout<<"(Svm) Training Model"<<endl;
				svm_model *model=svm_train(&prob,&param);
				if (verboseLevel > 1) {
					for (unsigned long i=0;i<(unsigned long)model->l;i++) 
						cout << "Sv["<<i<<"] ="<<model->sv_coef[0][i]<<endl;	
				}
				struct Svm_Hyperplane h;
				getHyperplane(h,model,labels,config); // only save w!
				saveHyperplane(h,outputFilename,config);				
				if (verbose) cout<<"(Svm) Model saved in "<<outputFilename<<endl;						
				svm_destroy_model(model);
			}
			/* Destroy for each training */
			svm_destroy_param(&param);
			delete prob.y;
			delete prob.x;
			delete[] x_space;
		}
	}	
	catch (Exception& e){cout << e.toString().c_str() << endl;}
return 0;
}

/*******************************************************************************
Test function
********************************************************************************/
int svmPredict(Config &config)
{
	String inputFilename = config.getParam("inputFilename");
	String outputFilename = config.getParam("outputFilename");
	String vectorFilesPath = config.getParam("vectorFilesPath");
	String modelFilesPath = config.getParam("modelFilesPath");	
	String vectorFilesExt = config.getParam("vectorFilesExtension");
	String modelFilesExt = config.getParam("modelFilesExtension");	

	try{
		String gender="M";
		int decision=0;
		ofstream outNist(outputFilename.c_str(),ios::out | ios::trunc); 
		XList inputNdx(inputFilename,config);
		XLine *pLine;
		while((pLine=inputNdx.getLine())!=NULL) {		
			String modelFilename=modelFilesPath+pLine->getElement(0)+modelFilesExt;
			Matrix <double> testInstances;
			testInstances.setDimensions(pLine->getElementCount()-1,config.getParam("vsize").toLong());
			if (verbose) cout << "(Svm) For "<<pLine->getElement(0) << ", add " <<pLine->getElementCount()-1 << " tests ";
			for (unsigned long i=1;i<pLine->getElementCount();i++) {
				Matrix <double> tmp;
				String testFile=vectorFilesPath+pLine->getElement(i)+vectorFilesExt;
				if (verboseLevel > 1) {cout << " [" << pLine->getElement(i)<<"] ";}	
				tmp.load(testFile,config);				
				if (tmp.cols()!=(unsigned long)config.getParam("vsize").toLong()) throw Exception("Param vectSize and vector size does not match",__FILE__,__LINE__);
				copyLine(testInstances,tmp,i-1);
			}	
			if (verbose) cout <<endl<<"(Svm) Test dimensions: ["<<testInstances.rows()<<","<<testInstances.cols()<<"]"<<endl;
			RealVector <double> scores;
			predict(testInstances,modelFilename,scores,config);		
			for (unsigned long i=1;i<pLine->getElementCount();i++) {
				if (verboseLevel >1) outputResultLine(scores[i-1], pLine->getElement(0),pLine->getElement(i) ,gender ,decision,cout);
				outputResultLine(scores[i-1], pLine->getElement(0),pLine->getElement(i) ,gender ,decision,outNist);
			}
		}
		outNist.close();
	}	
	catch (Exception& e){cout << e.toString().c_str() << endl;exit(1);}
return 0;
}
/********************************************************************************
Test Tnorm function (do not reload examples) : This should go away if we have a TabClientLine for SVM
********************************************************************************/
int svmPredictTnorm(Config &config)
{
	String inputFilename = config.getParam("inputFilename");
	String outputFilename = config.getParam("outputFilename");
	String vectorFilesPath = config.getParam("vectorFilesPath");
	String modelFilesPath = config.getParam("modelFilesPath");	
	String vectorFilesExt = config.getParam("vectorFilesExtension");
	String modelFilesExt = config.getParam("modelFilesExtension");	

	try{
		String gender="M";
		int decision=0;
		ofstream outNist(outputFilename.c_str(),ios::out | ios::trunc); 
		XList inputNdx(inputFilename,config);
		Matrix <double> testInstances;
		testInstances.setDimensions(inputNdx.getLine(0).getElementCount()-1,config.getParam("vsize").toLong());
		for (unsigned long i=1;i<inputNdx.getLine(0).getElementCount();i++) {
			Matrix <double> tmp;
			String testFile=vectorFilesPath+inputNdx.getLine(0).getElement(i)+vectorFilesExt;
			if (verboseLevel > 1) {cout << " [" << inputNdx.getLine(0).getElement(i)<<"] ";}			
			tmp.load(testFile,config);
			if (tmp.cols()!=(unsigned long)config.getParam("vsize").toLong()) throw Exception("Param vectSize and vector size does not match",__FILE__,__LINE__);
			copyLine(testInstances,tmp,i-1);
		}		
		XLine *pLine;		
		while((pLine=inputNdx.getLine())!=NULL) {	
			if (verbose) cout << "(Svm) For "<<pLine->getElement(0) << ", add " <<pLine->getElementCount()-1 << " tests ";			
			String modelFilename=modelFilesPath+pLine->getElement(0)+modelFilesExt;	
			if (verbose) cout <<endl<<"(Svm) Test dimensions: ["<<testInstances.rows()<<","<<testInstances.cols()<<"]"<<endl;
			RealVector <double> scores;
			predict(testInstances,modelFilename,scores,config);		
			for (unsigned long i=1;i<pLine->getElementCount();i++) {
				if (verboseLevel >1) outputResultLine(scores[i-1], pLine->getElement(0),pLine->getElement(i) ,gender ,decision,cout);
				outputResultLine(scores[i-1], pLine->getElement(0),pLine->getElement(i) ,gender ,decision,outNist);
			}
		}
		outNist.close();
	}	
	catch (Exception& e){cout << e.toString().c_str() << endl;exit(1);}
return 0;
}

double dotProduct(DoubleVector & v1,DoubleVector & v2) {
	double sum=0;
	for(unsigned long i=0;i<v1.size();i++) 
		sum+=v1[i]*v2[i];
return sum;	
}

#endif //!defined(ALIZE_Svm_cpp)


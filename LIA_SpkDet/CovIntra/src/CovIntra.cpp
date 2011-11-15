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

#if !defined(ALIZE_CovIntra_cpp)
#define ALIZE_CovIntra_cpp

#include <iostream>
#include <fstream>  // pour outFile
#include <cstdio>   // pour printf()
#include <cassert> // pour le debug pratique
#include <cmath>
#include <liatools.h>
#include "CovIntra.h"
extern "C" {
#include "svdlib.h" 
}

#if defined(_WIN32)
  #include <cfloat> // for _isnan()
  #define ISNAN(x) _isnan(x)
  #define ISINF(x) (!_finite(x))
#elif defined(linux) || defined(__linux) || defined(__CYGWIN__) || defined(__APPLE__)
  #define ISNAN(x) isnan(x)
  #define ISINF(x) isinf(x)
#else
  #error "Unsupported OS\n"
#endif


#define myTINY 1e-11
using namespace std;
using namespace alize;

// Divides a Vector by an Integer
void dividesSvByInteger(RealVector <double> & v, unsigned long nb) {
	if (nb==0) throw Exception("Divides by 0!",__FILE__,__LINE__);	
	for (unsigned long j=0;j<v.size();j++) 
		v[j]/=nb;
}

// Divides a Vector by an Integer
void dividesSvByInteger(Matrix <double> & v, unsigned long nb) {
	if (nb==0) throw Exception("Divides by 0!",__FILE__,__LINE__);
	for (unsigned long j=0;j<v.cols();j++) 
		v(0,j)/=nb;
}

bool badSv(Matrix <double> & v) {
	for (unsigned long j=0;j<v.cols();j++) {	
		if (abs(v(0,j))>myTINY) return false;
		if (ISNAN(v(0,j)) || ISINF(v(0,j))) return true;	
	}
return true;
}

void loadMeanSv(String id,Matrix <double> & v,Config & config) {
	MixtureServer ms(config);
	MixtureGD & curr=ms.loadMixtureGD(id);
	RealVector <double> currSv;
	if (config.getParam("vsize").toLong()!= (long)(ms.getDistribCount()*ms.getVectSize())) throw Exception("VectSize mismatch, check vsize parameter",__FILE__,__LINE__);
	currSv.setSize(config.getParam("vsize").toLong());
	if (config.existsParam("klgmm")) getKLVector(curr,currSv,config); 
	else modelToSv(curr,currSv);
	v.setDimensions(1,currSv.size());
	for (unsigned long j=0;j<currSv.size();j++) 	
		v(0,j)=currSv[j];
}

template <class T> class RealVectorStat {
	private:		
		unsigned long _size;
		unsigned long _nb;
		RealVector <T> _xAcc;
		RealVector <T> _xxAcc;
		RealVector <double> _mean;
		RealVector <double> _covInv;	
	public:
		explicit RealVectorStat(unsigned long size) {
			_size=size;
			_nb=0;
			_xAcc.setSize(_size);_xxAcc.setSize(_size);_mean.setSize(_size);_covInv.setSize(_size);
			_xAcc.setAllValues(0);_xxAcc.setAllValues(0);_mean.setAllValues(0.0);_covInv.setAllValues(0.0);
		}
		~ RealVectorStat() {};
		void acc(RealVector <T> &v) {
			if (_size!=v.size()) throw Exception("(RealVectorStat) size does not match!",__FILE__,__LINE__);
			for (unsigned long i=0;i<v.size();i++) {
				_xAcc[i]+=v[i];
				_xxAcc[i]+=v[i]*v[i];
			}
			_nb++;
		};
		RealVector <double> & covInv() {
			if (_nb==0)  throw Exception("(RealVectorStat) Nothing accumulated!",__FILE__,__LINE__);
			for (unsigned long i=0;i<_size;i++) {
				_mean[i]=((double)_xAcc[i])/_nb;
				_covInv[i]=1/((((double)_xxAcc[i])/_nb)-(_mean[i]*_mean[i]));
			}			
		return _covInv;
		};
		RealVector <double> & mean() {
			if (_nb==0)  throw Exception("(RealVectorStat) Nothing accumulated!",__FILE__,__LINE__);			
			for (unsigned long i=0;i<_size;i++)				
				_mean[i]=((double)_xAcc[i])/_nb;	
		return _mean;
		};
};

int CovIntra(Config & config) {
  try {
	bool gmm=false;
	if (config.existsParam("gmm")) {if(verbose) cout << "(CovIntra) Gmm mode: reading SV from gmm models"<<endl;gmm=true;}
	bool svd=false;
	if (config.existsParam("svd")) {if(verbose) cout << "(CovIntra) Computing the SVD"<<endl;svd=true;}	
	if (config.existsParam("klgmm")) {if(verbose) cout << "(CovIntra) Gmm KL mode: computing KL distance from gmm models"<<endl;gmm=true;}
	SVDVerbosity=verboseLevel;
	XList fileList(config.getParam("ndx"),config);
	XLine *pline;
	String *pModel;
	unsigned long svSize=config.getParam("vsize").toLong();
	unsigned long nbSpeakers=fileList.getLineCount();
	if (config.existsParam("nbSpeakers") && config.getParam("nbSpeakers").toLong() <= (long) nbSpeakers) nbSpeakers=config.getParam("nbSpeakers").toLong();
	
	// Compute desired Nb examples
	unsigned long nbEx=0;
	for (unsigned long s=0;s<nbSpeakers;s++)
		nbEx+=fileList.getLine(s).getElementCount();
	fileList.rewind();
	Matrix <double> UBMsv;
	loadMeanSv(config.getParam("inputWorldFilename"),UBMsv,config);
	DMat CCt=svdNewDMat(svSize,nbEx); //svdlib format	
	RealVectorStat <double> acc(svSize);
  	if (verbose) cout<<"Channel Matrix dimension: "<<svSize<<","<<nbEx <<endl << "Speaker Matrix dimension: "<<svSize<<","<<nbSpeakers<<endl;
	unsigned long idx=0;
	for (unsigned long r=0;r<nbSpeakers;r++) {
		pline=fileList.getLine();
		if (verbose) cout << "Sp ["<<r<<"] " <<pline->getElementCount() << " sessions, " <<endl;
		// Compute True Mean Statistic for Speaker in different sessions
		Matrix <double> meanSv;	
		meanSv.setDimensions(1,svSize);
		while((pModel=pline->getElement())!=NULL) {
			if (verboseLevel > 1) cout << "Load " <<*pModel<<endl;
			String filename;
			Matrix <double> currSv;
			if (!gmm) currSv.load(config.getParam("vectorFilesPath")+*pModel+config.getParam("vectorFilesExtension"),config);
			else loadMeanSv(*pModel,currSv,config);
			meanSv+=currSv;
		}
		
		/************************ This compute InterSpeaker variability in a R vector **************************************************/
		RealVector <double> vm;
		for (unsigned long i=0;i<meanSv.cols();i++) {
			double val=meanSv(0,i)/(pline->getElementCount());
			val-=UBMsv(0,i);
			vm.addValue(val);
		}
		acc.acc(vm);
		
		/************************* 2nd pass: Remove mean sv to get only channel (drop one) *******************************************/
		for (unsigned long c=0;c<pline->getElementCount();c++) {
			Matrix <double> currSv;					
			if (!gmm) currSv.load(config.getParam("vectorFilesPath")+pline->getElement(c)+config.getParam("vectorFilesExtension"),config);
			else loadMeanSv(pline->getElement(c),currSv,config);		
			Matrix <double> meanMinusCurrent(meanSv); // if drop one to compute mean on the line
			meanMinusCurrent-=currSv;
			dividesSvByInteger(meanMinusCurrent,(pline->getElementCount()-1));	
			currSv-=meanMinusCurrent;
			if (debug && badSv(currSv)) cout <<pline->getElement(c)<<" gives a bad channel supervector (zero|inf|nan value)"<<endl; // is faster as badSv is not evaluated if debug=0, at least that's what I think

			// Copy SV into matrix
			for(unsigned long val=0;val<currSv.cols(); val++) {
				CCt->value[val][idx]=currSv(0,val);
				if(abs(CCt->value[val][idx])<myTINY) CCt->value[val][idx]=0.0; // sparsify matrix
			}
			idx++;
		}
	}
	// Save interSp vector
	RealVector <double>& interSp=acc.mean();
	((Matrix <double>)interSp).saveDT((String)"mean.inter",config);

	String filename=config.getParam("channelMatrix");
	char *out=strdup(filename.c_str());
	if (verbose && !svd) cout << "(CovIntra) Covariance matrix saved in " << config.getParam("channelMatrix")<< ", need to perform eigenvalue analysis yourself or add --svd flag to the config file" <<endl;
	svdWriteDenseMatrix(CCt,out, SVD_F_DB);
	free(out);
	if (!svd) exit(1);
		
	/****** SVDLIBC part: EigenValue decomposition with SVD ********************/
	unsigned long nbIt=0;
	double kappa=1e-6;
	double las2end[2] = {-1.0e-30, 1.0e-30};	
	if (config.existsParam("nbIt")) nbIt=config.getParam("nbIt").toLong();
	if (config.existsParam("kappa")) kappa=config.getParam("kappa").toDouble();		
	if (config.existsParam("bound")) {
		las2end[1] = config.getParam("bound").toDouble();
		las2end[0] = -las2end[1];
		}
		
	SMat SCCt=svdConvertDtoS(CCt);	
	svdFreeDMat(CCt);		
	for(int i=0;i<SCCt->vals;i++)
		if (ISNAN(SCCt->value[i]) || ISINF(SCCt->value[i])) throw Exception("Error:nan or inf in sparse matrix file",__FILE__,__LINE__);
	unsigned long nbEv=config.getParam("nbEigenVectors").toLong();
	if (verbose) cout << "(CovIntra) Begin eigenValue problem resolution (using Doug Rohde's SVD C Library)" << endl;
	SVDRec res=svdLAS2(SCCt,nbEv,nbIt,las2end,kappa);
	svdFreeSMat(SCCt);
	/*Save EigValues*/
	Matrix <double> ev(1,res->d);
	for (int i=0;i<res->d;i++) ev(0,i)=res->S[i];
	String evFile=config.getParam("channelMatrix")+".S";
	if (verbose) cout << "(CovIntra) Singluar Values (square root) in " <<evFile<< endl;	
	ev.saveDT(evFile,config);
	/*Save EigVectors*/
	Matrix <double> NAP(res->Ut->rows,res->Ut->cols);
	if (verbose) cout<< "(CovIntra) Dimensions of U: ("<<NAP.rows()<<","<<NAP.cols()<<")"<< endl;
	for (int i=0;i<res->Ut->rows;i++)
		for (int j=0;j<res->Ut->cols;j++)
			NAP(i,j)=res->Ut->value[i][j];	
	svdFreeSVDRec(res);			
	if (verbose) cout << "(CovIntra) U saved in " << config.getParam("channelMatrix")<<endl;		
	NAP.save(config.getParam("channelMatrix"),config);
}
  catch (Exception& e) {cout << e.toString() << endl;}
return 0;
}



#endif

#if !defined(ALIZE_SuperVectors_cpp)
#define ALIZE_SuperVectors_cpp

#include "SuperVectors.h"
#include<iostream>
#include<fstream>
#include<cstdio>
#include<cassert>
#include<cmath>

using namespace alize;
using namespace std;

// comments in .h

void modelToSv(const MixtureGD &M,RealVector <double> &v){
	unsigned long modelSize=M.getDistribCount();
	unsigned long vectSize=M.getVectSize();	
	v.setSize(modelSize*vectSize);
	for (unsigned long i=0;i<modelSize;i++)
		for (unsigned long j=0;j<vectSize;j++)
			v[i*vectSize+j]=M.getDistrib(i).getMean(j);
}

void svToModel(RealVector  <double> &v, MixtureGD &M) {
	unsigned long modelSize=M.getDistribCount();
	unsigned long vectSize=M.getVectSize();	
	for (unsigned long i=0;i<modelSize;i++)
		for (unsigned long j=0;j<vectSize;j++)
			M.getDistrib(i).setMean(v[i*vectSize+j],j);
}

/*void projectOnSubSpace(RealVector<double> &x, Matrix<double> &U, RealVector <double>&proj) {
	// trick is to transform x'=S'Sx in x'=(U(U'x))' so that no need to compute U explicitely
	Matrix<double> Utmp=U;
	if (U.rows() > U.cols())
		Utmp.transpose();	
	Matrix<double> xt=((DoubleMatrix)x).transpose();
	Matrix<double> spkFactors=Utmp*xt;		
	if (verboseLevel > 1) {
		cout <<  "(SuperVectors) SpkFactors: (";
		for (unsigned long i=0;i<spkFactors.rows();i++)
			cout <<i<<":"<<spkFactors(i,0)<<","; 			
	}
	if (verboseLevel > 2) cout << "(SuperVectors) SpkFactors: "<<spkFactors.rows() <<" " <<spkFactors.cols() << "...";
	Matrix<double> Ut=Utmp.transpose();
	Matrix<double>  offset=Ut*spkFactors;
	if (verboseLevel > 2) cout << "(SuperVectors) offset: "<<offset.rows() <<" " <<offset.cols() <<"...";
	for (unsigned long i=0;i<offset.rows();i++) {
		proj[i]=offset(i,0);
	}
}*/

void projectOnSubSpace(RealVector<double> &x, Matrix<double> &U, RealVector <double>&proj) { //fast version
	// trick is to transform x'=S'Sx in x'=(U(U'x))' so that no need to compute U explicitely
	double *_u=U.getArray();
	RealVector <double> tmp;
	tmp.setSize(U.rows());
	tmp.setAllValues(0.0);
	proj.setAllValues(0.0);
	if (verboseLevel > 1) cout<<"x: ("<<x.size()<<") - U; ("<<U.rows()<<","<<U.cols()<<")"<<endl;
	for (unsigned long i=0;i<U.rows();i++) 
		for (unsigned long j=0;j<U.cols();j++) 
			tmp[i]+=_u[i*U.cols()+j]*x[j];
	U.transpose();
	_u=U.getArray();	
	for (unsigned long i=0;i<U.rows();i++) 
		for (unsigned long j=0;j<U.cols();j++) 
			proj[i]+=_u[i*U.cols()+j]*tmp[j];		
	U.transpose();
}


void computeNap(MixtureGD &M,Matrix <double> &U) {
	unsigned long svSize=M.getDistribCount()*M.getVectSize();	
	RealVector <double> v(svSize,svSize);
	v.setAllValues(0.0);
	RealVector <double> proj(svSize,svSize);		
	proj.setAllValues(0.0);	
	modelToSv(M,v);
	projectOnSubSpace(v,U,proj);
	v-=proj;
	svToModel(v,M);
}


// Estimate NAP channel effect with a fixed client supervector and put it into model
// 1. Project client vector on test vector
// 2. Project result on channel subspace
// 3. Compute NAP on client model
// 4. Add projection to client NAPed model
// Three same function are coded but none of them work right now!
void computeNAPChannelEffect(MixtureGD &T,MixtureGD &M,Matrix<double>&U) {
	unsigned long svSize=max(U.rows(),U.cols());
	RealVector <double> t(svSize,svSize);
	RealVector <double> m(svSize,svSize);
	RealVector <double> proj(svSize,svSize);
	RealVector <double> tmp(svSize,svSize);	
	RealVector <double> channel(svSize,svSize);	
	modelToSv(T,t);
	modelToSv(M,m);	
	double normT=0.0;
	double angle=0.0;
	for (unsigned long i=0;i<t.size();i++) // norm du test
		normT+=t[i]*t[i];
	for (unsigned long i=0;i<t.size();i++) // produit scalaire
		angle+=t[i]*m[i];
	if (verboseLevel > 1) cout << "# <t,m>/||t||2=" << angle/normT << endl;
	for (unsigned long i=0;i<t.size();i++) // norm du test
		proj[i]=(angle/normT)*t[i];	
	projectOnSubSpace(m,U,tmp); // Channel effect for test
	RealVector <double> mNap=m;
	mNap-=tmp;	
	projectOnSubSpace(proj,U,channel);
	channel+=mNap;
	if (debug) {
		double n=0.0;
		double a=0.0;	
		for (unsigned long i=0;i<t.size();i++) // norm du test
			n+=channel[i]*channel[i];
		for (unsigned long i=0;i<t.size();i++) // produit scalaire
			a+=t[i]*channel[i];
		cout << "# <model,test>/||model||2=" << a/n << endl;	
	}
	svToModel(channel,M);
}

// Get an estimate of channel effect for a model. 
//The image of projection of the SV client on the test vector is computed and the channel effect 
// estimate is done on this result
// This is computed with the Thales theorem
/*void computeNAPChannelEffect(MixtureGD &T,MixtureGD &M,Matrix<double>&U) {
	unsigned long svSize=max(U.rows(),U.cols());
	RealVector <double> t(svSize,svSize);
	RealVector <double> m(svSize,svSize);
	RealVector <double> proj(svSize,svSize);
	RealVector <double> tmp(svSize,svSize);
	RealVector <double> channel(svSize,svSize);	

	modelToSv(T,t);
	modelToSv(M,m);	
	projectOnSubSpace(t,U,proj); // Channel effect for test
	RealVector <double> tNap=t;
	tNap-=proj;
	projectOnSubSpace(m,U,tmp); // Channel effect for test
	RealVector <double> mNap=m;
	mNap-=tmp;

	double normT=0.0;
	double normM=0.0;
	for (unsigned long i=0;i<t.size();i++) // norm du test
		normT+=tNap[i]*tNap[i];
	for (unsigned long i=0;i<m.size();i++) // produit scalaire
		normM+=mNap[i]*mNap[i];
	cout << "# Expanding factor [" << sqrt(normM/normT) <<"]"<< endl;
	for (unsigned long i=0;i<channel.size();i++)
		channel[i]=sqrt(normM/normT)*proj[i];	

	channel+=mNap;
	double n=0.0;
	double a=0.0;	
	for (unsigned long i=0;i<t.size();i++) // norm du test
		n+=channel[i]*channel[i];
	for (unsigned long i=0;i<t.size();i++) // produit scalaire
		a+=t[i]*channel[i];
	cout << "# <model,test>/||model||2=" << a/n << endl;		
	svToModel(channel,M);
}*/

// Add UUt Mt to the napped SV for each model (same channel for every model)
/*void computeNAPChannelEffect(MixtureGD &T,MixtureGD &M,Matrix<double>&U) {
	unsigned long svSize=max(U.rows(),U.cols());
	RealVector <double> t(svSize,svSize);
	RealVector <double> m(svSize,svSize);
	RealVector <double> channel(svSize,svSize);
	modelToSv(T,t);
	projectOnSubSpace(t,U,channel);	
	computeNap(M,U);	
	modelToSv(M,m);		
	channel+=m;
	svToModel(channel,M);
}*/

//------------------------------------------------------------------------
// Weight part of fisher kernel -- Nicolas SCHEFFER (deprecated)
void getFisherWeightVector(const MixtureGD& world, const MixtureGD& clientMixture, RealVector<double> & v,Config& config) {
  double w,c;
  for (unsigned long i=0;i<world.getDistribCount();i++) {
    w=world.weight(i);
    c=clientMixture.weight(i);
    //v[i]=c/sqrt(w);
	v[i]=c/w;
  }
}


//------------------------------------------------------------------------
// Weight part of fisher kernel -- Nicolas SCHEFFER (deprecated)
void getKLVector(MixtureGD& model, RealVector<double> & v,Config& config) {
	unsigned long mSize=model.getDistribCount();
	unsigned long vSize=model.getVectSize();
	for (unsigned long i=0;i<mSize;i++) {
		DistribGD & d=model.getDistrib(i);
		double w=model.weight(i);
		for (unsigned long j=0;j<vSize;j++) {
			//v[i*vSize+j]=d.getMean(j)*(sqrt(w*d.getCovInv(j)));
			v[i*vSize+j]=d.getMean(j)*(sqrt(w*d.getCovInv(j)));
		}
	}
}

void getSuperVector(RealVector<double> &v,MixtureGD &aprioriModel,MixtureGD &clientMixture,Config &config) {
	if (config.getParam("superVector")=="SVMUBM"){ //weight supervector
	        v.setSize(clientMixture.getDistribCount(),true);
                getFisherWeightVector(aprioriModel,clientMixture,v,config); // For an SVM/UBM system (have to force client Mixture to be ML estimate MAPConst, alpha=1);
            }
	else if (config.getParam("superVector")=="KL") {//kl supervector
	        v.setSize(clientMixture.getDistribCount()*clientMixture.getVectSize(),true);		
		getKLVector(clientMixture,v,config);
	}
	else throw Exception("Cannot find supervector mode [kl|svmubm]",__FILE__,__LINE__);
}

#endif

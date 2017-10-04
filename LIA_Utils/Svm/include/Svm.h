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

#if !defined(ALIZE_Svm_h)
#define ALIZE_Svm_h

#include "alize.h"
#include "libsvm.h" 
using namespace alize;
using namespace std;

int svmTrain(alize:: Config &);
int svmPredict(alize:: Config &);
int svmPredictTnorm(alize:: Config &);	

//
// svm_model redeclaration as it is not in svm.h (bad)
//
struct svm_model
{
	svm_parameter param;	// parameter
	int nr_class;		// number of classes, = 2 in regression/one class svm
	int l;			// total #SV
	svm_node **SV;		// SVs (SV[l])
	double **sv_coef;	// coefficients for SVs in decision functions (sv_coef[k-1][l])
	double *rho;		// constants in decision functions (rho[k*(k-1)/2])
	double *probA;          // pariwise probability information
	double *probB;

	// for classification only

	int *label;		// label of each class (label[k])
	int *nSV;		// number of SVs for each class (nSV[k])
				// nSV[0] + nSV[1] + ... + nSV[k-1] = l
	// XXX
	int free_sv;		// 1 if svm_model is created by svm_load_model
				// 0 if svm_model is created by svm_train
};

class Instance { // to integrate, this is clean!
	private:
	Matrix <double> _ex;
	RealVector <unsigned long> _labels;
	XLine _id;
	public:
	Instance(unsigned long nbEx,unsigned long dimension) {
		_ex.setDimensions(nbEx,dimension);}
	double value(unsigned long i,unsigned long j) {
		return _ex(i,j);}
	unsigned long dimension() {
		return _ex.cols();}
	unsigned long nbEx() {
		return _ex.rows();}
	String& getId(unsigned long idx) {
		return _id.getElement(idx);}
	unsigned long getIdx(String & name) {
                for (unsigned long i=0;i<_id.getElementCount();i++)
                        if (_id.getElement(i)==name) return i;
                        else return -1;
				return -1;
        }
	void setId(String& name,unsigned long idx){
		_id.getElement(idx)=name;}
	void setLabels(RealVector <unsigned long>& labels) {
		_labels=labels;}
	void setLabel(unsigned long t,unsigned long idx) {
		_labels[idx]=t;}
	unsigned long getLabel(unsigned long idx){
		return _labels[idx];}
	unsigned long getLabel(String& name) {
		return this->getLabel(getIdx(name));}
	void setLabel(String & name,unsigned long t) {
		this->setLabel(t,getIdx(name));}
};


#endif // !defined(Svm)

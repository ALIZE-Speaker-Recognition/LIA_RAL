#if !defined(ALIZE_Svm_h)
#define ALIZE_Svm_h

#include "alize.h"
#include "svm.h" 
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

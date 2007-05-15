/****************************************
decembre 2004
*****************************************/
#if !defined(ALIZE_Hmm_h)
#define ALIZE_Hmm_h

#include <alize.h>


using namespace alize;
using namespace std;

class hmm{
	MixtureServer *ms;
	Config *conf;
	ObjectRefVector tabState;
	XLine tabStateName;
	DoubleVector transitions;
	
	
	public:
	unsigned long getNbState();
	unsigned long LoadState(String);
	unsigned long LoadState(String,String);
	unsigned long LoadState(MixtureGD&,String);	
	unsigned long addState(String);
	unsigned long deleteState(unsigned long indice);
	double getTransition(int i,int j);
	void getTransition(DoubleVector &);
	void setTransition(double,int , int );
	void reset();
	const hmm& operator=(const hmm& hmmc);
	MixtureGD& getDensity(unsigned long);
	void setDensity(MixtureGD& m, unsigned long nModel);
	void setDensity(String, unsigned long);
	const String &getStateName(unsigned long);
	void setStateName(unsigned long, String);
	hmm(MixtureServer &, Config &);
	hmm(const hmm&);
	~hmm();
	
	private:
	void assign(const hmm& hmmc);
	void save();
	void load();
};

// Return a copy of the transition matrix - should be included into hmm class - TODO
//DoubleVector &copyTransition(hmm& cHmm);
#endif

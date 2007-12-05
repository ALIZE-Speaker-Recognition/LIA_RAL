#if !defined(ALIZE_TopGauss_h)
#define ALIZE_TopGauss_h

#include <alize.h>
#include "liatools.h"

using namespace alize;
using namespace std;

class TopGauss {
	private:	
		RealVector <unsigned long> _idx; // gaussian idxs
		RealVector <unsigned long>_nbg; // nb gauss per frame
		RealVector <double>_snsw; // resting weight
		RealVector <double>_snsl; // resting lk
	
		unsigned long _nt; //nb frames
		unsigned long _nbgcnt; //nb gauss per file
	
	public:
		TopGauss() {
		};
		~TopGauss() {
		};
		double compute(MixtureGD &,FeatureServer &,String &,Config &);
		double get(MixtureGD &,FeatureServer &,String &,Config &);
		
		double compute(MixtureGD &,String &,Config &);
		RealVector <double> compute(MixtureGD &,FeatureServer &,Config &);
		
		double get(MixtureGD &,String &,Config &);
		RealVector <double> get(MixtureGD &,FeatureServer &,Config &);		
		
		void read(String &,Config &);
		void write(String &,Config &);
		
		RealVector <unsigned long> & nbg() {
			return _nbg;}
		RealVector <unsigned long> & idx() {
			return _idx;}			
		unsigned long frameToIdx(unsigned long &);
		unsigned long nbgcnt() {
			return _nbgcnt;}			
		unsigned long nt() {
			return _nt;}			
};
	
#endif
	

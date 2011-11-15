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

#if !defined(ALIZE_FactorAnalysis_h)
#define ALIZE_FactorAnalysis_h

#if defined(_WIN32)
#if defined(LIA_SPKTOOLS_EXPORTS)
#define LIA_SPKTOOLS_API __declspec(dllexport)
#else
#define LIA_SPKTOOLS_API __declspec(dllimport)
#endif
#else
#define LIA_SPKTOOLS_API
#endif

#include <liatools.h>

class TopGauss;

/// Translate filename into loc and sent
class LIA_SPKTOOLS_API Translate {
	XLine _id;
	ULongVector _locidx;
  	ULongVector _sessionidx;
	
	/// idx in XLine
  	long _idxOfID(const String &file) {
	_id.rewind();
  	for(unsigned long i=0;i<_id.getElementCount();i++)
		if(file==_id.getElement(i)) return i;
  	return -1;
  	}
	
	public:
	Translate() {};
		
	/// Fill structure with id and idxs		
	Translate(XList &fileList) {
		fileList.rewind(); XLine *pline;String *pFile;
		unsigned long sent=0;
		unsigned long loc=0;
		while((pline=fileList.getLine())!=NULL) { 
			while((pFile=pline->getElement())!=NULL) {
				_id.addElement(*pFile);
				_locidx.addValue(loc);
				_sessionidx.addValue(sent);
				sent++;
			}
			loc++;
		}
	};
	~Translate() {};
		
	/// Get speaker nb
	unsigned long locNb(const String &file) {
		long idx=_idxOfID(file);
		if (idx==-1) throw Exception("File not known",__FILE__,__LINE__);
	return _locidx[idx];
	}
	
	/// Get seesion nb
	unsigned long sessionNb(const String &file) {
		long idx=_idxOfID(file);
		if (idx==-1) throw Exception("File not known",__FILE__,__LINE__);
	return _sessionidx[idx];    
	}
};

class LIA_SPKTOOLS_API FactorAnalysisStat {
  private :
  	XList fileList; 
  	unsigned long _vsize;
	unsigned long _mixsize;
	unsigned long _supervsize;
	unsigned long _rang;
	unsigned long  _nb_speakers;
	unsigned long  _nb_sent;
	/// Regulation factor
	double _tau; 

	bool _topGauss;
  
    MixtureServer _ms;
	StatServer _ss;
  
	Matrix <double> _matX;
	Matrix <double> _matY;
	RealVector <double> _D;
	Matrix <double> _matU;
  
 	/// Accumulators 
	Matrix<double> _matS_X;
	Matrix<double> _matS_X_h;
	Matrix<double> _matN;
	Matrix<double> _matN_h;
	
	/// Copy of accumulators
  	Matrix<double> matcS_X;
	Matrix<double> matcS_X_h;
	Matrix<double> matcN;
	Matrix<double> matcN_h;	
	
	/// Inverse of L matrices
	RefVector <DoubleSquareMatrix> _l_h_inv;
	
	/// UBM SVs
	RealVector<double> _super_mean;
	RealVector<double> _super_invvar;
	
	/// Main functions to initialize values
	void _init(XList &,FeatureServer &,Config &);
	
  public :
	/// Association between filename and indexes of speaker and sessions
	Translate _ndxTable;
  
  	/// Constructors with XList of a string
  	FactorAnalysisStat(XList &,FeatureServer&,Config & config); 
  	FactorAnalysisStat(String &,FeatureServer&,Config & config);
	~FactorAnalysisStat() {
		_l_h_inv.deleteAllObjects();
	}
	/// Estimate General FA statistics N,S,F
	void computeAndAccumulateGeneralFAStats(FeatureServer&,Config & config);
	void computeAndAccumulateGeneralFAStats(SegCluster&,FeatureServer&,Config & config);	
	
	/// Estimate Channel factors loading
	void getXEstimate(); 
	
	/// Estimate Channel factors loading	
	void getYEstimate();
	
	/// Estimate SubSpace
	void getUEstimate(Config &); 
	void getUEstimateUnThreaded(); 
	#if defined(THREAD)
	void getUEstimateThreaded(unsigned long);// threaded version of estimate U	
	#endif	
	
	/// Cholesky inversion of L matrices
	void estimateAndInverseL(Config &);
	void estimateAndInverseLUnThreaded();
	#if defined(THREAD)	
	void estimateAndInverseLThreaded(unsigned long);// threaded version of inverse L	
	#endif
	
	/// Log Likelihood computations
	//void getLLK(Config & config);
	double getLLK(SegCluster &,MixtureGD &,FeatureServer &,Config &); 
	
	/// Substract stats to get X and Y estimate
	void substractSpeakerStats(); 
	void substractChannelStats();
	
	/// Feature normalization through smooth mixture transformation
	void normalizeFeatures(FeatureServer&,Config &);	
	void normalizeFeatures(SegCluster &selectedSegments,FeatureServer &fs,Config & config);	
	//void normalizeFeatures(TopGauss &tg,SegCluster &selectedSegments,FeatureServer &fs,MixtureGD &model,RealVector<double> m_xh_1,Config & config);	

	/// High Level functions, estimate X and Y and normalize features on a segcluster of featureServer
	void estimateXYAndNorm(FeatureServer &,Config &);	
	void estimateXYAndNorm(SegCluster&,FeatureServer &,Config &);

	/// Compute Supervectors for a file, either Ux or M+DY
	void getUX(RealVector <double> &,String&); 	
	void getMplusDY(RealVector <double> &,String&); 		

	/// Return models 
	void setYFromModel(MixtureGD&,String&); 
	void estimateXForKnownSpeaker(MixtureGD &,String&,Config &);
	void getTrueSpeakerModel(MixtureGD &,String&); 
	void getSpeakerModel(MixtureGD &,String&);
	void getSessionModel(MixtureGD &,String&);	
	void getFactorAnalysisModel(MixtureGD&,String&);
	
	/// Reset Accumulators
	void resetAcc() {
		_matS_X.setAllValues(0.0);	
		_matS_X_h.setAllValues(0.0);	
		_matN_h.setAllValues(0.0);
		_matN.setAllValues(0.0);	
		if (verbose) cout << "# FA Accumulators reset" << endl;
	}
	void resetXY() {
		_matY.setAllValues(0.0);
		_matX.setAllValues(0.0);			
	}
	
	/// Store Accumulators in temporary variables
	void storeAccs() {
		matcS_X=_matS_X;
		matcS_X_h=_matS_X_h;
		matcN_h=_matN_h; 
		matcN=_matN;
		if (verbose) cout << "(FactorAnalysisStat) FA Accs states stored" << endl;		
	}	
	
	/// Restore Accumulators in temporary variables
	void restoreAccs() {	
		_matS_X=matcS_X;
		_matS_X_h=matcS_X_h;
		_matN_h=matcN_h;
		_matN=matcN;	
		if (verbose) cout << "(FactorAnalysisStat) FA Accs states restored" << endl;			
	}	
	
	/// Save Accumulators on disk
	void saveAccs(Config &config) {
		_matS_X.save("S_X.mat",config);
		_matS_X_h.save("S_X_h.mat",config);
		_matN_h.save("N_h.mat",config); 
		_matN.save("N.mat",config);
		if (verbose) cout << "(FactorAnalysisStat) FA Accs states saved" << endl;		
	}		
	
	/// Load Accumulators from disk
	void loadAccs(Config &config) {	
		_matS_X.load("S_X.mat",config);
		_matS_X_h.load("S_X_h.mat",config);
		_matN_h.load("N_h.mat",config); 
		_matN.load("N.mat",config);
		if (verbose) cout << "(FactorAnalysisStat) FA Accs states loaded" << endl;		
	}
	
	/// Accessors
  	MixtureServer & getMixtureServer() {
		return this->_ms;}
  	StatServer & getStatServer() {
		return this->_ss;}
		XList & getFileList() {
		return this->fileList;}	
	unsigned long &getNbEx() {
		return _nb_sent;}	
	Matrix <double>& getU() {
		return _matU;}		
	Matrix <double>& getY() {
		return _matY;}	
	Matrix <double>& getX() {
		return _matX;}	
	RealVector <double>& getD() {
		return _D;}			
	unsigned long getRank() {
		return _rang;}			
	void setU(Matrix <double>& U) {
		_matU=U;}	
	void setX(Matrix <double>& X) {
		_matX=X;}	
	void setY(Matrix <double>& Y) {
		_matY=Y;}	
};

#endif

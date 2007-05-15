#if !defined(ALIZE_FactorAnalysis_h)
#define ALIZE_FactorAnalysis_h

#include <liatools.h>

class TopGauss;

/// Translate filename into loc and sent
class Translate {
	XLine _id;
	RealVector <unsigned long> _locidx;
  	RealVector <unsigned long> _sessionidx;
	
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

class FactorAnalysisStat {
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
  
	Matrix <double> _X;
	Matrix <double> _Y;
	RealVector <double> _D;
	Matrix <double> _U;
  
 	/// Accumulators 
	Matrix<double> _S_X;
	Matrix<double> _S_X_h;
	Matrix<double> _N;
	Matrix<double> _N_h;
	
	/// Copy of accumulators
  	Matrix<double> cS_X;
	Matrix<double> cS_X_h;
	Matrix<double> cN;
	Matrix<double> cN_h;	
	
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
	#ifdef THREAD
	void getUEstimateThreaded(unsigned long);// threaded version of estimate U	
	#endif	
	
	/// Cholesky inversion of L matrices
	void estimateAndInverseL(Config &);
	void estimateAndInverseLUnThreaded();
	#ifdef THREAD	
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
		_S_X.setAllValues(0.0);	
		_S_X_h.setAllValues(0.0);	
		_N_h.setAllValues(0.0);
		_N.setAllValues(0.0);	
		if (verbose) cout << "# FA Accumulators reset" << endl;
	}
	void resetXY() {
		_Y.setAllValues(0.0);
		_X.setAllValues(0.0);			
	}
	
	/// Store Accumulators in temporary variables
	void storeAccs() {
		cS_X=_S_X;
		cS_X_h=_S_X_h;
		cN_h=_N_h; 
		cN=_N;
		if (verbose) cout << "(FactorAnalysisStat) FA Accs states stored" << endl;		
	}	
	
	/// Restore Accumulators in temporary variables
	void restoreAccs() {	
		_S_X=cS_X;
		_S_X_h=cS_X_h;
		_N_h=cN_h;
		_N=cN;	
		if (verbose) cout << "(FactorAnalysisStat) FA Accs states restored" << endl;			
	}	
	
	/// Save Accumulators on disk
	void saveAccs(Config &config) {
		_S_X.save("S_X.mat",config);
		_S_X_h.save("S_X_h.mat",config);
		_N_h.save("N_h.mat",config); 
		_N.save("N.mat",config);
		if (verbose) cout << "(FactorAnalysisStat) FA Accs states saved" << endl;		
	}		
	
	/// Load Accumulators from disk
	void loadAccs(Config &config) {	
		_S_X.load("S_X.mat",config);
		_S_X_h.load("S_X_h.mat",config);
		_N_h.load("N_h.mat",config); 
		_N.load("N.mat",config);
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
		return _U;}		
	Matrix <double>& getY() {
		return _Y;}	
	Matrix <double>& getX() {
		return _X;}	
	RealVector <double>& getD() {
		return _D;}			
	unsigned long getRank() {
		return _rang;}			
	void setU(Matrix <double>& U) {
		_U=U;}	
	void setX(Matrix <double>& X) {
		_X=X;}	
	void setY(Matrix <double>& Y) {
		_Y=Y;}	
};

#endif

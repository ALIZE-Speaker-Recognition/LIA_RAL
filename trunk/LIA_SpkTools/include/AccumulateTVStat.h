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

#if !defined(ALIZE_AccumulateTVStat_h)
#define ALIZE_AccumulateTVStat_h

#if defined(_WIN32)
#if defined(LIA_SPKTOOLS_EXPORTS)
#define LIA_SPKTOOLS_API __declspec(dllexport)
#else
#define LIA_SPKTOOLS_API __declspec(dllimport)
#endif
#else
#define LIA_SPKTOOLS_API
#endif

#include <alize.h>
#include "liatools.h"

  /// This class represents an index of speakers and sessions
  /// A Transate contains the indexes of speakers and sessions
  /// 
  /// @author Anthony Larcher  alarcher@i2r.a-star.edu.sg
  /// @version 1.0
  /// @date 2012

class LIA_SPKTOOLS_API TVTranslate{
	XLine _id;
	ULongVector _locidx;
	ULongVector _sessionidx;
	RefVector<ULongVector> _locindices;
	RefVector<ULongVector> _sessindices;
	
	// idx in XLine
  	long _idxOfID(const String &file) {
		_id.rewind();
  		for(unsigned long i=0;i<_id.getElementCount();i++)
			if(file==_id.getElement(i)) return i;
  		return -1;
  	}
	
	
	public:
	TVTranslate() {};

	/// Fill structure with id and idxs	
	///
	TVTranslate(XList &fileList) {
		fileList.rewind(); XLine *pline;String *pFile;
		unsigned long sent=0;
		unsigned long loc=0;
		while((pline=fileList.getLine())!=NULL) {		// for each line of the NDX file
			while((pFile=pline->getElement())!=NULL) {	// for each file of the current line

				// if the file does not exist yet
				if(_id.getIndex(*pFile) == -1){
					_id.addElement(*pFile);
					
					ULongVector tmp(0,0);
					_locindices.addObject(*new ULongVector(0,0));
					_sessindices.addObject(*new ULongVector(0,0));
				}
				_sessionidx.addValue(sent);
				_sessindices.getObject(_id.getIndex(*pFile)).addValue(sent);

				_locidx.addValue(loc);
				_locindices.getObject(_id.getIndex(*pFile)).addValue(loc);

				sent++;
			}
			loc++;
		}
	};

	~TVTranslate() {};
		
	/// Return indices of speakers for a file
	ULongVector& locIndices(const String &file){
		unsigned long idx=_idxOfID(file);
		if (idx==-1) throw Exception("File not known",__FILE__,__LINE__);
		return _locindices[idx];
	}

	/// Return indices of sessions for a file
	ULongVector& sessIndices(const String &file){
		unsigned long idx=_idxOfID(file);
		if (idx==-1) throw Exception("File not known",__FILE__,__LINE__);
		return _sessindices[idx];
	}

	/// Return the idx of speakers
	/// @param file name of the feature file
	/// @return the idx of the given speaker
	///
	unsigned long locNb(const String &file) {
		unsigned long idx=_idxOfID(file);
		if (idx==-1) throw Exception("File not known",__FILE__,__LINE__);
		return _locidx[idx];
	}

	/// Return the idx of the first session of a speaker
	/// @param spk idx of the targeted speaker
	/// @return the idx of the first session of a speaker spk
	///
	unsigned long firstSession(unsigned long spk) {
		unsigned long i;
		for(i=0;;i++) if (_locidx[i]==spk) break;
		return (_sessionidx[i]);
	}
		
	/// Return the idx of a given session of a speaker
	/// @param spk idx of the speaker
	/// @param session number of the session for the given speaker
	/// @return the idx of the session session of the speaker spk
	///
	unsigned long numSession(unsigned long spk, unsigned long session) {
		unsigned long i, cnt;
		if(spk==0) return(session);
		for(i=session-1,cnt=0;;i--,cnt++) if (_locidx[i]!=spk) break;
		return (cnt);
		}

	/// Return the idx of speakers
	/// @param index of the feature file
	/// @return the idx of the given speaker
	///
	unsigned long locNb(unsigned long idx) {
		return _locidx[idx];
	}
	
	/// Return the idx of a sessions
	/// @param file feature file
	/// @return the idx of the session corresponding to the file
	///
	unsigned long sessionNb(const String &file) {
		long idx=_idxOfID(file);
		if (idx==-1) throw Exception("File not known",__FILE__,__LINE__);
	return _sessionidx[idx];
	}
};

  /// This class represents a accumulator of statistics. 
  /// A TVAcc contains the accumulators needed for TotalVariability
  /// estimation
  /// 
  /// @author Anthony Larcher  alarcher@i2r.a-star.edu.sg
  /// @version 1.0
  /// @date 2012

class LIA_SPKTOOLS_API TVAcc{

	private :

		MixtureServer _ms;
		StatServer _ss;
	
		XList _fileList; 
		unsigned long _vectSize;
		unsigned long _n_distrib;
		unsigned long _svSize;
		unsigned long _rankT;
		unsigned long _n_speakers;
		unsigned long _n_sessions;
	
		/// UBM SVs
		RealVector<double> _ubm_means;
		RealVector<double> _ubm_invvar;

		RealVector<double> _meanW;

		/// Accumulators 
		Matrix<double> _statN;
		Matrix<double> _statF;

		/// TotalVariability Objects
		Matrix<double> _T;
		Matrix<double> _W;
		
		///Accumulators for VSigmaVT
		RefVector<DoubleSquareMatrix> _TETt;
		
		///Accumulators
		Matrix<double> _A;
		Matrix<double> _C;

		///Minimum Divergence Accumulators
		DoubleSquareMatrix _R;
		DoubleVector _r;
		
		/// main fonction to initialise the accumulator
		/// @param list the XList containing speakers and sessions filenames
		/// @param config config filename
		///
		void _init(XList & list,Config & config);

		/// fonction to initialise the accumulator without loading data
		/// @param list the XList containing speakers and sessions filenames
		/// @param config config filename
		///
		void _init(Config & config);


	public :

		///Association between filename and indexes of speaker and sessions
		TVTranslate _ndxTable;
	
 		TVAcc(String &,Config &);
		TVAcc(XList &,Config &);
		TVAcc(Config &);

		//TVAcc(String &,Config &,String);
		//TVAcc(XList &,Config &, String);

		virtual ~TVAcc();
		virtual String getClassName() const;
	
		/// main fonction to accumulate statistics
		/// @param config config filename
		///
		void computeAndAccumulateTVStat(Config& config);
	
		/// Unthreaded fonction to accumulate statistics
		/// @param config config filename
		///
		void computeAndAccumulateTVStatUnThreaded(Config& config);
	
		/// Threaded fonction to accumulate statistics
		/// @param config config filename
		///
		#ifdef THREAD
			/// Threaded fonction to accumulate statistics
			/// @param config config filename
			///
			void computeAndAccumulateTVStatThreaded(unsigned long numThread, Config& config);
		#endif
	
		/// Fonction to accumulate statistics
		/// @param fs featureserver containing the required data
		/// @param config config filename
		///	
		void computeAndAccumulateTVStat(FeatureServer & fs,Config & config);
		
		/// Fonction to accumulate statistics
		/// @param sc SegCluster containing the selected segments
		/// @param fs featureserver containing the required data
		/// @param config config filename
		///
		void computeAndAccumulateTVStat(SegCluster & sc,FeatureServer & fs,Config & config);

		/// Function to reset the sufficient statistics accumulators
		///
		void resetAcc();

		/// Reset the tempporary variables used for matrix estimations
		///
		void resetTmpAcc();

		/// Load the Total Variability Matrix from a file
		/// @param file name of the matrix file
		/// @param config config filename
		///
		void loadT(const String& file, Config& config);

		/// Load the Total Variability Matrix from an existing matrix
		/// @param mat matrix object
		///
		void loadT(Matrix<double> & mat, Config& config);

		/// Load the Mean Estimate from Minimum Divergence
		/// @param meanEstimate DoubleVector object
		///
		void loadMeanEstimate(DoubleVector& meanEstimate);

		/// Load the null order statistics matrix of speakers
		/// @param config config filename
		///
		void loadN(Config&);

		/// Load the first order statistics matrix of speakers
		/// @param config config filename
		///
		void loadF_X(Config&);

		/// Initialise the Total Variability matrix by a Box-Muller random process
		///
		void initT(Config&);

		/// Save the Total Variability matrix
		/// @param file name of the matrix file to save
		///
		void saveT(const String& file, Config& config);
	
		/// Compute the VEVt matrices
		/// @param config config filename
		///
		void estimateTETt(Config & config);

		/// Compute the VEVt matrices without threads
		///
		void estimateTETtUnThreaded();

		#ifdef THREAD
			/// Compute the VEVt matrices using multi-threading
			/// @param threads number of threads to launch
			///
			void estimateTETtThreaded(unsigned long threads);
		#endif
		
		/// Get the VY supervector for a feature file
		/// @param vy DoubleVector to fill with the supervector
		/// @param name of the feature file
		///
		void getTW(DoubleVector & vy, String & file);

		/// Get the supervector for a given speaker
		/// @param vy DoubleVector to fill with the supervector
		/// @param spk index of the speaker
		///
		void getTW(DoubleVector &file, unsigned long spk);

		/// Get the MplusVY supervector for a given featurefile
		/// @param MplusVY DoubleVector to fill with the supervector
		/// @param file name of the feature file
		///
		void getMplusTW(DoubleVector &MplusVY, String& file);

		/// Get the MplusVY supervector for a given speaker
		/// @param MplusVY DoubleVector to fill with the supervector
		/// @param spk index of the speaker
		///
		void getMplusTW(DoubleVector &MplusVY, unsigned long spk);
		
		/// Compute the inverse of L matrices for TotalVariability estimation, estimate TotalVariability matrix and Y components
		/// @param config config filename
		///
		void estimateAandC(Config &config);

		/// Compute the inverse of L matrices for TotalVariability estimation, estimate TotalVariability matrix and Y components without threads
		/// @param config config filename
		///
		void estimateAandCUnthreaded(Config &config);

		#ifdef THREAD
			/// Compute the inverse of L matrices for TotalVariability estimation
			/// estimate TotalVariability matrix and Y components using multi-threads
			///
			void estimateAandCThreaded(unsigned long threads);
		#endif

		/// Estimate W components
		/// @param config config filename
		/// 
		void estimateW(Config &);

		/// Estimate W components without threads
		///
		void estimateWUnThreaded(Config&);

		#ifdef THREAD
			/// Estimate Y components using multi-threads
			/// @param threads number of threads to launch
			/// 
			void estimateWThreaded(unsigned long threads);
		#endif

		/// Approximate W components by using UBM weight
		/// @param config config filename
		/// 
		void estimateWUbmWeight(DoubleSquareMatrix &W, Config &);

		/// Approximate W components without threads
		///
		void estimateWUbmWeightUnThreaded(DoubleSquareMatrix &W, Config&);

		#ifdef THREAD
			/// Approximate W components using multi-threads
			/// @param threads number of threads to launch
			/// 
			void estimateWUbmWeightThreaded(DoubleSquareMatrix &W, unsigned long threads);
		#endif

		void estimateWEigenDecomposition(Matrix<double> D, Matrix<double> Q, Config& config);

		void estimateWEigenDecompositionUnThreaded(Matrix<double> D, Matrix<double> Q, Config& config);

		#ifdef THREAD
			/// Approximate W components using multi-threads
			/// @param threads number of threads to launch
			/// 
			void estimateWEigenDecompositionThreaded(Matrix<double> &D, Matrix<double> &Q, unsigned long threads);
		#endif

		/// Update the Total Variability Matrix
		///
		void updateTestimate();

		/// Return the XList
		///
		XList& getXList();

		/// Return the number of speakers of the TVAcc
		///
		unsigned long getNSpeakers();

		/// Return the number of distribution of the Universal Backgroud Model
		///
		unsigned long getNDistrib();

		/// Return the dimension of feature vectors
		///
		unsigned long getVectSize();

		/// Return the dimension of the supervectors
		///
		unsigned long getSvSize();

		/// Return the rank of the Total Variability Matrix
		///
		unsigned long getRankT();

		/// Return the Supervector of means of the UBM
		///
		DoubleVector& getUbmMeans();

		/// Return the Supervector of variances of the UBM
		///
		DoubleVector& getUbmInvVar();
		
		/// Return the Total Variability Matrix
		///
		Matrix<double> getT();

		/// Return the W matrix
		///
		Matrix<double> getW();

		/// Save matrix W
		/// @param yFilename filename to save the matrix
		///
		void saveW(String yFilename,Config &);

		/// Store the sufficient statistics accumulators for backup
		///
		void storeAccs();

		/// Restore the sufficient statistics accumulators
		///
		void restoreAccs();
		
		/// Subtract the mean component M from the statistics
		/// @param config config filename
		///
		void substractM(Config & config);

		/// Subtract the mean component M from the statistics without threads
		///
		void substractMUnThreaded();

		/// Subtract the mean component M from the statistics and whiten
		/// @param config config filename
		void normStatistics(Config &config);

		/// Subtract the mean component M from the statistics and whiten without threads
		///
		void normStatisticsUnThreaded();

		#ifdef THREAD
		/// Substract the mean component M from the speaker statistics and normalize using covariance using multi-threading
		/// @param threads number of threads to launch
		///
		void normStatisticsThreaded(unsigned long threads);
		#endif

		/// Substract the mean component M from the speaker statistics
		/// @param config config filename
		///
		void substractMplusTW(Config &);

		/// Substract the mean component M from the speaker statistics without threads
		///
		void substractMplusTWUnThreaded();

		#ifdef THREAD
		/// Substract the mean component M from the speaker statistics using multi-threading
		/// @param threads number of threads to launch
		///
		void substractMThreaded(unsigned long threads);

		/// Substract the speaker component MplusTV from the speaker statistics using multi-threading
		/// @param threads number of threads to launch
		///
		void substractMplusTWThreaded(unsigned long threads);
		#endif

		/// Return the speaker model ( M + TW ) for a given feature file
		/// @param mixture mixture to fill with the speaker model
		/// @param file name of the feature file
		///
		void getSpeakerModel(MixtureGD &mixture, String& file);

		/// Orthonormalize the Total Variability matrix
		///
		void orthonormalizeT();

		// Normalize the Total Variability matrix by using the square root of UBM inverse co-varariance
		void normTMatrix();

		///Save Accumulators on disk
		/// @param config config filename
		///
		void saveAccs(Config &config);
		
		///Compute LLK
		/// @param selectedSegments segment cluster to process
		/// @param model GMM model
		/// @param config config filename
		///
		double getLLK(SegCluster &selectedSegments,MixtureGD &model,FeatureServer&fs,Config & config);

		///Compute LLK
		/// @param selectedSegments segment cluster to process
		/// @param model GMM model
		/// @param config config filename
		///
		void verifyEMLK(Config& config);

		/// Return the Null order speaker statistic Matrix
		///
		Matrix <double>& getN();

		/// Return the First order speaker statistic Matrix
		///
		Matrix <double>& getF();


		/// Return a temporary DoubleSquareMatrix TETt
		/// @param idx index of the matrix to return
		///
		DoubleSquareMatrix& getTETt(unsigned long idx);

		///Update the Matrix according to the Minimum Divergence Criteria
		///
		void minDivergence();

		/// Save i-vector to a file
		/// @param config config filename
		///
		void saveWbyFile(Config &config);

		/// Compute Weighted covariance matrix from the TotalVariability matrix
		/// @param W the output weighted covariance matrix
		///
		void getWeightedCov(DoubleSquareMatrix &W, DoubleVector& weight, Config &config);

		/// Compute Weighted covariance matrix from the TotalVariability matrix without multi-threading
		/// @param W the output weighted covariance matrix
		///
		void getWeightedCovUnThreaded(DoubleSquareMatrix &W, DoubleVector& weight);

		#ifdef THREAD
		/// Compute Weighted covariance matrix from the TotalVariability matrix using multi-threading
		/// @param W the output weighted covariance matrix
		/// @param threads number of threads to launch
		///
		void getWeightedCovThreaded(DoubleSquareMatrix &W, DoubleVector& weight, unsigned long threads);
		#endif

		void computeEigenProblem(Matrix<double> &EP,Matrix<double> &eigenVect, long rank, Config& config);
		void computeEigenProblem(Matrix<double> &EP,Matrix<double> &eigenVect,Matrix<double> &eigenVal, long rank, Config& config);

		void approximateTcTc(Matrix<double> &D, Matrix<double> &Q, Config &config);
		void approximateTcTcUnThreaded(Matrix<double> &D, Matrix<double> &Q);

		#ifdef THREAD

		///
		void approximateTcTcThreaded(Matrix<double> &D, Matrix<double> &Q, unsigned long threads);
		#endif
};
#endif

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
	~TVTranslate() {};
		
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
		unsigned long _rankEV;
		unsigned long _n_speakers;
		unsigned long _n_sessions;
	
		/// UBM SVs
		RealVector<double> _ubm_means;
		RealVector<double> _ubm_invvar;

		RealVector<double> _meanY;

		/// Accumulators 
		Matrix<double> _matN;
		Matrix<double> _F_X;
		
		/// Copy of accumulators
		Matrix<double> _cN;
		Matrix<double> _cF_X;

		/// JFA Objects
		Matrix<double> _V;

		Matrix<double> _Y;
		
		///Accumulators for VSigmaVT
		RefVector<DoubleSquareMatrix> _vEvT;
		
		///Accumulators
		RefVector<DoubleSquareMatrix> _Aev;
		Matrix<double> _Cev;

		///Minimum Divergence Accumulators
		DoubleSquareMatrix _R;
		DoubleVector _r;
		
		/// main fonction to initialise the accumulator
		/// @param list the XList containing speakers and sessions filenames
		/// @param config config filename
		///
		void _init(XList & list,Config & config);


	public :

		///Association between filename and indexes of speaker and sessions
		TVTranslate _ndxTable;
	
 		TVAcc(String &,Config &);
		TVAcc(XList &,Config & );

		TVAcc(String &,Config &,String);
		TVAcc(XList &,Config &, String);


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
			/// Threaded fonction to accumulate statistics (STILL TO IMPLEMENTE)
			/// @param config config filename
			///
			void computeAndAccumulateTVStatThreaded(Config&);
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

		/// Load the Eigenvoice Matrix from a file
		/// @param file name of the matrix file
		/// @param config config filename
		///
		void loadEV(const String& file, Config& config);

		/// Load the Eigenvoice Matrix from an existing matrix
		/// @param mat matrix object
		///
		void loadEV(Matrix<double> & mat, Config& config);

		/// Load the Mean Estimate from Minimum Divergence
		/// @param meanEstimate DoubleVector object
		///
		void loadMeanEstimate(DoubleVector& meanEstimate);

		/// Load the null order statistics matrix of speakers
		/// @param config config filename
		///
		void loadN(Config&);

		/// Load the null order statistics matrix of sessions
		/// @param config config filename
		///
		void loadN_h(Config&);

		/// Load the first order statistics matrix of speakers
		/// @param config config filename
		///
		void loadF_X(Config&);

		/// Load the first order statistics matrix of sessions
		/// @param config config filename
		///
		void loadF_X_h(Config&);

		/// Initialise the EigenVoice matrix by a Box-Muller random process
		///
		void initEV(Config&);

		/// Save the EigenVoice matrix
		/// @param file name of the matrix file to save
		///
		void saveV(const String& file, Config& config);
	
		/// Compute the VEVt matrices
		/// @param config config filename
		///
		void estimateVEVT(Config & config);

		/// Compute the VEVt matrices without threads
		///
		void estimateVEVTUnThreaded();

		#ifdef THREAD
			/// Compute the VEVt matrices using multi-threading
			/// @param threads number of threads to launch
			///
			void estimateVEVTThreaded(unsigned long threads);
		#endif
		
		/// Get the VY supervector for a feature file
		/// @param vy DoubleVector to fill with the supervector
		/// @param name of the feature file
		///
		void getVY(DoubleVector & vy, String & file);

		/// Get the supervector for a given speaker
		/// @param vy DoubleVector to fill with the supervector
		/// @param spk index of the speaker
		///
		void getVY(DoubleVector &file, unsigned long spk);

		/// Get the VYplusDZ supervector for a given featurefile
		/// @param vyplusdz DoubleVector to fill with the supervector
		/// @param file name of the feature file
		///
		void getVYplusDZ(DoubleVector &vyplusdz, String &file);

		/// Get the VYplusDZ supervector for a given speaker
		/// @param vyplusdz DoubleVector to fill with the supervector
		/// @param spk index of the speaker
		///
		void getVYplusDZ(DoubleVector &, unsigned long spk);

		/// Get the MplusVY supervector for a given featurefile
		/// @param MplusVY DoubleVector to fill with the supervector
		/// @param file name of the feature file
		///
		void getMplusVY(DoubleVector &MplusVY, String& file);

		/// Get the MplusVY supervector for a given speaker
		/// @param MplusVY DoubleVector to fill with the supervector
		/// @param spk index of the speaker
		///
		void getMplusVY(DoubleVector &MplusVY, unsigned long spk);

		/// Compute the inverse of L matrices for EigenVoice estimation
		/// @param config config filename
		///
		void estimateAndInverseL_EV(Config&);

		/// Compute the inverse of L matrices for EigenVoice estimation without threads
		///
		void estimateAndInverseLUnThreaded_EV(Config& config);
	
		/// Estimate EigenVoice matrix and Y components
		/// @param config config filename
		/// 
		void estimateYandV(Config &);

		/// Estimate EigenVoice matrix and Y components without threads
		/// 
		void estimateYandVUnThreaded(Config& config);
		
		#ifdef THREAD
			/// Compute the inverse of L matrices for EigenVoice estimation using multi-threads
			/// @param threads number of threads to launch
			///
			void estimateAndInverseLThreaded_EV(unsigned long threads,Config& config);

			/// Estimate EigenVoice matrix and Y components using multi-threads
			/// @param threads number of threads to launch
			/// 
			void estimateYandVThreaded(unsigned long threads, Config& config);
		#endif
		
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

		/// Estimate Y components
		/// @param config config filename
		/// 
		void estimateY(Config &);

		/// Estimate Y components without threads
		///
		void estimateYUnThreaded(Config&);

		#ifdef THREAD
			/// Estimate Y components using multi-threads
			/// @param threads number of threads to launch
			/// 
			void estimateYThreaded(unsigned long threads);
		#endif

		/// Update the EigenVoice Matrix
		///
		void updateVestimate();

		/// Return the XList
		///
		XList& getXList();

		/// Return the number of speakers of the JFAAcc
		///
		unsigned long getNSpeakers();

		/// Return the number of sessions of the JFAAcc
		///
		unsigned long getNSessions();

		/// Return the number of distribution of the Universal Backgroud Model
		///
		unsigned long getNDistrib();

		/// Return the dimension of feature vectors
		///
		unsigned long getVectSize();

		/// Return the dimension of the supervectors
		///
		unsigned long getSvSize();

		/// Return the rank of the EigenVoice Matrix
		///
		unsigned long getRankEV();

		/// Return the Supervector of means of the UBM
		///
		DoubleVector& getUbmMeans();

		/// Return the Supervector of variances of the UBM
		///
		DoubleVector& getUbmInvVar();
		
		/// Return the EigenVoice Matrix
		///
		Matrix<double> getV();

		/// Return the Y matrix
		///
		Matrix<double> getY();

		/// Save matrix Y
		/// @param yFilename filename to save the matrix
		///
		void saveY(String yFilename,Config &);

		/// Store the sufficient statistics accumulators for backup
		///
		void storeAccs();

		/// Restore the sufficient statistics accumulators
		///
		void restoreAccs();
		
		/// Substract the speaker component MplusDZ from the speaker statistics
		/// @param config config filename
		///
		void substractM(Config & config);

		/// Substract the speaker component MplusDZ from the speaker statistics without threads
		///
		void substractMUnThreaded();

		/// Substract the speaker component MplusVYplusDZ from the speaker statistics
		/// @param config config filename
		///
		void substractMplusVY(Config &);

		/// Substract the speaker component MplusVYplusDZ from the speaker statistics without threads
		///
		void substractMplusVYUnThreaded();

		#ifdef THREAD
		/// Substract the speaker component MplusDZ from the speaker statistics using multi-threading
		/// @param threads number of threads to launch
		///
		void substractMThreaded(unsigned long threads);

		/// Substract the speaker component MplusVY from the speaker statistics using multi-threading
		/// @param threads number of threads to launch
		///
		void substractMplusVYThreaded(unsigned long threads);
		#endif

		/// Return the speaker model ( M + VY + DZ ) for a given feature file
		/// @param mixture mixture to fill with the speaker model
		/// @param file name of the feature file
		///
		void getSpeakerModel(MixtureGD &mixture, String& file);

		/// Orthonormalize the EigenVoice matrix
		///
		void orthonormalizeV();

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

		/// Return the Null order sessions statistic Matrix
		///
		Matrix <double>& getN_h();

		/// Return the First order speaker statistic Matrix
		///
		Matrix <double>& getF_X();

		/// Return the First order sessions statistic Matrix
		///
		Matrix <double>& getF_X_h();

		///Update the Matrix according to the Minimum Divergence Criteria
		///
		void minDivergence();


		void saveYbyFile(Config &config);

};
#endif

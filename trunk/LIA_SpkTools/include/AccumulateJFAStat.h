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

#if !defined(ALIZE_AccumulateJFAStat_h)
#define ALIZE_AccumulateJFAStat_h

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
  /// @date 2010

class LIA_SPKTOOLS_API JFATranslate{
	XLine _id;
	ULongVector _locidx;
	ULongVector _sessionidx;
//	RealVector<unsigned long > _locidx;
//  RealVector<unsigned long> _sessionidx;
	
	/// idx in XLine
  	long _idxOfID(const String &file) {
	_id.rewind();
  	for(unsigned long i=0;i<_id.getElementCount();i++)
		if(file==_id.getElement(i)) return i;
  	return -1;
  	}
	
	
	public:
	JFATranslate() {};

	/// Fill structure with id and idxs	
	///
	JFATranslate(XList &fileList) {
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
	~JFATranslate() {};
		
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
  /// A JFAAcc contains the accumulators needed for Joint Factor Analysis
  /// estimation
  /// 
  /// @author Anthony Larcher  alarcher@i2r.a-star.edu.sg
  /// @version 1.0
  /// @date 2010

class LIA_SPKTOOLS_API JFAAcc{

	private :

		MixtureServer _ms;
		StatServer _ss;
	
		XList _fileList; 
		unsigned long _vectSize;
		unsigned long _n_distrib;
		unsigned long _svSize;
		unsigned long _rankEC;
		unsigned long _rankEV;
		unsigned long _n_speakers;
		unsigned long _n_sessions;
	
		/// UBM SVs
		RealVector<double> _ubm_means;
		RealVector<double> _ubm_invvar;

		/// Accumulators 
		Matrix<double> _matN;
		Matrix<double> _N_h;
		Matrix<double> _F_X;
		Matrix<double> _F_X_h;
		
		/// Copy of accumulators
		Matrix<double> _cN;
		Matrix<double> _cN_h;
		Matrix<double> _cF_X;
		Matrix<double> _cF_X_h;

		/// JFA Objects
		RealVector<double> _D;
		Matrix<double> _V;
		Matrix<double> _matU;
		Matrix<double> _VU;

		Matrix<double> _Z;
		Matrix<double> _Y;
		Matrix<double> _matX;
		Matrix<double> _YX;
		
		///Accumulators for VSigmaVT
		RefVector<DoubleSquareMatrix> _vEvT;
		///Accumulators for USigmaUT
		RefVector<DoubleSquareMatrix> _uEuT;
		///Accumulators for VUSigmaVUT
		RefVector<DoubleSquareMatrix> _vuEvuT;
		
		///Accumulators for L EigenVoices and EigenChannel matrices
		RefVector<DoubleSquareMatrix> _l_spk_inv;
		RefVector<DoubleSquareMatrix> _l_sess_inv;
		RefVector<DoubleSquareMatrix> _l_yx_inv;
		
		///Accumulators
		RefVector<DoubleSquareMatrix> _Aev;
		Matrix<double> _Cev;
		RefVector<DoubleSquareMatrix> _Aec;
		Matrix<double> _Cec;
		
		/// main fonction to initialise the accumulator
		/// @param list the XList containing speakers and sessions filenames
		/// @param config config filename
		///
		void _init(XList & list,Config & config);

		/// main fonction to initialise the accumulator for I-Vector processing
		/// @param list the XList containing speakers and sessions filenames
		/// @param config config filename
		///
		void _init(XList & list,Config & config, String task);



	public :

		///Association between filename and indexes of speaker and sessions
		JFATranslate _ndxTable;
	
		JFAAcc(String &,Config &);
		JFAAcc(XList &,Config & );

		JFAAcc(String &,Config &,String);
		JFAAcc(XList &,Config &, String);


		virtual ~JFAAcc();
		virtual String getClassName() const;
	
		/// main fonction to accumulate statistics
		/// @param config config filename
		///
		void computeAndAccumulateJFAStat(Config& config);
	
		/// Unthreaded fonction to accumulate statistics
		/// @param config config filename
		///
		void computeAndAccumulateJFAStatUnThreaded(Config& config);
	
		/// Threaded fonction to accumulate statistics
		/// @param config config filename
		///
		#ifdef THREAD
			/// Threaded fonction to accumulate statistics (STILL TO IMPLEMENTE)
			/// @param config config filename
			///
			void computeAndAccumulateJFAStatThreaded(unsigned long numThread, Config& config);
		#endif
	
		/// Fonction to accumulate statistics
		/// @param fs featureserver containing the required data
		/// @param config config filename
		///	
		void computeAndAccumulateJFAStat(FeatureServer & fs,Config & config);
		
		/// Fonction to accumulate statistics
		/// @param sc SegCluster containing the selected segments
		/// @param fs featureserver containing the required data
		/// @param config config filename
		///
		void computeAndAccumulateJFAStat(SegCluster & sc,FeatureServer & fs,Config & config);

		/// Function to reset the sufficient statistics accumulators
		///
		void resetAcc();

		/// Reset the tempporary variables used for matrix estimations
		///
		void resetTmpAcc();

		/// Reset the tempporary variables used for matrix estimations
		///
		void resetTmpAcc(String task);

		/// Load the Eigenvoice Matrix from a file
		/// @param file name of the matrix file
		/// @param config config filename
		///
		void loadEV(const String& file, Config& config);

		/// Load the Eigenvoice Matrix from an existing matrix
		/// @param mat matrix object
		///
		void loadEV(Matrix<double> & mat, Config& config);

		/// Load the EigenChannel Matrix from a file
		/// @param file name of the matrix file
		/// @param config config filename
		///
		void loadEC(const String&, Config&);

		/// Load the EigenChannel Matrix from an existing matrix
		/// @param mat matrix object
		///
		void loadEC(Matrix<double> & mat, Config& config);

		/// Load the D Matrix from a file
		/// @param file name of the matrix file
		/// @param config config filename
		///
		void loadD(const String& file, Config& config);

		/// Load the D Matrix from an existing DoubleVector
		/// @param vec matrix object
		///
		void loadD(DoubleVector & vec);

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

		/// Initialise the EigenChannel matrix by a Box-Muller random process
		///
		void initEC(Config&);

		// Initialise the D matrix using the MAP paradigm
		///
		void initD(Config&);

		/// Initialise the VU matrix, concatenate the V and U matrices
		///
		void initVU();

		/// Save the EigenVoice matrix
		/// @param file name of the matrix file to save
		///
		void saveV(const String& file, Config& config);

		/// Save the EigenChannel matrix
		/// @param file name of the matrix file to save
		///
		void saveU(const String& file, Config& config);

		/// Save the D matrix
		/// @param file name of the matrix file to save
		///
		void saveD(const String& file, Config& config);
	
		/// Compute the VEVt matrices
		/// @param config config filename
		///
		void estimateVEVT(Config & config);

		/// Compute the VEVt matrices without threads
		///
		void estimateVEVTUnThreaded();

		/// Compute the UEUt matrices
		/// @param config config filename
		///
		void estimateUEUT(Config & config);

		/// Compute the UEUt matrices without threads
		///
		void estimateUEUTUnThreaded();

		/// Compute the VUEVUt matrices
		/// @param config config filename
		///
		void estimateVUEVUT(Config &);

		/// Compute the VUEVUt matrices without threads
		///
		void estimateVUEVUTUnThreaded();

		#ifdef THREAD
			/// Compute the VEVt matrices using multi-threading
			/// @param threads number of threads to launch
			///
			void estimateVEVTThreaded(unsigned long threads);

			/// Compute the UEUt matrices using multi-threading
			/// @param threads number of threads to launch
			///
			void estimateUEUTThreaded(unsigned long threads);

			/// Compute the VUEVUt matrices using multi-threading
			/// @param threads number of threads to launch
			///
			void estimateVUEVUTThreaded(unsigned long threads);
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

		/// Get the VUYX supervector for a feature file
		/// @param vy DoubleVector to fill with the supervector
		/// @param name of the feature file
		///
		void getVUYX(DoubleVector &vuyx, String & file);

		/// Get the VUYX supervector for a given speaker
		/// @param vuyx DoubleVector to fill with the supervector
		/// @param spk index of the speaker
		///
		void getVUYX(DoubleVector &vuyx, unsigned long spk);

		/// Get the DZ supervector for a given speaker
		/// @param dz DoubleVector to fill with the supervector
		/// @param spk index of the speaker
		///
		void getDZ(DoubleVector & dz, unsigned long spk);

		/// Get the UX supervector for a given session
		/// @param ux DoubleVector to fill with the supervector
		/// @param session index of the session
		///
		void getUX(DoubleVector & ux, unsigned long session);

		/// Get the UX supervector for a given session
		/// @param ux DoubleVector to fill with the supervector
		/// @param file name of the feature file
		///
		void getUX(DoubleVector &ux, String& file);

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
		
		/// Get the MplusDZ supervector for a given featurefile
		/// @param MplusDZ DoubleVector to fill with the supervector
		/// @param file name of the feature file
		///
		void getMplusDZ(DoubleVector &MplusDZ, String& file);

		/// Get the MplusDZ supervector for a given featurefile
		/// @param MplusDZ DoubleVector to fill with the supervector
		/// @param spk index of the speaker
		///
		void getMplusDZ(DoubleVector &MplusDZ, unsigned long spk);

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

		/// Get the MplusUX supervector for a given featurefile
		/// @param MplusUX DoubleVector to fill with the supervector
		/// @param file name of the feature file
		///
		void getMplusUX(DoubleVector &MplusUX, String& file);

		/// Get the MplusUX supervector for a given session
		/// @param MplusUX DoubleVector to fill with the supervector
		/// @param spk index of the session
		///
		void getMplusUX(DoubleVector &MplusUX, unsigned long session);

		/// Get the MplusVYplusDZ supervector for a given featurefile
		/// @param MplusVYplusDZ DoubleVector to fill with the supervector
		/// @param file name of the feature file
		///
		void getMplusVYplusDZ(DoubleVector &, String& file);

		/// Get the MplusVYplusDZ supervector for a given speaker
		/// @param MplusVYplusDZ DoubleVector to fill with the supervector
		/// @param index of the speaker
		///
		void getMplusVYplusDZ(DoubleVector &, unsigned long spk);

		/// Get the MplusVUYX supervector for a given featurefile
		/// @param MplusVUYX DoubleVector to fill with the supervector
		/// @param file name of the feature file
		///
		void getMplusVUYX(DoubleVector &MplusVUYX, String& file);

		/// Get the MplusVUYX supervector for a given speaker
		/// @param MplusVUYX DoubleVector to fill with the supervector
		/// @param index of the speaker
		///
		void getMplusVUYX(DoubleVector &MplusVUYX, unsigned long spk);	
		
		/// Compute the inverse of L matrices for EigenVoice estimation
		/// @param config config filename
		///
		void estimateAndInverseL_EV(Config&);

		/// Compute the inverse of L matrices for EigenVoice estimation without threads
		///
		void estimateAndInverseLUnThreaded_EV();

		/// Compute the inverse of L matrices for EigenChannel estimation
		/// @param config config filename
		///
		void estimateAndInverseL_EC(Config&);

		/// Compute the inverse of L matrices for EigenChannel estimation without threads
		///
		void estimateAndInverseLUnThreaded_EC();
		
		/// Compute the inverse of L matrices for EigenChannel and EigenVoice estimation
		/// @param config config filename
		///
		void estimateAndInverseL_VU(Config&);

		/// Compute the inverse of L matrices for EigenChannel and EigenVoice estimation without threads
		///
		void estimateAndInverseLUnThreaded_VU();
		
		/// Estimate EigenVoice matrix and Y components
		/// @param config config filename
		/// 
		void estimateYandV(Config &);

		/// Estimate EigenVoice matrix and Y components without threads
		/// 
		void estimateYandVUnThreaded();
		
		#ifdef THREAD
			/// Compute the inverse of L matrices for EigenVoice estimation using multi-threads
			/// @param threads number of threads to launch
			///
			void estimateAndInverseLThreaded_EV(unsigned long threads);

			/// Compute the inverse of L matrices for EigenChannel estimation using multi-threads
			/// @param threads number of threads to launch
			///
			void estimateAndInverseLThreaded_EC(unsigned long threads);

			/// Compute the inverse of L matrices for EigenChannel and EigenVoice estimation using multi-threads
			/// @param threads number of threads to launch
			///
			void estimateAndInverseLThreaded_VU(unsigned long threads);

			/// Estimate EigenVoice matrix and Y components using multi-threads
			/// @param threads number of threads to launch
			/// 
			void estimateYandVThreaded(unsigned long threads);
		#endif
		
		/// Estimate Y components
		/// @param config config filename
		/// 
		void estimateY(Config &);

		/// Estimate Y components without threads
		/// 
		void estimateYUnThreaded();

		/// Estimate EigenChannel matrix and X components
		/// @param config config filename
		/// 
		void estimateXandU(Config &);

		/// Estimate EigenChannel matrix and X components without threads
		/// 
		void estimateXandUUnThreaded();

		/// Estimate X components
		/// @param config config filename
		/// 
		void estimateX(Config &);

		/// Estimate X components without threads
		/// 
		void estimateXUnThreaded();

		#ifdef THREAD
			/// Estimate EigenChannel matrix and X components using multi-threads
			/// @param threads number of threads to launch
			/// 
			void estimateXandUThreaded(unsigned long threads);

			/// Estimate X components using multi-threads
			/// @param threads number of threads to launch
			/// 
			void estimateXThreaded(unsigned long threads);

			/// Estimate Y components using multi-threads
			/// @param threads number of threads to launch
			/// 
			void estimateYThreaded(unsigned long threads);
		#endif

		/// Estimate the EigenChannel Matrix
		/// 
		void estimateU();

		/// Estimate the D Matrix and Z components
		/// 
		void estimateZandD();

		/// Estimate the X and Y components
		/// 
		void estimateYX();

		/// Estimate the Z components Matrix
		/// 
		void estimateZ();

		/// Estimate the Z components Matrix using the MAP paradigm
		/// @param tau Maximum A Posteriori regulation factor
		/// 
		void estimateZMAP(double tau);
		
		/// Update the EigenVoice Matrix
		///
		void updateVestimate();

		/// Update the EigenChannel Matrix
		///
		void updateUestimate();
		
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

		/// Return the rank of the EigenChannel Matrix
		///
		unsigned long getRankEC();
		
		/// Return the Supervector of means of the UBM
		///
		DoubleVector& getUbmMeans();

		/// Return the Supervector of variances of the UBM
		///
		DoubleVector& getUbmInvVar();
		
		/// Return the EigenVoice Matrix
		///
		Matrix<double> getV();

		/// Return the EigenChannel Matrix
		///
		Matrix<double> getU();

		/// Return the D Matrix
		///
		DoubleVector getD();

		// Return the EigenVoice and EigenChannel matrices concatenated
		///
		Matrix<double> getVU();
		
		/// Return the Y matrix
		///
		Matrix<double> getY();

		/// Return the X matrix
		///
		Matrix<double> getX();

		/// Return the Z matrix
		///
		Matrix<double> getZ();

		/// Return the XY matrix
		///
		Matrix<double> getYX();

		/// Save matrix X
		/// @param xFilename filename to save the matrix
		///
		void saveX(String xFilename,Config &);

		/// Save matrix Y
		/// @param yFilename filename to save the matrix
		///
		void saveY(String yFilename,Config &);
		
		/// Save matrix Z
		/// @param zFilename filename to save the matrix
		///
		void saveZ(String zFilename,Config &);

		/// Cut the XY matrix to obtain two matrices X and Y
		///
		void splitYX();

		/// Store the sufficient statistics accumulators for backup
		///
		void storeAccs();

		/// Restore the sufficient statistics accumulators
		///
		void restoreAccs();
		
		/// Substract the speaker component MplusDZ from the speaker statistics
		/// @param config config filename
		///
		void substractMplusDZ(Config & config);

		/// Substract the speaker component MplusDZ from the speaker statistics without threads
		///
		void substractMplusDZUnThreaded();

		/// Substract the speaker component MplusVYplusDZ from the speaker statistics
		/// @param config config filename
		///
		void substractMplusVY(Config &);

		/// Substract the speaker component MplusVYplusDZ from the speaker statistics without threads
		///
		void substractMplusVYUnThreaded();


		/// Substract the channel component UX from the speaker statistics
		/// @param config config filename
		///
		void substractUX(Config & config);

		/// Substract the channel component UX from the speaker statistics without threads
		///
		void substractUXUnThreaded();
		
		/// Substract the channel component MplusUX from the speaker statistics
		/// @param config config filename
		///
		void substractMplusUX();				//creer la fonction multithread

		/// Substract the channel component MplusUX from the speaker statistics without threads
		///
		void substractMplusVUYX();				//creer la fonction multithread
		
		/// Substract the speaker component MplusVYplusDZ from the session statistics
		/// @param config config filename
		///
		void substractMplusVYplusDZ(Config &);

		/// Substract the speaker component MplusVYplusDZ from the session statistics without threads
		///
		void substractMplusVYplusDZUnThreaded();

		/// Substract the speaker component MplusDZ from the session statistics
		///
		void substractMplusDZByChannel();		//creer la version multithread
		
		#ifdef THREAD
		/// Substract the speaker component MplusDZ from the speaker statistics using multi-threading
		/// @param threads number of threads to launch
		///
		void substractMplusDZThreaded(unsigned long threads);

		/// Substract the channel component UX from the speaker statistics using multi-threading
		/// @param threads number of threads to launch
		///
		void substractUXThreaded(unsigned long threads);

		/// Substract the speaker component MplusVYplusDZ from the speaker statistics using multi-threading
		/// @param threads number of threads to launch
		///
		void substractMplusVYplusDZThreaded(unsigned long threads);

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

		/// Normalise features of a cluster by substracting the channel component
		/// @param segCluster segment cluster to process
		/// @param featureServer feature server to use
		/// @param config config filename
		///
		void normalizeFeatures(alize::SegCluster& segCluster, alize::FeatureServer& featureServer, alize::Config& config);

		/// Substract the Channel component from the features
		/// @param fs feature Server to use
		/// @param config config filename
		///
		void substractUXfromFeatures(FeatureServer &fs,Config &config);

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

};
#endif

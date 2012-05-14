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

#if !defined(ALIZE_PldaTools_h)
#define ALIZE_PldaTools_h

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

  /// This class represents a accumulator of statistics. 
  /// A TVAcc contains the accumulators needed for TotalVariability
  /// estimation
  /// 
  /// @author Anthony Larcher  alarcher@i2r.a-star.edu.sg
  /// @version 1.0
  /// @date 2012

class LIA_SPKTOOLS_API PldaDev{

	private :
	
		XList _fileList; 
		unsigned long _vectSize;
		unsigned long _n_speakers;
		unsigned long _n_sessions;
	
		/// Mean Vectors
		RealVector<double> _mean;
		Matrix<double> _speaker_means;

		/// Class definitions
		ULongVector _class;
		ULongVector _style;
		ULongVector _session_per_speaker;

		/// Data Matrix
		Matrix<double> _data;

	public :

		/// Constructors, load vectors from ASCII file or XList
 		PldaDev(String &,Config &);
		PldaDev(XList &,Config & );

		virtual ~PldaDev();
		virtual String getClassName() const;

		/// Compute global mean and mean per speaker
		///
		void computeAll();

		/// Get the size of data vectors
		///
		unsigned long getVectSize();

		/// Get the number of speakers
		///
		unsigned long getSpeakerNumber();

		/// Get the total number of sessions
		///
		unsigned long getSessionNumber();

		/// Get the matrix of vectors
		/// 
		Matrix<double> getData();

		/// Return the global mean of development data
		///
		RealVector<double>& getMean();

		/// Return the mean of a specific speaker
		/// @param spk the index of the speaker
		///
		void getSpeakerMean(unsigned long spk,RealVector<double>&);

		/// Normalize the length of all vectors to one (Euclidian Norm)
		///
		/// Global and speaker means are re-computed
		///
		void lengthNorm();

		/// Substrac a vector to all vectors in the PldaDev object
		/// @param mu mean vector to remove from all vectors
		///
		/// Global and speaker means are re-computed
		/// 
		void center(RealVector<double> &mu);

		/// Compute three co-variance matrices
		/// @param Sigma the total covariance matrix
		/// @param W the within class covariance matrix
		/// @param B the between class covariance matrix
		/// @param config configuration object
		///
		/// Data are centered within this function before computation of co-variance matrix
		///
		void computeCovMat(DoubleSquareMatrix &Sigma, DoubleSquareMatrix &W, DoubleSquareMatrix &B,Config &config);

		/// UnThreaded fonction to compute the three co-variance matrices
		/// @param Sigma the total covariance matrix
		/// @param W the within class covariance matrix
		/// @param B the between class covariance matrix
		/// @param threads number of threads to run
		///
		void computeCovMatUnThreaded(DoubleSquareMatrix &Sigma, DoubleSquareMatrix &W, DoubleSquareMatrix &B);


		/// Threaded fonction to compute the Choleski decomposition of Within Class Covariance matrix
		///
		#ifdef THREAD
			/// Threaded fonction to compute the three co-variance matrices
			/// @param Sigma the total covariance matrix
			/// @param W the within class covariance matrix
			/// @param B the between class covariance matrix
			/// @param threads number of threads to run
			///
			void computeCovMatThreaded(DoubleSquareMatrix &Sigma, DoubleSquareMatrix &W, DoubleSquareMatrix &B, unsigned long threads);
		#endif

		/// Compute the Choleski decomposition of Within Class Covariance matrix
		/// @param WCCN Choleski decomposition of the WCCN matrix
		///
		void computeWccnChol(DoubleSquareMatrix &WCCN, Config& config);

		/// Threaded fonction to compute the Choleski decomposition of Within Class Covariance matrix
		///
		#ifdef THREAD
			/// Threaded fonction to compute the Choleski decomposition of Within Class Covariance matrix
			/// @param config config filename
			/// @param threads number of threads to run
			///
			void computeWccnCholThreaded(DoubleSquareMatrix &WCCN, unsigned long threads);
		#endif

			/// Unthreaded fonction to compute the Choleski decomposition of Within Class Covariance matrix
			/// @param config config filename
			///
			void computeWccnCholUnThreaded(DoubleSquareMatrix &WCCN);


		/// Compute the Mahalanobis matrix for scoring (i.e. the inverse of the within-class co-variance matrix)
		/// @param M Mahalanobis matrix computed on the development set
		/// @param config configuration object
		///
		void computeMahalanobis(DoubleSquareMatrix &M, Config &config);

		/// Split the list of speaker in a given number of sublists and return the index of the first session of each sub-list 
		/// @param nbThread number of sub-list to split
		/// @param startIndex index of the first session of each sub-list
		///
		void splitPerSpeaker(unsigned long nbThread, RealVector<unsigned long> &startIndex);

		/// Compute the scatter matrices
		/// @param SW the within class covariance matrix
		/// @param SB the between class covariance matrix
		/// @params config configuration object
		///
		/// Data are centered within this function before computation of co-variance matrix
		///
		void computeScatterMat(DoubleSquareMatrix &SW, DoubleSquareMatrix &SB,Config &config);

		/// UnThreaded fonction to compute the scatter matrices
		/// @param SW the within class covariance matrix
		/// @param SB the between class covariance matrix
		/// @param threads number of threads to run
		///
		void computeScatterMatUnThreaded(DoubleSquareMatrix &W, DoubleSquareMatrix &B);

		/// Threaded fonction to compute the Choleski decomposition of Within Class Covariance matrix
		///
		#ifdef THREAD
			/// Threaded fonction to compute the the scatter matrices
			/// @param W the within class covariance matrix
			/// @param B the between class covariance matrix
			/// @param threads number of threads to run
			///
			void computeScatterMatThreaded(DoubleSquareMatrix &W, DoubleSquareMatrix &B, unsigned long threads);
		#endif

		/// Estimate parameters for Probabilistic Linear Discriminant Analysis model
		/// @param L LDA matrix computed on the development set
		/// @param config configuration object
		///
		void trainPLDA(Config &config);		// modifier les parametres...


};























  /// This class represents a accumulator of statistics. 
  /// A TVAcc contains the accumulators needed for TotalVariability
  /// estimation
  /// 
  /// @author Anthony Larcher  alarcher@i2r.a-star.edu.sg
  /// @version 1.0
  /// @date 2012

class LIA_SPKTOOLS_API PldaTest{

	private :
	
		XList _fileList; 
		unsigned long _vectSize;
		unsigned long _n_models;
		unsigned long _n_segments;
	
		XLine _modelLine;
		XLine _segLine;

		/// Data Matrix
		Matrix<double> _models;
		Matrix<double> _segments;
		Matrix<double> _scores;
		Matrix<unsigned long> _trials;

	public :

		/// Constructors, load vectors from ASCII file or XList
 		PldaTest(String &,Config &);
		PldaTest(XList &,Config & );

		virtual ~PldaTest();
		virtual String getClassName() const;


		/// Get the size of data vectors
		///
		unsigned long getVectSize();

		/// Get the number of models
		///
		unsigned long getModelsNumber();

		/// Get the number of test segments
		///
		unsigned long getSegmentsNumber();

		/// Get the matrix of models
		/// 
		Matrix<double> getModels();

		/// Get the matrix of segments
		/// 
		Matrix<double> getSegments();

		/// Get the matrix of scores
		/// 
		Matrix<double> getScores();

		/// Get the matrix of trials
		/// 
		Matrix<unsigned long> getTrials();

		/// Get the name of a model
		/// @param index is the index of the model in _modelLine
		///
		String getModelName(unsigned long index);

		/// Get the name of a segment
		/// @param index is the index of the segment in _segLine
		///
		String getSegmentName(unsigned long index);



		/// Normalize the length of all vectors to one (Euclidian Norm)
		///
		void lengthNorm();

		/// Substrac a vector to all vectors in the PldaTest object
		/// @param mu mean vector to remove from all vectors
		/// 
		void center(RealVector<double> &mu);

		/// Rotate the model and test vectors by multiplying on the left by M
		/// @param M the matrix to multiply
		///
		void rotateLeft(Matrix<double> &M);

		/// Comppute test using cosine distance
		/// 
		void cosineDistance(Config &config);

};
















#endif

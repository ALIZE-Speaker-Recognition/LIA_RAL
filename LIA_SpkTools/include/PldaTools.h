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

LIA_RAL is distributed in the hope that it will be useful
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
#include <Core>

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
		PldaDev();
 		PldaDev(String &,Config &);
		PldaDev(XList &,Config & );

		virtual ~PldaDev();
		virtual String getClassName() const;
	
		/// Load data and compute all
		/// @param ndxFilename the index file to load
		/// @param config the config
		///
		void load(String & ndxFilename,Config & config);

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
		
		/// Get the total number of sessions for a give speaker
		///
		unsigned long getSpeakerSessionNumber(unsigned long);

		/// Get the vector of the total number of sessions per speaker
		///
		ULongVector& getSpeakerSessionNumber();


		/// Get the matrix of vectors
		/// 
		Matrix<double> getData();
		
		double getData(unsigned long, unsigned long);

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

		/// Substrac a vector to all vectors in the PldaDev object
		/// @param mu mean vector to remove from all vectors (Eigen::VectorXd)
		///
		/// Global and speaker means are re-computed
		/// 
		void center(Eigen::VectorXd &mu);
		
		/// Rotate the development vectors by multiplying on the left by M
		/// @param M the matrix to multiply
		///
		void rotateLeft(Matrix<double> &M);

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

		/// Compute the total co-variance matrices in Eigen format
		/// Beware, the matrix is not normalized by the number of sessions
		/// @param Sigma the total covariance matrix
		/// @param config configuration object
		///
		/// Data are centered within this function before computation of co-variance matrix
		///
		void computeCovMatEigen(Eigen::MatrixXd &Sigma,Config &config);

		/// UnThreaded fonction to compute the total co-variance matrix
		/// @param Sigma the total covariance matrix
		///
		void computeCovMatEigenUnThreaded(Eigen::MatrixXd &Sigma);

		#ifdef THREAD
			/// Threaded fonction to compute the total co-variance matrix
			/// @param Sigma the total covariance matrix
			/// @param threads number of threads to run
			///
			void computeCovMatEigenThreaded(Eigen::MatrixXd &Sigma, unsigned long threads);
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
			void computeWccnCholThreaded(DoubleSquareMatrix &WCCN, unsigned long threads,Config & config);
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

		/// Split the list of speaker in a given number of sublists and return the index of the first session of each sub-list 
		/// @param nbThread number of sub-list to split
		/// @param startIndex index of the first session of each sub-list
		/// @param spkStartIndex index of the first speaker of each sub-list
		///
		void splitPerSpeaker(unsigned long nbThread, RealVector<unsigned long> &startIndex, RealVector<unsigned long> &spkStartIndex);

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

		/// Estimate normalization parameters or Spherical Nuisance Normalization
		/// @param config the configuration file
		///
		void sphericalNuisanceNormalization(Config& config);


		/// Estimate parameters for Probabilistic Linear Discriminant Analysis model
		/// @param config configuration object
		///
		void trainPLDA(Config &config);		// modifier les parametres...

		/// Compute Linear Discriminant Analysis matrix
		/// @param L LDA matrix computed on the development set
		/// @param ldaRank rank of the LDA matrix
		/// @param config configuration object
		///
		void computeLDA(Matrix<double> &L, long ldaRank, Config &config);

		/// Compute Eigen Value decomposition
		/// @param EP matrix to decompose
		/// @param eigenVect matrix of eigen vectors
		/// @param rank number of eigen vectors to keep
		/// @param config configuration object
		///
		void computeEigenProblem(Matrix<double> &EP,Matrix<double> &eigenVect, long rank, Config& config);

		/// Compute Eigen Value decomposition
		/// @param EP matrix to decompose
		/// @param eigenVect matrix of eigen vectors
		/// @param eigenVal diagonal matrix of eigen values
		/// @param rank number of eigen vectors to keep
		/// @param config configuration object
		///
		void computeEigenProblem(Matrix<double> &EP, Matrix<double> &eigenVect , Matrix<double> &eigenVal , long rank,Config &config);

		unsigned long getClass(unsigned long);

		ULongVector& getClass();
};






/// This class contains a PLDA model 
  /// 
  /// @author Anthony Larcher  alarcher@i2r.a-star.edu.sg
  /// @version 1.0
  /// @date 2012

class LIA_SPKTOOLS_API PldaModel{

	private :
	
		//Development data
		PldaDev _Dev;				// Development data

		// Model parameteres
		unsigned long _rankF;
		unsigned long _rankG;
		unsigned long _vectSize;

		// Model Matrices
		Eigen::VectorXd _originalMean;
		Eigen::MatrixXd _F;				// Speaker subspace 
		Eigen::MatrixXd _G;				// Channel subspace
		Eigen::MatrixXd _Sigma;			// Precision matrix
		Eigen::VectorXd _Delta;			// new Mean computed with Minimum Divergence

		// Temporary accumulators
		Eigen::MatrixXd _sigmaObs;
		Eigen::MatrixXd _invSigma;
		Eigen::MatrixXd _Ftweight;
		Eigen::MatrixXd _Gtweight;
		Eigen::MatrixXd _GtweightG;
		Eigen::MatrixXd _FtweightG;
		Eigen::MatrixXd _invGtweightGplusEye;

		// Accumulateurs pour l'apprentissage
//		Eigen::MatrixXd _Eh;	// TO DO utilise ou ?
		Eigen::MatrixXd _EhhSum;
		Eigen::MatrixXd _xhSum;
		Eigen::MatrixXd _U;


	public :

		/// Constructors
		PldaModel();

		PldaModel(String mode, Config &config);	
		// le PldaModel peut etre initialise de 2 facons differentes (train et test)

		void initTrain(PldaDev, Config &);
		// initialise l'objet avec les accumulateurs de statistiques
		// recupere les donnees de dev qui ont eventuellement deja ete traitees en parallele de donnees de test
		// initialise le model avec des matrices existantes ou init random
		// initialise les accumulateurs Eh, Ehh et u
		
		void initTest(Config &);
		//initialise l'objet avec seulement les donnees necessaires au test
		// charge depuis un objet PldaModel ou les matrices separemment

		virtual ~PldaModel();


		virtual String getClassName() const;

		void initModel(Config &config);

		void initF(Config& config);

		void splitPerSpeaker(unsigned long nbThread, RealVector<unsigned long> &startIndex);

		void initG(Config& config);

		void updateMean();

		void centerData();

		void updateModel(Config&);

		void em_iteration(Config &config,unsigned long it);

		void getExpectedValues(Config&, unsigned long it);
		void getExpectedValuesUnThreaded(Config&);

		#ifdef THREAD
			/// Threaded fonction to estimate PLDA model
			/// @param config config filename
			///
			void getExpectedValuesThreaded(unsigned long numThread, Config& config, unsigned long it);
		#endif
	
		void saveModel(Config &config);


		void mStep(unsigned long it,Config& config);

		/// Load EigenVoice matrix
		/// @param filename name of the file to load
		/// @param config the configuration
		/// 
		void loadF(String filename, Config& config);

		/// Load EigenChannel matrix
		/// @param filename name of the file to load
		/// @param config the configuration
		/// 
		void loadG(String filename, Config& config);

		/// Load Noise matrix
		/// @param filename name of the file to load
		/// @param config the configuration
		/// 
		void loadNoise(String filename, Config& config);

		/// Load Original mean vector
		/// @param filename name of the file to load
		/// @param config the configuration
		/// 
		void loadOriginalMean(String filename, Config& config);

		/// Load Mean vector
		/// @param filename name of the file to load
		/// @param config the configuration
		/// 
		void loadDelta(String filename, Config& config);

		/// Precomputation of temporary variables
		/// 
		void preComputation();

		void save(String filename, Config& config);		// a ecrire

		void load(String filename, Config& config);		// a ecrire

		Eigen::MatrixXd getFtweight();

		Eigen::MatrixXd getFtweightG();

		Eigen::MatrixXd getInvGtweightGplusEye();

		Eigen::MatrixXd getGtweight();

		Eigen::MatrixXd getF();

		unsigned long getRankF();

		unsigned long getRankG();

		PldaDev& getDev();

};
















  /// This class represents a accumulator of statistics. 
  /// A PldaTest contains the enrollment and test i-vectors for testing
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
		unsigned long _n_enrollment_segments;
		unsigned long _n_test_segments;
		unsigned long _n_enrollSessions_max;
	
		XLine _modelIDLine;					// size is the number of model, each element of the list is unique
		XLine _modelSessionslLine;			// size is the total number of enrollment sessions over all models, contains the name of the file to load
		XLine _modelIndexLine;				// size is the total number of enrollment sessions over all models, contains the corresponding modelID
		XLine _segLine;

		/// Data Matrix
		Matrix<double> _models;
		Matrix<double> _segments;
		Matrix<double> _scores;
		BoolMatrix _trials;

	public :

		/// Constructors, load vectors from ASCII file or XList
		PldaTest(String &, String & ,Config &);
 		PldaTest(String &,Config &);
		PldaTest(XList &,Config & );
		PldaTest(Config &);
		PldaTest();

		virtual ~PldaTest();
		virtual String getClassName() const;

		/// Load PldaTest data from different sources
		/// @param config the configuration
		///
		void load(Config &);

		void splitPerModel(unsigned long nbThread, RealVector<unsigned long> &startIndex, RealVector<unsigned long> &startModel);

		/// Get the size of data vectors
		///
		unsigned long getVectSize();

		/// Get the number of models
		///
		unsigned long getModelsNumber();

		/// Get the maximum number of enrollment segments
		///
		unsigned long getMaxEnrollmentSession();

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
		BoolMatrix getTrials();

		/// Get the name of a model
		/// @param index is the index of the model in _modelIDLine
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

		/// Apply Spherical Nuisance Normalization
		/// @param config the configuration
		/// 
		void sphericalNuisanceNormalization(Config& config);

		/// Comppute test using cosine distance
		/// @param config the configuration
		/// 
		void cosineDistance(Config &config);
		
		/// Compute test using cosine distance
		/// 
		void mahalanobisDistance(Matrix<double>& Mah, Config &config);


		void twoCovScoringMixPart(DoubleSquareMatrix &G, Config &config);
		void twoCovScoringMixPartUnThreaded(DoubleSquareMatrix &G, Config& config);
		void twoCovScoringMixPartThreaded(DoubleSquareMatrix &G, unsigned long NUM_THREADS);


		/// Compute test using two covariance model
		/// @param W the within class covariance matrix
		/// @param B the between class covariance matrix
		/// @param config the configuration
		///
		void twoCovScoring(DoubleSquareMatrix& W, DoubleSquareMatrix& B, Config &config);


		void pldaScoring(Config &config, Eigen::MatrixXd FTJF, double alpha_one, Eigen::MatrixXd K_one, unsigned long rankF);
		void pldaScoringUnThreaded(Config& config, Eigen::MatrixXd FTJF, double alpha_one, Eigen::MatrixXd K_one, unsigned long rankF);
		void pldaScoringThreaded(Config &config, Eigen::MatrixXd FTJF, double alpha_one, Eigen::MatrixXd K_one, unsigned long rankF, unsigned long NUM_THREADS);

		/// Compute test using Probabilistic Linear Discriminant Analysis native scoring
		/// @param pldaModel the PLDA generative model
		/// @param config the configuration
		///
		void pldaNativeScoring(PldaModel& pldaModel, Config &config);

		/// Compute test using Probabilistic Linear Discriminant Analysis and mean of enrollment sesisojns for each model
		/// @param pldaModel the PLDA generative model
		/// @param config the configuration
		///
		void pldaMeanScoring(PldaModel& pldaModel, Config &config);

		/// Save test segments vector after normalization
		/// @param outputDir output directory to store normalized vectors
		/// @param config the configuration
		///
		void saveSegments(String outputDir, Config& config);

		/// Save enrollment and test segments vector after normalization
		/// @param outputDir output directory to store normalized vectors
		/// @param config the configuration
		///
		void saveVectors(String outputDir, Config& config);

};
















#endif

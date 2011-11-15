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

#if !defined(ALIZE_Tools_h)
#define ALIZE_Tools_h

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
#include "ClusteringCriterion.h"

using namespace alize;
using namespace std;


// initMAPConfig: Initialize MAP adaptation parameters according to the context of Seg
LIA_SPKTOOLS_API MAPCfg& initMAPConfig(Config &, String);

// modifyMAPConfig: Modify MAP adaptation parameters according to the context
LIA_SPKTOOLS_API void modifyMAPConfig(Config &, MAPCfg &, String);

// learnModelFromSplit: building of a model by EM and split initialisation
LIA_SPKTOOLS_API MixtureGD& learnModelFromSplit(Config &, SegCluster &, StatServer &, FeatureServer &, MixtureServer &, unsigned long, unsigned long, unsigned long);
//LIA_SPKTOOLS_API MixtureGD& learnModelFromSplit(Config &, SegCluster &, StatServer &, FeatureServer &, MixtureServer &, unsigned long, unsigned long, unsigned long, RealVector <double> &);

// learnModel: building of a model by EM and model already initialised
LIA_SPKTOOLS_API void learnModel(Config &, SegCluster &, StatServer &, FeatureServer &, MixtureServer &, MixtureGD &);

// findAlreadyUsed: verify if the selected segment is suitable or not
LIA_SPKTOOLS_API bool findAlreadyUsed(SegCluster &, unsigned long, unsigned long, String);

// CreateInitHMM: building of an HMM with 1 state, GMM initialized by model adaptation
LIA_SPKTOOLS_API void CreateInitHMM(Config &, MAPCfg &, SegCluster &, MixtureGD &, hmm &, StatServer &, FeatureServer &, MixtureServer &);

// CreateInitHMM: building of an HMM with 1 state, GMM initialized by EM
LIA_SPKTOOLS_API void CreateInitHMM(Config &, SegCluster &, hmm &, StatServer &, FeatureServer &, MixtureServer &);

// findNameClass: Extract the acoustic class from the speaker label
LIA_SPKTOOLS_API void findNameClass(String &, String &);

// isNewSpeaker: verify wether the speaker label looks like L0x (new speaker)
LIA_SPKTOOLS_API bool isNewSpeaker(String);

// A function for initializing an hmm with several states
LIA_SPKTOOLS_API void InitHMM(Config &, SegServer &, MixtureGD &, hmm &);
LIA_SPKTOOLS_API void InitHMM(Config &, SegServer &, MixtureGD &, MixtureGD &, hmm &);
LIA_SPKTOOLS_API void InitHMMOverlap(Config &, SegServer &, MixtureGD &, hmm &);

// computeTransitions: HMM transition computation: either equiprobable or with a constant pii and pij ratio 
LIA_SPKTOOLS_API void computeTransitions(Config &, hmm &, String, double, SegServer &);
LIA_SPKTOOLS_API void computeTransitions(Config &, hmm &, SegServer &);

// comparison between two segmentations according to different criteria
LIA_SPKTOOLS_API bool isComparable(Config &, SegServer &, SegServer &, DoubleVector &, DoubleVector &);

// comparison between two segmentations according to different criteria
LIA_SPKTOOLS_API bool isComparableAndVerifSpeaker(Config &, SegServer &, SegServer &, DoubleVector &, DoubleVector &);

// comparison between two segmentations according to the boundaries
LIA_SPKTOOLS_API bool isDifferentSegmentation(Config &, SegCluster &, SegCluster &);

// comparison between two viterbi paths according to likelihoods
LIA_SPKTOOLS_API bool isLessEpsilon(Config &, DoubleVector &, DoubleVector &, double &);

// Verify the validity of the last speaker according to the corresponding duration
LIA_SPKTOOLS_API bool isValidLengthOfLastSpeaker(Config &, SegServer &);

// Compute the likelihood betwwen a segmentation and an HMM
LIA_SPKTOOLS_API double ComputePath(Config &, hmm &, DoubleVector, SegServer &, XList &, StatServer &);

// Compute the likelihood betwwen a state of an HMM and associated set of segments
LIA_SPKTOOLS_API double ComputeLLROnOneCluster(Config &, hmm &, SegServer &, unsigned long, StatServer &, FeatureServer &, DoubleVector &);

// Saving a segmentation
LIA_SPKTOOLS_API void saveSegmentation(Config &, SegServer &, FeatureServer &, String &, int);
LIA_SPKTOOLS_API void saveAcousticSegmentation(Config &, SegServer &, FeatureServer &, String &, int);

// clean the last speaker cluster for a best segment selection
LIA_SPKTOOLS_API void CleanSpeaker(Config &, SegCluster &, SegCluster &, unsigned long, LabelServer &);

// look for empty speaker (without any segments and delete them
LIA_SPKTOOLS_API void NoDataSpeakerVerification(Config &, hmm &, SegServer &);

LIA_SPKTOOLS_API int compare (const void  *, const void *);

LIA_SPKTOOLS_API double meanLikelihoodExt(StatServer &, FeatureServer &, MixtureGD &, unsigned long, unsigned long, Config &);
LIA_SPKTOOLS_API void accumulateStatLLKExt(StatServer &, FeatureServer &, MixtureStat &, unsigned long, unsigned long, Config &);

LIA_SPKTOOLS_API void viterbiDecoding(Config &, hmm &, SegCluster &, SegServer &, StatServer &, FeatureServer &, LabelServer &, DoubleVector &);
LIA_SPKTOOLS_API void copyPathInCluster(Config &, SegServer &, const ULongVector &, hmm &, LabelServer &, unsigned long, String);

LIA_SPKTOOLS_API void segAdaptation(Config &, MAPCfg &, hmm &, MixtureGD &, SegServer &, SegServer &, StatServer &, FeatureServer &, MixtureServer &);
LIA_SPKTOOLS_API void segEM(Config &, hmm &, SegServer &, SegServer &, StatServer &, FeatureServer &, MixtureServer &, unsigned long);

LIA_SPKTOOLS_API MixtureGD& createWorld(Config &, SegCluster &, StatServer &, FeatureServer &, MixtureServer &);
LIA_SPKTOOLS_API MixtureGD& createWorld(Config &, SegCluster &, StatServer &, FeatureServer &, MixtureServer &, unsigned long);

LIA_SPKTOOLS_API MixtureGD &mixtureInitBySplit(Config &, MixtureServer &, StatServer &, FeatureServer &, SegCluster &, DoubleVector &, DoubleVector &, unsigned long, TrainCfg &, unsigned long);

#endif

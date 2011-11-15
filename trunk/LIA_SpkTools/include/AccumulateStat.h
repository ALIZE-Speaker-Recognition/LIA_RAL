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

#if !defined(ALIZE_AccumulateStat_h)
#define ALIZE_AccumulateStat_h

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
using namespace alize;
using namespace std;


//-------------------------------------------------------------------------
//-- Accumulate the log likelihood for the selected frames and a given llk statistic accumulator
LIA_SPKTOOLS_API void accumulateStatLLK(StatServer &ss,FeatureServer &fs,MixtureStat &llkAcc, unsigned long idxBeginFrame,unsigned long nbFrames,
		       Config &config);  
// one a Segment
LIA_SPKTOOLS_API void accumulateStatLLK(StatServer &ss,FeatureServer &fs,MixtureStat &llkAcc,Seg* seg,Config &config);
// One on Cluster
LIA_SPKTOOLS_API void accumulateStatLLK(StatServer &ss,FeatureServer &fs,MixtureStat &llkAcc,SegCluster &selectedSegments,Config &config);

//-------------------------------------------------------------------------
//-- Accumulate the log likelihood for the selected frames and a given llk statistic accumulator
LIA_SPKTOOLS_API void accumulateStatLLK(StatServer &ss,FeatureServer &fs,MixtureStat &llkAcc, unsigned long idxBeginFrame,unsigned long nbFrames,double &weight,
		       Config &config);  
// one a Segment
LIA_SPKTOOLS_API void accumulateStatLLK(StatServer &ss,FeatureServer &fs,MixtureStat &llkAcc,Seg* seg,double &weight,Config &config);
// One on Cluster
LIA_SPKTOOLS_API void accumulateStatLLK(StatServer &ss,FeatureServer &fs,MixtureStat &llkAcc,SegCluster &selectedSegments,double &weight,Config &config);

//-------------------------------------------------------------------------
//-- accumulate the statistic for EM, using a current accumulator 
//-- CAUTION: THE ACCUMULATOR SHOULD BE INITIALIZED (resetEM) BEFORE THE CALL
//--          A GET CALL SHOULD BE DONE AFTER THE CALL  
// One a part of the feature stream
LIA_SPKTOOLS_API double accumulateStatEM(StatServer &ss,FeatureServer &fs,MixtureStat &emAcc,
			       unsigned long idxBeginFrame,unsigned long nbFrames,Config &config);
// one a Segment
LIA_SPKTOOLS_API double accumulateStatEM(StatServer &ss,FeatureServer &fs,MixtureStat &emAcc,Seg* seg,Config &config);
// One on Cluster
LIA_SPKTOOLS_API double accumulateStatEM(StatServer &ss,FeatureServer &fs,MixtureStat &emAcc,SegCluster &selectedSegments,Config &config);

// This is internal, you should only call accumulateStatEM
// Threaded version
#if defined(THREAD)
LIA_SPKTOOLS_API double accumulateStatEMThreaded(StatServer &ss,FeatureServer &fs,MixtureStat &emAcc,SegCluster &selectedSegments,Config &config);
#endif
// UnThreaded version
LIA_SPKTOOLS_API double accumulateStatEMUnThreaded(StatServer &ss,FeatureServer &fs,MixtureStat &emAcc,SegCluster &selectedSegments,Config &config);
// Weighted Version
LIA_SPKTOOLS_API double accumulateStatEM(StatServer &ss,FeatureServer &fs,MixtureStat &emAcc,SegCluster &selectedSegments,double & weight, Config &config);

/* Alex Preti things
double accumulateStatEM(StatServer &ss,FeatureServer &fs,MixtureStat &emAcc,
			       unsigned long idxBeginFrame,unsigned long nbFrames,double &weight,Config &config);
// one a Segment
double accumulateStatEM(StatServer &ss,FeatureServer &fs,MixtureStat &emAcc,Seg* seg,double &weight,Config &config);
// One on Cluster
double accumulateStatEM(StatServer &ss,FeatureServer &fs,MixtureStat &emAcc,SegCluster &selectedSegments,double &weight,Config &config);*/

//-------------------------------------------------------------------------
//-- accumulate the statistic on the frames (mean and cov), using a current 
//-- FrameAcc
//-- CAUTION: THE ACCUMULATOR SHOULD BE INITIALIZED BEFORE THE FIRST CALL
//--          A COMPUTE_ALL AND A GET CALL SHOULD BE DONE AFTER THE CALLS  
LIA_SPKTOOLS_API void accumulateStatFrame(FrameAcc & frameAcc,FeatureServer &fs,
			 unsigned long idxBeginFrame,unsigned long nbFrames,Config &config);
// one a Segment
LIA_SPKTOOLS_API void accumulateStatFrame(FrameAcc & frameAcc,FeatureServer &fs,Seg* seg,Config &config);
// One on Cluster
LIA_SPKTOOLS_API void accumulateStatFrame(FrameAcc &frameAcc,FeatureServer &fs,SegCluster &selectedSegments,Config &config);


//-------------------------------------------------------------------------
//-- accumulate the statistic on the frames raw distribution of each coefficient using a current 
//-- CAUTION: 
//         *THE ACCUMULATOR SHOULD BE INITIALIZED BEFORE THE FIRST CALL
//          initHistoTab()
//         *THE HISTO SHOULD BE COMPUTED BEFORE TO USE THE STAT
//          computeHistoTab()
//         *The histoTab should be freezen after use
//          freezeHistoTab();
//      
// Init the Histo array (one by coeff)
LIA_SPKTOOLS_API double areaHisto(const Histo & histo,unsigned long bin);
LIA_SPKTOOLS_API double areaHisto(const Histo & histo,unsigned long bin, double nonObserved);
LIA_SPKTOOLS_API double linearInterpolation(double val,double lower,double higher);
LIA_SPKTOOLS_API void freezeHistoTab(Histo* &histoT);
LIA_SPKTOOLS_API void initHistoTab(Histo* &histoT,unsigned long size, unsigned long nbBins);
//void resetHistoTab(Histo* histoT,unsigned long size);
LIA_SPKTOOLS_API void computeHistoTab(Histo* histoT,unsigned long size);
LIA_SPKTOOLS_API void accumulateHistoFrame(Histo  *histoT,FeatureServer &fs,
			  unsigned long idxBeginFrame,unsigned long nbFrames,Config &config);
// one a Segment
LIA_SPKTOOLS_API void accumulateHistoFrame(Histo *histoT,FeatureServer &fs,Seg* seg,Config &config);
// One on Cluster
LIA_SPKTOOLS_API void accumulateHistoFrame(Histo *histoT,FeatureServer &fs,SegCluster &selectedSegments,Config &config);


#endif //!defined(ALIZE_AccumulateStat_h)
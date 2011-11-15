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

#if defined(_WIN32)
#if defined(LIA_SPKTOOLS_EXPORTS)
#define LIA_SPKTOOLS_API __declspec(dllexport)
#else
#define LIA_SPKTOOLS_API __declspec(dllimport)
#endif
#else
#define LIA_SPKTOOLS_API
#endif

#include<iostream>
#include<fstream>
#include<cstdio>
#include<cassert>
#include<cmath>
#include "TrainTools.h"

using namespace alize;
using namespace std;


/**********************************************************
* mergeCluster: merge two clusters
***********************************************************/
LIA_SPKTOOLS_API SegCluster& mergeCluster(SegCluster& c1, SegCluster& c2, SegServer& segTemp,String merge="NULL");


/**********************************************************
* clrCrit: Application of CLR criterion for clustering
***********************************************************/
LIA_SPKTOOLS_API double clrCrit(Config& config, SegCluster& c1, SegCluster &c2, StatServer& ss, FeatureServer& fs, MixtureGD& m1, MixtureGD& m2, MixtureGD& world);

/**********************************************************
* glrCrit: Application of GLR criterion for clustering
***********************************************************/
LIA_SPKTOOLS_API double gllrCrit(Config& config, SegCluster& c1, SegCluster &c2, SegCluster &c12, StatServer& ss, FeatureServer& fs, MixtureGD& m1, MixtureGD& m2, MixtureGD& m12);

/**********************************************************
* bicCrit: Application of GLR criterion for clustering
***********************************************************/
LIA_SPKTOOLS_API double bicCrit(Config& config, SegCluster& c1, SegCluster &c2, SegCluster &c12, StatServer& ss, FeatureServer& fs, MixtureGD& m1, MixtureGD& m2, MixtureGD& m12);

/**********************************************************
* bicCrit: Application of Delta BIC criterion for clustering (IDIAP)
***********************************************************/
LIA_SPKTOOLS_API double deltabicCrit(Config& config, SegCluster& c1, SegCluster &c2, SegCluster &c12, StatServer& ss, FeatureServer& fs, MixtureGD& m1, MixtureGD& m2, MixtureGD&
m12);

/**********************************************************
* clusteringCriterionByAdapt: Application of different criterion for clustering (models are trained by adaptation)
***********************************************************/
LIA_SPKTOOLS_API double clusteringCriterionByAdapt(Config& config, Seg *segment1, Seg *segment2, StatServer& ss, FeatureServer& fs, MixtureGD& world, String crit);


/**********************************************************
* clusteringCriterion: Application of different criterion for clustering (models are trained by EM/ML - Initialization by external model)
***********************************************************/
LIA_SPKTOOLS_API double clusteringCriterion(Config& config, Seg *segment1, Seg *segment2, StatServer& ss, FeatureServer& fs, MixtureGD& world, String crit);

/**********************************************************
* clusteringCriterionWithoutWorldInit: Application of different criterion for clustering (models are trained by EM/ML - Initialization by external model)
***********************************************************/
LIA_SPKTOOLS_API double clusteringCriterionWithoutWorldInit(Config& config, SegCluster& c1, SegCluster& c2, StatServer& ss, FeatureServer& fs,MixtureGD& world,String crit);

/**********************************************************
* clusteringCriterionWithoutWorldInit: Application of different criterion for clustering (models are trained by EM/ML - Initialization by external model)
***********************************************************/
LIA_SPKTOOLS_API double clusteringCriterionWithoutWorldInit(Config& config, SegCluster& c1, MixtureGD& m1, SegCluster& c2, MixtureGD& m2, StatServer& ss, 
		FeatureServer& fs,MixtureGD& world,String crit);

/**********************************************************
* clusteringCriterion: Application of different criterion for clustering (models are trained by EM/ML - Initialization by external model)
***********************************************************/
LIA_SPKTOOLS_API double clusteringCriterionWithoutWorldInit(Config& config, Seg *segment1, Seg *segment2, StatServer& ss, FeatureServer& fs,MixtureGD& world,String crit);

/**********************************************************
* isSimilarSegment: compare two segments with BIC, GLR, ..., criterions
***********************************************************/
LIA_SPKTOOLS_API bool isSimilarSegment(Config& config, Seg *segment1, Seg *segment2, StatServer& ss, FeatureServer& fs, MixtureGD& world, String crit);

/**********************************************************
* bestFittingSegment: search the best segments of a cluster (in terms of normalized or not likelihood) and return it
***********************************************************/
LIA_SPKTOOLS_API Seg *bestFittingSegment(Config& config, SegCluster& cluster, MixtureGD& m, StatServer& ss, FeatureServer& fs);
LIA_SPKTOOLS_API Seg *bestFittingSegment(Config& config, SegCluster& cluster, MixtureGD& m, StatServer& ss, FeatureServer& fs, MixtureGD& world);
LIA_SPKTOOLS_API Seg *bestFittingSegment(Config& config, SegCluster& cluster, MixtureGD& m, StatServer& ss, FeatureServer& fs, hmm& actualHMM,unsigned long except);

/**********************************************************
* bestFittingCLuster: search the best cluster for a segment (in terms of normalized likelihood) and return it
***********************************************************/
LIA_SPKTOOLS_API unsigned long bestFittingCluster(Config& config, hmm& actualHMM, SegServer& actualSeg, Seg *segment, StatServer& ss, FeatureServer& fs, unsigned long exceptInd=200);

/**********************************************************
* intraCluster: evaluate purity intra cluster
***********************************************************/
LIA_SPKTOOLS_API void intraCluster(Config& config, hmm& actualHMM, SegServer& actualSeg, StatServer& ss, FeatureServer& fs, MixtureGD& world, String crit);

/**********************************************************
* interCluster: Evaluate purity inter cluster
***********************************************************/
LIA_SPKTOOLS_API void interCluster(Config& config, hmm& actualHMM, SegServer& actualSeg, StatServer& ss, FeatureServer& fs, MixtureGD& world, String crit);

/**********************************************************
 * clusteringCriterionWithoutWorldInit: Application of different criterion for clustering (models are trained by EM/ML - Initialization by external model)
 ***********************************************************/
LIA_SPKTOOLS_API double clusteringCriterionWithoutWorldInitOneGaus(Config& config, SegCluster& c1, SegCluster& c2, StatServer& ss, FeatureServer& fs,String crit);

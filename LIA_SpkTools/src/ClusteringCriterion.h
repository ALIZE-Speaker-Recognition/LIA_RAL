//****************************************
// Audio segmentation functions
// A part of LIA_SpkSeg system
// First version February 2006
// Authors:
//
// Corinne Fredouille (corinne.fredouille@lia.univ-avignon.fr)
//*****************************************

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
* 
* Author C. Fredouille February 2006
***********************************************************/
SegCluster& mergeCluster(SegCluster& c1, SegCluster& c2, SegServer& segTemp,String merge="NULL");


/**********************************************************
* clrCrit: Application of CLR criterion for clustering
* 
* Author C. Fredouille February 2006
***********************************************************/
double clrCrit(Config& config, SegCluster& c1, SegCluster &c2, StatServer& ss, FeatureServer& fs, MixtureGD& m1, MixtureGD& m2, MixtureGD& world);

/**********************************************************
* glrCrit: Application of GLR criterion for clustering
* 
* Author C. Fredouille February 2006
***********************************************************/
double gllrCrit(Config& config, SegCluster& c1, SegCluster &c2, SegCluster &c12, StatServer& ss, FeatureServer& fs, MixtureGD& m1, MixtureGD& m2, MixtureGD& m12);

/**********************************************************
* bicCrit: Application of GLR criterion for clustering
* 
* Author C. Fredouille February 2006
***********************************************************/
double bicCrit(Config& config, SegCluster& c1, SegCluster &c2, SegCluster &c12, StatServer& ss, FeatureServer& fs, MixtureGD& m1, MixtureGD& m2, MixtureGD& m12);

/**********************************************************
* bicCrit: Application of Delta BIC criterion for clustering (IDIAP)
* 
* Author C. Fredouille February 2006
***********************************************************/
double deltabicCrit(Config& config, SegCluster& c1, SegCluster &c2, SegCluster &c12, StatServer& ss, FeatureServer& fs, MixtureGD& m1, MixtureGD& m2, MixtureGD&
m12);

/**********************************************************
* clusteringCriterionByAdapt: Application of different criterion for clustering (models are trained by adaptation)
* 
* Author C. Fredouille February 2006
***********************************************************/
double clusteringCriterionByAdapt(Config& config, Seg *segment1, Seg *segment2, StatServer& ss, FeatureServer& fs, MixtureGD& world, String crit);


/**********************************************************
* clusteringCriterion: Application of different criterion for clustering (models are trained by EM/ML - Initialization by external model)
* 
* Author C. Fredouille February 2006
***********************************************************/
double clusteringCriterion(Config& config, Seg *segment1, Seg *segment2, StatServer& ss, FeatureServer& fs, MixtureGD& world, String crit);

/**********************************************************
* clusteringCriterionWithoutWorldInit: Application of different criterion for clustering (models are trained by EM/ML - Initialization by external model)
* 
* Author C. Fredouille February 2006
***********************************************************/
double clusteringCriterionWithoutWorldInit(Config& config, SegCluster& c1, SegCluster& c2, StatServer& ss, FeatureServer& fs,MixtureGD& world,String crit);

/**********************************************************
* clusteringCriterionWithoutWorldInit: Application of different criterion for clustering (models are trained by EM/ML - Initialization by external model)
* 
* Author C. Fredouille February 2006
***********************************************************/
double clusteringCriterionWithoutWorldInit(Config& config, SegCluster& c1, MixtureGD& m1, SegCluster& c2, MixtureGD& m2, StatServer& ss, 
		FeatureServer& fs,MixtureGD& world,String crit);

/**********************************************************
* clusteringCriterion: Application of different criterion for clustering (models are trained by EM/ML - Initialization by external model)
* 
* Author C. Fredouille February 2006
***********************************************************/
double clusteringCriterionWithoutWorldInit(Config& config, Seg *segment1, Seg *segment2, StatServer& ss, FeatureServer& fs,MixtureGD& world,String crit);

/**********************************************************
* isSimilarSegment: compare two segments with BIC, GLR, ..., criterions
* 
* Author C. Fredouille February 2006
***********************************************************/
bool isSimilarSegment(Config& config, Seg *segment1, Seg *segment2, StatServer& ss, FeatureServer& fs, MixtureGD& world, String crit);

/**********************************************************
* bestFittingSegment: search the best segments of a cluster (in terms of normalized or not likelihood) and return it
* 
* Author C. Fredouille February 2006
***********************************************************/
Seg *bestFittingSegment(Config& config, SegCluster& cluster, MixtureGD& m, StatServer& ss, FeatureServer& fs);
Seg *bestFittingSegment(Config& config, SegCluster& cluster, MixtureGD& m, StatServer& ss, FeatureServer& fs, MixtureGD& world);
Seg *bestFittingSegment(Config& config, SegCluster& cluster, MixtureGD& m, StatServer& ss, FeatureServer& fs, hmm& actualHMM,unsigned long except);



/**********************************************************
* bestFittingCLuster: search the best cluster for a segment (in terms of normalized likelihood) and return it
* 
* Author C. Fredouille February 2006
***********************************************************/
unsigned long bestFittingCluster(Config& config, hmm& actualHMM, SegServer& actualSeg, Seg *segment, StatServer& ss, FeatureServer& fs, unsigned long exceptInd=200);


/**********************************************************
* intraCluster: evaluate purity intra cluster
* 
* Author C. Fredouille February 2006
***********************************************************/
void intraCluster(Config& config, hmm& actualHMM, SegServer& actualSeg, StatServer& ss, FeatureServer& fs, MixtureGD& world, String crit);

/**********************************************************
* interCluster: Evaluate purity inter cluster
* 
* Author C. Fredouille February 2006
***********************************************************/
void interCluster(Config& config, hmm& actualHMM, SegServer& actualSeg, StatServer& ss, FeatureServer& fs, MixtureGD& world, String crit);

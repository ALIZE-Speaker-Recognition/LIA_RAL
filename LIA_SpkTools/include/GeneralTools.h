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

#if !defined(ALIZE_GeneralTools_h)
#define ALIZE_GeneralTools_h

#if defined(_WIN32)
#if defined(LIA_SPKTOOLS_EXPORTS)
#define LIA_SPKTOOLS_API __declspec(dllexport)
#else
#define LIA_SPKTOOLS_API __declspec(dllimport)
#endif
#else
#define LIA_SPKTOOLS_API
#endif

#include "alize.h"
#include "liatools.h"

using namespace alize;
using namespace std;

//--------------------------------------------------------------------------------
//----------------------         ComputeTest stuff    ----------------------------
const unsigned long CST_MAX_CLIENT_LINE=500;

class LIA_SPKTOOLS_API TabClientLine{
  unsigned long nbMaxModelLine;
  MixtureServer *ms;
  Config *conf;
  MixtureGD ** tabModel;
  XLine *pline;
  unsigned long nbModelsLine;
  
public:
  unsigned long nbClientLine();
  String& getClientName(unsigned long);
  MixtureGD& getClientModel(unsigned long);
  unsigned long loadLine(XLine*,String="",bool=true,bool=false);  
  TabClientLine(MixtureServer &, Config &, unsigned long=CST_MAX_CLIENT_LINE);
  ~TabClientLine();
};

//---------------------------------------------------------------------------------
//-----------------------  TabHisto definition ----------------------------------
class LIA_SPKTOOLS_API TabHisto{
  unsigned long _nb;                        // Nb of histo
  ObjectRefVector _tabHisto;				// histo tabs
  XLine _id;                                // the client ID;
public:
  TabHisto(unsigned long nbBins, unsigned long maxTargetLine){_nb=0;for (unsigned int i=0;i<maxTargetLine;i++) _tabHisto.addObject(*(new Histo(nbBins)));}// build a default histo numbers and default nbBins
  ~TabHisto(){_tabHisto.deleteAllObjects();}                           // Verify if the memory is freezen for ...Vector and the XLine - TODO
  unsigned long getSize(){return _nb;}
  void accumulateValueInTab(const String &id,double score);
  Histo & getHisto(const String &id);
  void computeHistoInTab(const String &id);
  String & getId(unsigned long n);
  unsigned long getIndex(const String &id);
  Histo & getHistoFromVect(unsigned long n);
};

 
//-----------------------------------------------------------------------------------
// Compute the entropy from an Histo
LIA_SPKTOOLS_API double computeEntropy(Histo & hist);

//-----------------------------------------------------------------------------------
// Compute mean from an Histo
LIA_SPKTOOLS_API double computeMean(Histo & hist);


LIA_SPKTOOLS_API long setDecision(double LLRClient, double decisionThreshold);

//---------------------------------------------------------------------------------
//-----------------------  ScoreAccum definition ----------------------------------
// JFB
class LIA_SPKTOOLS_API ScoreAccum{
  unsigned long _nb;                        // Nb of  clients
  DoubleVector _score;                      // The scores
  ULongVector _nbFrame;                     // The accumulated value
  XLine _id;                                // the client ID;
public:
  ScoreAccum(){_nb=0;}
  ~ScoreAccum(){}                           // Verify if the memory is freezen for ...Vector and the XLine - TODO
  unsigned long getSize(){return _nb;}
  void addAndAccumulate(const String &id,double score, unsigned long nbFrames);
  double getScore(const String &id);
  double getScore(unsigned long n);
  String & getId(unsigned long n);
  unsigned long getIndex(const String &id);
};

// Class for component selection
struct TabWeightElem{
  double weight;
  unsigned long distrib; 
  DistribGD *distribP;
};

LIA_SPKTOOLS_API int _compF(const void * op1, const void *op2);

class LIA_SPKTOOLS_API TabWeight{
  unsigned long _size;
  TabWeightElem *_tab;
  unsigned long _nbTop;
  void _sortByWeight(const MixtureGD &data){
    for(unsigned i=0;i<_size;i++){
      _tab[i].weight=data.weight(i);
      _tab[i].distrib=i;
      _tab[i].distribP=&(data.getDistrib(i));
    }
    qsort ( (void *) _tab,_size, sizeof(TabWeightElem),_compF);
  }
  void _nbTopDyn(double threshold){
    _nbTop=0;
    for (double tot=0;(_nbTop<_size)&&(tot<threshold);_nbTop++) tot+=_tab[_nbTop].weight;
    if (verbose) cout <<"Dyn nb top["<<_nbTop<<"]"<<endl;
  }
public:
  TabWeight(){};
  TabWeight(const MixtureGD &);
  TabWeight(const MixtureGD &,unsigned long);
  TabWeight(const MixtureGD &,double);
  void init(const MixtureGD &,unsigned long);
  void init(const MixtureGD &,double);
  unsigned long getNbTop(){return _nbTop;}
  unsigned long getDistrib(unsigned long id){return _tab[id].distrib;}
  double getWeight(unsigned long id){return _tab[id].weight;}
  void setDistrib(unsigned long id,unsigned long distrib,double weight){_tab[id].distrib=distrib;_tab[id].weight=weight;}
  TabWeightElem *getArray(){return _tab;}
  ~TabWeight(){delete []_tab;}
};


// Random picking of frames (bagging) functions, based on segment/cluster processing
// It is independent of the segment length
// Try to decrease the number of segments for fasting the world model training
// Used mainly in TrainTools.cpp 
// Author: JFB
LIA_SPKTOOLS_API bool baggedFrame(double baggedFrameProbability);
LIA_SPKTOOLS_API void baggedSegments(SegCluster &,SegCluster &,double,unsigned long=3,unsigned long=7);
LIA_SPKTOOLS_API void baggedSegmentsConstraint(SegCluster &,SegCluster &,double,unsigned long=3,unsigned long=7);
// Works on a set of bagged clusters - only one reading of the
// segments and multiple selections, one by bagged cluster
LIA_SPKTOOLS_API void baggedSegments(SegCluster &selectedSegments,SegCluster &baggedSeg,unsigned long nbBagged,double & baggedProbability,
		    unsigned long minimumLength,unsigned long maximumLength);
// Same but returns both selected and unselected clusters
// Take care, both clusters should be created in the same server - NOT TESTED
LIA_SPKTOOLS_API void baggedSegments(SegCluster &selectedSegments,SegCluster &baggedSelected,SegCluster &baggedUnselected,double baggedProbability,unsigned long minimumLength,unsigned long maximumLength);
//-------------------------------------------------------------------------
//-- Compute the mean log likelihood for the Selected frames and a given model
LIA_SPKTOOLS_API double meanLikelihood(StatServer &ss,FeatureServer &fs,MixtureGD &model,unsigned long idxBeginFrame,unsigned long nbFrames,Config &config);
// one a cluster
LIA_SPKTOOLS_API double meanLikelihood(StatServer &ss,FeatureServer &fs,MixtureGD &model,SegCluster &selectedSegments,Config &config);
// on a set of input streams - OLD and NEw, TODO FIRST SHOULD BE DELETED
LIA_SPKTOOLS_API double meanLikelihood(StatServer &ss,FeatureServer **fsTab,SegCluster **segTab,unsigned long nbStream,MixtureGD &model,Config &config);
LIA_SPKTOOLS_API double meanLikelihood(StatServer &ss,ObjectRefVector &FeatServ,ObjectRefVector &ClusterSeg,MixtureGD &model,DoubleVector &decision,Config &config);

//-------------------------------------------------------------------------
//-- Compute the mean and cov of selected the data using segmental mode
LIA_SPKTOOLS_API void globalMeanCov (FeatureServer &fs,SegCluster &selectedSegments,FrameAcc & globalFrameAcc,Config &config);  
// On a complete feature stream
LIA_SPKTOOLS_API void globalMeanCov (FeatureServer &fs,FrameAcc & globalFrameAcc,Config &config);  

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
LIA_SPKTOOLS_API void computeHistoTab(Histo* histoT,unsigned long size);
LIA_SPKTOOLS_API void accumulateHistoFrame(Histo  *histoT,FeatureServer &fs,
			  unsigned long idxBeginFrame,unsigned long nbFrames,Config &config);
// one a Segment
LIA_SPKTOOLS_API void accumulateHistoFrame(Histo *histoT,FeatureServer &fs,Seg* seg,Config &config);
// One on Cluster
LIA_SPKTOOLS_API void accumulateHistoFrame(Histo *histoT,FeatureServer &fs,SegCluster &selectedSegments,Config &config);



// ----------------------------------------------------------------------------------------------------------
// Feature Warping giving a source(tab of histo, one by coeff) and a target distribution 
// for a segment and cluster (segment is considerred to be the minimum time unit to perform this)
LIA_SPKTOOLS_API void computeWarp(Histo *histoT,Histo &destH,FeatureServer & fs,unsigned long begin, unsigned long length,Config &config);
// on a segment
LIA_SPKTOOLS_API void computeWarp(Histo *histoT,Histo &destH, FeatureServer & fs, Seg* seg, Config & config);
// on a cluster
LIA_SPKTOOLS_API void computeWarp(Histo *histoT,Histo &destH, FeatureServer & fs, SegCluster & selectedSegments, Config & config);
 
// ----------------------------------------------------------------------------------------------------------
// Feature Mean subtraction and Cov reduction for a segment; cluster

LIA_SPKTOOLS_API void cms(String &,FeatureServer &,Config &);

LIA_SPKTOOLS_API void computeZeroOne(FrameAccGD &frameAccu,FeatureServer & fs,unsigned long begin, unsigned long length,Config &config) ;
// on a segment
LIA_SPKTOOLS_API void computeZeroOne(FrameAccGD & frameAccu, FeatureServer & fs, Seg* seg, Config & config);
// on a cluster
LIA_SPKTOOLS_API void computeZeroOne(FrameAccGD & frameAccu, FeatureServer & fs, SegCluster & selectedSegments, Config & config);
// 
LIA_SPKTOOLS_API void computeZeroOne(const DoubleVector &mean,const DoubleVector& cov,FeatureServer & fs,unsigned long begin, unsigned long length,Config &config) ;
// on a segment
LIA_SPKTOOLS_API void computeZeroOne(const DoubleVector &mean,const DoubleVector& cov, FeatureServer & fs, Seg* seg, Config & config);
// on a cluster
LIA_SPKTOOLS_API void computeZeroOne(const DoubleVector &mean,const DoubleVector& cov, FeatureServer & fs, SegCluster & selectedSegments, Config & config);

// Feature writing in an output stream w - could be used for multiple segmen,ts from multiple files to one file
LIA_SPKTOOLS_API void outputFeatureFile(Config &config, FeatureServer &fs, Feature & f, FeatureFileWriter &w);
// on a segment
LIA_SPKTOOLS_API void outputFeatureFile(Config &config, FeatureServer &fs, Seg * seg, FeatureFileWriter &w);
// on a cluster
LIA_SPKTOOLS_API void outputFeatureFile(Config &config, FeatureServer &fs, SegCluster & selectedSegments, FeatureFileWriter &w);
// on a part of a file
LIA_SPKTOOLS_API void outputFeatureFile(Config &config, FeatureServer &fs, unsigned long begin,unsigned long length, FeatureFileWriter &w);

// Feature Mapping related
LIA_SPKTOOLS_API unsigned long getBestGaussian(Mixture & M, Feature & f);
LIA_SPKTOOLS_API void mapDataToDistrib(double & data, const double meanData, const double covData, const double meanMap, const double covMap);

LIA_SPKTOOLS_API void featureMapping(MixtureServer & ms, FeatureServer & fs,unsigned long begin, unsigned long length,Config &config);
// on a segment
LIA_SPKTOOLS_API void featureMapping(MixtureServer & ms, FeatureServer & fs,Seg * seg,Config &config);

// on a cluster
LIA_SPKTOOLS_API void featureMapping(MixtureServer & ms, FeatureServer & fs,SegCluster & selectedSegments,Config &config);


// Model based likelihood computation 
LIA_SPKTOOLS_API double likelihoodGD(const DistribGD&,const DistribGD &);           
LIA_SPKTOOLS_API double likelihoodGD(const MixtureGD&,const MixtureGD &);           
LIA_SPKTOOLS_API double likelihoodGD(const MixtureGD&,const MixtureGD &,TabWeight&);
LIA_SPKTOOLS_API double likelihoodGD(const MixtureGD& ,const MixtureGD &,TabWeight &,TabWeight &);

#endif //!defined(ALIZE_GeneralTools_h)

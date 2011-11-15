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

#if !defined(ALIZE_SegTools_h)
#define ALIZE_SegTools_h

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
#include <liatools.h>


using namespace alize;
using namespace std;

// Display a cluster
LIA_SPKTOOLS_API void showCluster(SegCluster& cluster);
// return in clusterOutput the clusterOne without the data in clusterTwothe number of frame in the cluster
LIA_SPKTOOLS_API unsigned long totalFrame(SegCluster& cluster);
// return the last frame of a segment
LIA_SPKTOOLS_API unsigned long endSeg(Seg *seg);


//------------------------------------------------------------------------------------------------------------------------------
// Output a cluster in a  label file
// Author: JFB
LIA_SPKTOOLS_API void outputLabelFile(SegCluster &selectedSeg,String FileName,Config &config);


/*********************************************************
* fonction d'analyse par segments d'une liste de fichiers "listSeg" 
* Paramètres entrées :
* pointeur vers le serveur de label "labelServer"
* le nombre de modèle pour apprentissage "nbrModels"
* tableau avec les labels d'apprentissage de chaque modele "labelaAnaliser"
* un pointeur vers une fonction "pfunc" qui est la fonction qui represente l'action à appliquer (EM,MAP,etc...)
* deux pointeurs vers void "pv1" et "pv2" utilisé par le pointeur vers la fonction 
*
* Auteur : Dan Istrate novembre 2004
**********************************************************/
LIA_SPKTOOLS_API void analiseParSegments(Config &config, XList *listSeg, LabelServer *labelServer, int nbrModels, XLine *labelaAnaliser,SegServer &segmentsServer,void (*pfunc)(int,Feature&,Config &config,void *pv1,void *pv2),void *pv1,void *pv2);



//******* Create a default cluster with the label labelToAnalyse, including all the input file frames
//** From N Scheffer
//** Modifications JFB*/
LIA_SPKTOOLS_API void createDefSeg (SegServer & ss, FeatureServer & fs,Config &config);

//************************************** 
// Verify the labels - if a segment is finishing after the last frame, trunc it*/
LIA_SPKTOOLS_API void verifyClusterFile(SegServer& segmentsServer,FeatureServer& fs,Config& config);

// From a time to a frame index, taking into account possible bugs on boundaries
LIA_SPKTOOLS_API unsigned long timeToFrameIdx(real_t time,real_t frameLength);
LIA_SPKTOOLS_API real_t frameIdxToTime(unsigned long idx,real_t frameLength);

// Functions for reading the label files in a segment server
// Load the segments of a file in one cluster
LIA_SPKTOOLS_API void  loadLabelFile(SegCluster &cluster,String fileName,String path, String extension,Config &config);
// the main function for loading the labels - do the job for a given file
LIA_SPKTOOLS_API void loadClusterFile(String &fileName,SegServer& segmentsServer,LabelServer& labelServer,Config& config);
// for a list of file
LIA_SPKTOOLS_API void initializeClusters(const XLine& listFiles,SegServer& segmentsServer,LabelServer& labelServer,Config& config);
// For a list of input files, stored into a XList
LIA_SPKTOOLS_API void initializeClusters(const XList& listXFiles,SegServer& segmentsServer,LabelServer& labelServer,Config& config);
// For a filename parameter, which could be a feature file or a feature file list with an .lst extention
LIA_SPKTOOLS_API void initializeClusters(String &file,SegServer& segmentsServer,LabelServer& labelServer,Config& config);

LIA_SPKTOOLS_API void initializeClusters(String &file,SegServer& segmentsServer,LabelServer& labelServer,Config& config,bool nomComplet);


// A function for create and initialize a Viterbi accumulator
LIA_SPKTOOLS_API ViterbiAccum& createAndInitializeViterbiAccum(StatServer &ss, hmm &cHmm);

// accumulateStatViterbi() is used for computing the viterbi path
// TAKE CARE: THE accumulator should be reseted before the call
LIA_SPKTOOLS_API void   accumulateStatViterbi(FeatureServer &fs,ViterbiAccum &va,unsigned long beginIdx,unsigned long length,Config &config);

// On a segment
LIA_SPKTOOLS_API void accumulateStatViterbi(FeatureServer &fs,ViterbiAccum &va,Seg *seg,Config &config);

// On a segment with PreSeg
// Dan Istrate juillet 2005
LIA_SPKTOOLS_API void accumulateStatViterbi(FeatureServer &fs,ViterbiAccum &va,Seg *seg,Config &config,SegServer& TransitSegServer,DoubleVector&
TransitionsFort,DoubleVector& TransitionsFaible,unsigned long NbState);
    
// Remove all the segments of a cluster
LIA_SPKTOOLS_API void removeAllSeg(SegServer &segServer,SegCluster &cluster);

// initialize an existing segmentation by deleting all the segments and creating an empty cluster by state/speaker
LIA_SPKTOOLS_API void initializeCluster(SegServer &currentSeg,hmm& cHmm,LabelServer &labelServer);

// removePartOfSeg() removes a part of a segment from a cluster
LIA_SPKTOOLS_API bool removePartOfSeg(SegCluster &clusterSeg,SegServer& segServer,LabelServer &labelServer,
		     String segSourceName,unsigned long segStart,unsigned long segLength);

LIA_SPKTOOLS_API void displayAllClusters(Config& config, SegServer& seg);

LIA_SPKTOOLS_API void displayOneCluster(Config& config, SegCluster &clusterT);

LIA_SPKTOOLS_API void displayAllSegments(Config& config, SegServer& seg);

LIA_SPKTOOLS_API void displayAllSegmentsFromRef(Config& config, String &fileInit, unsigned long fileSize);

//return the duration of all the segments in a cluster
LIA_SPKTOOLS_API unsigned long computeClusterTime(SegCluster& cluster);

// look for a specific label in a list
LIA_SPKTOOLS_API long findLabel(XLine classToAnalyse,String labelToFind);

LIA_SPKTOOLS_API void moveSegmentFromOneClusterToAnother(LabelServer& labelServer,Seg *segment, SegCluster &currentCluster, SegCluster &newCluster);

/**********************************************************
* longerSegment: search the longer  of a cluster (in terms of time) and return it
* 
* Author C. Fredouille February 2006
***********************************************************/
LIA_SPKTOOLS_API Seg *longerSegment(Config& config, SegCluster& cluster);

/**********************************************************
* findClusterIndex: search the index of a cluster in a segment server
* 
* Author C. Fredouille February 2006
***********************************************************/
LIA_SPKTOOLS_API unsigned long findClusterIndex(String name,SegServer& segTmp);

#endif

// SegTools.h
// This file is a part of LIA Softwares LIA_SpkDet and LIA_SpkSeg, based on ALIZE toolkit 
// LIA_SpkDet is a free, open tool for speaker recognition
// LIA_SpkSeg is a free, open tool for speaker segmentation
// This project is a development project initiated and funded by the LIA lab.
// See www.lia.univ-avignon.fr
// 
// ALIZE is needed for LIA_SpkDet and LIA_SpkSeg
// for more information about ALIZE, see http://www.lia.univ-avignon.fr/heberges/ALIZE/
// first version Noveember 2004
// New version February 2005
//
// Copyright (C) 2004
//  Laboratoire d'informatique d'Avignon [www.lia.univ-avignon.fr]
// Main authors
//  Dan Istrate [dan.istrate@lia.univ-avignon.fr]]
//  Jean-Francois Bonastre [jean-francois.bonastre@lia.univ-avignon.fr]
//      
// LIA_SpkDet and LIA_SpkSeg are free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
// You should have received a copy of the GNU General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// The LIA team as well as the ALIZE project want to highlight the limits of voice authentication
// in a forensic context. 
// The following paper proposes a good overview of this point:
// [Bonastre J.F., Bimbot F., Boe L.J., Campbell J.P., Douglas D.A., Magrin-chagnolleau I.,
//  Person  Authentification by Voice: A Need of Caution,
//  Eurospeech 2003, Genova]
// The conclusion of the paper of the paper is proposed bellow:
// [Currently, it is not possible to completely determine whether the
//  similarity between two recordings is due to the speaker or to other
//  factors, especially when: (a) the speaker does not cooperate, (b) there
//  is no control over recording equipment, (c) recording conditions are not 
//  known, (d) one does not know whether the voice was disguised and, to a
//  lesser extent, (e) the linguistic content of the message is not
//  controlled. Caution and judgment must be exercised when applying speaker
//  recognition techniques, whether human or automatic, to account for these
//  uncontrolled factors. Under more constrained or calibrated situations,
//  or as an aid for investigative purposes, judicious application of these
//  techniques may be suitable, provided they are not considered as infallible.
//  At the present time, there is no scientific process that enables one to
//  uniquely characterize a person=92s voice or to identify with absolute
//  certainty an individual from his or her voice.]
//
// Contact Jean-Francois Bonastre (jean-francois.bonastre@lia.univ-avignon.fr) for
// more information about the licence or the use of LIA_SpkDet or LIA_SpkSeg


#if !defined(ALIZE_SegTools_h)
#define ALIZE_SegTools_h

#include <alize.h>
#include <liatools.h>


using namespace alize;
using namespace std;

// Display a cluster
void showCluster(SegCluster& cluster);
// return in clusterOutput the clusterOne without the data in clusterTwothe number of frame in the cluster
unsigned long totalFrame(SegCluster& cluster);
// return the last frame of a segment
unsigned long endSeg(Seg *seg);


//------------------------------------------------------------------------------------------------------------------------------
// Output a cluster in a  label file
// Author: JFB
void outputLabelFile(SegCluster &selectedSeg,String FileName,Config &config);


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
void analiseParSegments(Config &config, XList *listSeg, LabelServer *labelServer, int nbrModels, XLine *labelaAnaliser,SegServer &segmentsServer,void (*pfunc)(int,Feature&,Config &config,void *pv1,void *pv2),void *pv1,void *pv2);



//******* Create a default cluster with the label labelToAnalyse, including all the input file frames
//** From N Scheffer
//** Modifications JFB*/
void createDefSeg (SegServer & ss, FeatureServer & fs,Config &config);

//************************************** 
// Verify the labels - if a segment is finishing after the last frame, trunc it*/
void verifyClusterFile(SegServer& segmentsServer,FeatureServer& fs,Config& config);

// From a time to a frame index, taking into account possible bugs on boundaries
unsigned long timeToFrameIdx(real_t time,real_t frameLength);
real_t frameIdxToTime(unsigned long idx,real_t frameLength);

// Functions for reading the label files in a segment server
// Load the segments of a file in one cluster
void  loadLabelFile(SegCluster &cluster,String fileName,String path, String extension,Config &config);
// the main function for loading the labels - do the job for a given file
void loadClusterFile(String &fileName,SegServer& segmentsServer,LabelServer& labelServer,Config& config);
// for a list of file
void initializeClusters(const XLine& listFiles,SegServer& segmentsServer,LabelServer& labelServer,Config& config);
// For a list of input files, stored into a XList
void initializeClusters(const XList& listXFiles,SegServer& segmentsServer,LabelServer& labelServer,Config& config);
// For a filename parameter, which could be a feature file or a feature file list with an .lst extention
void initializeClusters(String &file,SegServer& segmentsServer,LabelServer& labelServer,Config& config);

void initializeClusters(String &file,SegServer& segmentsServer,LabelServer& labelServer,Config& config,bool nomComplet);


// A function for create and initialize a Viterbi accumulator
ViterbiAccum& createAndInitializeViterbiAccum(StatServer &ss, hmm &cHmm);

// accumulateStatViterbi() is used for computing the viterbi path
// TAKE CARE: THE accumulator should be reseted before the call
void   accumulateStatViterbi(FeatureServer &fs,ViterbiAccum &va,unsigned long beginIdx,unsigned long length,Config &config);

// On a segment
void accumulateStatViterbi(FeatureServer &fs,ViterbiAccum &va,Seg *seg,Config &config);

// On a segment with PreSeg
// Dan Istrate juillet 2005
void accumulateStatViterbi(FeatureServer &fs,ViterbiAccum &va,Seg *seg,Config &config,SegServer& TransitSegServer,DoubleVector&
TransitionsFort,DoubleVector& TransitionsFaible,unsigned long NbState);
    
// Remove all the segments of a cluster
void removeAllSeg(SegServer &segServer,SegCluster &cluster);

// initialize an existing segmentation by deleting all the segments and creating an empty cluster by state/speaker
void initializeCluster(SegServer &currentSeg,hmm& cHmm,LabelServer &labelServer);

// removePartOfSeg() removes a part of a segment from a cluster
bool removePartOfSeg(SegCluster &clusterSeg,SegServer& segServer,LabelServer &labelServer,
		     String segSourceName,unsigned long segStart,unsigned long segLength);

void displayAllClusters(Config& config, SegServer& seg);

void displayOneCluster(Config& config, SegCluster &clusterT);

void displayAllSegments(Config& config, SegServer& seg);

void displayAllSegmentsFromRef(Config& config, String &fileInit, unsigned long fileSize);

//return the duration of all the segments in a cluster
unsigned long computeClusterTime(SegCluster& cluster);

// look for a specific label in a list
long findLabel(XLine classToAnalyse,String labelToFind);

void moveSegmentFromOneClusterToAnother(LabelServer& labelServer,Seg *segment, SegCluster &currentCluster, SegCluster &newCluster);

/**********************************************************
* longerSegment: search the longer  of a cluster (in terms of time) and return it
* 
* Author C. Fredouille February 2006
***********************************************************/
Seg *longerSegment(Config& config, SegCluster& cluster);

/**********************************************************
* findClusterIndex: search the index of a cluster in a segment server
* 
* Author C. Fredouille February 2006
***********************************************************/
unsigned long findClusterIndex(String name,SegServer& segTmp);

#endif

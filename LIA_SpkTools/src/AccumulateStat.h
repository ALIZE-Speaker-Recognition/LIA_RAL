// AccumulateStat.h
// This file is a part of LIA Softwares LIA_SpkDet and LIA_SpkSeg, based on ALIZE toolkit 
// LIA_SpkDet is a free, open tool for speaker recognition
// LIA_SpkSeg is a free, open tool for speaker segmentation
// This project is a development project initiated and funded by the LIA lab.
// See www.lia.univ-avignon.fr
// 
// ALIZE is needed for LIA_SpkDet and LIA_SpkSeg
// for more information about ALIZE, see http://www.lia.univ-avignon.fr/heberges/ALIZE/
// First Version 07/15/2004
// New version February 2005
//
// Copyright (C) 2004
//  Laboratoire d'informatique d'Avignon [www.lia.univ-avignon.fr]
// Main author:
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

#if !defined(ALIZE_AccumulateStat_h)
#define ALIZE_AccumulateStat_h
#include <alize.h>
#include "liatools.h"
using namespace alize;
using namespace std;


//-------------------------------------------------------------------------
//-- Accumulate the log likelihood for the selected frames and a given llk statistic accumulator
void accumulateStatLLK(StatServer &ss,FeatureServer &fs,MixtureStat &llkAcc, unsigned long idxBeginFrame,unsigned long nbFrames,
		       Config &config);  
// one a Segment
void accumulateStatLLK(StatServer &ss,FeatureServer &fs,MixtureStat &llkAcc,Seg* seg,Config &config);
// One on Cluster
void accumulateStatLLK(StatServer &ss,FeatureServer &fs,MixtureStat &llkAcc,SegCluster &selectedSegments,Config &config);

//-------------------------------------------------------------------------
//-- Accumulate the log likelihood for the selected frames and a given llk statistic accumulator
void accumulateStatLLK(StatServer &ss,FeatureServer &fs,MixtureStat &llkAcc, unsigned long idxBeginFrame,unsigned long nbFrames,double &weight,
		       Config &config);  
// one a Segment
void accumulateStatLLK(StatServer &ss,FeatureServer &fs,MixtureStat &llkAcc,Seg* seg,double &weight,Config &config);
// One on Cluster
void accumulateStatLLK(StatServer &ss,FeatureServer &fs,MixtureStat &llkAcc,SegCluster &selectedSegments,double &weight,Config &config);

//-------------------------------------------------------------------------
//-- accumulate the statistic for EM, using a current accumulator 
//-- CAUTION: THE ACCUMULATOR SHOULD BE INITIALIZED (resetEM) BEFORE THE CALL
//--          A GET CALL SHOULD BE DONE AFTER THE CALL  
// One a part of the feature stream
double accumulateStatEM(StatServer &ss,FeatureServer &fs,MixtureStat &emAcc,
			       unsigned long idxBeginFrame,unsigned long nbFrames,Config &config);
// one a Segment
double accumulateStatEM(StatServer &ss,FeatureServer &fs,MixtureStat &emAcc,Seg* seg,Config &config);
// One on Cluster
double accumulateStatEM(StatServer &ss,FeatureServer &fs,MixtureStat &emAcc,SegCluster &selectedSegments,Config &config);

// This is internal, you should only call accumulateStatEM
// Threaded version
#ifdef THREAD
double accumulateStatEMThreaded(StatServer &ss,FeatureServer &fs,MixtureStat &emAcc,SegCluster &selectedSegments,Config &config);
#endif
// UnThreaded version
double accumulateStatEMUnThreaded(StatServer &ss,FeatureServer &fs,MixtureStat &emAcc,SegCluster &selectedSegments,Config &config);
// Weighted Version
double accumulateStatEM(StatServer &ss,FeatureServer &fs,MixtureStat &emAcc,SegCluster &selectedSegments,double & weight, Config &config);

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
void accumulateStatFrame(FrameAcc & frameAcc,FeatureServer &fs,
			 unsigned long idxBeginFrame,unsigned long nbFrames,Config &config);
// one a Segment
void accumulateStatFrame(FrameAcc & frameAcc,FeatureServer &fs,Seg* seg,Config &config);
// One on Cluster
void accumulateStatFrame(FrameAcc &frameAcc,FeatureServer &fs,SegCluster &selectedSegments,Config &config);


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
double areaHisto(const Histo & histo,unsigned long bin);
double areaHisto(const Histo & histo,unsigned long bin, double nonObserved);
double linearInterpolation(double val,double lower,double higher);
void freezeHistoTab(Histo* &histoT);
void initHistoTab(Histo* &histoT,unsigned long size, unsigned long nbBins);
//void resetHistoTab(Histo* histoT,unsigned long size);
void computeHistoTab(Histo* histoT,unsigned long size);
void accumulateHistoFrame(Histo  *histoT,FeatureServer &fs,
			  unsigned long idxBeginFrame,unsigned long nbFrames,Config &config);
// one a Segment
void accumulateHistoFrame(Histo *histoT,FeatureServer &fs,Seg* seg,Config &config);
// One on Cluster
void accumulateHistoFrame(Histo *histoT,FeatureServer &fs,SegCluster &selectedSegments,Config &config);


#endif //!defined(ALIZE_AccumulateStat_h)

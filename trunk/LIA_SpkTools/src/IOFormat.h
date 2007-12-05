// IOFormat.h
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

#if !defined(ALIZE_IOFormat_h)
#define ALIZE_IOFormat_h

#include <alize.h>

using namespace alize;
using namespace std;

//******************************************************
// These functions are used to ouput Objects into a vector
//******************************************************
//------------------------------------------------------------------------
// Output Weight of a  model in a vector -- Nicolas SCHEFFER (depracated)
void outputWeightVector(MixtureGD& world,MixtureGD& clientMixture, double frameCount, Config& config); 
//------------------------------------------------------------------------
// Output vector in SVMLight Format -- Nicolas SCHEFFER
void outputSVMLightVector(RealVector<double	>&,const String &,Config&); 
//------------------------------------------------------------------------
// Output Matrix in SVMLight Format
void outputDSM(DoubleSquareMatrix &,MixtureGD&,Config&);

//******************************************************
// These functions are related to file format
//******************************************************

//-------------------------------------------------------------------------
// LIA Internal Format
//Output a result line in a file in LIA format: F client 1 test 0 20 -0.02
void outputResultLine(double LLR, const String & clientName, const String & seg, const String &gender, int decision, ostream& strm);

void outputResultLine(double LLR, const String & clientName, const String & seg,real_t begin,real_t end, const String &gender, int decision, ostream& strm);

void outputResultLIARALLine(const String &gender, const String & clientName, const String & channel, const String & seg, const String & start, const String & duration, double LLR, ostream& strm);

//-------------------------------------------------------------------------
//Output a result line in NIST 04 format
void outputResultNIST04Line(const String &trainTypeTest, const String &adaptationMode, const String &segTypeTest, const String &gender, const String & clientName, const String & seg, const String &decision, double LLR, ostream& strm);

//-------------------------------------------------------------------------
//Output a result line in ETF format
void outputResultETFLine(const String &source, const String &channel, const String &start, double duration, const String &type, const String &sub, const String &event, double score, const String &decision, ostream& strm);

//-------------------------------------------------------------------------
//Output a result line in MDTM format
void outputResultMDTMLine(const String &source, const String &channel, const String &start, double duration, const String &type, double conf, const String &sub, ostream& strm);
//-------------------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------------------
// Affect the fields values depending on the format
//-------------------------------------------------------------------------------------------------------

// Return LIA info fields
void setLIAInfoFields(unsigned long & genderField, unsigned long & nameField, unsigned long & decisionField, unsigned long & segField, unsigned long & scoreField);

// Return NIST info fields
void setNIST04InfoFields(unsigned long & genderField, unsigned long & nameField, unsigned long & decisionField, unsigned long & segField, unsigned long & scoreField);

#endif //!defined(ALIZE_IOFormat_h)

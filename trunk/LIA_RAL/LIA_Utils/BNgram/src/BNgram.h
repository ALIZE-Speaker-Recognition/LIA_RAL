// BNgram.h
// This file is a part of LIA Softwares LIA_SpkDet and LIA_SpkSeg, based on ALIZE toolkit 
// LIA_SpkDet is a free, open tool for speaker recognition
// LIA_SpkSeg is a free, open tool for speaker segmentation
// This project is a development project initiated and funded by the LIA lab.
// See www.lia.univ-avignon.fr
//   
// ALIZE is needed for LIA_SpkDet and LIA_SpkSeg
// for more information about ALIZE, see http://www.lia.univ-avignon.fr/heberges/ALIZE/
// First version 
//
// Copyright (C) 2005
//  Laboratoire d'informatique d'Avignon [www.lia.univ-avignon.fr]
// Main authors
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

// Author:
// Jean-Francois Bonastre (jean-francois.bonastre@univ-avignon.fr)
// First version November 12/30/2005

#if !defined(ALIZE_BNgram_h)
#define ALIZE_BNgram_h

#include <iostream>  
#include <fstream>  // pour outFile
#include <cstdio>   // pour printf()
#include <cassert> // pour le debug pratique
#include <cmath>
#include <liatools.h>

// NGRAM FUNCTIONS
struct BNgramNode{
  short int idx;         // First symb of the node
  unsigned long *count;  // Count array for the symbs of the node
  BNgramNode**child;     // child array for the next level
  BNgramNode *brother;   // link for the next node of this level
};
class BNgram{
  // unsigned long _nb;
  unsigned long _totalCount;
  short int _maxLength;
  unsigned long *_circArray;
  unsigned long _readIdx;
  short int  _nbElt;
  BNgramNode *_seed;
  void _freeTree(BNgramNode*);
  BNgramNode*_new();
  void _add(short int,bool,short int);
  BNgramNode *_newNode(short int,BNgramNode*);
  BNgramNode *_findInsert(BNgramNode*,short int);
  void _show(BNgramNode*,String,ostream &);
public:
  BNgram(unsigned long,unsigned long);
  ~BNgram();
  void beginSeq(short int,short int=1 );
  void addSymb(short int,short int=1 );
  void endSeq(short int=1);
  void show(ostream &);
}; 

void testBNgram(Config &config);
#endif

// BNgram.cpp
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

#if !defined(ALIZE_BNgram_cpp)
#define ALIZE_BNgram_cpp

#include <iostream>  
#include <fstream>  // pour outFile
#include <cstdio>   // pour printf()
#include <cassert> // pour le debug pratique
#include <cmath>
#include <liatools.h>
#include "BNgram.h"


// BTREE based NGRAM FUNCTIONS

// Private
BNgramNode * BNgram::_newNode(short int symb,BNgramNode*br=NULL){
  BNgramNode * ptr=new BNgramNode;
  ptr->child=new BNgramNode* [_nbElt];
  ptr->count=new unsigned long [_nbElt];
  for (short int i=0;i<_nbElt;i++){ptr->child[i]=NULL;ptr->count[i]=0;}
  ptr->idx=(symb/_nbElt)*_nbElt;
  ptr->brother=br;
  return ptr;
}
BNgramNode * BNgram::_findInsert(BNgramNode* ptr,short int idx){
  // TAKE CARE, A LIST ALWAYS BEGIN BY A NODE WITH IDX=0
  if (idx<(ptr->idx+_nbElt))return ptr;                   // current node is the good one
  while ((ptr->brother!=NULL) && (idx>=(ptr->brother->idx+_nbElt))) ptr=ptr->brother;
  if (ptr->brother==NULL) {ptr->brother=_newNode(idx); return ptr->brother;}          // new node at the end
  if ((idx>=ptr->brother->idx)&& (idx<ptr->brother->idx+_nbElt)) return ptr->brother; // the node exists
  ptr->brother=_newNode(idx,ptr->brother); return ptr->brother;                       // new node in the midle;
}
void BNgram::_add(short int symb,bool first,short int nb=1){
  // TAKE CARE, A LIST ALWAYS BEGIN BY A NODE WITH IDX=0
  // 1: seed memorizes the good seed for the n level. It is set and created if necessary by the level n-1
  static BNgramNode* seed; 
  static short int level; // level is the current ngram order set by the level n-1 (i.e. at the end of this function)
  if (first) {
    seed=_seed;           // If first level (unigram), begin at the Btree seed, else searchSeed and level are ready
    level=0;
    _totalCount++;
  }
  if (debug) cout <<"add symb["<<symb<<"] Level["<<level<<"]"<<endl;
  // 2: find the good node and add the count 
  seed=_findInsert(seed,symb); // find the good node, create it if necessary (3 takes care of the begin of the list)
  seed->count[symb%_nbElt]+=nb;
  
  //3:  Now, prepare the work for the next level
  if (seed->child[symb%_nbElt]==NULL) seed->child[symb%_nbElt]=_newNode(0);
  seed=seed->child[symb%_nbElt];
  level++;
}
void BNgram::_freeTree(BNgramNode *seed){
  if (seed==NULL) return;
  for (short int i=0;i<_nbElt;i++)_freeTree(seed->child[i]);
  delete [] seed->child;
  delete [] seed->count;
  delete seed;
}
void BNgram::_show(BNgramNode * seed,String ngramS="",ostream & str=cout){
  if (seed==NULL) return;
  for (BNgramNode* ptr=seed;(ptr!=NULL);ptr=ptr->brother)
    for (short int i=0;i<_maxLength;i++)
      if (ptr->count[i]!=0){
	String newNgramS=ngramS+" "+String::valueOf(ptr->idx+i);
	str<<newNgramS<<" "<<ptr->count[i]<<endl;
	_show(ptr->child[i],newNgramS,str);
      } 
}
// Public
BNgram::BNgram(unsigned long maxLength,unsigned long nbElementByNode=100){
  _maxLength=maxLength;
  _totalCount=0;
  //  _nb=0;
  _circArray=new unsigned long [_maxLength];
  _readIdx=0;
  _nbElt=nbElementByNode;
  _seed=_newNode(0);
}
BNgram::~BNgram(){
  delete [] _circArray;
  _freeTree(_seed);
}
void BNgram::beginSeq(short int symb,short int nb){
  _readIdx=0;
  _circArray[_readIdx++]=symb;
}
void BNgram::addSymb(short int symb,short int nb){
  _circArray[_readIdx%_maxLength]=symb;
  if (_readIdx>=(unsigned long) (_maxLength-1)){
    _add(_circArray[(_readIdx-(_maxLength-1))%_maxLength],true,nb);
    for (unsigned long i=(_readIdx-(_maxLength-1))+1;i<=_readIdx;i++)
      _add(_circArray[i%_maxLength],false,nb);
  }
  _readIdx++;
}
void BNgram::endSeq(short int nb){ // Take care, readIdx was wrongly incremented 
  for (unsigned long length=_maxLength-1;length>0;length--)
    if (_readIdx>=length){
      _add(_circArray[(_readIdx-length)%_maxLength],true,nb);
      for (unsigned long i=(_readIdx-length)+1;i<_readIdx;i++)
	_add(_circArray[i%_maxLength],false,nb);
    }
}

void BNgram::show(ostream & str=cout){
  String ngramS="";
  str<<"symbol total count["<<_totalCount<<"]"<<endl;
  _show(_seed,ngramS,str);
}

// TEST MAIN FUNCTON
void testBNgram(Config &config){
  short int order=config.getParam("maxOrder").toLong();
  short int symb=-1;
  BNgram ngram(order,2);
  cin>>symb;
 
  bool begin=true;
  while (symb!=-1){
    if (begin){
       ngram.beginSeq(symb);
       begin=false;
    }
    else ngram.addSymb(symb);
    cin>>symb;
  }
  ngram.endSeq();
  ngram.show();
}
#endif

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

#if !defined(ALIZE_SequenceExtractor_cpp)
#define ALIZE_SequenceExtractor_cpp

#include <iostream>  
#include <fstream>  // pour outFile
#include <cstdio>   // pour printf()
#include <cassert> // pour le debug pratique
#include <cmath>
#include <string>

#include "liatools.h"
#include "BNGram.h"
#include "SequenceExtractor.h"

 
//*************************************************************************************************************************
//************ Common part tree functions
//************ common part tree is the main class for the ssequence exctraction. It represents the  ngrams in the tree
//************* and implements the sequence extraction functions using this structure
//*************************************************************************************************************************
CommonPartTree::CommonPartTree()
{
  _totalCount=0;
  _totalChildCount=0;
  _seed=NULL;
}
CommonPartTree:: CommonPartTree(const CommonPartTree &obj){
throw Exception("recopy constructor not available"
        , __FILE__, __LINE__);
}
CommonPartTree CommonPartTree:: operator=(const CommonPartTree & obj){
throw Exception("= constructor not available"
        , __FILE__, __LINE__);
}

CommonPartTree::~CommonPartTree(){
  _freeTree(_seed);
}
void CommonPartTree::_freeTree(CommonPartTreeNode *ptr){
  if (ptr==NULL) return;
  _freeTree(ptr->br);
  _freeTree(ptr->ch);
  delete ptr;
}
CommonPartTreeNode* CommonPartTree::_newNode(const short int symb,const unsigned long count,
					    CommonPartTreeNode *ch,CommonPartTreeNode *br){
  CommonPartTreeNode *ptr=new CommonPartTreeNode;
  ptr->symb=symb;
  ptr->count=count;
  ptr->totalChildCount=0;
  ptr->ch=ch;
  ptr->br=br;
  return ptr;
}

CommonPartTreeNode* CommonPartTree::_findInsert(const short int symb,unsigned long count, CommonPartTreeNode * ptr){
  if (ptr==NULL)       // We arrived on a leaf
    return _newNode(symb,count,NULL,NULL);
  if (ptr->symb!=symb) {
    CommonPartTreeNode *pTmp=_findInsert(symb,count,ptr->br);
    if (ptr->br==NULL) ptr->br=pTmp;
    return pTmp;
  }
  else return ptr;
}

void CommonPartTree::addNGram(NGram & nGram){
  for (unsigned long idx=0; idx<nGram.getSize();idx++){ // For each ngram
    unsigned long count=nGram.getCount(idx);
    CommonPartTreeNode *currentP=_findInsert(nGram.getSymbol(idx,0),count,_seed);
    if (_seed==NULL) _seed=currentP;  // special case for first insertion
    CommonPartTreeNode *tmpP=NULL;
    for (unsigned long order=1;order<nGram.getOrder();order++){
      tmpP=currentP;
      currentP=_findInsert(nGram.getSymbol(idx,order),count,currentP->ch);
      if (tmpP->ch==NULL) tmpP->ch=currentP;
    } // end of ngram loop  
    if (nGram.getOrder()==1)  _totalChildCount+=count; // special case for unigrams
    else tmpP->totalChildCount+=count;
  } // End ngram loop
}

// findPartSeq finds the path corresponding to the seq in the tree and return the node corresponding to 
// the requested order. Retunr null if no node is available
CommonPartTreeNode * CommonPartTree::_findPartSeq(Seq &seq,short int order,CommonPartTreeNode *ptr){
  if (ptr==NULL) return NULL;
  if (seq.getLength()==0) return _seed;      
  if (order>=seq.getLength()) return NULL;                // 
  if (seq[order].isIn(ptr->symb)){ 
    if (order==seq.getLength()-1) return ptr;               // We are on the last part of the sequence
    else return _findPartSeq(seq,order+1,ptr->ch);
  }
  else return _findPartSeq(seq,order,ptr->br);
}

// findMaxSeq search the sequence with the maximum count and length and return the corresponding count
// Could work from scratch or from a current seg and an order where to work
// It launches recursive private functions

unsigned long CommonPartTree::_findMaxSeq(CommonPartTreeNode *ptr,Seq &seq, short int & order){
  if (ptr==NULL){ order=0; return 0;}
  short int orderCh=order+1;
  short int orderBr=order;
  Seq seqBr=seq;
  seq.setLength(order+1);seq[order].setSymb(ptr->symb);
  unsigned long brCount=_findMaxSeq(ptr->br,seqBr,orderBr);
  unsigned long chCount=_findMaxSeq(ptr->ch,seq,orderCh);
 
  if (orderBr<orderCh){seq[order].setSymb(ptr->symb); order=orderCh;return chCount;}  // child is the longuest sequence and win
  // child and/or borther == NULL
  if (orderCh==0){
    if (orderBr==0) return ptr->count;                                   // Brother==NULL, current is the winner
    if (orderBr>order) {order=orderBr; seq=seqBr;return brCount;}       // Brother find the longest sequence
    // equal length and child==NULL
    if (ptr->count>brCount)  return ptr->count;                              // current is the winner
    else {order=orderBr; seq=seqBr;return brCount;}                     // Brother is the winner
  }
  // child and Brother !=NULL
  if (orderBr>orderCh){order=orderBr; seq=seqBr;return brCount;}        // Brother find the longest sequence and win
  // Equal length (and child and brother !=NULL)     
  if (brCount>(chCount)){order=orderBr;seq=seqBr;return brCount;}       // Brother has the max count and win
  else{seq[order].setSymb(ptr->symb);order=orderCh;return chCount;}     // child has the max count and win
}
unsigned long CommonPartTree::findMaxSeq(Seq &seq){
  short int order=0; 
  if (_seed!=NULL){
    return _findMaxSeq(_seed,seq,order);
  }
  return 0;
}
unsigned long CommonPartTree::findMaxEndSeq(Seq &seq){
  if (seq.getLength()==0) return findMaxSeq(seq);
  if (_seed!=NULL){
    CommonPartTreeNode * ptr=_findPartSeq(seq,0,_seed);
    if (ptr==NULL) return 0;
    short int order=seq.getLength();
    if (ptr->ch != NULL) return _findMaxSeq(ptr->ch,seq,order);
    else return ptr->count;
  }
  return 0;
}

// It suppress the nodes corresponding to the sequence in the tree - interface function and recursive function
// Take care, a sequence could contain more than one symbol by node...
CommonPartTreeNode* CommonPartTree::_suppressSeq(CommonPartTreeNode * ptr,Seq &seq,
						short int order, unsigned long &childCountDelta){
  if (order>=seq.getLength())
    throw Exception("sequence longer than the tree" , __FILE__, __LINE__);
  CommonPartTreeNode* ptrRet=ptr;
  CommonPartTreeNode* ptrTmp=ptr;
  while ((ptr!=NULL)&& (!seq[order].isIn(ptr->symb))) {ptrTmp=ptr;ptr=ptr->br;}
  if (seq[order].isIn(ptr->symb)){      // We are in the good tree branch
    if (order==seq.getLength()-1) {     // and at the end of the sequence
      childCountDelta=ptr->count;
      ptrTmp->br=ptr->br;
      if (ptrRet==ptr) ptrRet=ptr->br;
      delete ptr; 
    }
    else{ // We are in the good branch but not at the end of the sequence
      ptr->ch=_suppressSeq(ptr->ch,seq,order+1,childCountDelta);  
      if (ptr->count<childCountDelta) throw Exception("count problem in the tree, childcount < delta" , __FILE__, __LINE__);
      ptr->totalChildCount-=childCountDelta;
      ptr->count-=childCountDelta;
      if (ptr->count==0) {               // We are in a node with a new count=0, to be destroyed
	ptrTmp->br=ptr->br;
	if (ptrRet==ptr) ptrRet=ptr->br;
      }
    }
  }
  else{ // Branch nor find
    childCountDelta=0;
    ptrRet= NULL;
  }
  return ptrRet;
}

void CommonPartTree::suppressSeq(Seq &seq){
  short int order=0;
  if (seq.getLength()==0) return;
  unsigned long childCountDelta=0;
  _seed=_suppressSeq(_seed,seq,order,childCountDelta);
  _totalChildCount-=childCountDelta;
}
void CommonPartTree::_show(CommonPartTreeNode *ptr,unsigned long order){
  if (ptr==NULL) return;
  cout <<"Order["<<order<<"] Symbol["<<ptr->symb<<"] count["<<ptr->count<<"]Childcount["<<ptr->totalChildCount<<"]"<<endl;
  _show(ptr->ch,order+1);
  _show(ptr->br,order);
}
void CommonPartTree::show(){
  cout <<"TotalCount["<<_totalCount<<"]totalChildCount["<<_totalChildCount<<"]"<<endl;
  _show(_seed,0);
}
//*************************************************************************************************************************
//************ SymbTab 
//************ SymTab is simply a code for the symbols. It allows to group several input symbols in only one node
//************ It is ssimply an array of booleans. The boolean i says if the input symbol i is present or not in the node
//************ (this multiple symbols functionality is not used currently, implemented for future use)
//*************************************************************************************************************************
void SymbTab::reserveMem(){
  _symbTab= new bool[_nbSymb];
}
SymbTab:: SymbTab(const SymbTab &obj){
  _nbSymb=obj._nbSymb;
  reserveMem();
  for (int i=0;i<_nbSymb;i++) _symbTab[i]=obj._symbTab[i];
} 
SymbTab& SymbTab:: operator=(const SymbTab &obj){
  if (this!=(&obj)){
    if (_symbTab!=NULL) delete [] _symbTab;
    _nbSymb=obj._nbSymb;
    reserveMem();
    for (int i=0;i<_nbSymb;i++) _symbTab[i]=obj._symbTab[i];
  }
  return (*this);
}
bool SymbTab:: operator==(const SymbTab & obj){
  if ((_symbTab==NULL)|| (obj._symbTab==NULL)) throw Exception("symbol array not allowed" , __FILE__, __LINE__);
  if (_nbSymb!=obj._nbSymb) throw Exception("symbol arrays differ in size" , __FILE__, __LINE__);
  for (int i=0;i<_nbSymb;i++) if (_symbTab[i]!=obj._symbTab[i]) return false;
  return true;
}
bool SymbTab:: operator!=(const SymbTab & obj){
  if ((_symbTab==NULL)|| (obj._symbTab==NULL)) throw Exception("symbol array not allowed" , __FILE__, __LINE__);
  if (_nbSymb!=obj._nbSymb) throw Exception("symbol arrays differ in size" , __FILE__, __LINE__);
  for (int i=0;i<_nbSymb;i++) if (_symbTab[i]!=obj._symbTab[i]) return true;
  return false;
}
bool SymbTab::isIn(short int symb){
  if (_symbTab==NULL) throw Exception("symbol array not allowed" , __FILE__, __LINE__);
  if ((symb<0)||(symb>=_nbSymb)){
    cout << "symb=["<<symb<<"]"<<endl;
    throw Exception("symbol out of boundaries" , __FILE__, __LINE__);
  }
  return _symbTab[symb];
}
void SymbTab::setSymb(short int symb){
  if ((symb<0) || (symb>_nbSymb))throw Exception("index out of bound"
						 , __FILE__, __LINE__);
  _symbTab[symb]=true;
}
void SymbTab::show(){
  cout <<"[";
  for (int i=0;i<_nbSymb;i++) 
    if (_symbTab[i]) cout <<i<<" ";
  cout <<"]"<<endl;
}

//*****************************************************************************************************************************************
//*** ReadMemory groups all the functionalities for reading the input files, in a segmental mode
//*** It reads the symbol one by one in the segment and says when the end of file or the end of segment is reached
//*** It also allows to comeback in the history for the deconding, for syntax analysys using non LL1 grammary (used here
//*** It will allow binary file support asap
//*****************************************************************************************************************************************

ReadMemory::ReadMemory(String inputFilename,int bufSize=100,unsigned long begin=0,unsigned long length=0){
  _inputFile.open(inputFilename.c_str());
  _bufSize=bufSize;
  _buf=new short int[_bufSize];
  _idx=0;
  _realIdx=0;
  _begin=begin;
  _length=length;
  _segmental=(length>0);
}
ReadMemory::~ReadMemory(){
  delete [] _buf;
  _inputFile.close();
}
bool ReadMemory::notEof(){
  if (!_segmental) return (!_inputFile.eof());
  else return ((_idx<_begin+_length)&&(!_inputFile.eof()));
}
bool ReadMemory::eof(){
  return !(notEof());
}
bool ReadMemory::readSymb(short int &symb,string & oov){ 
  if (notEof()){
    char c;
    _inputFile>>c;
    while (notEof() && (c<=32)) _inputFile>>c; 
    if (notEof()){
      _inputFile.putback(c);
      if ((c >= '0') && (c <= '9')) _inputFile >> symb;
      else{
	_inputFile>>oov;
	symb=-1;
      }
      if (debug) cout <<"symbol["<<_idx<<"]=["<<symb<<"]"<<endl;
      return true;
    }

  }
  return false;
}
bool ReadMemory::lecture(bool init=false){
  short int tmpS=-1;
  string oov;
  if ((!init) && (_idx<_realIdx)){
    _idx++; 
    return true;
  }
  if ((init)&&(_segmental))
    while (_idx<_begin){
      if (readSymb(tmpS,oov)){
	_idx++;
	_realIdx++;
	if (tmpS==-1) cerr << "Warning, oov["<<oov<<"] detected and ignored, idx["<<_realIdx<<"]"<<endl;
      } else return false;
    }
  
  if (notEof()){
    if (!init) {_realIdx++;_idx++;}
    while (readSymb(tmpS,oov)){
      if (tmpS==-1) cerr << "Warning, oov["<<oov<<"] detected and ignored, idx["<<_realIdx<<"]"<<endl;
      else {_buf[_realIdx%_bufSize]=tmpS;return true ;}
      _realIdx++;_idx++;
    }
  }
  return false;
}
unsigned long ReadMemory::getIdx(){
  return _idx;
}
short int ReadMemory::getCurrentSymb(){
  return _buf[_idx%_bufSize];
}
void ReadMemory::setIdx(unsigned long idx){
  _idx=idx;
}

void Seq::_reserveMem(){
  _array=new SymbTab[_maxLength]; 
  for (int i=0;i<_maxLength;i++){  
    _array[i].setNbSymb(_nbInputSymb);
    _array[i].reserveMem();
  }
}
Seq::Seq(const Seq& obj){
  _maxLength=obj._maxLength;
  _length=obj._length;
  _nbInputSymb=obj._nbInputSymb;   
  _reserveMem();
  for (int i=0;i<_maxLength;i++) _array[i]=obj._array[i];   
}
Seq::Seq(short int maxLength, short int nbInputSymb){
  _maxLength=maxLength;
  _nbInputSymb=nbInputSymb;
  _length=0;
  _reserveMem();    
}
Seq::~Seq(){
  delete [] _array;
}
Seq & Seq::operator=(const Seq& obj){
  if (this !=(&obj)){
    delete [] _array;
    _maxLength=obj._maxLength;
    _nbInputSymb=obj._nbInputSymb;   
    _length=obj._length;
    _reserveMem();
    for (int i=0;i<_maxLength;i++) _array[i]=obj._array[i];   
  }
  return *(this);
}
SymbTab& Seq::operator[](short int idx){
  if (idx>=_maxLength) throw Exception("index > sequence length"
				      , __FILE__, __LINE__);
  return _array[idx];
}
void Seq::setLength(short int length){
  if (length>_maxLength) throw Exception("index > sequence length", __FILE__, __LINE__);
  _length=length;
}
short int Seq::getLength(){
  return _length;
}
void Seq::init(short int length=0){
  if ((length>_maxLength)|| (length<0)) throw Exception("index > max sequence length" , __FILE__, __LINE__);
  for (int i=length;i<_maxLength;i++)
    _array[i].init(false);
  _length=length;
}
void Seq::show(){
  cout <<"sequence"<<endl;
  for (int i=0;i<_length;i++){
    cout <<"["<<i<<"]";_array[i].show();}
}
// Add a symbole at the end of the sequence
void Seq:: add(SymbTab *ptr){   
  if (_length==_maxLength) throw Exception("out of bound" , __FILE__, __LINE__);
  _array[_length]=(*ptr);
  _length++;
}


// SEQUENCE DECODER
SequenceDecoder::SequenceDecoder(short int nbInputSymb){          
  _nbOutputSeq=0;
  _nbOutputSeqPart=0;
  _nbInputSymb=nbInputSymb;
  _seed=NULL;
}
SequenceDecoder:: SequenceDecoder(const SequenceDecoder &obj){
throw Exception("recopy constructor not available"
        , __FILE__, __LINE__);
}
SequenceDecoder SequenceDecoder:: operator=(const SequenceDecoder & obj){
throw Exception("= constructor not available"
        , __FILE__, __LINE__);
}
SequenceDecoder::~SequenceDecoder(){
  _freeTree(_seed);
}
void SequenceDecoder::_freeTree(SequenceDecoderNode *seed){
  if (seed==NULL) return;
  _freeTree(seed->br);
  _freeTree(seed->ch);
  delete seed->symbTab;
  delete seed;
}
SequenceDecoderNode* SequenceDecoder::_newNode(const SymbTab &symbTab,const short int outputSymb,
					      SequenceDecoderNode *ch,SequenceDecoderNode *br){
  SequenceDecoderNode *ptr=new SequenceDecoderNode;
  ptr->symbTab=new SymbTab(symbTab);
  ptr->outputSymb=outputSymb;
  ptr->ch=ch;
  ptr->br=br;
  return ptr;
}
SequenceDecoderNode* SequenceDecoder::_findInsert(const SymbTab &symbTab,SequenceDecoderNode *ptr){
  if (ptr==NULL) return _newNode(symbTab,-1,NULL,NULL); // We arrived on a leaf
  while ((ptr->br!=NULL) && (*(ptr->symbTab)!=symbTab)) ptr=ptr->br;
  if (*(ptr->symbTab)==symbTab) return ptr;
  else{
    ptr->br=_newNode(symbTab,-1,NULL,NULL);
    return ptr->br;
  }
}

void SequenceDecoder::addSequence(Seq &sequence, const short int outputSymb){
  if (sequence.getLength()==0) throw Exception(" Try to add a null length sequence" , __FILE__, __LINE__); 
  SequenceDecoderNode *currentP=_findInsert(sequence[0],_seed);
  if (_seed==NULL) _seed=currentP;
  for (int step=1;step<sequence.getLength();step++){  // For each of the other steps in the sequence
    SequenceDecoderNode* tmpP=_findInsert(sequence[step],currentP->ch);
    if (currentP->ch==NULL)currentP->ch=tmpP;
    currentP=tmpP;
  } 
  if (currentP->outputSymb != -1) throw Exception("Try to add a sequence on an existing one" , __FILE__, __LINE__); 
  currentP->outputSymb=outputSymb;
  _nbOutputSeqPart++;
}
void SequenceDecoder::show(){
  _show(_seed,0);
}
void SequenceDecoder::_show(SequenceDecoderNode *ptr,short int order){
  if (ptr==NULL) return;
  cout <<"Order["<<order<<"]"<<endl;
  (*ptr->symbTab).show();
  if (ptr->outputSymb!=-1) cout<<"outputSymb["<<ptr->outputSymb<<"]"<<endl;
  _show(ptr->ch,order+1);
  _show(ptr->br,order);
}

void SequenceDecoder::_toFile(SequenceDecoderNode *ptr,short int symb,Seq & currentSeq,ostream & outputFile){
  if (ptr==NULL)return;
  currentSeq.add(ptr->symbTab);
  if (ptr->outputSymb==symb){ // We find a sequence corresponding to the current output symbol
    outputFile << "Symb["<<symb<<"]=";
    for (int i=0; i<currentSeq.getLength();i++){
      outputFile<<"{";
      for (int s=0;s<currentSeq[i].getNbSymb();s++)
	if (currentSeq[i].isIn(s)) outputFile <<s<<" ";
      outputFile<<"}";
    }
    outputFile <<endl;
  }
  _toFile(ptr->ch,symb,currentSeq,outputFile);
  currentSeq.init(currentSeq.getLength()-1);
  _toFile(ptr->br,symb,currentSeq,outputFile);
}
void SequenceDecoder::toFile(ostream & outputFile){
  for(unsigned long int i=0; i<_nbOutputSeq;i++){
    Seq currentSeq(100,_nbInputSymb); 
    currentSeq.init(0); 
    _toFile(_seed,i,currentSeq,outputFile);
  }
}
SequenceDecoderNode * SequenceDecoder::_load(istream &inputFile){
  string tmp;
  if (inputFile.eof()) return NULL;
  inputFile>>tmp;
  if (tmp=="nil") return NULL;
  SequenceDecoderNode *ptr=NULL;
  SequenceDecoderNode* retPtr=NULL;
  bool init=true;
  if (tmp=="begin"){
    while (tmp!="nil"){
      SymbTab sym(_nbInputSymb);
      SequenceDecoderNode * newPtr=_newNode(sym,-1,NULL,NULL);
      newPtr->ch=_load(inputFile);
      inputFile>>newPtr->outputSymb;
      short int s;
      do{
	inputFile>>s;
	if (s!=-1)newPtr->symbTab->setSymb(s);
      } while (s!=-1);
      if (init) {retPtr=newPtr;ptr=newPtr;init=false;}
      else{ptr->br=newPtr;ptr=ptr->br;}
      inputFile>>tmp;
    }
    ptr->br=NULL;
    return retPtr;
  }
  else throw Exception("nil or eof is missing" , __FILE__, __LINE__); 
}
void SequenceDecoder::load(istream &inputFile){
  inputFile>>_nbInputSymb;
  inputFile>>_nbOutputSeqPart;
  inputFile>> _nbOutputSeq;
  _seed=_load(inputFile);
}
void SequenceDecoder::_save(SequenceDecoderNode *ptr,ostream & outputFile){
  if (ptr==NULL) outputFile<<"nil"<<endl;
  else{
    while (ptr!=NULL){
      outputFile<<"begin"<<endl;
      _save(ptr->ch,outputFile);
      outputFile<<ptr->outputSymb<<" ";
      for (int s=0;s<ptr->symbTab->getNbSymb();s++)
	if (ptr->symbTab->isIn(s)) outputFile<<s<<" ";
      outputFile<<"-1"<<endl;
      ptr=ptr->br;
    }
    outputFile<<"nil"<<endl;
  }
}
void SequenceDecoder::save(ostream & outputFile){
  outputFile<<_nbInputSymb<<endl;
  outputFile<<_nbOutputSeqPart<<endl;;
  outputFile<< _nbOutputSeq<<endl;
  _save(_seed,outputFile);
}
// _decode() finds one sequence. currentIdx is on the first symbol to decode. At the end, the input stream is on the next symbol to decode
bool SequenceDecoder::_decode(SequenceDecoderNode *ptr,
			      ReadMemory & inp,ostream & outputFile,unsigned long &begin,short int &decodedSeq){
  if (debug) cout << "_decode["<<inp.getIdx()<<"]=["<<inp.getCurrentSymb()<<"]"<<endl;
  while ((ptr!=NULL)&&(!ptr->symbTab->isIn(inp.getCurrentSymb()))) ptr=ptr->br;
  if (ptr==NULL) return false;
  
  if (ptr->ch==NULL){ //end of the sequence
    if (debug) cout <<"end of sequence detected idx["<<inp.getIdx()<<"]"<<endl;
    if (ptr->outputSymb==-1) cerr<< "Warning: not ended sequence at idx["<<inp.getIdx()<<"]"<<endl;
    else {outputFile <<begin<<" "<<inp.getIdx()<<" "<<ptr->outputSymb<<endl;decodedSeq=ptr->outputSymb;};
    inp.lecture();
    begin=inp.getIdx();
    return true;
  }
  
  if (ptr->outputSymb==-1){ // We are in the sequence but no output at this level: go at the next level
    if (debug) cout <<"no output at this level idx["<<inp.getIdx()<<"], go to the next level"<<endl;
    if (!inp.lecture()) return false;
    return _decode(ptr->ch,inp,outputFile,begin,decodedSeq);
  }
  // If eof, we output the current solution
  if (inp.eof()) {
    outputFile <<begin<<" "<<inp.getIdx()<<" "<<ptr->outputSymb<<endl;decodedSeq=ptr->outputSymb;
    //  inp.lecture();
    // begin=inp.getIdx();
    return true;
  }
  // Now, we have at least one solution, maybe several..., we have to select the longuest one
  unsigned long idx=inp.getIdx();
  if (debug) cout <<"save idx["<<idx<<"]symb="<<inp.getCurrentSymb()<<endl;
  inp.lecture();
  if (inp.eof() || (!_decode(ptr->ch,inp,outputFile,begin,decodedSeq))){   // We didn't find a more longer sequence
    outputFile <<begin<<" "<<idx<<" "<<ptr->outputSymb<<endl;decodedSeq=ptr->outputSymb;
    inp.setIdx(idx);                             // Reseting the input stream at the last good value
    if (debug) cout << "comeback to idx["<<idx<<"]symb="<<inp.getCurrentSymb()<<endl;
    inp.lecture();
    begin=inp.getIdx();
  }
  return true; 
}

void SequenceDecoder::decode(String inputFilename,ostream & outputFile,unsigned long segBegin,unsigned long segLength,
			     bool overlapMode,bool ngramMode,BNgram &ngram) {
  ReadMemory inp(inputFilename,100,segBegin,segLength); // 100 is the buffer size, should be >= the max sequence length. 
                                                        // if length==0 -> segmental mode not used
  if (verbose) cout <<"Begin decoding for seg["<<segBegin<<","<<segLength<<"] overlap:"<<overlapMode<<endl;
  if (!inp.lecture(true)) return;                       // true is used only for the first call, false by default
  unsigned long begin=inp.getIdx();
  bool flag=true;
  while (inp.notEof()){
    unsigned long saveIdx=inp.getIdx();
    short int seq;            // decoded sequence
    bool ret=_decode(_seed,inp,outputFile,begin,seq);
    if (!ret){ 
      cout <<"WARNING, Seq unknown beginning by symb["<<inp.getCurrentSymb()<<"]idx["<<inp.getIdx()<<"]"<<endl;
      if (!inp.lecture()) return;
      begin=inp.getIdx();
    }
    if (ngramMode) if (flag){ ngram.beginSeq(seq); flag=false;} else ngram.addSymb(seq);
    if (overlapMode){ // We want sequences with a max overlap i.e. a sequence is researched for each input symbol
      inp.setIdx(saveIdx);
      if (verbose) cout <<"overlapMode, comeback to symb idx:"<<saveIdx<<endl;
      if (!inp.lecture()) return;
      begin=inp.getIdx();
    }
    if (verbose) cout <<"decode new sequence at idx["<<inp.getIdx()<<"]"<<endl;
  }  
  ngram.endSeq();
}


//*************************************************************************************************************************
//*****           Main functions (interface)
//**************************************************************************************************************************

// sequenceDecoder - Using a sequence decoder tree, this function reads a sequence of input symbols and transcodes it, 
// producing a list of labels begin, end (in frames) sequence number
int sequenceDecoder(Config &config){
  String inputTreeFilename=config.getParam("decoderFilename");  // The complete filename for the decoder tree
  String outputFilename=config.getParam("outputFilename");      // The complete output filename
  String inputFilename=config.getParam("inputFilename");        // The complete input filename 
  String labelSelectedFrames;
  String labelFilename;
  bool label=false;
  if (config.existsParam("labelFilename")){
    labelFilename=config.getParam("labelFilename");
    labelSelectedFrames=config.getParam("labelSelectedFrames");              // label for selected frames - Only the frames from segments 
    label=true;
  }
  bool overlapMode=false;
  if (config.existsParam("overlapMode")) overlapMode=config.getParam("overlapMode").toBool();
  bool ngramMode=false;
  String nOutputF;
  unsigned long oNgramOrder=0;
  if (config.existsParam("ngramOutputFilename")){
    ngramMode=true;nOutputF=config.getParam("ngramOutputFilename");
    oNgramOrder=config.getParam("ngramOutputOrder").toLong();
  }
  BNgram outNgram(oNgramOrder,100);
  SequenceDecoder decTree(0);
  ifstream inputTreeFile(inputTreeFilename.c_str());
  decTree.load(inputTreeFile);
  ofstream outputFile(outputFilename.c_str(),ios::out| ios::trunc);
  if (label){
    SegServer segmentsServer;                                     // Create the segment server for managing the segments/clusters
    LabelServer labelServer;                                      // Create the label server for managing the segment/cluster names
    loadClusterFile(labelFilename,segmentsServer,labelServer,config); // load the label file
    long codeSelectedFrame=labelServer.getLabelIndexByString(labelSelectedFrames);   // Get the index of the cluster with in interest audio segments
    if (codeSelectedFrame==-1) throw Exception("No segment with ["+labelSelectedFrames+"] label" , __FILE__, __LINE__);
    SegCluster& selectedSegments=segmentsServer.getCluster(codeSelectedFrame); // Gives the cluster of the selected/used segments   
    Seg* seg;                                                     // reset the reader at the begin of the input stream
    selectedSegments.rewind(); 
    while((seg=selectedSegments.getSeg())!=NULL){                 // For each of the selected segments
      decTree.decode(inputFilename,outputFile,seg->begin(),seg->length(),overlapMode,ngramMode,outNgram); // decode a given segment
    }
  }
  else decTree.decode(inputFilename,outputFile,0,0,overlapMode,ngramMode,outNgram); // all the frames of the file
  if (ngramMode){
    ofstream outputFileNgram(nOutputF.c_str(),ios::out| ios::trunc);
    outNgram.show(outputFileNgram);
  }
  
  return 0;
}
// sequence extractor builds a common part tree defining a set of variable length symbols with the as same 
// as possible probabilities to occur.
// the input data are the grams, from 1 to maxOrder. Each contains the observed ngram with the corresponding probability
// the required number of symbols (nbSym) should also be defined a priori. It gives pSym=1/nbSym, the target probability for each symbol
// The algo consists in selecting the maximum length sequence (S0,S1,..,Sn) with the maximum probability and to aglomerate it with common 
// part sequences until pSym is reached. I.E, making a symbol by combining (a,b,c,d) with (a,b,c,e) (a,b,c,f)...
// The variable length symbol is then inserted in a common part tree: the result !
// the algo is repeated until the number of symbol is reached
int sequenceExtractor(Config& config){

  int maxOrder=config.getParam("maxOrder").toLong();              // Getting the order of ngrams. Ngrams from 1 to maxOrder should be present
  unsigned long maxNgram=config.getParam("maxNgram").toLong();    // Max ngram by file
  int nbSymb=config.getParam("nbInputSymb").toLong();             // Number of input symbols
  const int nbSeq=config.getParam("nbOutputSymb").toLong();       // Number of sequence to detect
  bool equalInputInfo=false;
  if (config.existsParam("equalInputInfo")) equalInputInfo=config.getParam("equalInputInfo").toBool(); // Work with count*size of the sequence
  String ngramFilename=config.getParam("ngramFilename");          // Base filename for the ngram. No order filename have the form ngramFilenameN.ext 
  String ngramExt=config.getParam("ngramExt");                    // Ext for the ngram files, should include the .
  String outputFilenameInfo=config.getParam("outputInfoFilename");        // The complete output filename for info file
  ofstream outputFileInfo(outputFilenameInfo.c_str(),ios::out| ios::trunc);
  String outputFilename=config.getParam("outputFilename");        // The complete output filename for the decoder tree
  ofstream outputFile(outputFilename.c_str(),ios::out| ios::trunc);
 
  CommonPartTree initialTree;                             // The tree with all the info  
  if (verbose) cout << "reading the ngram files"<<endl;
  for (int order=1;order<=maxOrder;order++){              // Load the ngrams, put it in a common part tree
    NGram nGram(order,maxNgram);
    nGram.load(ngramFilename+String::valueOf(order)+ngramExt,config);
    if (debug) nGram.showTable();
    initialTree.addNGram(nGram);                        // insert all the ngram in the tree 
  }
  initialTree.setTotalCount(initialTree.getTotalChildCount());  // All the symbols are seen in the unigrams
  if (verbose) cout <<"total input symbols:"<<initialTree.getTotalChildCount()<<endl;
  if (debug) {cout<<"initial tree"<<endl; initialTree.show();}

  // Initialize the tree 
  SequenceDecoder seqTree(nbSymb);
  Seq currentSeq(maxOrder,nbSymb); 
  // Begin the main loop
  unsigned long remainingInputSymb=initialTree.getTotalChildCount();
  bool fEnd=false;
  for (int seq=0;(seq<nbSeq) && (!fEnd);seq++){ // Main sequence loop - a sequence is decided at each iteration
    if (verbose) cout << "begin of the main sequence detection loop, seq:"<<seq<<endl;
    unsigned long targetSeqCount=remainingInputSymb/(nbSeq-seq);                  // Count for each variable length symbol
    if (verbose) cout << "totalchildcount"<<initialTree.getTotalChildCount()<<" Remaining input symb:"<<remainingInputSymb<<endl;
    if (verbose) cout <<"target Seq count["<<targetSeqCount<<"]"<<endl;
    currentSeq.init(); // reset the current seq
    unsigned long seqCount=initialTree.findMaxSeq(currentSeq); // Finding the maximum length sequence with the max count 
    if (equalInputInfo) seqCount*=currentSeq.getLength();      // If you want to take into account the number of input symbols and not only the count
    if (currentSeq.getLength()==0) fEnd=true;
    else{ 
      if (verbose){
	cout <<"First Seq["<<seq<<"] Length["<<currentSeq.getLength();
	if (equalInputInfo)cout<<"] Seq input symbols count[";
	else cout <<"] SeqCount[";
	cout<<seqCount<<"]"<<endl;currentSeq.show();
      }
      initialTree.suppressSeq(currentSeq); // suppress the sequence in the initialTree
      if (debug){ cout <<"tree after suppress"<<endl; initialTree.show();}
      seqTree.addSequence(currentSeq,seq);
      if (debug) {cout <<"seq tree after adding initial sequence"<<endl;seqTree.show();}
      short int length=currentSeq.getLength()-1;// -1 because we want to suppress the last level
      while ((seqCount<targetSeqCount)&&(length>=0)){
	bool end=false;
	while ((!end)&&(length>=0)&&(seqCount<targetSeqCount)){
	  currentSeq.init(length);                             // withdrath the info of the end of the sequence
	  if (debug) {cout <<"current path"<<endl;currentSeq.show();}
	  unsigned long deltaSeqCount=initialTree.findMaxEndSeq(currentSeq);
	  if (equalInputInfo) deltaSeqCount*=currentSeq.getLength();                       // If you want to take into account the number of input symbols and not only the count
	  end=(deltaSeqCount==0)||(currentSeq.getLength()==0);
	  if (!end){
	    seqCount+=deltaSeqCount;
	    if (verbose){
	      cout <<"Seq["<<seq<<"] add new seq Length["<<currentSeq.getLength();
	      if (equalInputInfo)cout<<"] Seq input symbols count[";
	      else cout <<"] SeqCount[";
	      cout<<seqCount<<"]"<<endl;currentSeq.show();
	    }
	    initialTree.suppressSeq(currentSeq); // suppress the sequence in the initialTree
	    if (debug){ cout <<"tree after suppress seq"<<endl; initialTree.show();}
	    seqTree.addSequence(currentSeq,seq);
	    if (debug){cout <<"seq tree after adding seq"<<endl;seqTree.show();}
	    length=currentSeq.getLength()-1;
	  }
	  else length--;
	}
      }
      remainingInputSymb-=seqCount;
      if (seqCount==0) fEnd=true;
      else{
	outputFileInfo<<"Sequence["<<seq<<"] Count["<<seqCount<<"]"<<endl;
	seqTree.setNbOutputSeq(seq+1);
      }
    }
  } 
  
  if (verboseLevel>1){cout <<"FINAL SEQUENCE TREE"<<endl; seqTree.show();}
  if (verbose) seqTree.toFile(outputFileInfo);
  seqTree.save(outputFile);
 
  return 0;
}

#endif

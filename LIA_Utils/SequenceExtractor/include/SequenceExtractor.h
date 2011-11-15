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

#if !defined(ALIZE_SequenceExtractor_h)
#define ALIZE_SequenceExtractor_h

#include "alize.h"
#include "LabelNGram.h"
#include "BNGram.h"

class SymbTab{
  short int _nbSymb;
  bool *_symbTab;
 public:
  void init(bool val=false){for (int i=0;i<_nbSymb;_symbTab[i++]=val);}
  void reserveMem();
  void setNbSymb(short int nbSymb){_nbSymb=nbSymb;}
  SymbTab(short int nbSymb){_nbSymb=nbSymb;reserveMem();init(false);}
  SymbTab(){_nbSymb=0;_symbTab=NULL;}
  SymbTab(const SymbTab &);
  ~SymbTab(){if (_symbTab!=NULL) delete [] _symbTab;};
  SymbTab& operator =(const SymbTab & symbTab);
  bool operator==(const SymbTab &);
  bool operator!=(const SymbTab &);
  bool isIn(short int symb);
  void setSymb(short int symb);
  short int getNbSymb(){return _nbSymb;}
  void show();
};

class ReadMemory{
  ifstream _inputFile;
  unsigned long _idx;
  unsigned long _realIdx;
  short int *_buf;
  int _bufSize;
  unsigned long _begin;
  unsigned long _length;
  bool _segmental;
  bool readSymb(short int &, string &);
public:
  ReadMemory(String,int,unsigned long,unsigned long);
  ~ReadMemory();
  bool notEof();
  bool eof();
  bool lecture(bool);
  unsigned long getIdx();
  short int getCurrentSymb();
  void setIdx(unsigned long idx);
};

class Seq{
  short int _maxLength;
  short int _length;
  short int _nbInputSymb;
  SymbTab *_array;
  void _reserveMem();
 public:
  Seq(short int,short int);
  Seq(const Seq&);
  ~Seq();
  Seq & operator=(const Seq &);
  SymbTab &operator[](short int);
  void setLength(short int order);
  short int getLength();
  void init(short int);
  void show();
  void add(SymbTab *ptr);
 
};

struct CommonPartTreeNode{
  short int symb; 
  unsigned long count;
  unsigned long totalChildCount;
  CommonPartTreeNode *ch;
  CommonPartTreeNode *br;
};
class CommonPartTree{
  CommonPartTreeNode * _seed;
  unsigned long _totalCount;
  unsigned long _totalChildCount;
  void _freeTree(CommonPartTreeNode* seed);
  CommonPartTreeNode* _newNode(const short int symb,const unsigned long count,CommonPartTreeNode *ch,CommonPartTreeNode *br);
  CommonPartTreeNode* _findInsert(const short int symb,unsigned long count,CommonPartTreeNode * ptr);
  unsigned long _findMaxSeq(CommonPartTreeNode *ptr,Seq &seq, short int &);
  CommonPartTreeNode* _suppressSeq(CommonPartTreeNode * ptr,Seq &seq,
				  short int order, unsigned long &childCountDelta);
  void _show(CommonPartTreeNode *ptr,unsigned long order);
  CommonPartTreeNode* _findPartSeq(Seq &,short int,CommonPartTreeNode *);
  public:
  CommonPartTree();
  CommonPartTree(const CommonPartTree &);
  CommonPartTree operator =(const CommonPartTree &);
  ~CommonPartTree();
  void addNGram(NGram& nGram);
  unsigned long findMaxSeq(Seq&);
  unsigned long findMaxEndSeq(Seq&);
  void suppressSeq(Seq &seq);
  unsigned long getTotalChildCount(){return _totalChildCount;}
  void setTotalCount(unsigned long total){_totalCount=total;}
  unsigned long getTotalCount(){return _totalCount;}
  void show();
};



struct SequenceDecoderNode{
  SymbTab *symbTab;
  short int outputSymb;
  SequenceDecoderNode *ch;
  SequenceDecoderNode *br;
};

class SequenceDecoder{
  short int _nbInputSymb;          // Number of input symbols, i.e. number of gaussian index
  unsigned long _nbOutputSeqPart;  // Current number of detected sequences (part of output seq)
  unsigned long _nbOutputSeq;      // Nb of different output symbols in the tree
  SequenceDecoderNode *_seed;
  void _freeTree(SequenceDecoderNode* seed);
  SequenceDecoderNode* _newNode(const SymbTab &symbTab,const short int outputSymb,SequenceDecoderNode *ch,SequenceDecoderNode *br);
  SequenceDecoderNode* _findInsert(const SymbTab &symbTab,SequenceDecoderNode *ptr);
  void _show(SequenceDecoderNode *,short int);
  void _toFile(SequenceDecoderNode *,short int,Seq &,ostream &);
  void _save(SequenceDecoderNode *,ostream & );
  SequenceDecoderNode *_load(istream &);
  bool _decode(SequenceDecoderNode *ptr,ReadMemory & inp,ostream & outputFile,unsigned long &begin,short int &);
 public:
  SequenceDecoder(short int);
  SequenceDecoder(const SequenceDecoder &);
  SequenceDecoder operator =(const SequenceDecoder &);
  ~SequenceDecoder();
  short int getNbInputSymb(){return _nbInputSymb;}
  short int getNbOutputSeq(){return (short)_nbOutputSeq;}
  void setNbOutputSeq(unsigned long nb){_nbOutputSeq=nb;}
  unsigned long getNbOutputSeqPart(){return _nbOutputSeqPart;}
  void addSequence(Seq &, const short int);
  void show();
  void toFile(ostream &);
  void save(ostream &);
  void load(istream &);
  void decode(String,ostream &,unsigned long,unsigned long,bool,bool,BNgram &); 
};  

// main functions
int sequenceExtractor(alize::Config&);
int sequenceDecoder(alize::Config &);

#endif // 

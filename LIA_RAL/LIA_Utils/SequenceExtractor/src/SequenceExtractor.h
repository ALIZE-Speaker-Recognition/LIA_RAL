#if !defined(ALIZE_SequenceExtractor_h)
#define ALIZE_SequenceExtractor_h

#include "alize.h"
#include "LabelNGram.h"
#include "BNgram.h"


extern bool debug;
extern bool verbose;

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
  short int getNbOutputSeq(){return _nbOutputSeq;}
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

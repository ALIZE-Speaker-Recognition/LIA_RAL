#if !defined(ALIZE_LabelNGram_h)
#define ALIZE_LabelNGram_h
#include <iostream>
#include <fstream>  // pour outFile
#include "alize.h"
using namespace alize;
using namespace std;
extern bool debug;
extern bool verbose;

const short int OOV=-1;
//-------------------------------------------------------------------------
class NGram{
  unsigned long _order; // Order of the Ngram
  unsigned long _nb;    // number of Ngram 
  short int * _sym;     // Description of the NGram, one after one
  unsigned long *_count; 
  unsigned long _totalCount;
public:
  
  NGram(unsigned long order, unsigned long nb);
  NGram(const NGram &ng);
  ~NGram();
  const NGram& operator=(const NGram& ng);

  unsigned long getOrder(){return _order;}
  unsigned long getSize(){return _nb;}
  void setSize(const unsigned size);
  short int getSymbol(const unsigned idx,const unsigned long o);
  unsigned long getCount(const unsigned idx);
  void setCount(const unsigned idx, const unsigned long& count);
  unsigned long getTotalCount();
  void setTotalCount(const unsigned long& count);
  void setSymbol(const unsigned,const unsigned long, const short int, unsigned long);
  void showTable(ostream &out=cout);
  void load(const String filename,Config &config);

};
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
int labelNGram(alize::Config&);

#endif // !defined(Scoring)

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

#if !defined(ALIZE_NGram_cpp)
#define ALIZE_LabelNGram_cpp
#include <cstdio>   // pour printf()
#include <cassert> // pour le debug pratique
#include <cmath>
#include <liatools.h>
#include "LabelNGram.h"


using namespace alize;
using namespace std;
NGram::NGram(unsigned long order, unsigned long nb){
    _order=order;
    _nb=nb;
    _sym=new short int [_order*_nb];
    _count=new unsigned long [_nb];
    _totalCount=0;
  }
NGram:: ~NGram(){
  _order=0;
  _nb=0;
  delete [] _sym;
  delete [] _count;
  }

NGram::NGram(const NGram &ng){ 
  _order=ng._order;
  _nb=ng._nb;
  _totalCount=ng._totalCount;
  _sym=new short int [_order*_nb];
  _count=new unsigned long [_nb];
  memcpy(_sym, ng._sym, _nb*_order*sizeof(short int)); 
  memcpy(_count, ng._count, _nb*sizeof(unsigned long)); 
}
const NGram& NGram::operator=(const NGram& ng)
{
  if (this==&ng)
    return ng;
  if ((ng._nb*ng._order)!=(_nb*_order)){
    delete _sym;
    delete _count;
    _sym=new short int [ng._order*ng._nb];
    _count=new unsigned long [ng._nb];
  }
  _nb=ng._nb;
  _order=ng._order;
  _totalCount=ng._totalCount;
  memcpy(_sym, ng._sym, _nb*_order*sizeof(short int)); 
  memcpy(_count, ng._count, _nb*sizeof(unsigned long)); 
  return *this;
}
void NGram::setSize(const unsigned size){
  if (size<=_nb) _nb=size;
  else   throw Exception("Resize is allowad only for reducing the size"
			 , __FILE__, __LINE__);
  short int * sym=new short int [_order*_nb];
  unsigned long *count=new unsigned long [_nb];
  memcpy(sym,_sym, _nb*_order*sizeof(short int)); 
  memcpy(count,_count, _nb*sizeof(unsigned long)); 
  delete _sym;
  delete _count;
  _sym=sym;
  _count=count;
}
 
short int NGram::getSymbol(const unsigned idx,const unsigned long o){
  if ((idx<0) || (idx>=_nb))
      throw Exception("out of array"
        , __FILE__, __LINE__);
  return _sym[(idx*_order)+o];
}
unsigned long NGram::getCount(const unsigned idx){
  if ((idx<0) || (idx>=_nb))
      throw Exception("out of array"
        , __FILE__, __LINE__);
  return _count[idx];
}
void NGram::setCount(const unsigned idx, const unsigned long &count){
  if ((idx<0) || (idx>=_nb))
      throw Exception("out of array"
        , __FILE__, __LINE__);
  _count[idx]=count;
}
unsigned long NGram::getTotalCount(){
  return _totalCount;
}
void NGram::setTotalCount(const unsigned long &count){
  _totalCount=count;
}
void NGram::setSymbol(const unsigned idx,const unsigned long o, const short int sym, unsigned long count=0){
  if ((idx<0) || (idx>=_nb))
      throw Exception("out of array"
        , __FILE__, __LINE__);
 _sym[(idx*_order)+o]=sym; 
}
void NGram::showTable(ostream &out){
  for(unsigned idx=0;idx<_nb;idx++){
    out<<"Sym["<<idx<<"]=";
    for (unsigned long s=0;s<_order;s++)
      out<<"["<<getSymbol(idx,s)<<"]"; 
    out<<"count["<<getCount(idx)<<"]"<<endl;
  }
}
    

//-------------------------------------------------------------------------
// Load the NGRAM table, selecting the nbSelected first
void NGram::load(const String filename,Config &config){
  XList input(filename,config);
  XLine *linep;                                                          // Pointor on the current test line
  input.getLine(0);
 
  unsigned long idx=0;
  while (((linep=input.getLine()) != NULL)&&(idx<getSize())){
    for (unsigned long i=0;i<getOrder();i++){
      short int a=linep->getElement(i).toLong();
      setSymbol(idx,i,a);
    }
     if (linep->getElementCount()==(getOrder()+1)){
	unsigned long count=linep->getElement(getOrder()).toLong();
	setCount(idx,count);_totalCount+=count;}
     else setCount(idx,0);
    idx++;
  }
  if (idx!=getSize()){
    cout << "WARNING ! Number of ngram in the file["<<idx<<"] < to the number requested ["<<getSize()<<"]"<<endl;
    setSize(idx);
  }
  if (verboseLevel>1){
    cout <<"load symbol table from ["<<filename <<"]"<<endl;
    showTable();
  }
}
bool isNGram(short int *sym,NGram &tabS,unsigned long & tag){      
  bool find=false;
  unsigned long idx;
  for (idx=0;(!find) && (idx<tabS.getSize());idx++){
    find=true;
    for (unsigned long s=0;(find) && (s<tabS.getOrder());s++)
      find=(sym[s]==tabS.getSymbol(idx,s));
  }
  if (find){
    tag=idx;
    return true;
  }
  else return false;
}

short int recognizeSymbol(unsigned long &idxFrame,unsigned long end,ULongVector &tabS){
  unsigned long sym=tabS[idxFrame];
  while ((idxFrame<end)&&(tabS[idxFrame]==sym))idxFrame++;
  return sym;
}
//-------------------------------------------------------------------------
// Should changed for  a circulate array
void moveTab(unsigned long *begin,short int *sym,unsigned long *end,unsigned long order){
  for (unsigned long i=0;i<order-1;i++){
    begin[i]=begin[i+1];
    end[i]=end[i+1];
    sym[i]=sym[i+1];
  }      
}
//-------------------------------------------------------------------------
void computeLabelNGram(NGram & NG,SegCluster &cluster,SegCluster &clusterOut,ULongVector &tabS,unsigned long nbSym){
  unsigned long begin[100]; // max order 100...
  short int sym[100];
  unsigned long end[100];
  SegServer & segServerOut=clusterOut.getServer();                           // Get the clusterserver reelated to the output
  cluster.rewind();
  Seg* seg;                                                                  // Reset the reader at the begin of the input stream
  while((seg=cluster.getSeg())!=NULL){                                       // For each of the selected segments
    unsigned long idxFrame=seg->begin(); 
    unsigned long endS=endSeg(seg);
    if (endS>=nbSym) endS=nbSym;                                          // Just if there is less symbol in the file than in the label 
    if (idxFrame>endS) idxFrame=endS;
    unsigned long beginOOV=idxFrame;
    bool oov=true;
    if (debug) cout <<"begin Seg["<<idxFrame<<"]"<<endl;

    for (unsigned long n=0;(idxFrame<endS) &&(n<NG.getOrder()-1);n++){                      // Recognize the (n-1) first symbols
      begin[n]=idxFrame;
      sym[n]=recognizeSymbol(idxFrame,endS,tabS);
      end[n]=idxFrame-1;
      if (debug) cout <<"sym ["<<sym[n]<<"] begin["<<begin[n]<<"] end["<<end[n]<<"] idxframe["<<idxFrame<<"]"<<endl;
    }
    while(idxFrame<endS){
      begin[NG.getOrder()-1]=idxFrame;
      sym[NG.getOrder()-1]=recognizeSymbol(idxFrame,endS,tabS);
      end[NG.getOrder()-1]=idxFrame-1;
      if (debug) cout <<"sym ["<<sym[NG.getOrder()-1]<<"] begin["<<begin[NG.getOrder()-1]
		      <<"] end["<<end[NG.getOrder()-1]<<"] idxframe["<<idxFrame<<"]"<<endl;
      unsigned long tag;
      if (isNGram(sym,NG,tag)){
	if ((oov)&&(beginOOV<begin[0])){
	  if (debug) cout <<"OOV1  begin["<<beginOOV <<"] end["<<begin[0]-1<<"]"<<endl;
	  Seg &segTmp=segServerOut.createSeg(beginOOV,begin[0]-beginOOV,0,"oov",seg->sourceName());
	  clusterOut.add(segTmp);
	}
	if (debug) cout <<"NGRAM ["<<tag<<"] begin["<<begin[0] <<"] end["<<end[NG.getOrder()-1]<<"]"<<endl;
	Seg &segTmp=segServerOut.createSeg(begin[0],end[NG.getOrder()-1]-begin[0]+1,0,String::valueOf(tag),seg->sourceName());
	clusterOut.add(segTmp);
	beginOOV=idxFrame;
	oov=false;
      }
      else oov=true;
      moveTab(begin,sym,end,NG.getOrder());
    }
    if (oov){
      Seg &segTmp=segServerOut.createSeg(beginOOV,idxFrame-beginOOV,0,"oov",seg->sourceName());
      clusterOut.add(segTmp);
      if (debug) cout <<"OOV2  begin["<<beginOOV <<"] end["<<idxFrame-1<<"]"<<endl;  
    }
  }
}


//-------------------------------------------------------------------------
unsigned long loadSymbol(const String &filename,const String &type,ULongVector & ret,Config &config){
 
  if (type=="ascii"){
    unsigned long nbSym=0;
    XList infile(filename,config);
    XLine list=infile.getAllElements();
    ret.setSize(list.getElementCount());
    nbSym=list.getElementCount();
    for (unsigned long i=0;i<nbSym;i++){
      String *tmp=list.getElement();
      if ((*tmp)!="oov")
	ret[i]=tmp->toLong();
      else ret[i]=OOV;
    }
    if (debug) cout << "nb sym:"<<nbSym<<endl;
    return nbSym;
  }
  else throw Exception(type+" file type non recognized for a symbol file"
        , __FILE__, __LINE__);
}  

//-------------------------------------------------------------------------
int labelNGram(Config& config)
{
  if (config.existsParam("debug"))debug=true; else debug=false;  
  if (config.existsParam("verbose"))verbose=true; else verbose=false;
  String extOutputLabel=".sym.lbl";                                               // the extension of the output files    
  if (config.existsParam("saveLabelFileExtension")) extOutputLabel=config.getParam("saveLabelFileExtension");   
  String pathOutput="./";                                                // the path of the output files    
  if (config.existsParam("labelOutputPath")) pathOutput=config.getParam("labelOutputPath");    
  String extSymbol=".sym";                                               // the extension of the symbol files   
  if (config.existsParam("symbolFileExtension")) extSymbol=config.getParam("symbolFileExtension");   
  String pathSymbol="./";   
  if (config.existsParam("symbolPath")) pathSymbol=config.getParam("symbolPath");   
  String formatSymbol="ascii";   
  if (config.existsParam("symbolFormat")) pathSymbol=config.getParam("symbolFormat");  
 
  String NGramFilename=config.getParam("NGramFilename");                                        
  unsigned long NGramOrder=3;
  if (config.existsParam("NGramOrder")) NGramOrder=config.getParam("NGramOrder").toLong();  
  unsigned long NGramSelected=16;
  if (config.existsParam("NGramSelected")) NGramSelected=config.getParam("NGramSelected").toLong();  
  NGram NGramTable(NGramOrder,NGramSelected);
  NGramTable.load(NGramFilename,config);                       // Load the NGRAM table, selecting the NGramSelected first
 
  String inputFilename=config.getParam("inputFilename");
  String labelSelectedFrames=config.getParam("labelSelectedFrames");
  XLine inputFileList;
  try{
    if (inputFilename.endsWith(".lst")){ // input is file containing a list of filenames
      XList tmp(inputFilename,config);
      inputFileList=tmp.getAllElements();
    }
    else inputFileList.addElement(inputFilename); // a single filename
    String *p;
    while ((p=inputFileList.getElement())){
      String& filename=*p;
      if (verbose)
	cout <<"labelNGram file["<<filename<<"] Table["<<NGramFilename<<"] Order["<<NGramOrder<<"] Selected["<<NGramSelected<<"]"<<endl;
      SegServer segServer;                
      LabelServer labelServer;
      loadClusterFile(filename,segServer,labelServer,config);
      long codeSelectedFrame=labelServer.getLabelIndexByString(labelSelectedFrames);       // Get the index of the selected cluster
      if (codeSelectedFrame==-1){                                                             // No data for this model !!!!!!!!!!!!!!
      cout << " WARNING - NO DATA with the label["<<labelSelectedFrames<<"] in file ["<<filename<<"]"<<endl;
      exit(0);
      }
      SegCluster& cluster=segServer.getCluster(codeSelectedFrame);                                   // Gives the cluster of the selected/used segments
      ULongVector tabS;
      unsigned long nbSym=loadSymbol(pathSymbol+filename+extSymbol,formatSymbol,tabS,config);        // Read the stream of symbols
      SegServer segServerOutput;
      SegCluster& clusterOut=segServerOutput.createCluster(0,labelSelectedFrames,cluster.sourceName());  
      //  
      computeLabelNGram(NGramTable,cluster,clusterOut,tabS,nbSym);
      //
    
      if (verbose){
	cout <<"File["<<filename<<"]" <<endl;
	cout << "Output the new label file in ["<<pathOutput+filename+extOutputLabel <<"]"<<endl;
      }
      outputLabelFile(clusterOut,pathOutput+filename+extOutputLabel,config);
    } // end file loop
  } // fin try
  
  
  catch (Exception& e)
    { 
      cout << e.toString().c_str() << endl;
    }
  return 0;
}



#endif

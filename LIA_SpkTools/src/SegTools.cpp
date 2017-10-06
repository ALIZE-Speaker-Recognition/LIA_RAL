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

#if !defined(ALIZE_SegTools_cpp)
#define ALIZE_SegTools_cpp

#include "SegTools.h"
#include<iostream>
#include<fstream>
#include<cstdio>
#include<cassert>
#include<cmath>

using namespace alize;
using namespace std;


// Display the cluster
void showCluster(SegCluster& cluster){
  cluster.rewind();
  Seg *seg;                                                                              
  while((seg=cluster.getSeg())!=NULL)
    cout<<"Seg File["<<seg->sourceName()<<"]["<<seg->begin()<<","<<seg->length()<<"]"<<endl;
}
// return the number of frame in the cluster
unsigned long totalFrame(SegCluster& cluster)
{
  unsigned long time=0;
  cluster.rewind();
  Seg *seg;                                                                              
  while((seg=cluster.getSeg())!=NULL)
    time+=seg->length();
  return time;
}

// return the number of frame in the server of segment
unsigned long totalFrame(SegServer& s)
{
unsigned long time=0;
for(unsigned long icluster=0; icluster<s.getClusterCount(); icluster++){
  SegCluster& cluster=s.getCluster(icluster);
  time+=totalFrame(cluster);
}
return time;
}

// return the last frame of a segment
unsigned long endSeg(Seg *seg)
{
  return seg->begin()+seg->length()-1;
}
//
//------------------------------------------------------------------------------------------------------------------------------
// Output a cluster in a  label file
void outputLabelFile(SegCluster &selectedSeg,String FileName,Config &config)
{
  ofstream outFile(FileName.c_str(),ios::out | ios::trunc);                      // Initialise the output file
  if (!outFile.good()) 
	  throw Exception("Can't create label file",__FILE__,__LINE__);
  real_t frameLength=config.getParam("frameLength").toDouble();
  selectedSeg.rewind();                                                      // Reset the cluster of selected segments
  Seg * seg;
  while((seg=selectedSeg.getSeg())!=NULL){                                       // For all the segments
    outFile <<frameIdxToTime(seg->begin(),frameLength)<<" "<<frameIdxToTime(seg->begin()+seg->length()-1,frameLength)<<" "<<seg->string()<<endl;
    if (debug)cout <<frameIdxToTime(seg->begin(),frameLength)<<" "<<frameIdxToTime(seg->begin()+seg->length()-1,frameLength)<<" "<<seg->string()<<endl;
  }
}



//******* Create a default cluster with the label labelToAnalyse, including all the input file frames
//** From N Scheffer & JFB
void createDefSeg (SegServer & segServer, FeatureServer & fs, Config & config)
{
  unsigned long segLength=fs.getFeatureCount();
  unsigned long begin=0;   
  // labelCode for selected frame
  unsigned long codeSelectedFrame=1;
  segServer.createSeg(begin,segLength,codeSelectedFrame);
}

// From a time to a frame index, taking into account possible bugs on boundaries
// JFB
unsigned long timeToFrameIdx(real_t time,real_t frameLength)
{
  double tmpd,res;
  res=modf(time/frameLength,&tmpd);
  if (res>0.99999)
    return (unsigned long)tmpd+1;
  else return (unsigned long)tmpd;
}
real_t frameIdxToTime(unsigned long idx,real_t frameLength)
{
  unsigned long t= (unsigned long)(idx*1000*frameLength);
  real_t res=(real_t)t/1000.0;
  return res;
}
//**************************************
// Verify the labels - if a segment is finishing after the last frame, trunc it
void verifyClusterFile(SegServer& segmentsServer,FeatureServer& fs,Config& config)
{
  for (unsigned long clusterIdx=0;clusterIdx<segmentsServer.getClusterCount();clusterIdx++){  // For each cluster
    if (debug) cout << "(SegTools) verifying label for ["<<clusterIdx<<"]"<<endl; 
    SegCluster& cluster=segmentsServer.getCluster(clusterIdx);                                // Get the cluster
    if (debug) cout << "(SegTools) verifying label for ["<<cluster.string()<<"]"<<endl; 
    Seg* seg;
    cluster.rewind();  
    while((seg=cluster.getSeg())!=NULL)                                                       // For each segment
      if ((seg->begin()+seg->length())>fs.getFeatureCountOfASource(seg->sourceName())){ // The end of the segment is after the last frame
		unsigned long newLength=fs.getFeatureCountOfASource(seg->sourceName())-seg->begin();
		if (debug|| verbose)cout<<"Warning File["<< seg->sourceName()<<"], Truncate Segment begin:"<<seg->begin()<<" length:"<<seg->length()<<
		      " file size:"<<fs.getFeatureCountOfASource(seg->sourceName())<<" new length:"<<newLength<<endl;
		seg->setLength(newLength);		      
      }
  } //End cluster loop
}

//***************************************
//* Reading a label file and adding the segments in the server. Each label/name is in
//* a separate cluster
//* Dan Istrate, JFB,  november 2004
//***************************************


// Add a segment to the server, in an existing cluster - same name or in a new cluster - new name
void addSegment(const Config & config, 
		const String &file,const String &name, unsigned long segFrameBegin, unsigned long segFrameLength,
		LabelServer& labelServer,SegServer& segmentsServer)
{
  unsigned long codeLabel,labelCount;     
  if (debug ||(verboseLevel > 2)) cout << "add a segment, name["<<name<<"] begin frame["<<segFrameBegin<<"] length in frame["<<segFrameLength<<"]"<<endl;
  labelCount=labelServer.size();
  Label labelTemp(name);
  codeLabel=labelServer.addLabel(labelTemp);
  if(codeLabel==labelCount){ // Create a new cluster Cluster
      SegCluster& cluster=segmentsServer.createCluster(codeLabel,name," ");
      cluster.add(segmentsServer.createSeg(segFrameBegin,segFrameLength,labelServer.addLabel(labelTemp),name,file));
    }
  else{// the cluster exists, it is first extract and the segment is added
    SegCluster& cluster=segmentsServer.getCluster(codeLabel);
    cluster.add(segmentsServer.createSeg(segFrameBegin,segFrameLength,labelServer.addLabel(labelTemp),name,file));
  }
}

// Load the segments of a file in one cluster
void  loadLabelFile(SegCluster &cluster,String fileName,String path, String extension,Config &config)
{
  if (verboseLevel>1) cout << ">> Proceeding label file reading in one cluster for ["<<fileName<<"] <<"<<endl;
  double frameLength = config.getParam("frameLength").toDouble();         // length in s of a frame
  String fileLabel=path+fileName+extension;         // build the complete filename
  XList labelSet(fileLabel,config);
  XLine *segmentp;
  while ((segmentp=labelSet.getLine()) != NULL){                             // Get the lines/segment one by one
    const String& beginString=segmentp->getElement(0);                     // Get the begin time
    const String& endString=segmentp->getElement(1);                       // Get the end time
    const String& name=segmentp->getElement(2);                            // Get the label
    unsigned long segFrameBegin=timeToFrameIdx(beginString.toDouble(),frameLength); // begin in frame
    unsigned long segFrameLength=timeToFrameIdx(endString.toDouble(),frameLength)-segFrameBegin+1;
    cluster.add(cluster.getServer().createSeg(segFrameBegin,segFrameLength,0,name,fileName));
  }  
}


// the main function for loading the labels - do the job for a given file
void loadClusterFile(String &fileName,SegServer& segmentsServer,LabelServer&
labelServer,Config& config)
{
   String labelPath="./";
  if (config.existsParam("labelFilesPath")) labelPath= config.getParam("labelFilesPath");                     // The path of the lable files
  String labelFileExtension=".lbl";	
  if(config.existsParam("labelFilesExtension"))	 // The extension of the lable files
    labelFileExtension=config.getParam("labelFilesExtension");	 
  double frameLength = config.getParam("frameLength").toDouble();         // length in s of a frame
  String inputLabelFormat="LIARAL";
  if (config.existsParam("inputLabelFormat"))                            // the label file format is decided by the user (not LIARAL std format)
    inputLabelFormat=config.getParam("inputLabelFormat");             
  if(debug) cout<<"(SegTools) The label format is "<<inputLabelFormat<<endl;
  if (verboseLevel>1)  cout << "(SegTools) >> Proceeding label file reading for ["<<fileName<<"] <<"<<endl;
  //  double tmpt;	
  bool addDefault=false;
  String defaultLabel;
  if  (config.existsParam("addDefaultLabel")){
    addDefault=true;
    defaultLabel=config.getParam("defaultLabel");
  }
  String fileLabel;
  fileLabel=labelPath+fileName+labelFileExtension;         // build the complete filename - to be included in ALIZE default - TODO
  

  XList labelSet;
  bool noLabelFile=false;
  try{
    labelSet.load(fileLabel,config);                 // if the file exists, read it in a Xlist - one segment by line
  }
  catch(FileNotFoundException& e){
    noLabelFile=true;
  }
  if ((noLabelFile)&& (addDefault==false)){
    cout<<"(SegTools)  Error: For "<< fileName << "["<<fileLabel<<"]"<<endl;
    throw Exception("There is no correspondent label file",__FILE__,__LINE__);		
  }
  if (noLabelFile){ // no label file and the default option is on - creating a default segment
    FeatureServer fsFake(config,fileName);
    unsigned long segFrameLength=fsFake.getFeatureCount();
    unsigned long segFrameBegin=0;   
    if (verbose) cout << "(SegTools) Default label file for ["<<fileName<<"], nb features["<<segFrameLength<<"]"<<endl;
    addSegment(config,fileName,defaultLabel,segFrameBegin,segFrameLength,labelServer,segmentsServer); 
  } 
  else{ // a Label file is existing
    bool visit=false;
    XLine *segmentp;
    while ((segmentp=labelSet.getLine()) != NULL){                             // Get the lines/segment one by one
      visit=true;
      if (inputLabelFormat=="LIARAL"){                                         // A line is composed by [begin_time end_time name]
	const String& beginString=segmentp->getElement(0);                     // Get the Client ID (id)
	const String& endString=segmentp->getElement(1);                       // Get the Client ID (id)
	const String& name=segmentp->getElement(2);                            // Get the Client ID (id)
	unsigned long segFrameBegin=timeToFrameIdx(beginString.toDouble(),frameLength); // begin in frame
	unsigned long segFrameLength=timeToFrameIdx(endString.toDouble(),frameLength)-segFrameBegin+1;
	addSegment(config,fileName,name,segFrameBegin,segFrameLength,labelServer,segmentsServer); 
      }
      else if (inputLabelFormat=="mdtm"){// mdtm format. An example[20030418_0700_0800_FRANCEINTER_DGA 1 0.000 5.928 speaker NA male Jean-Jacques_Bernard]
	const String& name=segmentp->getElement(7);  
	const String& beginString=segmentp->getElement(2); 
	const String& lengthString=segmentp->getElement(3);
	unsigned long segFrameBegin=timeToFrameIdx(beginString.toDouble(),frameLength);                               // begin in frame
	unsigned long segFrameLength=timeToFrameIdx(lengthString.toDouble(),frameLength);                             // length in frame
	addSegment(config,fileName,name,segFrameBegin,segFrameLength,labelServer,segmentsServer); 	
      }
      else cout << "label type["<<inputLabelFormat<<"] unknown" <<endl;
    }
    if (!visit) cout <<"WARNING, label file ["<<fileName <<"] empty"<<endl;
  }
}


// for a list of file
void initializeClusters(const XLine& listFiles,SegServer& segmentsServer,LabelServer& labelServer,Config& config)
{
  listFiles.rewind();
  String *filep;
  while ((filep=listFiles.getElement()) != NULL){   // for each label file
    String file=(*filep);                           // get the filename 
    loadClusterFile(file,segmentsServer,labelServer,config);
  }//while seg 
}
// For a list of input files, stored into a XList
void initializeClusters(const XList& listXFiles,SegServer& segmentsServer,LabelServer& labelServer,Config& config)
{
  XLine &listFiles=listXFiles.getAllElements(); 
  initializeClusters(listFiles,segmentsServer,labelServer,config);
}
// For a filename parameter, which could be a feature file or a feature file list with an .lst extention
void initializeClusters(String &file,SegServer& segmentsServer,LabelServer& labelServer,Config&
config)
{
   if (file.endsWith(".lst")){                                                                        // If the file parameter is the name of a XList file
    XList listXFiles(file,config);                                 // (a text file with a list of filenames)
    initializeClusters(listXFiles,segmentsServer,labelServer,config);
  }
  else                                                                                               // It is a single filename/basename for feature and label files
    loadClusterFile(file,segmentsServer,labelServer,config); 
}


/************
* Functions mainly used by the LIA_SpkSeg package
************/


// A function for create and initialize a Viterbi accumulator
ViterbiAccum& createAndInitializeViterbiAccum(StatServer &ss, hmm &cHmm){
  ViterbiAccum& va=ss.createViterbiAccum(); // CREATE 
  for(unsigned int i=0;i<cHmm.getNbState();i++) // Add the states
    va.addState(cHmm.getDensity(i));
  for(unsigned int i=0;i<cHmm.getNbState();i++) // Set the transitions
    for(unsigned int j=0;j<cHmm.getNbState();j++){
      va.logTransition(i,j)=log(cHmm.getTransition(i,j));
    }
  return va;
}

// accumulateStatViterbi() is used for computing the viterbi path
// TAKE CARE: THE accumulator should be reseted before the call
void   accumulateStatViterbi(FeatureServer &fs,ViterbiAccum &va,unsigned long beginIdx,unsigned long length,Config &config)
{
  double fudge=1.0;
  if(config.existsParam("fudge")) fudge=config.getParam("fudge").toDouble();
  unsigned long viterbiBufferLength=config.getParam("viterbiBufferLength").toLong();
  for(unsigned long ifeature=beginIdx;ifeature<=beginIdx+length-viterbiBufferLength;ifeature+=viterbiBufferLength)
    va.computeAndAccumulate(fs, ifeature, viterbiBufferLength,fudge);
  if (length%viterbiBufferLength!=0)
    va.computeAndAccumulate(fs,beginIdx+(length/viterbiBufferLength*viterbiBufferLength), length%viterbiBufferLength,fudge); // End of the segment
}
// On a segment
void accumulateStatViterbi(FeatureServer &fs,ViterbiAccum &va,Seg *seg,Config &config)
{  
  unsigned long begin=seg->begin()+fs.getFirstFeatureIndexOfASource(seg->sourceName());              // Find the index of the first frame of the file in the buffer
  accumulateStatViterbi(fs,va,begin,seg->length(),config);
}

// On a segment with PreSeg
// Dan Istrate juillet 2005
void accumulateStatViterbi(FeatureServer &fs,ViterbiAccum &va,Seg *seg,Config &config,SegServer& TransitSegServer,DoubleVector&
TransitionsFort,DoubleVector& TransitionsFaible,unsigned long NbState)
{  
  unsigned long begin=seg->begin()+fs.getFirstFeatureIndexOfASource(seg->sourceName());              // Find the index of the first frame of the file in the buffer
  double fudge=config.getParam("fudge").toDouble();
  unsigned long viterbiBufferLength=config.getParam("viterbiBufferLength").toLong();
  bool drapeauFort=false; 
  SegCluster& clusterFort=TransitSegServer.getCluster(0);
  unsigned long lc;
 
  for(unsigned long ifeature=begin;ifeature<=begin+seg->length()-viterbiBufferLength;ifeature+=viterbiBufferLength)
  {
    if(clusterFort.getFeatureLabelCode(ifeature+(viterbiBufferLength/2),lc)&& !drapeauFort) 
    {//passage en transitions fortes
    	for(unsigned long i=0;i<NbState;i++) // Set the transitions
		    for(unsigned long j=0;j<NbState;j++){
		      va.logTransition(i,j)=log(TransitionsFort[i+j*NbState]);
    		}	
	drapeauFort=true;
    }
    if(!clusterFort.getFeatureLabelCode(ifeature+(viterbiBufferLength/2),lc)&& drapeauFort) 
    {//passage en transitions faible
    	for(unsigned long i=0;i<NbState;i++) // Set the transitions
		    for(unsigned long j=0;j<NbState;j++){
		      va.logTransition(i,j)=log(TransitionsFaible[i+j*NbState]);
    		}	
	drapeauFort=false;
    }
    va.computeAndAccumulate(fs, ifeature, viterbiBufferLength,fudge);  
  }
  if (seg->length()%viterbiBufferLength!=0)
    va.computeAndAccumulate(fs,begin+(seg->length()/viterbiBufferLength*viterbiBufferLength), seg->length()%viterbiBufferLength,fudge); // End of the segment
}
    
// Remove all the segments of a cluster
void removeAllSeg(SegServer &segServer,SegCluster &cluster)
{
  while(cluster.getCount()!=0){
    segServer.remove(cluster.get(0));		
  }
}


// initialize an existing segmentation by deleting all the segments and creating an empty cluster by state/speaker
void initializeCluster(SegServer &currentSeg,hmm& cHmm,LabelServer &labelServer){

unsigned long cl;
Label l;
currentSeg.removeAllClusters(); // suppress the clusters (one by speaker/state)
currentSeg.removeAllSegs();     // suppress the segments 
for(unsigned long i=0;i<cHmm.getNbState();i++){
	l.setString(cHmm.getStateName(i));
	cl=labelServer.addLabel(l);	
	currentSeg.createCluster(cl,cHmm.getStateName(i)," ");	
}
}


// removePartOfSeg() removes a part of a segment from a cluster
bool removePartOfSeg(SegCluster &clusterSeg,SegServer& segServer,LabelServer &labelServer,
		     String segSourceName,unsigned long segStart,unsigned long segLength){

clusterSeg.rewind();
  
for(unsigned long i=0;i<segServer.getSegCount();i++){
	Seg& seg=segServer.getSeg(i);

	if((seg.sourceName()==segSourceName)&&(segStart>=seg.begin())&&(segStart+segLength<=seg.begin()+seg.length())){ // Find !!
		unsigned long start=seg.begin();
		unsigned long longueur=seg.length();
		if(start==segStart){                                          // Same begin                      
			if((longueur-segLength)!=0){
				seg.setBegin(start+segLength);               
				seg.setLength(longueur-segLength);
			}else{
				segServer.remove(seg);
				return true;
			}//else
      		}
		else 
			if((segStart+segLength)==(start+longueur)){              // Same end
				seg.setLength(longueur-segLength);
			}
	      		else{
				seg.setLength(segStart-start+1);
				clusterSeg.add(segServer.createSeg((int)segStart+segLength-1,start+longueur-segStart-segLength+1,
					   labelServer.addLabel(seg.string()),seg.string(),seg.sourceName()));
	      		}
	      return true;
	}//if
    
}//for
return false;
}


void displayAllClusters(Config& config, SegServer& seg)
{

if(seg.getClusterCount()!=0)
	for(unsigned long icluster=0;icluster<seg.getClusterCount();icluster++)
	{
		SegCluster& clusterT=seg.getCluster(icluster);
		displayOneCluster(config, clusterT);
	}
else cout<<"No cluster in Segments Server"<<endl;
}//display


void displayOneCluster(Config& config, SegCluster &clusterT)
{
  Seg* segment;
  double frameLength=config.getParam("frameLength").toDouble();
  clusterT.rewind();
  if(clusterT.string() != ""){
	  cout << "Name: " << clusterT.string() << endl;
  }
  else{
	  cout << "Name: empty" << endl;  
  }
  if(clusterT.getCount()!=0)
  	while(((segment=clusterT.getSeg())!=NULL))
	{
		cout << segment->begin()*frameLength << " " << (segment->begin()+segment->length()-1)*frameLength << " " << segment->string() << endl;
	}
  else cout<<"No segments in cluster "<<clusterT.string()<<endl;
}


void displayAllSegments(Config& config, SegServer& seg)
{
double frameLength=config.getParam("frameLength").toDouble();
unsigned long frameCnt=0;
unsigned long iSeg=0;

if(seg.getSegCount()==0){
	cout<<"No segments in Segments Server"<<endl;
	return;
}

// Display of current segmentation in color graphics
unsigned long fileSize = config.getParam("fileSize").toLong();

while(iSeg < seg.getSegCount()){
	Seg& segment=seg.getSeg(iSeg);
	unsigned long b=segment.begin();
	
	// gap between two segments
	while((frameCnt < fileSize) && (frameCnt < b)){
		cout << "\033[1m\033[24m\033[0mX";
		frameCnt += 30;
	}
	while((frameCnt < fileSize) &&(frameCnt < (b+segment.length()-1))){
		String lab = segment.string();
		String s = lab[lab.length()-1];
		try{
			unsigned long c=33+s.toLong();
			cout << "\033[1m\033[24m\033[" << c << "m" << s;
			frameCnt += 30;		
		}
		catch(Exception& e){
			cout << "\033[1m\033[24m\033[0m" << s;
			frameCnt += 30;		
		}
	}
	iSeg++;
}
while(frameCnt < fileSize){
	cout << "\033[1m\033[24m\033[0mX";
	frameCnt += 30;
}

cout << "\033[0m\n\n" << endl << endl;

// Display of current segmentation in text	
if(verboseLevel == 2){
	for(unsigned long iseg=0;iseg<seg.getSegCount();iseg++){
		Seg& segment=seg.getSeg(iseg);
		cout << segment.begin()*frameLength << " " << (segment.begin()+segment.length()-1)*frameLength << " " << segment.string() << endl;
	}
}			

// if desired, display of reference segmentation
if(config.existsParam("refSegDisplay")){
	cout << config.getParam("refSegDisplay") << endl << endl;
}	
}//display


void displayAllSegmentsFromRef(Config& config, String &fileInit, unsigned long fileSize)
{
double frameLength=config.getParam("frameLength").toDouble();
unsigned long frameCnt=0;
unsigned long iSeg=0;
String refSegDisplay = "";

SegServer refSegServer;
LabelServer refLabel;

// Reference file reading
if(verbose) cout << ">> Reference Segmentation for comparison purpose <<" << endl;
String fileRefPath=config.getParam("fileRefPath");
String lblPath=config.getParam("labelFilesPath");
config.setParam("labelFilesPath",fileRefPath);
initializeClusters(fileInit,refSegServer,refLabel,config);                   // Reading the segmentation files for each feature input file
config.setParam("labelFilesPath",lblPath);

if(refSegServer.getSegCount()==0){
	cout<<"No segments in Segments Server"<<endl;
	return;
}


while(iSeg < refSegServer.getSegCount()){
	Seg& segment=refSegServer.getSeg(iSeg);
	unsigned long b=segment.begin();
	
	// gap between two segments
	while((frameCnt < fileSize) && (frameCnt < b)){
		refSegDisplay += "\033[1m\033[24m\033[0mX";
		frameCnt += 30;
	}
	while((frameCnt < fileSize) &&(frameCnt < (b+segment.length()-1))){
		unsigned long s = segment.labelCode();
		unsigned long c=33+s;
		refSegDisplay += "\033[1m\033[24m\033["+ refSegDisplay.valueOf(c)+"m"+refSegDisplay.valueOf(s);
		frameCnt += 30;		
	}
	iSeg++;
}
while(frameCnt < fileSize){
	refSegDisplay += "\033[1m\033[24m\033[0mX";
	frameCnt += 30;
}

refSegDisplay += "\033[0m\n\n";

if(verboseLevel == 3){
	for(unsigned long iseg=0;iseg<refSegServer.getSegCount();iseg++){
		Seg& segment=refSegServer.getSeg(iseg);
		refSegDisplay += refSegDisplay.valueOf(segment.begin()*frameLength)+" "+refSegDisplay.valueOf((segment.begin()+segment.length()-1)*frameLength)+" "+segment.string()+"\n";
	}
	refSegDisplay += "\n";
}
config.setParam("refSegDisplay", refSegDisplay);	

} // display


// look for a specific label in a list
long findLabel(XLine classToAnalyse,String labelToFind){

for(unsigned long i=0;i<classToAnalyse.getElementCount();i++){
	if(classToAnalyse.getElement(i)==labelToFind) return(i);
}

return(-1);//not found

}//findLabel

void moveSegmentFromOneClusterToAnother(LabelServer& labelServer, Seg *segment, SegCluster &currentCluster, SegCluster &newCluster){
	currentCluster.remove(*segment);
	segment->setLabelCode(newCluster.labelCode());	
	segment->setString((labelServer.getLabel(newCluster.labelCode())).getString());	
	
	newCluster.add(*segment);
}

/**********************************************************
* longerSegment: search the longer  of a cluster (in terms of time) and return it
* 
* Author C. Fredouille February 2006
***********************************************************/
Seg *longerSegment(Config& config, SegCluster& cluster){

DoubleVector time;                                        // The LLR value and the startindex of each of the potential trial
cluster.rewind();
Seg *segment;	
while((segment=cluster.getSeg())!=NULL){
		time.addValue((double)segment->length());
}//while segment

unsigned long ind = time.getIndexOfLargestValue();

cluster.rewind();
for (unsigned long i=0; i<=ind; i++){
	segment=cluster.getSeg();
	
}
	
cluster.rewind();
return segment;
}

/**********************************************************
* findClusterIndex: search the index of a cluster in a segment server
* 
* Author C. Fredouille February 2006
***********************************************************/
long findClusterIndex(String name,SegServer& segTmp){
  for(unsigned long icluster=0;icluster<segTmp.getClusterCount();icluster++){
    SegCluster& cluster=segTmp.getCluster(icluster);	
    if(cluster.string()==name)
      return icluster;
    
  }//for i
  return -1;
}//findClusterIndex


#endif


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

#if !defined(ALIZE_NormFeat_cpp)
#define ALIZE_NormFeat_cpp

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <liatools.h>
#include "NormFeat.h"
#include "AccumulateJFAStat.h"

using namespace alize;
using namespace std;

// Display Stats
void displayStats(FrameAccGD &frameAccu,FeatureServer & fs,ostream & out,Config &config) {  
  const DoubleVector & featureMean = frameAccu.getMeanVect();                    // Get the mean vector
  const DoubleVector & featureStd = frameAccu.getStdVect();                      // Get the std vector (sqrt(cov))
  unsigned long vectsize=fs.getVectSize();                                       // Get the vect size (number of coeff)
	for (unsigned int i = 0; i < vectsize; i++) 
		{out<<featureMean[i]<<" "<<featureStd[i]<< endl;}
}

void loadStatsFromFile(String filename,DoubleVector & extMean, DoubleVector & extCov) {
	Config tmp; // why is that ???? N.S.
	XList extStats(filename,tmp);
	for(unsigned long i=0;i<extStats.getLineCount();i++) {
		if(extStats.getLine(i).getElementCount()!=2) cerr << "W: Stat file does not match standard format (mean cov)" << endl;
		extMean[i]=extStats.getLine(i).getElement(0).toDouble();
		extCov[i]=extStats.getLine(i).getElement(1).toDouble();
	}
}
// Window mode stuff
bool emptyWindow(SegCluster &window){
  window.rewind();
  return (window.getSeg()==NULL);
}
// return the first seg of a cluster
Seg * firstSeg(SegCluster &window){
  window.rewind();
  Seg *seg= window.getSeg();
  return seg;
}
// return the last seg of a cluster
Seg * lastSeg(SegCluster &window){
  window.rewind();
  Seg *seg=window.getSeg();
  if (seg==NULL) throw Exception("empty window" , __FILE__, __LINE__);
  Seg *retSeg=seg;
  while((seg=window.getSeg())!=NULL) retSeg=seg;
  return retSeg;
}
// Compute the time interval (in s )between the first and last frames of the window
double windowLength(SegCluster &window,double frameLength){
  double tBegin=frameIdxToTime(firstSeg(window)->begin(),frameLength);
  Seg *lastS=lastSeg(window);
  double tEnd=frameIdxToTime(lastS->begin()+lastS->length(),frameLength);
  return (tEnd-tBegin);
}

// Used in order to fill the window
// Add segments and return true if the length is not reached
// Do nothing and return false if the length is reached
bool inWindow(Seg &seg,SegCluster &window,double windowDuration, double frameLength){     
  if (emptyWindow(window)) { window.addCopy(seg);return true;}
  if (windowLength(window,frameLength)<=windowDuration){
    window.addCopy(seg);
    return true;
  }
  return false;                             // the window is full . Do nothing
}

// suppress the first segment of the window
void suppressFirst(SegCluster &window){
  window.rewind();
  Seg *seg=window.getSeg();
  if (seg!=NULL) window.remove(*seg);
}
// remove the first segment, add the current segment 
void moveWindow(SegCluster &window,Seg &seg,Seg *&currentSeg,double duration, double frameLength){
  Seg *fSeg;
  window.addCopy(seg);
  if (debug) cout << "suppress seg[";
  do {
    if (debug) cout <<firstSeg(window)->begin()<<","<<firstSeg(window)->length()<<"][";
    suppressFirst(window);fSeg=firstSeg(window);
  } while ((currentSeg!=NULL)&&(fSeg->begin()<currentSeg->begin()) &&(windowLength(window,frameLength)>duration));
  if (debug) cout <<"] from the window and add seg["<<seg.begin()<<","<<seg.length()<<"]"<<endl;
  if (currentSeg==NULL) currentSeg=&seg;        // We deal with all the seg in the window before the move, next one is the added one
}
//compute the mean and cov for the current window
void computeWindowParam(Config &config,const DoubleVector &gMean,const DoubleVector &gCov,FeatureServer & fs,
			SegCluster &window,double windowDuration,double windowComp,double frameLength,
			DoubleVector &mean,DoubleVector &cov){
  FrameAccGD frameAccu;                                                      // Defines a frame accumulator for mean and cov computation 
  frameAccu.reset();
  accumulateStatFrame(frameAccu,fs, window, config);
  mean = frameAccu.getMeanVect();     // Get the mean vector
  cov = frameAccu.getStdVect();       // Get the std vector
  if (debug) {cout <<"Compute param on"<<endl;showCluster(window);}
  double duration=windowLength(window,frameLength);
  if (verboseLevel>1) cout <<"Window  real time["<<duration<<"]";
  //if (duration>windowDuration) duration=windowDuration; ????????????????????????????????
  double alpha=frameIdxToTime(totalFrame(window),frameLength);
  if (verboseLevel>1) cout <<" data time["<<alpha <<"]";
  alpha/=duration;
  if (alpha>1) alpha =1.0;
  if (verboseLevel>1) cout <<" Alpha["<<alpha<<"]";
  alpha*=windowComp;
  if (verboseLevel>1) cout <<" after Compensation["<<alpha<<"]"<<endl;;
  for (unsigned long c=0;c<fs.getVectSize();c++){
    mean[c]=(alpha*mean[c])+((1-alpha)*gMean[c]);
    cov[c]=(alpha*cov[c]*cov[c]) +((1-alpha)*gCov[c]*gCov[c])
      +(((1-alpha)*alpha)*pow(mean[c]-gMean[c],2));
    cov[c]=sqrt(cov[c]);
  }
     
}
bool inMiddle(Seg * &currentSeg,SegCluster &window,double windowDuration,double frameLength){
  double beginW=frameIdxToTime(firstSeg(window)->begin(),frameLength); 
  return (frameIdxToTime(currentSeg->begin(),frameLength)<=beginW+(windowDuration/2.0)); 
}
bool nextSeg(Seg * &currentSeg,SegCluster &window){
  window.rewind();
  Seg *seg=NULL;
  while(((seg=window.getSeg())!=NULL) &&(seg->begin()!=currentSeg->begin())); // TODO ADD OPERATOR != == between seg
  if (seg!=NULL) currentSeg=window.getSeg();
  else throw Exception("currentSeg is not in the window" , __FILE__, __LINE__);
  if (currentSeg==NULL) return false;
  else return true;
}

// get the glocal ubm offset
void getUbmOffset(RealVector <double>& ubm_offset,Matrix<double>& U,MixtureGD & world,Config & config) {
	if (verbose) cout << "NAP mode, get global UBM offset" << endl;	
	unsigned long svSize=world.getDistribCount()*world.getVectSize();
	ubm_offset.setSize(svSize);
	RealVector <double> ubm_sv(svSize,svSize);		
	if (verbose) cout << "Size of SVs: [" << ubm_sv.size() << "]...";
	modelToSv(world,ubm_sv);
	projectOnSubSpace(ubm_sv,U,ubm_offset);
}

// remove ubm offset on a frame
void removeOffsetOnFrame(RealVector <double> & occ,Feature & f, RealVector <double> & ubm_offset,unsigned long nbDistrib) {
	for (unsigned long i=0;i<f.getVectSize();i++) {
		if (verboseLevel > 4) cout << "dim " << i << ":" ;		
		double _offset=0.0;// get offset for this cepstrum coeff
		for (unsigned long j=0;j<nbDistrib;j++) { // get corresponding values (TODO: on top ten)	
			if (verboseLevel > 4) cout << j << "["<<j*f.getVectSize()+i<<"] ";		
			_offset+=ubm_offset[j*f.getVectSize()+i]*occ[j];		// to check ... propal from Anthony
		}
		f[i]-=_offset;
	}
}

// Compute Occ and Remove offset
void featureChannelCompNAP(RealVector <double> &ubm_offset,MixtureGD & world,FeatureServer &fs,StatServer & ss,SegCluster & selectedSegments,Config &config) {
	MixtureGDStat &occAcc=ss.createAndStoreMixtureGDStat(world);
	Seg* seg;                                                                         // reset the reader at the begin of the input stream
	selectedSegments.rewind();      
	while((seg=selectedSegments.getSeg())!=NULL){                                     // For each of the selected segments
	  unsigned long idxBeginFrame=seg->begin()+fs.getFirstFeatureIndexOfASource(seg->sourceName()); 
	  fs.seekFeature(idxBeginFrame); 
	  Feature f;
	  for (unsigned long idxFrame=0;idxFrame<seg->length();idxFrame++){               // For each frame of the segment
		fs.readFeature(f);
		occAcc.computeAndAccumulateOcc(f);
		RealVector <double> occ=occAcc.getOccVect();
		removeOffsetOnFrame(occ,f,ubm_offset,world.getDistribCount());
	  }
	}
	if (verbose) cout << "done" << endl;
}
	
int normFeat (Config & config) {
	
  String inputFeatureFileName =config.getParam("inputFeatureFilename");          // input feature - could be a simple feature file or a list of filename
  String labelSelectedFrames  =config.getParam("labelSelectedFrames");           // Only the frames from segments with this label  will be used
  bool segmentalMode=false;                                                      // the norm is done segment by segment
  bool fileMode=false;                                                           // the norm is done on the entire file (default)
  bool windowMode=false;                                                         // the norm is done time window by time window
  if (config.existsParam("segmentalMode"))segmentalMode=config.getParam("segmentalMode").toBool();
  if (config.existsParam("fileMode")) fileMode=config.getParam("fileMode").toBool(); // DEFAULT MODE
  if (config.existsParam("windowMode")) windowMode=config.getParam("windowMode").toBool();
  if (segmentalMode && windowMode) // Only one of this modes should be active
    throw Exception (" Incompatible modes (segmental/windowMode)" , __FILE__, __LINE__);
  if ((!fileMode)&&(!segmentalMode) && (!windowMode)) fileMode=true;
  double windowDuration=5; // Default length for window=5s
  double windowComp=1;     // Default compensation for window vs global means and cov (=1, no compensation)
  if (windowMode){
    if (config.existsParam("windowDuration")) windowDuration=config.getParam("windowDuration").toDouble();
    if (config.existsParam("windowComp")) windowComp=config.getParam("windowComp").toDouble();
  }
  bool writeAllFeature=true; // Output a vector for all input vectors (selected and not selected vectors) - DEFAULT=on
  if (config.existsParam("writeAllFeatures")) writeAllFeature=config.getParam("writeAllFeatures").toBool();    // Define if all the feature (selected or not) should be written 
  bool cms=false;
  bool extStat=false;
  bool var=false;
  bool warp=false;
  if (config.existsParam("cmsOnly")) cms=config.getParam("cmsOnly").toBool();
  if (config.existsParam("varOnly")) cms=config.getParam("varOnly").toBool();
  if (config.existsParam("warp"))    warp=config.getParam("warp").toBool();

  if (config.existsParam("externalStatsFilename")) extStat=true;
  if (cms && var) throw Exception (" cmsOnly and varOnly are not compatible)" , __FILE__, __LINE__);
  if (cms && warp)throw Exception (" cmsOnly and warp are not compatible)" , __FILE__, __LINE__);
  if (var && warp)throw Exception (" varOnly and warp are not compatible)" , __FILE__, __LINE__);
  if (extStat && warp)throw Exception (" externalStats and warp are not compatible)" , __FILE__, __LINE__);

  if (verbose){
    cout << "NormFeat";
    if (!segmentalMode && !windowMode && !fileMode) cout << "NO NORM IS SELECTED"<<endl;
    if (segmentalMode) cout << " Segmental Mode - Normalisation  by segment"<<endl;
    if (windowMode) cout << "Window mode, Window Duration["<<windowDuration<<"], Window compensation ["<<windowComp<<"]"<<endl;
    if (fileMode) cout << "File Mode- Normalisation by file by file" <<endl;
    if (fileMode && (windowMode || segmentalMode)) cout <<" a global 0,1 norm is applied after segmental/window mode"<<endl;
    if (extStat) cout << "Loading statistics from external file: ["<<config.getParam("externalStatsFilename")<<"]"<<endl;
    if (cms) cout <<"CMS: mean only mode"<<endl;
    if (var) cout <<"VAR: variance only mode"<<endl;
    if (warp) cout <<"Warping mode"<<endl;
    if (!(cms) && (!var) && (!warp)) cout << "var and mean normalisation"<< endl; 

  }
  double frameLength=0.01; // length of a frame ins, by default 10 ms
  if (config.existsParam("frameLength")) frameLength=config.getParam("frameLength").toDouble();
  // Warping mode global variables and initialisation
  Histo gaussHisto; // Targeted distribution histogram
  Histo *rawHistoT;  // Tab of histo. Used to compute the raw histo
		   // coeff by coeff
  unsigned long nbBin=40;
  bool warp01=false;
  if (config.existsParam("warpBinCount")) nbBin=config.getParam("warpBinCount").toLong();
  if (warp){
      if (config.existsParam("warpGaussHistoFilename")) gaussHisto.load(config.getParam("warpGaussHistoFilename"));
      else {
	  unsigned long gaussHistoSampleCount=50000;
	  unsigned long gaussHistoBinCount=100;
	  gaussHisto=makeGausHisto(gaussHistoSampleCount,0, 1,gaussHistoBinCount);      // Make the targeted distribution	   
      }
      initHistoTab(rawHistoT,1,1);  // Make an empty histo by coeff
				    // TODO(suppress)
      if (config.existsParam("warpAdd01Norm")) warp01=config.getParam("warpAdd01Norm").toBool();
}
 
  // Begin !!
  XLine inputFeatureFileNameList;                                                // The (feature) input filename list
  if ( inputFeatureFileName.endsWith(".lst")){                                   // If the file parameter is the name of a XList file
	XList inputFileNameXList(inputFeatureFileName,config);                   // Read the filename list file
	inputFeatureFileNameList=inputFileNameXList.getAllElements();            // And put the filename in a list if the file is a list of feature filenames
  }
  else {                                                                         // It was a simple feature file and not a filename list
    inputFeatureFileNameList.addElement(inputFeatureFileName);                   // add the filename in the list
  }
  try{
	String *file;
	while ((file=inputFeatureFileNameList.getElement())!= NULL){                         // Loop on each feature file
		String & featureFileName=(*file);                                            // Current file basename
		FeatureServer fs(config,featureFileName);
		DoubleVector extMean(fs.getVectSize(),fs.getVectSize()); // Init double vectors, does it have to be so harmful!!!
		DoubleVector extCov(fs.getVectSize(),fs.getVectSize());
		if (extStat) {
		  loadStatsFromFile(config.getParam("externalStatsFilename"),extMean, extCov);
		}
		// Begin the real work, file by file
		SegServer segmentsServer;                                              // Create the segment server for managing the segments/clusters
		LabelServer labelServer;                                               // Create the label server, for indexing the segments/clusters
		initializeClusters(featureFileName,segmentsServer,labelServer,config); // Reading the segmentation files for each feature input file
		verifyClusterFile(segmentsServer,fs,config);                           // Verify if the segments ending before the end of the feature files...
		unsigned long codeSelectedFrame=labelServer.getLabelIndexByString(labelSelectedFrames);// Get the index of the cluster with in interest audio segments
		if (verbose) cout << "NormFeat file normalisation["<<featureFileName<<"]"<<endl;
		SegCluster& selectedSegments=segmentsServer.getCluster(codeSelectedFrame); // The general cluster of the selected/used segments   
		SegServer segServerWindow;                                                 // specific segment server for window mode
		SegCluster& window=segServerWindow.createCluster(0,labelSelectedFrames,selectedSegments.sourceName());  // the window
		//
		if (debug || (verboseLevel>2)){
		  DoubleVector globalMean;
		  DoubleVector globalCov;
		  FrameAccGD globalFrameAcc;
		  unsigned long nbFrame=computeMeanCov(config,fs,selectedSegments,globalMean,globalCov); // Compute the global mean and covariance
		  cout <<"global mean and cov before norm, number of frame= ["<<nbFrame<<"]"<<endl;
		  for (unsigned i=0; i < fs.getVectSize(); i++)cout << "mean[" << i << "=" << globalMean[i] << "]\tcov[" << globalCov[i] << "]" << endl;
		}
		// SEGMENTAL MODE
		if (segmentalMode) {                                                       // Segmental mode (the norm is done seg by seg)
		  Seg *seg;
		  selectedSegments.rewind(); 
		  while((seg=selectedSegments.getSeg())!=NULL) {                     // For each of the selected segments
		    if (debug) cout << "current Seg["<<seg->begin()<<","<<seg->length()<<"]"<<endl;
		    if (!warp){
		      DoubleVector mean,cov;
		      if(extStat) {                                                    // load external stats from file (NS)
			mean=extMean;
			cov=extCov;
		      }
		      else {
			FrameAccGD frameAccu;                                                      // Defines a frame accumulator for mean and cov computation
			frameAccu.reset();                                                         // reset the frame accumulator
			accumulateStatFrame(frameAccu,fs, seg, config);
			mean = frameAccu.getMeanVect();     // Get the mean vector
			cov = frameAccu.getStdVect();       // Get the std vector (sqrt(cov)) (remove const and ref to modify)
		      }
		      if (cms) cov.setAllValues(1.0);	// remove cov effect if cms only (N.S) 
		      if (var) mean.setAllValues(0.0);      // remove effects if var only (jfb)
		      computeZeroOne(mean,cov,fs,seg,config);
		    } // end of mean and/or cov normalisation
		    if (warp){ // Gaussian warping mode
			freezeHistoTab(rawHistoT); //TODO replace by
			initHistoTab(rawHistoT,fs.getVectSize(),nbBin);//Reset
			accumulateHistoFrame(rawHistoT,fs,seg,config);
			computeHistoTab(rawHistoT,fs.getVectSize());
			computeWarp(rawHistoT,gaussHisto,fs,seg,config);
		    } // end of warping mode
		  } // end of a given segment
		} // end of classical segmental mode
		// WINDOW MODE
		if (windowMode){ // JFB March 20 2006
		  // compute the global mean and cov
		  FrameAccGD frameAccu;                                                      // Defines a frame accumulator for mean and cov computation
		  frameAccu.reset();    
		  accumulateStatFrame(frameAccu, fs, selectedSegments, config); 
		  DoubleVector gMean,gCov;             // global (on all semected segments) mean and cov vectors
		  gMean= frameAccu.getMeanVect();      // Get the mean vector
		  gCov = frameAccu.getStdVect();       // Get the std vector (sqrt(cov)) (remove const and ref to modify)
		  DoubleVector mean,cov;               // mean and cov estimations for a given window
		  Seg *currentSeg=NULL;                // the central segment, where the transformation should be applied
		  Seg *seg;
		  selectedSegments.rewind(); 
		  while((seg=selectedSegments.getSeg())!=NULL) {      // For each of the selected segments
		    if (debug){
		      cout <<"selected Seg["<<seg->begin()<<","<<seg->length()<<"]"<<endl;
		      cout <<"current seg";
		      if (currentSeg==NULL) cout <<" NULL"<<endl;
		      else cout <<"["<<currentSeg->begin()<<","<<currentSeg->length()<<"]"<<endl;
		      cout <<"window"<<endl;
		      showCluster(window);
		    }
		    if (inWindow(*seg,window,windowDuration,frameLength)){// add segments until the window length is reached
		      if (currentSeg==NULL) currentSeg=seg;
		    }
		    else{                                                 // Normal stuff: the window is full, compute the norm, suppress the first segment
		      bool empty=emptyWindow(window);
		      if (empty) {window.addCopy(*seg);currentSeg=seg;}   // special case where one seg is the window (is larger than the time constraint)
 		      if (!warp){                                         // compute the mean and cov for the window (mean and cov normalisation modes)
			  computeWindowParam(config,gMean,gCov,fs,window,windowDuration,windowComp,frameLength,mean,cov); 
			  if (cms) cov.setAllValues(1.0);	              // remove cov effect if cms only (N.S) 
			  if (var) mean.setAllValues(0.0);               // remove effects if var only (jfb)	  
		      }
		      else {                                             // warping mode, compute the histos
			  freezeHistoTab(rawHistoT); //TODO replace by
			  initHistoTab(rawHistoT,fs.getVectSize(),nbBin);//Reset
			  accumulateHistoFrame(rawHistoT,fs,window,config);
			  computeHistoTab(rawHistoT,fs.getVectSize());
			
		      }
	  
		      bool ok=true;
		      while ((ok) && (inMiddle(currentSeg,window,windowDuration,frameLength))){             // until we reach the middle of the window 
			if (debug) cout << "Apply norm on Seg ["<<currentSeg->begin()<<","<<currentSeg->length()<<"]"<<endl;
			if (!warp) // mean and cov normalisation modes
			    computeZeroOne(mean,cov,fs,currentSeg,config);// compute the norm for the current segment
			else       // warping mode
			    computeWarp(rawHistoT,gaussHisto,fs,currentSeg,config);
			ok=nextSeg(currentSeg,window);
		      }          
		      if (empty) suppressFirst(window);
		      else moveWindow(window,*seg,currentSeg,windowDuration,frameLength);        // remove the first segments, add the current segment 
		    }
		  } 
		  // Deal with the last window
		  if (debug) cout <<"Deal with last window"<<endl;
		  if (!emptyWindow(window)){
		      if (!warp){ // mean and norm normalisation modes
			  computeWindowParam(config,gMean,gCov,fs,window,windowDuration,windowComp,frameLength,mean,cov); // compute the mean and cov for the current window
			  if (cms) cov.setAllValues(1.0);	          // remove cov effect if cms only (N.S) 
			  if (var) mean.setAllValues(0.0);           // remove effects if var only (jfb)
		      }
		      else { // warping mode
			  initHistoTab(rawHistoT,fs.getVectSize(),nbBin);//Reset
			  accumulateHistoFrame(rawHistoT,fs,window,config);
			  computeHistoTab(rawHistoT,fs.getVectSize());
		      }
		      
		      do{
			  if (debug) cout << "Apply norm on Seg ["<<currentSeg->begin()<<","<<currentSeg->length()<<"]"<<endl;
			  if (!warp) computeZeroOne(mean,cov,fs,currentSeg,config);    // compute the norm for the current segment
			  else computeWarp(rawHistoT,gaussHisto,fs,currentSeg,config); // warping mode
		      }
		      while (nextSeg(currentSeg,window));        // until the end;
		      
		  }
		}
		if (fileMode){                                   // cluster mode (the norm is done on the complete cluster
		    if (!warp){ // mean and cov normalisation modes
			DoubleVector mean,cov;
			if(extStat) {// load external stats from file (NS)
			    mean=extMean;
			    cov=extCov;}
			else {
			    FrameAccGD frameAccu;                                                      // Defines a frame accumulator for mean and cov computation
			    frameAccu.reset();   
			    accumulateStatFrame(frameAccu,fs, selectedSegments, config); 
			    mean = frameAccu.getMeanVect();     // Get the mean vector
			    cov  = frameAccu.getStdVect();      // Get the std vector (sqrt(cov)) (remove const and ref to modify)
			}
			if (cms) cov.setAllValues(1.0);	// remove cov effect if cms only (N.S) 
			if (var) mean.setAllValues(0.0);      // remove effects if var only (jfb)
			computeZeroOne(mean,cov,fs, selectedSegments, config);
		    }
		    else{    // Warping mode
			freezeHistoTab(rawHistoT); //TODO replace by
			initHistoTab(rawHistoT,fs.getVectSize(),nbBin);//Reset
			accumulateHistoFrame(rawHistoT,fs,selectedSegments,config);
			computeHistoTab(rawHistoT,fs.getVectSize());
			rawHistoT[0].saveGnuplot("toto");
			computeWarp(rawHistoT,gaussHisto,fs,selectedSegments,config);
		    } 
		}                                   
		// Warp optionnal additional normalisation
		if (warp && warp01){
		    DoubleVector mean,cov;
		    FrameAccGD frameAccu;                                                      // Defines a frame accumulator for mean and cov computation
		    frameAccu.reset();   
		    accumulateStatFrame(frameAccu,fs, selectedSegments, config); 
		    mean = frameAccu.getMeanVect();     // Get the mean vector
		    cov  = frameAccu.getStdVect();      
		    computeZeroOne(mean,cov,fs, selectedSegments, config);  
		}
		// Output the normalized features - Take care: only the saveFeatureFileExtension parameter make a difference between the input and output files
		if (!(config.existsParam("featureFlags")))
		    config.setParam("featureFlags",fs.getFeatureFlags().getString());// Put the file flag in the config (could be different for each file   
		cout << "Writing to: " << featureFileName << endl;
		FeatureFileWriter w(featureFileName, config);                             // build a featurefile writer to output the features (real features)
		if (writeAllFeature) {                                                    // Output all the features- feature count id the same
		    SegServer fakeSegServer;                                          // Create a new fake segment server
		    fakeSegServer.createCluster(0);                                   // Create a new cluster
		    SegCluster& fakeSeg=fakeSegServer.getCluster(0);                  // Get the cluster               
		    fakeSeg.add(fakeSegServer.createSeg(0,fs.getFeatureCount(),codeSelectedFrame,
							labelSelectedFrames,featureFileName));                            // Add a segment with all the features
		    outputFeatureFile(config,fs,fakeSeg,w);       	          // output all the features - giving the same file length
		}
		else
		    outputFeatureFile(config,fs,selectedSegments, w);    // Output only the selected features - giving a shorter output 
		if (debug || (verboseLevel>2)){
		    DoubleVector globalMean;
		    DoubleVector globalCov;
		    FrameAccGD globalFrameAcc;
		    unsigned long nbFrame=computeMeanCov(config,fs,selectedSegments,globalMean,globalCov); // Compute the global mean and covariance
		    cout <<"global mean and cov after norm, number of frame= ["<<nbFrame<<"]"<<endl;
		    for (unsigned i=0; i < fs.getVectSize(); i++)cout << "mean[" << i << "=" << globalMean[i] << "]\tcov[" << globalCov[i] << "]" << endl;
		}
		// Freeze
		if (warp) freezeHistoTab(rawHistoT);
	}// End feature file loop

  }// end try
  
  catch (Exception & e)
      {
	  cout << e.toString ().c_str () << endl;
      }
  return 0;
}

int infoFeat (Config & config) {
	
  String inputFeatureFileName =config.getParam("inputFeatureFilename");          // input feature - could be a simple feature file or a list of filename
  String labelSelectedFrames  =config.getParam("labelSelectedFrames");           // Only the frames from segments with this label  will be used
  bool segmentalMode=false;
  if (config.existsParam("segmentalMode"))segmentalMode=config.getParam("segmentalMode").toBool(); // selected mode for segmental computation (1 norm by file or by segment)
  if (verbose){
    cout << "NormFeat -- Display Stats on Feature";
    if (segmentalMode) cout << " Segmental Mode - Normalisation  by segment"<<endl;
    else cout << "File Mode - Normalisation by file by file" <<endl;
  }
  XLine inputFeatureFileNameList;                                                // The (feature) input filename list
  if ( inputFeatureFileName.endsWith(".lst")){                                   // If the file parameter is the name of a XList file
	XList inputFileNameXList(inputFeatureFileName,config);                       // Read the filename list file
	inputFeatureFileNameList=inputFileNameXList.getAllElements();                       // And put the filename in a list if the file is a list of feature filenames
  }
  else {                                               // It was a simple feature file and not a filename list
    inputFeatureFileNameList.addElement(inputFeatureFileName);                   // add the filename in the list
  }
  try{
    MixtureServer ms(config); 
 
    String *file;
    while ((file=inputFeatureFileNameList.getElement())!= NULL){                   // Loop on each feature file
      String & featureFileName=(*file);                                            // Current file basename
      FeatureServer fs(config,featureFileName);                                    // Reading the feature file
      SegServer segmentsServer;                                                             // Create the segment server for managing the segments/clusters
      LabelServer labelServer;                                                              // Create the lable server, for indexing the segments/clusters
      initializeClusters(featureFileName,segmentsServer,labelServer,config);               // Reading the segmentation files for each feature input file
      verifyClusterFile(segmentsServer,fs,config);                                          // Verify if the segments ending before the end of the feature files...
      unsigned long codeSelectedFrame=labelServer.getLabelIndexByString(labelSelectedFrames);            // Get the index of the cluster with in interest audio segments
      if (verbose)cout << "Display Stats for file ["<<featureFileName<<"]"<< endl;
      SegCluster& selectedSegments=segmentsServer.getCluster(codeSelectedFrame);        // Gives the cluster of the selected/used segments   
      FrameAccGD frameAccu;                                                             // Defines a frame accumulator for mean and cov computation
      frameAccu.reset();                                                                // reset the frame accumulator
      Seg *seg;                                                                         // current selectd segment
      selectedSegments.rewind();                                                        // reset the reader at the begin of the input stream
      while((seg=selectedSegments.getSeg())!=NULL){                                     // For each of the selected segments
	accumulateStatFrame(frameAccu,fs,seg,config);                   // Accumulate the stat for the segment
	if (segmentalMode){                                                             // Normalize segment by segment mode is selected
	cout << "Warning: Display Stats only for seg: " << seg->begin() << " " << seg->begin()+seg->length() << endl;
	displayStats(frameAccu,fs,cout,config);                      // Display Stats
	frameAccu.reset();                                                            // Reset the accumulator for the next segment      
	}
      }
      if (!segmentalMode){                                                              // The mean/cov are computed on all the segments
	selectedSegments.rewind();                                                      // reset the reader at the begin of the input stream
	globalMeanCov (fs,selectedSegments,frameAccu,config);
	String filename=fs.getNameOfASource(0);
	ofstream out((filename+".stat").c_str());
	cout << "Warning: Saving Stats only in " << filename << ".stat" << endl;
	displayStats(frameAccu,fs,out,config);                      // Display Stats
	out.close();
	}
	}
 }// end try
  catch (Exception & e)
      {
	cout << e.toString ().c_str () << endl;
      }
  return 0;
}

int featMap (Config & config) {

String inputFeatureFileName =config.getParam("inputFeatureFilename");          // input feature - could be a simple feature file or a list of filename
String labelSelectedFrames  =config.getParam("labelSelectedFrames");           // Only the frames from segments with this label  will be used
bool writeAllFeature=(config.getParam("writeAllFeatures")=="true");            // Define if all the feature (selected or not) should be written
  if (config.getParam_debug()) debug=true;  else debug=false;
  if (config.existsParam("verbose")) verbose=true; else verbose=false; 
if (verbose){
	cout << "NormFeat: FeatureMapping Mode - Normalisation by file by file" <<endl;
}
XLine inputFeatureFileNameList;										// The (feature) input filename list
if ( inputFeatureFileName.endsWith(".lst")) {                                                    // If the file parameter is the name of a XList file
	XList inputFileNameXList(inputFeatureFileName,config);				// Read the filename list file
	inputFeatureFileNameList=inputFileNameXList.getAllElements();                       // And put the filename in a list if the file is a list of feature filenames
}
else {                                               // It was a simple feature file and not a filename list
	inputFeatureFileNameList.addElement(inputFeatureFileName);                   // add the filename in the list
}
try{
	MixtureServer ms(config); 
	ms.loadMixtureGD(config.getParam("inputRootModelFilename"));						// Load the model to map to
	ms.loadMixtureGD(config.getParam("inputSubFilename"));
	
	String *file;
	while ((file=inputFeatureFileNameList.getElement())!= NULL){                   // Loop on each feature file
		String & featureFileName=(*file);                                            // Current file basename
		FeatureServer fs(config,featureFileName);                                    // Reading the feature file
		SegServer segmentsServer;                                                             // Create the segment server for managing the segments/clusters
		LabelServer labelServer;                                                              // Create the lable server, for indexing the segments/clusters
		initializeClusters(featureFileName,segmentsServer,labelServer,config);               // Reading the segmentation files for each feature input file
		verifyClusterFile(segmentsServer,fs,config);                                          // Verify if the segments ending before the end of the feature files...
		unsigned long codeSelectedFrame=labelServer.getLabelIndexByString(labelSelectedFrames);            // Get the index of the cluster with in interest audio segments
		if (verbose) cout << "Feature Mapping file normalisation["<<featureFileName<<"]"<< endl;
		SegCluster& selectedSegments=segmentsServer.getCluster(codeSelectedFrame);        // Gives the cluster of the selected/used segments   
		FrameAccGD frameAccu;                                                             // Defines a frame accumulator for mean and cov computation
		frameAccu.reset();                                                                // reset the frame accumulator
		featureMapping(ms,fs,selectedSegments,config); // perform mapping on the whole segment
		
		// Output the normalized features - Take care: only the saveFeatureFileExtension parameter make a difference between the input and output files
		if (!(config.existsParam("featureFlags")))
			config.setParam("featureFlags",fs.getFeatureFlags().getString());                 // Put the file flag in the config (could be different for each file   
			//config.setParam("featureServerMask",config.getParam("outputFeatureServerMask")); // To be replaced 
		FeatureFileWriter w(featureFileName, config);                                        // build a featurefile writer to output the features
		
		if (writeAllFeature) {                                                            // Output all the features- feature count id the same
			SegServer fakeSegServer;                                                        // Create a new fake segment server
			fakeSegServer.createCluster(0);                                                 // Create a new cluster
			SegCluster& fakeSeg=fakeSegServer.getCluster(0);                                // Get the cluster               
			fakeSeg.add(fakeSegServer.createSeg(0,fs.getFeatureCount(),codeSelectedFrame,
							    labelSelectedFrames,featureFileName));       // Add a segment with all the features
			outputFeatureFile(config,fs,fakeSeg,w);       	                        // output all the features - giving the same file length
		}
		else
		outputFeatureFile(config,fs,selectedSegments, w);                       // Output only the selected features - giving a shorter output 
	}                                                                                   // End feature file loop
} // end try
catch (Exception & e)
{
	cout << e.toString ().c_str () << endl;
}
return 0;
}

void featureWarping(Matrix <double> table,FeatureServer & fs,SegCluster & selectedSegments,Config & config) {
  Seg *seg;                                                                         // current selectd segment
  selectedSegments.rewind();                                                      // reset the reader at the begin of the input stream
  while((seg=selectedSegments.getSeg())!=NULL){                                   // For each of the selected segments
	unsigned long begin=seg->begin()+fs.getFirstFeatureIndexOfASource(seg->sourceName()); // Idx of the first frame of the current file in the feature server
	fs.seekFeature(begin);
	Feature f;
	for (unsigned long idxFrame=0;idxFrame<seg->length();idxFrame++){                          // for all the features of the segment
	fs.readFeature(f,0); 		
	// do it
	fs.writeFeature(f);
	}
  }	
}

int featWarp (Config & config) {

String inputFeatureFileName =config.getParam("inputFeatureFilename");          // input feature - could be a simple feature file or a list of filename
String labelSelectedFrames  =config.getParam("labelSelectedFrames");           // Only the frames from segments with this label  will be used
bool writeAllFeature=(config.getParam("writeAllFeatures")=="true");            // Define if all the feature (selected or not) should be written
  if (config.getParam_debug()) debug=true;  else debug=false;
  if (config.existsParam("verbose")) verbose=true; else verbose=false; 
if (verbose){
	cout << "NormFeat: FeatureWarping Mode - Normalisation by file by file" <<endl;
}
XLine inputFeatureFileNameList;										// The (feature) input filename list
if ( inputFeatureFileName.endsWith(".lst")) {                                                    // If the file parameter is the name of a XList file
	XList inputFileNameXList(inputFeatureFileName,config);				// Read the filename list file
	inputFeatureFileNameList=inputFileNameXList.getAllElements();                       // And put the filename in a list if the file is a list of feature filenames
}
else {                                               // It was a simple feature file and not a filename list
	inputFeatureFileNameList.addElement(inputFeatureFileName);                   // add the filename in the list
}
try{
	MixtureServer ms(config); 
	ms.loadMixtureGD(config.getParam("inputRootModelFilename"));						// Load the model to map to
	ms.loadMixtureGD(config.getParam("inputSubFilename"));
	
	String *file;
	while ((file=inputFeatureFileNameList.getElement())!= NULL){                   // Loop on each feature file
		String & featureFileName=(*file);                                            // Current file basename
		FeatureServer fs(config,featureFileName);                                    // Reading the feature file
		SegServer segmentsServer;                                                             // Create the segment server for managing the segments/clusters
		LabelServer labelServer;                                                              // Create the lable server, for indexing the segments/clusters
		initializeClusters(featureFileName,segmentsServer,labelServer,config);               // Reading the segmentation files for each feature input file
		verifyClusterFile(segmentsServer,fs,config);                                          // Verify if the segments ending before the end of the feature files...
		unsigned long codeSelectedFrame=labelServer.getLabelIndexByString(labelSelectedFrames);            // Get the index of the cluster with in interest audio segments
		if (verbose) cout << "Feature Warping file normalisation["<<featureFileName<<"]"<< endl;
		SegCluster& selectedSegments=segmentsServer.getCluster(codeSelectedFrame);        // Gives the cluster of the selected/used segments   
		Matrix <double> table;
		table.load("table");
		featureWarping(table,fs,selectedSegments,config); // perform mapping on the whole segment
		
		// Output the normalized features - Take care: only the saveFeatureFileExtension parameter make a difference between the input and output files
		if (!(config.existsParam("featureFlags")))
			config.setParam("featureFlags",fs.getFeatureFlags().getString());                 // Put the file flag in the config (could be different for each file   
			//config.setParam("featureServerMask",config.getParam("outputFeatureServerMask")); // To be replaced 
		FeatureFileWriter w(featureFileName, config);                                        // build a featurefile writer to output the features
		
		if (writeAllFeature) {                                                            // Output all the features- feature count id the same
			SegServer fakeSegServer;                                                        // Create a new fake segment server
			fakeSegServer.createCluster(0);                                                 // Create a new cluster
			SegCluster& fakeSeg=fakeSegServer.getCluster(0);                                // Get the cluster               
			fakeSeg.add(fakeSegServer.createSeg(0,fs.getFeatureCount(),codeSelectedFrame,
							    labelSelectedFrames,featureFileName));       // Add a segment with all the features
			outputFeatureFile(config,fs,fakeSeg,w);       	                        // output all the features - giving the same file length
		}
		else
		outputFeatureFile(config,fs,selectedSegments, w);                       // Output only the selected features - giving a shorter output 
	}                                                                                   // End feature file loop
} // end try
catch (Exception & e)
{
	cout << e.toString ().c_str () << endl;
}
return 0;
}

int normFeatNAP (Config & config) {

String inputFeatureFileName =config.getParam("inputFeatureFilename");          // input feature - could be a simple feature file or a list of filename
String labelSelectedFrames  =config.getParam("labelSelectedFrames");           // Only the frames from segments with this label  will be used
bool writeAllFeature=(config.getParam("writeAllFeatures")=="true");            // Define if all the feature (selected or not) should be written
  if (config.getParam_debug()) debug=true;  else debug=false;
  if (config.existsParam("verbose")) verbose=true; else verbose=false; 
if (verbose){
	cout << "NormFeat: NAP Mode - Normalisation by file by file" <<endl;
}
XLine inputFeatureFileNameList;										// The (feature) input filename list
if ( inputFeatureFileName.endsWith(".lst")) {                                                    // If the file parameter is the name of a XList file
	XList inputFileNameXList(inputFeatureFileName,config);				// Read the filename list file
	inputFeatureFileNameList=inputFileNameXList.getAllElements();                       // And put the filename in a list if the file is a list of feature filenames
}
else {                                               // It was a simple feature file and not a filename list
	inputFeatureFileNameList.addElement(inputFeatureFileName);                   // add the filename in the list
}
try{
	MixtureServer ms(config); 
	StatServer ss(config);
	MixtureGD & world=ms.loadMixtureGD(config.getParam("inputWorldFilename"));						// Load the world model
	// Compute UBM channel effect
	RealVector <double> ubm_offset; // for nap
	Matrix <double> channelMatrix;
	channelMatrix.load(config.getParam("initChannelMatrix"),config);
	getUbmOffset(ubm_offset,channelMatrix,world,config);
	

	String *file;
	while ((file=inputFeatureFileNameList.getElement())!= NULL){                   // Loop on each feature file
		String & featureFileName=(*file);                                            // Current file basename
		FeatureServer fs(config,featureFileName);                                    // Reading the feature file
		SegServer segmentsServer;                                                             // Create the segment server for managing the segments/clusters
		LabelServer labelServer;                                                              // Create the lable server, for indexing the segments/clusters
		initializeClusters(featureFileName,segmentsServer,labelServer,config);               // Reading the segmentation files for each feature input file
		verifyClusterFile(segmentsServer,fs,config);                                          // Verify if the segments ending before the end of the feature files...
		unsigned long codeSelectedFrame=labelServer.getLabelIndexByString(labelSelectedFrames);            // Get the index of the cluster with in interest audio segments
		if (verbose) cout << "Eigenchannel file normalisation["<<featureFileName<<"]"<< endl;
		SegCluster& selectedSegments=segmentsServer.getCluster(codeSelectedFrame);        // Gives the cluster of the selected/used segments   
		FrameAccGD frameAccu;                                                             // Defines a frame accumulator for mean and cov computation
		frameAccu.reset();                                                                // reset the frame accumulator
		featureChannelCompNAP(ubm_offset,world,fs,ss,selectedSegments,config); // perform compensation on the whole segment
		
		// Output the normalized features - Take care: only the saveFeatureFileExtension parameter make a difference between the input and output files
		if (!(config.existsParam("featureFlags")))
			config.setParam("featureFlags",fs.getFeatureFlags().getString());                 // Put the file flag in the config (could be different for each file   
			//config.setParam("featureServerMask",config.getParam("outputFeatureServerMask")); // To be replaced 
		FeatureFileWriter w(featureFileName, config);                                        // build a featurefile writer to output the features
		
		if (writeAllFeature) {                                                            // Output all the features- feature count id the same
			SegServer fakeSegServer;                                                        // Create a new fake segment server
			fakeSegServer.createCluster(0);                                                 // Create a new cluster
			SegCluster& fakeSeg=fakeSegServer.getCluster(0);                                // Get the cluster               
			fakeSeg.add(fakeSegServer.createSeg(0,fs.getFeatureCount(),codeSelectedFrame,
							    labelSelectedFrames,featureFileName));       // Add a segment with all the features
			outputFeatureFile(config,fs,fakeSeg,w);       	                        // output all the features - giving the same file length
		}
		else
		outputFeatureFile(config,fs,selectedSegments, w);                       // Output only the selected features - giving a shorter output 
	}                                                                                   // End feature file loop
} // end try
catch (Exception & e)
{
	cout << e.toString ().c_str () << endl;
}
return 0;
}

int normFeatFA (Config & config) {

String inputFeatureFileName =config.getParam("inputFeatureFilename");          // input feature - could be a simple feature file or a list of filename
String labelSelectedFrames  =config.getParam("labelSelectedFrames");           // Only the frames from segments with this label  will be used
bool writeAllFeature=(config.getParam("writeAllFeatures")=="true");            // Define if all the feature (selected or not) should be written
if (verbose)	cout << "(NormFeat) EigenChannel Mode - Normalisation by file by file" <<endl;

XLine inputFeatureFileNameList;										// The (feature) input filename list
if ( inputFeatureFileName.endsWith(".lst")) {                                                    // If the file parameter is the name of a XList file
	XList inputFileNameXList(inputFeatureFileName,config);				// Read the filename list file
	inputFeatureFileNameList=inputFileNameXList.getAllElements();                       // And put the filename in a list if the file is a list of feature filenames
}
else {                                               // It was a simple feature file and not a filename list
	inputFeatureFileNameList.addElement(inputFeatureFileName);                   // add the filename in the list
}
try{
	String *file;
	while ((file=inputFeatureFileNameList.getElement())!= NULL){                   // Loop on each feature file
		String & featureFilename=(*file);                                            // Current file basename
		if (verbose) cout << "(NormFeat) Eigenchannel file normalisation["<<featureFilename<<"]"<< endl;
		FeatureServer fs(config,featureFilename);
		String labelSelectedFrames=config.getParam("labelSelectedFrames");
		
		SegServer segmentsServer;
		LabelServer labelServer;
		initializeClusters(featureFilename,segmentsServer,labelServer,config);
		verifyClusterFile(segmentsServer,fs,config);
		unsigned long codeSelectedFrame=labelServer.getLabelIndexByString(labelSelectedFrames);
		SegCluster& selectedSegments=segmentsServer.getCluster(codeSelectedFrame);  

		selectedSegments.rewind();
  
		FactorAnalysisStat FAAcc(featureFilename,fs,config);        
		FAAcc.estimateXYAndNorm(selectedSegments,fs,config);
		cms(featureFilename,fs,config);

		// Output the normalized features - Take care: only the saveFeatureFileExtension parameter make a difference between the input and output files
		if (!(config.existsParam("featureFlags")))
			config.setParam("featureFlags",fs.getFeatureFlags().getString());                 // Put the file flag in the config (could be different for each file   
		FeatureFileWriter w(featureFilename, config);                                        // build a featurefile writer to output the features
		
		if (writeAllFeature) {                                                            // Output all the features- feature count id the same
			SegServer fakeSegServer;                                                        // Create a new fake segment server
			fakeSegServer.createCluster(0);                                                 // Create a new cluster
			SegCluster& fakeSeg=fakeSegServer.getCluster(0);                                // Get the cluster               
			fakeSeg.add(fakeSegServer.createSeg(0,fs.getFeatureCount(),codeSelectedFrame,
							    labelSelectedFrames,featureFilename));       // Add a segment with all the features
			outputFeatureFile(config,fs,fakeSeg,w);       	                        // output all the features - giving the same file length
		}
		else
		outputFeatureFile(config,fs,selectedSegments, w);                       // Output only the selected features - giving a shorter output 
	}                                                                                   // End feature file loop
} // end try
catch (Exception & e)
{
	cout << e.toString ().c_str () << endl;
}
return 0;
}

//-------------------------------------------------------------------------------------------------------
//	Perform LFA normalisation of the features (Substract the Channel component in the feature space)
//-------------------------------------------------------------------------------------------------------
int normFeatLFA (Config & config) {

	String inputFeatureFileName =config.getParam("inputFeatureFilename");          // input feature - could be a simple feature file or a list of filename

	String labelSelectedFrames  =config.getParam("labelSelectedFrames");           // Only the frames from segments with this label  will be used
	bool writeAllFeature=(config.getParam("writeAllFeatures")=="true");            // Define if all the feature (selected or not) should be written
	if (verbose)	cout << "(NormFeat) Feature LFA Mode - Normalisation file by file" <<endl;

	XLine inputFeatureFileNameList;										// The (feature) input filename list
	if ( inputFeatureFileName.endsWith(".lst")) {						// If the file parameter is the name of a XList file
		XList inputFileNameXList(inputFeatureFileName,config);			// Read the filename list file
		inputFeatureFileNameList=inputFileNameXList.getAllElements();	// And put the filename in a list if the file is a list of feature filenames
	}
	else {																// It was a simple feature file and not a filename list
		inputFeatureFileNameList.addElement(inputFeatureFileName);		// add the filename in the list
	}

try{
	String *file;
		
	//Load Eigenchannel Matrix
	Matrix<double> U;
	String uName = config.getParam("matrixFilesPath") + config.getParam("eigenChannelMatrix") + config.getParam("loadMatrixFilesExtension");
	U.load (uName, config);

	while ((file=inputFeatureFileNameList.getElement())!= NULL){                   // Loop on each feature file

		String & featureFilename=(*file);                                            // Current file basename
		if (verbose) cout << "(NormFeat) Eigenchannel file normalisation["<<featureFilename<<"]"<< endl;

		FeatureServer fs(config,featureFilename);

		String labelSelectedFrames=config.getParam("labelSelectedFrames");
		
		SegServer segmentsServer;
		LabelServer labelServer;
		initializeClusters(featureFilename,segmentsServer,labelServer,config);

		verifyClusterFile(segmentsServer,fs,config);

		unsigned long codeSelectedFrame=labelServer.getLabelIndexByString(labelSelectedFrames);
		SegCluster& selectedSegments=segmentsServer.getCluster(codeSelectedFrame);  
		selectedSegments.rewind();

		JFAAcc jfaAcc(featureFilename,config);

		jfaAcc.loadEC(U, config);

		///Compute JFA stats
		jfaAcc.computeAndAccumulateJFAStat(config);

		jfaAcc.substractMplusDZByChannel();

		jfaAcc.substractMplusUX();

		///Estimate uEuT for the test
		jfaAcc.estimateUEUT(config);

		///Estimate and inverse L matrices
		jfaAcc.estimateAndInverseL_EC(config);

		///Estimate X for the test segment
		jfaAcc.estimateX(config);

		jfaAcc.estimateZMAP(config.getParam("regulationFactor").toLong());

		///Estimate XYZ on the test segment and normalise the features by substracting Ux
		jfaAcc.substractUXfromFeatures(fs,config);

		if(config.getParam("cms").toBool()){
			cms(featureFilename,fs,config);
		}

		// Output the normalized features - Take care: only the saveFeatureFileExtension parameter make a difference between the input and output files
		if (!(config.existsParam("featureFlags")))
			config.setParam("featureFlags",fs.getFeatureFlags().getString());				// Put the file flag in the config (could be different for each file   
		FeatureFileWriter w(featureFilename, config);										// build a featurefile writer to output the features

		if (writeAllFeature) {																// Output all the features- feature count id the same
			SegServer fakeSegServer;                                                        // Create a new fake segment server
			fakeSegServer.createCluster(0);                                                 // Create a new cluster
			SegCluster& fakeSeg=fakeSegServer.getCluster(0);                                // Get the cluster               
			fakeSeg.add(fakeSegServer.createSeg(0,fs.getFeatureCount(),codeSelectedFrame,
							    labelSelectedFrames,featureFilename));						// Add a segment with all the features
			outputFeatureFile(config,fs,fakeSeg,w);       									// output all the features - giving the same file length
		}
		else
		outputFeatureFile(config,fs,selectedSegments, w);									// Output only the selected features - giving a shorter output 
	}																						// End feature file loop
} // end try
catch (Exception & e)
{
	cout << e.toString ().c_str () << endl;
}
return 0;
}



#endif // !defined(ALIZE_NormFeat_cpp)
  

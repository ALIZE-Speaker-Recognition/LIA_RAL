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

#if !defined(ALIZE_EnergyDetector_cpp)
#define ALIZE_EnergyDetector_cpp

#include <iostream>
#include <fstream>  
#include <cassert> 
#include <cmath>
#include <liatools.h>
#include "SegTools.h"
#include "EnergyDetector.h"

#define DirectInterpolation 1

using namespace alize;
using namespace std;

// Plot the enrgy model
void plotEnergyDistrib(MixtureGD &mixt)
{
  unsigned long distribCount = mixt.getDistribCount();
  cout << "EnergyModel"<<endl;
  for (unsigned long c=0; c<distribCount; c++)
    cout << "Component["<<c<<"] Mean["<<mixt.getDistrib(c).getMean(0)<<
		 "] Cov["<<mixt.getDistrib(c).getCov(0)<<"] Weight["<<mixt.weight(c)<<"]"<<endl;
}	
//---------------------------------------------------------------------------------------------------------------
// find the highest energy component...
unsigned long findMaxEnergyDistrib(MixtureGD &mixt)
{
  unsigned long distribCount = mixt.getDistribCount();
  unsigned long cmpMax=0;
  for (unsigned long c=1; c<distribCount; c++)
    if (mixt.getDistrib(c).getMean(0)>mixt.getDistrib(cmpMax).getMean(0))
      cmpMax=c;
  if (verbose) cout << "Highest component["<<cmpMax<<"] Mean["<<mixt.getDistrib(cmpMax).getMean(0)<<
		 "] Cov["<<mixt.getDistrib(cmpMax).getCov(0)<<"] Weight["<<mixt.weight(cmpMax)<<"]"<<endl;
  return cmpMax;
}	
// find the lowest energy component...
unsigned long findMinEnergyDistrib(MixtureGD &mixt)
{
  unsigned long distribCount = mixt.getDistribCount();
  unsigned long cmpMin=0;
  for (unsigned long c=1; c<distribCount; c++)
    if (mixt.getDistrib(c).getMean(0)<mixt.getDistrib(cmpMin).getMean(0))
      cmpMin=c;
  if (verbose) cout << "Lowest component["<<cmpMin<<"] Mean["<<mixt.getDistrib(cmpMin).getMean(0)<<
		 "] Cov["<<mixt.getDistrib(cmpMin).getCov(0)<<"] Weight["<<mixt.weight(cmpMin)<<"]"<<endl;
  return cmpMin; 
}	

double computeEnergyThreshold(FeatureServer & fs,double pSelect,unsigned long nbBins=100)
{
  Histo histo(nbBins);                                             // Create an histo accumulator with 100 bins
  Feature f;                                                     // reset the reader at the begin of the input stream
  fs.reset();                                                     // feature server reset 
for (unsigned long ind=0;fs.readFeature(f); ind++) // feature loop
    histo.accumulateValue(f[0]);                               // Accumulate the energy in the histo Accumulator
  histo.computeHisto();                                           // Compute the histo
  long i=nbBins-1;  
  real_t count=0;
while((i>=0) && (count<=pSelect)){                               // Find the bin corresponding to the percentage of data wanted
    count+=histo.count(i)*(histo.higherBound(i)-histo.lowerBound(i));
    i--;
 }
  double threshold;
  if (i>=0) threshold=histo.higherBound(i);                        // Set the threshold to the higherBound of the next bin
  else threshold=histo.lowerBound(0);                              // if 100% of data should be selected
  if (verbose)  cout << "Percentage wanted["<<(int) (pSelect*100.0) <<"]Energy threshold["<<threshold<<"]"<<endl;
  return threshold;	
}

// Build the segments with the energized frames
unsigned long selectFrames(FeatureServer &fs,SegServer & segServer,double threshold,SegCluster &selectedSeg,SegCluster &outputSeg,String labelOutput,String fileName)
{
  unsigned long countFrames=0;
  fs.reset();                                                       // feature server reset
  unsigned long ind=0;
  unsigned long begin=0;
  bool in=false;
  Seg *seg;                                                         // current selectd segment
  selectedSeg.rewind();                                             // reset the reader at the begin of the input stream
  while((seg=selectedSeg.getSeg())!=NULL){                          // For each input segments
    for (unsigned long idx=seg->begin();idx<seg->begin()+seg->length();idx++){ // for each frame
      Feature f;
      fs.seekFeature(idx);
      fs.readFeature(f);
      if (f[0]>threshold){                                         // the frame is selected
	countFrames++;
	if (in==false){                                             // Begin of a new segment         
	  in=true;                                                  
	  begin=ind;
	}
      }
      else if (in){                                                // End of a segment
	in=false;
	Seg & segFake=segServer.createSeg(begin,ind-begin,0,       // Create a segment - Take care : length=end-begin+1 but ind =end+1 !!
					  labelOutput,fileName);
	outputSeg.add(segFake);                                  // Add a segment 	
      }
      ind++;                                                       // Increment the frame index
    }                                                              // end of one input segment
    if (in){                                                       // deal with the last energized segmeent inside the current input segment
      in=false;
      Seg & segFake=segServer.createSeg(begin,ind-begin+1,0,       // Create a segment 
					labelOutput,fileName);
      outputSeg.add(segFake);                                    // Add a segment  - Take care : length=end-begin+1 and ind=end in this case !!
    }                 
  }                                                              // end feature loop                   
  
  return countFrames;
}



//-------------------------------------------------------------------------
// Iniitialise the energy distrib to be as we want, we have an a priori on this distrib
//------------------------------------------------------------------------
MixtureGD &energyMixtureInit(MixtureServer &ms, StatServer &ss, FeatureServer &fs,MixtureGD 
&world,SegCluster &selectedSegments,const DoubleVector &globalCov, Config& config)
{
   unsigned long vectSize = world.getVectSize();	
  unsigned long distribCount = world.getDistribCount();
 
  double mean=-2.0;
  double meanIncrement;
  if (distribCount>1) meanIncrement=4.0/(double)(distribCount-1);
  else meanIncrement=1;
  for (unsigned long indg=0;indg<distribCount;indg++)  {        // For each component
    DistribGD& d = world.getDistrib(indg);                      // Get it
    for (unsigned long c=0;c<vectSize;c++,mean+=meanIncrement){                     // copy it
      d.setCov(1.0, c);
      d.setMean(mean, c);
    }
    d.computeAll();
  }  
  world.equalizeWeights();                                      // set weight = 1/distribCount for each distrib 
  return world;
}


//---------------------------------------------------------------------------------------------------------------
// Energy detector
// TAKE CARE !!!! INPUT DATA IS A PARAMETER FILE WITH ONLY THE ENERGY AS THE COEFF 0
// Use the featureMask ALIZE option for selecting the ernegy
SegCluster&  energyDetector(Config& config,SegServer &segServer,String &featureFileName)
{
  int nbTrainIt = config.getParam("nbTrainIt").toLong();	         // number of train it
  String thresholdMode;                                                  // Select the mode for the threshold (default=weight)
  if (config.existsParam("thresholdMode"))                               // weight -> select wh(+alpha*wm) % of the frames 
    thresholdMode= config.getParam("thresholdMode");                     //     (wh is the weight of the highest comp, wm, the weight of the midle component
  else thresholdMode="meanStd";                                          // meanStd -> the threshold is mean of the highest comp - alpha*Std of the same comp
  double alpha = config.getParam("alpha").toDouble();                    // alpha is the percentage of central gaussian framess taken into account
  double flooring = config.getParam("varianceFlooring").toDouble();      // Variance flooring and ceiling 
  double ceiling  = config.getParam("varianceCeiling").toDouble();  
  String labelSelectedFrames = config.getParam("labelSelectedFrames");   // Label of the selected frames/segments in the input file
  String labelOutput = config.getParam("labelOutputFrames");             // Label of the selected frames in the output file
  
  if (verbose) cout << "Proceeding Energy based silence detection for ["<<featureFileName<<"]"<<endl;
  FeatureServer fs(config,featureFileName);                            // Reading the feature file and create a feature server to deal with the data
  LabelServer labelServer;                                                              // Create the lable server, for indexing the segments/clusters
  initializeClusters(featureFileName,segServer,labelServer,config);                     // Reading the segmentation files for each feature input file
  unsigned long codeSelectedFrame=labelServer.getLabelIndexByString(labelSelectedFrames);            // Get the index of the cluster with in interest audio segments    
  SegCluster& selectedSegments=segServer.getCluster(codeSelectedFrame);   // Gives the cluster of the selected/used segments   
  //FrameAcc energyAccu(fs.getVectSize());
  //computeZeroOneOnACluster(energyAccu, fs, selectedSegments, config);
  unsigned long codeOutput=segServer.getClusterCount();
  segServer.createCluster(codeOutput);                                    // Create a new cluster for the output segments
  SegCluster& outputSeg=segServer.getCluster(codeOutput);                 // Get the cluster  
  MixtureServer ms(config);                                          // Create a mixture server in order to build an energy model
  StatServer ss(config,ms);                                          // Create a statistic server for model learning
  
  FrameAccGD globalFrameAcc;                                         // Create a frame accumulator for mean/cov computation
  globalMeanCov (fs,selectedSegments,globalFrameAcc,config);         // Compute the global mean and cov on the input segments
  
  DoubleVector globalMean=globalFrameAcc.getMeanVect();              // Get the Mean
  DoubleVector globalCov=globalFrameAcc.getCovVect();                // Get the Cov
  if (verboseLevel>1){
	cout <<"global mean and cov"<<endl;
	for (unsigned i=0; i < fs.getVectSize(); i++)cout << "mean[" << i << "=" << globalMean[i] << "]\tcov[" << globalCov[i] << "]" << endl;
  }
  MixtureGD & energyModel=ms.createMixtureGD();	               // Creating the energy model
  //  energyMixtureInit(ms,ss,fs,energyModel,selectedSegments,globalCov,config); // Init the energy model with randomly picked data 
  energyMixtureInit(ms, ss, fs,energyModel,selectedSegments,globalCov,config);  // Fixed init
  if (verboseLevel>1)plotEnergyDistrib(energyModel);
  MixtureStat &emAcc=ss.createAndStoreMixtureStat(energyModel);
  for (int trainIt=0; trainIt<nbTrainIt; trainIt++){                 // Training It loop
    emAcc.resetEM();                                                 // EM stuff
    double llkPreviousIt=accumulateStatEM(ss,fs,emAcc,selectedSegments,config); // Accumulate EM statistics 
    energyModel = emAcc.getEM();                                     // Get the EM estimate
    varianceControl(energyModel,flooring,ceiling,globalCov);         // Apply the variance normalisation		  	  
    if (verbose) cout << "Partial Train it["<<trainIt-1 <<"] (-1 means initial partial it) LLK="<<llkPreviousIt<<" Nb Frames="
		      <<emAcc.getEMFeatureCount()<<endl;
    if (verboseLevel>1)plotEnergyDistrib(energyModel);
  } 
  unsigned long nbInitialFrames= (unsigned long) emAcc.getEMFeatureCount(); // TODO, REPLACE THIS BY A COUNT OF selectedSegments duration - 
  double threshold=0;                                                  // The Energy threshold 
  unsigned long higher=findMaxEnergyDistrib(energyModel);	    
  if (thresholdMode=="weight"){                                      // The selection is based on the weight of the highest energy component
    double selectedWeight=energyModel.weight(higher);
    if  (energyModel.getDistribCount()== 3){
      unsigned long lower=findMinEnergyDistrib(energyModel);
      unsigned long middle=3-(higher+lower);
      double lossH=likelihoodLoss(energyModel.getDistrib(middle),energyModel.weight(middle),
				  energyModel.getDistrib(higher),energyModel.weight(higher));
      double lossL=likelihoodLoss(energyModel.getDistrib(middle),energyModel.weight(middle),
				  energyModel.getDistrib(lower),energyModel.weight(lower));
      if (verbose) cout  << "lossL["<<lossL<<"] lossH["<<lossH<<"]"<<endl;
      if (lossH<lossL)
	selectedWeight+=alpha*energyModel.weight(middle);		  
    }
    threshold=computeEnergyThreshold(fs,selectedWeight);		
  }
  else if (thresholdMode=="meanStd"){                                     // the Energy threshold is the mean of the highest energy component minus alpha*std 
    threshold=energyModel.getDistrib(higher).getMean(0) - (alpha*sqrt(energyModel.getDistrib(higher).getCov(0)));
  }
  else cout << thresholdMode<< "unknown"<<endl;
  unsigned long nbCurrentFrames=selectFrames(fs,segServer,threshold,selectedSegments,outputSeg,labelOutput,featureFileName); // Build the segments with the energized frames
  if (verbose){
    double selected=(float) nbCurrentFrames/ (float)nbInitialFrames;
    cout <<"File ["<<featureFileName<<"]	Number of initial frames ["<<nbInitialFrames<<"]	Number of selected frames ["
	 <<nbCurrentFrames << "]	percentage selected [" << (int) (selected*100) << "]" << endl;
  }
  return outputSeg;
}


#endif //!defined(ALIZE_EnergyDetector_cpp)

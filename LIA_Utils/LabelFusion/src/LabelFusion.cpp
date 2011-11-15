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

#if !defined(ALIZE_LabelFusion_cpp)
#define ALIZE_LabelFusion_cpp

#include <iostream>
#include <fstream>  // pour outFile
#include <cstdio>   // pour printf()
#include <cassert> // pour le debug pratique
#include <cmath>
#include <liatools.h>
#include "LabelFusion.h"


using namespace alize;
using namespace std;



// return in clusterOutput the clusterOne without the data in clusterTwo
void interSeg(SegCluster& clusterOne,SegCluster& clusterTwo,SegCluster& clusterOutput)
{
  unsigned long overlapTime=0;
  clusterOne.rewind();
  clusterTwo.rewind();
  SegServer & segServerOutput=clusterOutput.getServer();                                     // Get the clusterserver reelated to the output
  Seg *segOne;                                                                               // Will give the current segment in the hypothese
  Seg *segTwo=clusterTwo.getSeg();                                                           // Will give the current segment in the ref
  while((segOne=clusterOne.getSeg())!=NULL){                                                 // For each segment of the first file
    while ((segTwo) && (endSeg(segTwo)<segOne->begin()))
      segTwo=clusterTwo.getSeg();
    if ((segTwo) && (segTwo->begin()<=endSeg(segOne))){                                      // there is an  overlap between the two segments
      if (debug){
	cout << "overlap between :"<<endl;
	cout << segOne->sourceName()<<" "<<segOne->begin()<<" "<<segOne->length()<<endl;
	cout << segTwo->sourceName()<<" "<<segTwo->begin()<<" "<<segTwo->length()<<endl;   
      }
      if (segTwo->begin()<=segOne->begin()){                                                  // The seg2 begin before the seg one
	if (endSeg(segTwo)>=endSeg(segOne)){                                                  // the seg2 ending after the seg 1, seg1 is suppressed
	  if (debug) cout << "seg suppressed"<<endl;
	  overlapTime+=segOne->length();
	}
	else{                                                                                 // seg2 is ending before seg1, seg1 is shorter
	  Seg &seg=segServerOutput.createSeg(endSeg(segTwo)+1,endSeg(segOne)-endSeg(segTwo)+1,0,segOne->string(),segOne->sourceName());
	  if (debug) cout << "Modify1 the seg "<<seg.begin()<<" "<<seg.length()<<endl;   
	  clusterOutput.add(seg);
	  overlapTime+=segOne->length()-seg.length();
	}
      }
      else{                                                                                   // Seg2 is beginning after seg1
	Seg &seg=segServerOutput.createSeg(segOne->begin(),segTwo->begin()-segOne->begin(),0,segOne->string(),segOne->sourceName());
	if (debug) cout << "Modify2 the seg "<<seg.begin()<<" "<<seg.length()<<endl;   
	clusterOutput.add(seg);
	overlapTime+=segOne->length();
	overlapTime-=seg.length();
	if (endSeg(segTwo)<endSeg(segOne)){                                                  // Seg2 is ending before SegOne, adding a seg
	  Seg &seg=segServerOutput.createSeg(endSeg(segTwo)+1,endSeg(segOne)-endSeg(segTwo),0,segOne->string(),segOne->sourceName());
	  if (debug) cout << "Adding the seg "<<seg.begin()<<" "<<seg.length()<<endl;   
	  clusterOutput.add(seg);
	  overlapTime-=seg.length();
	} 
      } 
    }
    else{                                                                                     // No overlap, put the seg in the output cluster
      Seg & seg=segServerOutput.createSeg(segOne->begin(),segOne->length(),0,segOne->string(),segOne->sourceName());
      if (debug) cout << "saving the seg "<<seg.begin()<<" "<<seg.length()<<endl;   
      clusterOutput.add(seg);
    }
  }
  if (verbose) cout <<"InterSeg, total overlap time (suppressed frames) ["<<overlapTime<<"]"<<endl;
}

bool selFrame(double frm,double length, double threshold)
{
  return (frm/length>=threshold);  
}
// return in clusterOut the clusterIn after applkying morphological rules
// Apply a time window (length winLength)  and select/converts the label only if 
//    at least selectThreshold of the window frame are in the cluster. 
void morphologicalFilter(SegCluster& clusterIn,SegCluster& clusterOut,
			 unsigned long winLength,double selectThreshold)
{
  if (verbose) cout <<"beginning the morphological filter"<<endl;
  SegServer & segServerOut=clusterOut.getServer();                                       // Get the clusterserver reelated to the output
  // Step I
  clusterIn.rewind(); 
  Seg *segIn=clusterIn.getSeg();                                                         // Will give the current segment in the input
  if (segIn!=NULL){
    unsigned long suppressedTime=0;
    unsigned long addedTime=0;
    String label=segIn->string();
    String sourceName=segIn->sourceName();
    unsigned long winIdx=segIn->begin();                                                 // Pos of the window
    unsigned long winFrm=segIn->length();                                                // Number of frames (from the clusterIn) in the window
    unsigned long end=endSeg(segIn);
    while(segIn!=NULL){                                                                  // General loop for visiting all the seg
      while (((segIn=clusterIn.getSeg())!=NULL)&&(endSeg(segIn)<winIdx+winLength)){      // For one window
	winFrm+=segIn->length();
	end=endSeg(segIn);
      }
      if ((winFrm)&& selFrame(winFrm,winLength,selectThreshold)){       // A good window is selected
	unsigned long length=end-winIdx+1;
	Seg &seg=segServerOut.createSeg(winIdx,length,0,label,sourceName);
	if (debug) cout << "Adding the seg "<<seg.begin()<<" "<<seg.length()<<endl;   
	clusterOut.add(seg);
	addedTime+=(length-winFrm);
      }
      else suppressedTime+=winFrm;
      if (segIn!=NULL){                                                                  // It is not the end
	winIdx=segIn->begin();
	winFrm=segIn->length();
	end=endSeg(segIn);
      }
    }
    if (verbose){
      long diffTime=addedTime-suppressedTime;
      cout <<"morpho filter, suppressed frame["<<suppressedTime<<"] added frame["<<addedTime<<
		   "] total["<<diffTime<<"]"<<endl;
    }
  }
}
void morphologicalFilter(SegCluster& clusterIn,SegCluster& clusterOut,Config &config)
{
  unsigned long winLength=config.getParam("winLength").toLong();
  double selectThreshold=config.getParam("selectThreshold").toDouble();
  morphologicalFilter(clusterIn,clusterOut,winLength,selectThreshold);
}
//-------------------------------------------------------------------------
int labelFusion(Config& config)
{
  String extOutput=".lbl";                                               // the extension of the output files    
  if (config.existsParam("saveLabelFileExtension")) extOutput=config.getParam("saveLabelFileExtension");   
  String pathOutput;//="./";                                                // the path of the output files    
  String pathInput;//="./";                                                   // the path of the input files    
  if (config.existsParam("labelOutputPath")) pathOutput=config.getParam("labelOutputPath");

  if (config.existsParam("labelInputPath")){
	pathInput=config.getParam("labelInputPath");
	config.setParam("labelFilesPath", config.getParam("labelInputPath"));
  }
 
  String fileOut=config.getParam("labelOneFilename");
  String fileOne=config.getParam("labelOneFilename");
  String fileTwo=config.getParam("labelTwoFilename");
  if (config.existsParam("outputFilename"))
    fileOut=config.getParam("outputFilename");
  String labelSelectedFrames=config.getParam("labelSelectedFrames");

  try{
    SegServer segServerOne;                
    SegServer segServerTwo;	
    LabelServer labelServerOne;
    LabelServer labelServerTwo;

    loadClusterFile(fileOne,segServerOne,labelServerOne,config);
    if (debug) cout <<"label 1 loaded"<<endl;
    loadClusterFile(fileTwo,segServerTwo,labelServerTwo,config);
    if (debug) cout <<"label 2 loaded"<<endl;
    long codeSelectedFrameOne=labelServerOne.getLabelIndexByString(labelSelectedFrames);       // Get the index of the selected cluster
    if (codeSelectedFrameOne==-1){                                                             // No data for this model !!!!!!!!!!!!!!
      cout << " WARNING - NO DATA with the label["<<labelSelectedFrames<<"] in file ["<<fileOne<<"]"<<endl;
      exit(0);
    }
    long codeSelectedFrameTwo=labelServerTwo.getLabelIndexByString(labelSelectedFrames);       // Get the index of the selected cluster
    if (codeSelectedFrameTwo==-1){                                                             // No data for this model !!!!!!!!!!!!!!
      cout << " WARNING - NO DATA with the label["<<labelSelectedFrames<<"] in file ["<<fileTwo<<"]"<<endl;
    }
    SegCluster& clusterOne=segServerOne.getCluster(codeSelectedFrameOne);                 // Gives the cluster of the selected/used segments fil 1    
    SegCluster& clusterTwo=segServerTwo.getCluster(codeSelectedFrameTwo);                 // Gives the cluster of the selected/used segments file 2
    
    SegServer segServerOutput;
    SegCluster& clusterInt=segServerOutput.createCluster(0,labelSelectedFrames,clusterOne.sourceName());  
    interSeg(clusterOne,clusterTwo,clusterInt);
    SegCluster& clusterOutput=segServerOutput.createCluster(1,labelSelectedFrames,clusterOne.sourceName());    
    morphologicalFilter(clusterInt,clusterOutput,config);
    if (verbose){
      unsigned long init=totalFrame(clusterOne);
      unsigned long final=totalFrame(clusterOutput);
      long suppressed=init-final;
      cout <<"File["<<fileOne<<"] Initial number of Frame["<<init<<"] Final number of frame["<<final<<"] Suppressed ["<<suppressed<<"]"<<endl;
      cout << "Output the new label file in ["<<pathOutput+fileOut+extOutput <<"]"<<endl;
    }
    outputLabelFile(clusterOutput,pathOutput+fileOut+extOutput,config);
  } // fin try
  
  
  catch (Exception& e)
    { 
      cout << e.toString().c_str() << endl;
    }
  return 0;
}

//-------------------------------------------------------------------------
// Same than labelFusion but applies only the morphological filter one 1 file
int labelMorphing(Config& config)
{
  String extOutput=".lbl";                                               // the extension of the output files    
  if (config.existsParam("saveLabelFileExtension")) extOutput=config.getParam("saveLabelFileExtension");   
  String pathOutput;//="./";                                                // the path of the output files 
  String pathInput;//="./";                                                   // the path of the input files 
  if (config.existsParam("labelOutputPath")){ pathOutput=config.getParam("labelOutputPath");
  }
	
  if (config.existsParam("labelInputPath")){
	pathInput=config.getParam("labelInputPath");
	config.setParam("labelFilesPath", config.getParam("labelInputPath"));
  }
  
  String fileIn= config.getParam("labelFilename");
  String fileOut=fileIn;
  if (config.existsParam("outputFilename"))
    fileOut=config.getParam("labelFilename");
  String labelSelectedFrames=config.getParam("labelSelectedFrames");

  try{
    SegServer segServerIn;                
    LabelServer labelServerIn;
    loadClusterFile(fileIn,segServerIn,labelServerIn,config);
    if (debug) cout <<"label  loaded"<<endl;
    long codeSelectedFrameIn=labelServerIn.getLabelIndexByString(labelSelectedFrames);       // Get the index of the selected cluster
    if (codeSelectedFrameIn==-1){                                                             // No data for this model !!!!!!!!!!!!!!
      cout << " WARNING - NO DATA with the label["<<labelSelectedFrames<<"] in file ["<<fileIn<<"]"<<endl;
      exit(0);
    }
    SegCluster& clusterIn=segServerIn.getCluster(codeSelectedFrameIn);                 // Gives the cluster of the selected/used segments fil 1    
   
    SegServer segServerOutput;
    SegCluster& clusterOutput=segServerOutput.createCluster(1,labelSelectedFrames,clusterIn.sourceName());    
    morphologicalFilter(clusterIn,clusterOutput,config);
    if (verbose){
      unsigned long init=totalFrame(clusterIn);
      unsigned long final=totalFrame(clusterOutput);
      long suppressed=init-final;
      cout <<"File["<<fileIn<<"] Initial number of Frame["<<init<<"] Final number of frame["<<final<<"] Suppressed ["<<suppressed<<"]"<<endl;
      cout << "Output the new label file in ["<<fileOut+extOutput <<"]"<<endl;
    }

    outputLabelFile(clusterOutput,pathOutput+fileOut+extOutput,config);
  } // fin try
  
  
  catch (Exception& e)
    { 
      cout << e.toString().c_str() << endl;
    }
  return 0;
}


#endif //!defined(ALIZE_LabelFusion_cpp)

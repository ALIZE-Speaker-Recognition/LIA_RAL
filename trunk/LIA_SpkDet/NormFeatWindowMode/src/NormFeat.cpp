// NormFeat.cpp
// This file is a part of LIA Software LIA_SpkDet, based on ALIZE toolkit 
// LIA_SpkDet  is a free, open tool for speaker recognition
// LIA_SpkDet is a development project initiated and funded by the LIA lab.
// See www.lia.univ-avignon.fr
// 
// ALIZE is needed for LIA_SpkDet
// for more information about ALIZE, see http://www.lia.univ-avignon.fr/heberges/ALIZE/
//
// Copyright (C) 2004
//  Laboratoire d'informatique d'Avignon [www.lia.univ-avignon.fr]
//  Jean-Francois Bonastre [jean-francois.bonastre@lia.univ-avignon.fr]
//      
// LIA_SpkDet is free software; you can redistribute it and/or
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
// more information about the licence or the use of LIA_SpkDet
// Main Author : Alexandre PRETI
// version : 31 oct 2007

#if !defined(ALIZE_NormFeat_cpp)
#define ALIZE_NormFeat_cpp

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <liatools.h>
#include "NormFeat.h"

using namespace alize;
using namespace std;



//VERSION WITHOUT MEMORY
void computeCMVnorm(Feature &f, DoubleVector &cMean, DoubleVector &cCov ){
	
		
	for (unsigned long c = 0; c < cMean.size(); c++)
	{
	    f[c] = (f[c] - cMean[c])/cCov[c];
	}	
}



//INIT CMVNorm on a matrix of windowFrame Features (RefVector)

void loadMeanAndCovParam(RefVector <Feature> &windowFrameAcc, DoubleVector &cMean,DoubleVector &cCov, long &windowDuration )
{
    FrameAccGD local;
    //create local Feature
    Feature f(cMean.size());
    int i = 0;
	
    if(windowFrameAcc.size() == (unsigned long)windowDuration) 
    {
        while(i < windowDuration)
        {
	   f = windowFrameAcc[i];
	   local.accumulate(f);
	   i++;
        }
        cMean = local.getMeanVect();
        cCov  = local.getStdVect();
        local.reset();
	
    }
    else
    {
	cout << "Pb with buffer length" <<endl;
    }
   
}

//Update CMVNorm parameters when adding a new feature in the window (faster execution than loadMeanAndCovParam)

void updateMeanAndCovParam(Feature &f, DoubleVector &cMean,DoubleVector &cCov, long &windowDuration, long &frameCount, long &lookHead )
{	
    //Update window mean & var parameters
	
    //Forget factor Beta, need almost double precision
    double Beta;
    //Do not update mean & var parameters as current feature was already processed   
    if(frameCount < lookHead)
    {
	Beta = 1.0;    
    }
    //Update mean & var parameters for "not seen" feature  
    else
    {
        Beta = (((double)windowDuration - 1) / (double)windowDuration) ;
    }
    
    for (unsigned long i = 0; i < cMean.size(); i++)
    {
	cMean[i] = Beta * cMean[i]  + (1 - Beta)*f[i];
	cCov[i]  = sqrt( cCov[i] * cCov[i] * Beta + (1 - Beta) * (f[i]*f[i])); 
        
    }
    
   
}


//COMPUTE parameters on a windowDuration window, update the window by deleting the 1st frame and adding the current frame
void computeCMVparameters(Feature &f, long &windowDuration, RefVector <Feature> &windowFrameAcc, DoubleVector &cMean, DoubleVector &cCov, long &frameCount  )
{	
	// suppress one frame & add the current frame
	
	loadMeanAndCovParam(windowFrameAcc,cMean,cCov,windowDuration );
	
	windowFrameAcc.removeObject(0);
	windowFrameAcc.addObject(f);
	
	
}


/*
Variance and Mean cepstral normalization

Init : used zero padding + initDelay frames before norm OR init from a reference vector (TO DO)										     
For updating norm statistics a forget factor is used = windowFrame-1/windowFrame

TAKE CARE : works only on voiced feature, doesn't use labels

*/
int normFeatOnlineMode (Config & config) {
	
  String inputFeatureFileName =config.getParam("inputFeatureFilename");          // input feature - could be a simple feature file or a list of filename
  String labelSelectedFrames  =config.getParam("labelSelectedFrames");           // Only the frames from segments with this label  will be used 
 
  long windowDuration = 300, // Default length for window = 3s
       frameCount = 0,
       lookHead   = 0;
  bool 	initWithDelay = false;
  
  if (config.existsParam("initWithDelay"))
  { 
      lookHead = config.getParam("initWithDelay").toLong();
      initWithDelay = true;
	
  }

  if(lookHead > windowDuration)
  {
      cout << "BAD INIT DELAY, SET TO WINDOW DURATION"<<endl;
      lookHead = windowDuration;
  }
  
  bool writeAllFeature    = true; // Output a vector for all input vectors (selected and not selected vectors) - DEFAULT=on
  if (config.existsParam("writeAllFeatures"))
  {
      writeAllFeature=config.getParam("writeAllFeatures").toBool();    // Define if all the feature (selected or not) should be written 
  }
  
  if (verbose){
    cout << "NormFeat";
    cout << " Window mode, Window Duration["<<windowDuration<<"]"<<endl;
    cout << "var and mean normalisation"<< endl; 
  }
  double frameLength = 0.01; // length of a frame ins, by default 10 ms
  
  if (config.existsParam("frameLength")) frameLength=config.getParam("frameLength").toDouble();
 
  
	
  RefVector <Feature> memFrameAcc;
  RefVector <Feature> windowFrameAcc;
  
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
	
	  // Loop on each feature file
	while ((file=inputFeatureFileNameList.getElement())!= NULL)
	{                         
	    String & featureFileName=(*file);
	    FeatureFileReader fr(featureFileName, config);
		
	    if (!(config.existsParam("featureFlags")))
		config.setParam("featureFlags",fr.getFeatureFlags().getString());
	    
	    if (config.existsParam("windowDuration"))
            {	  
                windowDuration = config.getParam("windowDuration").toLong();
            }
  
	    cout << "Normalizing " << featureFileName << endl;
	    FeatureFileWriter w(featureFileName, config);                             // build a featurefile writer to output the features (real features)
		
	 
	    // current mean and cov vectors  
            DoubleVector cMean(fr.getVectSize(),fr.getVectSize()),cCov(fr.getVectSize(),fr.getVectSize());             
	    frameCount = 0;
	    // Loop on each feature
	    
	    // AT START ZERO PADDING + COPY OF LOOKHEAD DATA
	    if(initWithDelay)
	    {
		while(frameCount < windowDuration)
	        {   
		    Feature &f = *new Feature(fr.getVectSize());
		    
		    //ZERO PADDING
		    if ( frameCount < (windowDuration - lookHead))
		    {
		        f.reset(); 
		    }
		    else 
		    {
		        if(!fr.readFeature(f))
			{	
				//Case when lookhead > file size : window size = file size
				cout << "Reach end of file, new window size :  " << frameCount<<endl;
				windowDuration=frameCount;
				break;
			}
		    }
		
		    windowFrameAcc.addObject(f);
		    frameCount++;
	        }
		
		
		loadMeanAndCovParam(windowFrameAcc,cMean,cCov,windowDuration);
	        //END FILL 1st WINDOW (buffer)
		  
	    }
	    // AT START INIT WITH A MEAN VECTOR
	    else
	    {
		 // TO DO  
	    }
	    //START NORM	
	    fr.seekFeature(0);
	    frameCount = 0;
	   
	    Feature ftmp;
	    //NORMALISE FEATURES		    
	    while(fr.readFeature(ftmp))
	    {    
		Feature &f = *new Feature(ftmp);
		frameCount++;
		//UPDATE PARAMETERS 
		updateMeanAndCovParam(f, cMean,cCov, windowDuration, frameCount, lookHead );
		//computeCMVparameters(f, windowDuration, windowFrameAcc, cMean, cCov, frameCount  );
		//NORMALISE FEATURE		    
		computeCMVnorm(f, cMean, cCov);
		//cout <<"Frame Count :"<<frameCount<<endl;
		w.writeFeature(f);
		delete &f; //DELETE THIS IF USE OF computeCMVparameters
		
	    }
	    windowFrameAcc.deleteAllObjects();
	    
	    
        }//END LOOP ON FILES
 }	
		

  catch (Exception & e)
  {
      cout << e.toString ().c_str () << endl;
  }
  return 0;
}
		
		

#endif // !defined(ALIZE_NormFeat_cpp)

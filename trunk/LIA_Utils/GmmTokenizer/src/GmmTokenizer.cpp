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

#if !defined(ALIZE_GMMTokenizer_cpp)
#define ALIZE_GMMTokenizer_cpp

#include <iostream>
#include <fstream>  // pour outFile
#include <cstdio>   // pour printf()
#include <cassert> // pour le debug pratique
#include <cmath>
#include "GmmTokenizer.h"
#include "liatools.h"

using namespace alize;
using namespace std;

void computeConfusionMatrix(Feature & f,StatServer &ss, MixtureGDStat &acc, unsigned long & nBest, Matrix <unsigned long> & mce_matrix) {
	
	acc.computeAndAccumulateLLK(f,1.0,DETERMINE_TOP_DISTRIBS);     // Determine the winning components and compute world LLK
	const LKVector& v = ss.getTopDistribIndexVector();
        const unsigned long bestDistrib=v[0].idx;
        for (unsigned long i=0; i<nBest;i++)
		mce_matrix(bestDistrib,v[i].idx)++;
	acc.resetLLK();                    // Reset the world LLK accumulator
}

void computeConfusionMatrix(Seg * seg, FeatureServer & fs, StatServer & ss,  MixtureGDStat & acc, unsigned long & nBest,Matrix <unsigned long> &mce_matrix) {
	unsigned long idxBeginFrame=seg->begin()+fs.getFirstFeatureIndexOfASource(seg->sourceName()); 
	fs.seekFeature(idxBeginFrame); 
	Feature f;
	for (unsigned long idxFrame=0;idxFrame<seg->length();idxFrame++){    // For each frame of the segment
		fs.readFeature(f);
		computeConfusionMatrix(f,ss,acc,nBest,mce_matrix);
	} // frame loop
}

void computeConfusionMatrix(SegCluster & selectedSegments,  FeatureServer & fs,  StatServer & ss,  MixtureGD & world, unsigned long & nBest, Matrix <unsigned long> & mce_matrix) {
	MixtureGDStat &acc=ss.createAndStoreMixtureGDStat(world);
	Seg* seg;                                                                         // reset the reader at the begin of the input stream
	selectedSegments.rewind();                              
	while((seg=selectedSegments.getSeg())!=NULL){                                     // For each of the selected segments
		computeConfusionMatrix(seg, fs, ss,acc,nBest, mce_matrix) ;
	} // segments loop
}


void computeSymbols(Feature & f, MixtureGD & world, StatServer & ss, ULongVector & stream, Config & config) {
	ss.computeAndAccumulateLLK(world, f,DETERMINE_TOP_DISTRIBS);     // Determine the winning components and compute world LLK
	const LKVector& v = ss.getTopDistribIndexVector();
        const unsigned long bestDistrib=v[0].idx;
	stream.addValue(bestDistrib);
}

// on a segment
void computeSymbols(Seg* seg, FeatureServer & fs, MixtureGD & world, StatServer & ss, ULongVector & stream,Config & config) {
	unsigned long idxBeginFrame=seg->begin()+fs.getFirstFeatureIndexOfASource(seg->sourceName()); 
	fs.seekFeature(idxBeginFrame); 
	Feature f;
	for (unsigned long idxFrame=0;idxFrame<seg->length();idxFrame++){
		fs.readFeature(f);	
		computeSymbols(f,world,ss,stream,config);
	}
}
// on a cluster
void computeSymbols(SegCluster & selectedSegments, FeatureServer & fs, MixtureGD & world, StatServer & ss, ULongVector & stream, Config & config) {	
	Seg* seg;                       // reset the reader at the begin of the input stream
	selectedSegments.rewind();
	while((seg=selectedSegments.getSeg())!=NULL){ 
		computeSymbols(seg, fs, world, ss, stream, config); 	
	}
}

int GaussianConfusionMatrix(Config & config)
{
         String inputNDXFileName = config.getParam("inputFeatureFilename");                        // NDX inputfile filename - described the experience 
        String inputWorldFilename = config.getParam("inputWorldModelName");                   // World model file used for the LLR computation
        String labelSelectedFrames =config.getParam("labelSelectedFrames");              // label for selected frames - Only the frames from segment with this label  will be used
        //double frameLength = config.getParam("frameLength").toDouble();                  // length in s of a frame
        String matrixName=config.getParam("matrixOutputName");
        unsigned long nBest =config.getParam("topDistribsCount").toLong();  
	
        try{
                XList ndx(inputNDXFileName,config);                                    // Read the test definition file (ndx)
                XLine *linep;                                                          // Pointor on the current test line
                ndx.getLine(0);
                MixtureServer ms(config);
                StatServer ss(config, ms);
                MixtureGD& world = ms.loadMixtureGD(inputWorldFilename);               // Load the world model
		if (verbose) cout << "Confusion Matrix Memory Allocation" << endl;
                unsigned long model_size=world.getDistribCount(); 
                Matrix <unsigned long> mce_matrix;
		mce_matrix.setDimensions(model_size,model_size);
                
                while ((linep=ndx.getLine()) != NULL){                                 // Loop on each line of the ndx input file
                        String &featureFileName=linep->getElement(0);                        // Get the testfile basename
                        FeatureServer fs(config,featureFileName);                            // Reading the feature file
                        SegServer segmentsServer;                                                             // Create the segment server for managing the segments/clusters
                        LabelServer labelServer;                                                              // Create the lable server, for indexing the segments/clusters
                        initializeClusters(featureFileName,segmentsServer,labelServer,config);                // Reading the segmentation files for each feature input file
                        verifyClusterFile(segmentsServer,fs,config);                                          // Verify if the segments ending before the end of the feature files...
                        unsigned long codeSelectedFrame=labelServer.getLabelIndexByString(labelSelectedFrames);            // Get the index of the cluster with in interest audio segments
                        SegCluster& selectedSegments=segmentsServer.getCluster(codeSelectedFrame); // Gives the cluster of the selected/used segments   
                        if (verbose) cout << "test seg["<<featureFileName<<"]"<< endl;
                        computeConfusionMatrix(selectedSegments, fs, ss, world, nBest, mce_matrix);
		} // ndx loop
	mce_matrix.save(matrixName,config);
      
  } // fin try
 
  catch (Exception& e){ 
    cout << e.toString().c_str() << endl;
  }
  return 0;	     
        
}


int GMMTokenizer(Config & config)
{
        String inputNDXFileName = config.getParam("inputFeatureFilename");                        // NDX inputfile filename - described the experience 
        String inputWorldFilename = config.getParam("inputWorldModelName");                   // World model file used for the LLR computation
        String labelSelectedFrames =config.getParam("labelSelectedFrames");              // label for selected frames - Only the frames from segment with this label  will be used     
        String matrixName, symbolsFilesPath;
        symbolsFilesPath=config.getParam("symbolsFilesPath");
        
        try{
                XList ndx(inputNDXFileName,config);                                    // Read the test definition file (ndx)
                XLine *linep;                                                          // Pointor on the current test line
                ndx.getLine(0);
                MixtureServer ms(config);
                StatServer ss(config, ms);
                MixtureGD& world = ms.loadMixtureGD(inputWorldFilename);               // Load the world model
                
                while ((linep=ndx.getLine()) != NULL){                                 // Loop on each line of the ndx input file
                        String &featureFileName=linep->getElement(0);                        // Get the testfile basename
                        FeatureServer fs(config,featureFileName);                            // Reading the feature file
                        SegServer segmentsServer;                                                             // Create the segment server for managing the segments/clusters
                        LabelServer labelServer;                                                              // Create the lable server, for indexing the segments/clusters
                        initializeClusters(featureFileName,segmentsServer,labelServer,config);                // Reading the segmentation files for each feature input file
                        verifyClusterFile(segmentsServer,fs,config);                                          // Verify if the segments ending before the end of the feature files...
                        unsigned long codeSelectedFrame=labelServer.getLabelIndexByString(labelSelectedFrames);            // Get the index of the cluster with in interest audio segments
                        SegCluster& selectedSegments=segmentsServer.getCluster(codeSelectedFrame); // Gives the cluster of the selected/used segments   
                        if (verbose) cout << "test seg["<<featureFileName<<"]"<< endl;
                        ULongVector stream;
						computeSymbols(selectedSegments, fs, world, ss, stream,config);
                        String output="./"+symbolsFilesPath+"/"+featureFileName+".sym";
                        stream.save(output);
                } // ndx loop
 } // fin try
 
  catch (Exception& e){ 
    cout << e.toString().c_str() << endl;
  }
  return 0;	

}

#endif //!defined(ALIZE_GMMTokenizer_cpp)

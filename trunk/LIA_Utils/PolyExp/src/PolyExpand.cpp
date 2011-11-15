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

#include <iostream>
#include "alize.h"
#include "liatools.h"
#include <cmath>
using namespace std;
using namespace alize;

//**************************************************************************************************
// ---- Polynomial expension
// For a single frame 
void computeExpansion(Feature & f,Feature & expF) {
	expF[0]=1;// first value in expansion is one
	unsigned long size=f.getVectSize();
	for (unsigned long i=0;i<size;i++) // first [1 VectSize+1] values in the expension are the feature
		expF[i+1]=f[i];
	
	unsigned long idx=0;	// index of the expansion 
	//calculate expansion
	for (unsigned long i=0;i<size+1;i++){
		for (unsigned long j=i;j<size+1;j++){
			for (unsigned long k=j;k<size+1;k++){
				expF[idx] = expF[i]*expF[j]*expF[k] ;// all combinations with repetion
				idx++;
			}	
		}
	}
	
}

// Frame Level
void computeAndAccumulateExpansion(FeatureServer & fs,FrameAccGD & avgExp,unsigned long idxBeginFrame,unsigned long nbFrames,Config & config) {
	
	fs.seekFeature(idxBeginFrame);                                           // go to the frame in the buffer (and load it if needed)
	unsigned long vectSize=fs.getVectSize();	
	unsigned long expSize = (vectSize + 3) * (vectSize + 2) * (vectSize + 1) / 6 ;
	if (debug) cout << "Expension size: "<<expSize<<endl;
	Feature expF(expSize);
	for (unsigned long n=0;n<nbFrames;n++){
		Feature f;
		fs.readFeature(f);
		// Defines size of expansion
		computeExpansion(f,expF);
		avgExp.accumulate(expF);
	}
}

// Segment level
void computeAndAccumulateExpansion(FeatureServer & fs, FrameAccGD & avgExp, Seg * seg,Config & config) {
	unsigned long begin=seg->begin()+fs.getFirstFeatureIndexOfASource(seg->sourceName()); // Find the index of the first frame of the file in the buffer
	computeAndAccumulateExpansion(fs,avgExp,begin,seg->length(),config);
}
	
	
// Cluster level
void computeAndAccumulateExpansion(FeatureServer & fs, FrameAccGD & avgExp, SegCluster & selectedSegments,Config & config) {
	Seg* seg;                                                     // reset the reader at the begin of the input stream
	selectedSegments.rewind();      
	while((seg=selectedSegments.getSeg())!=NULL)                  // For each of the selected segments
		computeAndAccumulateExpansion(fs,avgExp,seg,config);
  }
// ---------------------- End of poly exp
//**************************************************************************************************

void multiplyByR(DoubleVector & avgExp, XList & rMat) {
	for (unsigned long i=0;i<avgExp.size();i++) {
		avgExp[i]*=rMat.getLine(i).getElement(0).toDouble();
	}
}

void zNorm(DoubleVector & avgExp, XList & rMat) {
	for (unsigned long i=0;i<avgExp.size();i++) {
		avgExp[i]-=rMat.getLine(i).getElement(1).toDouble(); // minus mean
		avgExp[i]/=rMat.getLine(i).getElement(0).toDouble(); // divides std
	}
}

void computeRSqrt(DoubleVector & R,unsigned long count) {
	for (unsigned long i=0;i<R.size();i++) {
		R[i]/=count;
		R[i]=1/sqrt(R[i]);
	}
}

void outputR(const DoubleVector & v1,const DoubleVector & v2,Config &config) {
	String outRFile=config.getParam("computeR");
	ofstream out(outRFile.c_str());
	for (unsigned long i=0;i<v1.size();i++) {
		out << v1[i] << " " <<  v2[i] << endl;
	}
	out << endl;
	out.close();
}
void outputInstanceSVMLight(const DoubleVector & v, String & filename,Config & config) {
	String exType=config.getParam("exType");
	ofstream out(filename.c_str());
	out << exType << " ";
	for (unsigned long i=0;i<v.size();i++) {
		out << i+1 << ":" << v[i] << " ";
	}
	out << endl;
	out.close();
}
void outputInstance(const DoubleVector & v, String & filename,Config & config) {
	String format=config.getParam("format");
	if (format=="SVMLight") outputInstanceSVMLight(v,filename,config);
	else {cerr << "E: Format unknown" << endl;}
}

// Compute speaker session var with full matrix
int PolyExpand(Config & config){
  try {
	XList inputList(config.getParam("inputFeatureFilename"),config);
	XLine * pLine;
	bool computeR=config.existsParam("computeR");
	bool normalize=config.existsParam("normalize");
	String labelSelectedFrames=config.getParam("labelSelectedFrames");
	bool verbose=false;
	if (config.existsParam("verbose")) verbose=config.getParam("verbose").toBool();
	XList xMat;
	if (normalize) {
		xMat.load(config.getParam("normalize"),config);
	}
// Accumulate avg expansion over multiple files
	FrameAccGD avgExp;
	while((pLine=inputList.getLine())!=NULL) {
		String filename=pLine->getElement(0);
		if (verbose) cout << "Processing file: ["<<filename<<"] ... ";
		FeatureServer fs(config,filename);	
                SegServer segmentsServer;                                              // Create the segment server for managing the segments/clusters
                LabelServer labelServer;                                               // Create the lable server, for indexing the segments/clusters
                initializeClusters(filename,segmentsServer,labelServer,config); // Reading the segmentation files for each feature input file
                verifyClusterFile(segmentsServer,fs,config);                    // Verify if the segments ending before the end of the feature files...
		unsigned long codeSelectedFrame=labelServer.getLabelIndexByString(labelSelectedFrames);// Get the index of the cluster with in interest audio segments
		SegCluster& selectedSegments=segmentsServer.getCluster(codeSelectedFrame); // Gives the cluster of the selected/used segments
		// Accumulates average poly expansion
		computeAndAccumulateExpansion(fs,avgExp,selectedSegments,config);
		if (verbose) cout << "Done" <<endl;
		// get mean exp
		if (!computeR) {
			DoubleVector avgExpVect=avgExp.getMeanVect();
			if (normalize) multiplyByR(avgExpVect,xMat);
			String outFile=config.getParam("vectorFilesPath")+filename+config.getParam("vectorFilesExtension");
			outputInstance(avgExpVect,outFile,config);
			avgExp.reset();
		}
	}
	if (computeR) {
		const DoubleVector & meanR=avgExp.getMeanVect();
		DoubleVector R=avgExp.getxAccVect();
		computeRSqrt(R,avgExp.getCount());
		// output R matrix;
		outputR(R,meanR,config);
	}
}
  catch (Exception& e) {cout << e.toString() << endl;}
return 0;
}

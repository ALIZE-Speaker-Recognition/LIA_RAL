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

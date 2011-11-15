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

/**
 * \file AcousticSegmentation.cpp
 * \author LIA
 * \version 2.0
 * \date june 2011
 *
 * \brief Short description of the speech activity detection (SAD) behavior
 *
 * The iterative SAD process is based on Viterbi decoding and model adaptation applied to a n-state HMM.\n
 * The n states represent the acoustic events expected in the signal (e.g: speech, music, speech+music, narrow/wide band, ...). They are each initialized with a GMM model trained on appropriate, separate data using an EM/ML algorithm (see LIA_SpkDet packages for the GMM model training).
 * State transition probabilities are fixed according to the parameters \a transitionMethod and \a gammaTransition.\n
 * Finally, minimum duration rules are applied in order to refine the acoustic segmentation yielded by the iterative process (possibility of replacing too short segments with adjacent segment class).
 *
**/

#if !defined(ALIZE_AcousticSegmentation_cpp)
#define ALIZE_AcousticSegmentation_cpp

#include "Tools.h"
#include "AcousticSegmentation.h"

using namespace alize;
using namespace std;

/*! \brief Adapts the gmm models of the n-state HMM according to the current segmentation
 */

void segAdaptationClass(Config& config,MAPCfg &mapCfg, hmm& actualHMM,ObjectRefVector &classModelTab,SegServer& previousSeg,SegServer& actualSeg,
		   StatServer& ss,FeatureServer &fs,MixtureServer& ms){
  
NoDataSpeakerVerification(config,actualHMM,actualSeg);

for(unsigned long icluster=0;icluster<actualSeg.getClusterCount();icluster++){ // For each cluster
    
	SegCluster& clusterA=actualSeg.getCluster(icluster);         // Get the current cluster of segments 
	
	if (verbose) cout<<" Adapt segments for speaker "<<clusterA.string()<<endl;

	if(previousSeg.getClusterCount() != 0){
			MixtureGD& classModel=(MixtureGD&)classModelTab.getObject(icluster);
			actualHMM.setDensity(classModel, icluster);		
			MixtureGD& m=(MixtureGD&)actualHMM.getDensity(icluster);
			adaptModel(config,ss,ms,fs,clusterA,classModel,m, mapCfg);
	}
}//for icluster
}//segAdaptation

/*! \brief Verifies whether the stop criteria (for adaptation/decoding step) are reached
 *
 * \return True if the segmentations (previous and current) between two consecutive SAD iterations are identical
 * otherwise false
 *
 * \remark The comparison between both segmentations is defined by :
 * -# the difference in the number of speakers
 * -# the difference in values of the Viterbi paths according to an epsilon value given in the configuration file
 *
 * \attention The parameter \a decodingEpsilon is required for step 2.\n
 * The value of the parameter \a decodingEpsilon defines the absolute margin of acceptable difference in values of Viterbi paths between both segmentations.
 */

bool isStopAcousticSegmentation(Config& config,SegServer& previousSeg,SegServer& segmentServer,DoubleVector& actualViterbiProb,DoubleVector& previousViterbiProb,hmm& actualHMM){

double epsilon=config.getParam("decodingEpsilon").toDouble();

unsigned long icluster;

//Comparison between the last two segmentations
if(verbose) cout<<">> Comparison between the 2 segmentations <<"<<endl;
if(previousSeg.getClusterCount()!=segmentServer.getClusterCount()){
	cout<<"ERROR : Problem with number of clusters! ("<<previousSeg.getClusterCount()<<" != "<<segmentServer.getClusterCount()<<endl;
	exit(-1);
}

for(icluster=0;icluster<previousSeg.getClusterCount();icluster++){
	if(isDifferentSegmentation(config, previousSeg.getCluster(icluster), segmentServer.getCluster(icluster))){
		if(!isLessEpsilon(config,actualViterbiProb,previousViterbiProb,epsilon))
			return true;
	}
}
if(debug) cout << "isComparable => true !!!"<< endl;
return false;

}//isStopAcousticSegmentation

/*! \brief Aggregate successive segments labeled in the same way
 *
 * \attention This procedure is run only if the parameter \a aggregateSegment is set to true.
 */

void aggregateSegment(SegServer& segmentServer){

unsigned long iseg=0;
if (verbose) cout << "Aggregate segments" << endl;
while(iseg<(segmentServer.getSegCount()-1)){
	Seg& segment1=segmentServer.getSeg(iseg);
	
	unsigned long ind=iseg+1;
	unsigned long l=segment1.length();
	unsigned long currentEnd=endSeg(&segment1);
	while((ind<segmentServer.getSegCount()) && (segment1.string() == segmentServer.getSeg(ind).string()) && ((currentEnd+1)==segmentServer.getSeg(ind).begin()) ){
		l+=segmentServer.getSeg(ind).length();
		currentEnd=endSeg(&(segmentServer.getSeg(ind)));
		ind++;
	}
	if(l>segment1.length()){ // aggregate !!!
		for(unsigned long i=iseg+1; i<ind; i++){
			if(verbose)
				cout << "segment: " << (segmentServer.getSeg(iseg+1)).begin() << " " << endSeg(&(segmentServer.getSeg(iseg+1))) << " " << (segmentServer.getSeg(iseg+1)).string() << " deleted" << endl;
			segmentServer.remove(segmentServer.getSeg(iseg+1));
		}
		segment1.setLength(l);
	}
	iseg++;
}
}

/*! \brief Apply length rules on the segments according to the config file
 *
 * \remark The successive segments with a same label for which the accumulated length is lower to the parameter \a class<n>Length are re-labelled with the previous label.\n
 * If the parameter \a aggregateSegment is set to true, the successive segments are aggregated within the same label.
 */

void applyLengthRules(Config config,SegServer& segmentServer, hmm &actualHMM){

bool change;

String musicLabel="M";
if (config.existsParam("musicLabel")) musicLabel=config.getParam("musicLabel");
bool aggregSegment=false;
if(config.existsParam("aggregateSegment")) aggregSegment=config.getParam("aggregateSegment").toBool();

for(unsigned long i=0; i<actualHMM.getNbState(); i++){
if(actualHMM.getLength(i) != 0){
do{
String previous=(segmentServer.getSeg(0)).string();
if (verbose) cout << "Apply Length Rules for " << actualHMM.getStateName(i) << endl;
unsigned long iseg=0;
change=false;
if(segmentServer.getSegCount() == 1) return;
while(iseg<(segmentServer.getSegCount()-1)){

	Seg& segment1=segmentServer.getSeg(iseg);
	if(segment1.string() == actualHMM.getStateName(i)){	
		unsigned long ind=iseg+1;
		unsigned long l=segment1.length();
		// accumulate length of successive segments labelling similarly
		while((segment1.string() == segmentServer.getSeg(ind).string()) && (ind<segmentServer.getSegCount()-1)){
			l+=segmentServer.getSeg(ind).length();
			ind++;
		}
	
		cout << "segment: " << segment1.string() << "length:" << l << " pour " << actualHMM.getLength(segment1.string()) << endl;
		if((actualHMM.getLength(segment1.string()) != 0) && (l<=actualHMM.getLength(segment1.string()))){
			// Transform all the segment strings to the next one
			for(unsigned long i=iseg; i<ind; i++){
				//String className=config.getParam("class"+String::valueOf(actualHMM.getReplacement(segment1.string())));

				//(segmentServer.getSeg(i)).setString(className); // here, the replacement label is fixed in the configuration
				if((segmentServer.getSeg(i)).string() != previous){ // in the case where the previous and current are the same !
					if(verbose)
						cout << "segment: " << (segmentServer.getSeg(i)).begin() << " " << endSeg(&(segmentServer.getSeg(i))) << " " << (segmentServer.getSeg(i)).string() << " becomes " << previous << endl;			
					(segmentServer.getSeg(i)).setString(previous);  // here, it is the previous label seen
					change=true;
				}
			}
		}
		iseg=ind;
	}
	else
		iseg++;
		
	if((segmentServer.getSeg(iseg-1)).string() != musicLabel)
		previous=(segmentServer.getSeg(iseg-1)).string();
}
}
while(change);

// for the last segment
unsigned long iseg=segmentServer.getSegCount()-1;
Seg& segment1=segmentServer.getSeg(iseg);
if((iseg-1) >= 0) {
	String previous=(segmentServer.getSeg(iseg-1)).string();
	if(segment1.string() == actualHMM.getStateName(i)){
		unsigned long l=segment1.length();
		if((actualHMM.getLength(segment1.string()) != 0) && (l<actualHMM.getLength(segment1.string()))){
			// Transform all the segment strings to the next one
//			String className=config.getParam("class"+String::valueOf(actualHMM.getReplacement(segment1.string())));
			if(verbose)
				cout << "segment: " << (segmentServer.getSeg(iseg)).begin() << " " << endSeg(&(segmentServer.getSeg(iseg))) << " " << (segmentServer.getSeg(iseg)).string() << " becomes " << previous << endl;			

//			(segmentServer.getSeg(iseg)).setString(className);
			(segmentServer.getSeg(iseg)).setString(previous);
		}
	}
}
if(aggregSegment) aggregateSegment(segmentServer); // if demanded, aggregate segments after each class 
}
}
}

/*! \brief The iterative SAD process (Speech Activity Detection)
 */

void AcousticSegmentationProcess(Config& config,SegServer& resegOutputServer,SegServer& segmentServer,
		    StatServer& ss,MixtureServer& ms,FeatureServer &fs,LabelServer& labelServer){


unsigned long iterationNb=0;
if (config.existsParam("iterationNb")) iterationNb=config.getParam("iterationNb").toLong();

resegOutputServer.removeAllClusters();
resegOutputServer.removeAllSegs();

MAPCfg mapCfg(config);
  
// Read the file and the current segmentation
if(verbose) cout<<">> Read the segmentation <<"<<endl;
  
labelServer.clear();

DoubleVector actualViterbiProb,previousViterbiProb;
SegServer previousSeg, decodingSeg;	
previousSeg.removeAllSegs();
previousSeg.removeAllClusters();

decodingSeg=segmentServer;
//SegCluster &cluster=decodingSeg.getCluster(0);
decodingSeg.removeAllClusters();
SegCluster& cluster=decodingSeg.createCluster();
for(unsigned long iseg=0; iseg<decodingSeg.getSegCount(); iseg++){
	cluster.add(decodingSeg.getSeg(iseg));
}

hmm actualHMM(ms,config); //Create the HMM

unsigned long classNb=config.getParam("classNb").toLong();
ObjectRefVector classModelTab; 
for(unsigned long i=1; i<=classNb; i++){
	String className=config.getParam("class"+String::valueOf(i));
	if(verbose) cout << "Add Model: " << className << endl;
	MixtureGD &m=ms.loadMixtureGD(className);
	classModelTab.addObject(m);
	// segment Length rules
	unsigned long length=0;
	if(config.existsParam("class"+String::valueOf(i)+"Length")){
		length=config.getParam("class"+String::valueOf(i)+"Length").toLong();
	}
	unsigned long replacement=1;
	if(config.existsParam("class"+String::valueOf(i)+"Replace")){
		replacement=config.getParam("class"+String::valueOf(i)+"Replace").toLong();
	}
	if(verbose) cout << "class: " << className << " replaced with " << replacement << " lenght " << length << endl;
	actualHMM.LoadState(m,className,length,replacement);
}

computeTransitions(config,actualHMM,segmentServer);


if((verbose) && (verboseLevel==2))  cout<<"HMM has "<<actualHMM.getNbState()<<" states"<<endl;
actualViterbiProb.clear();
previousViterbiProb.clear();
actualViterbiProb.setSize(cluster.getCount());
previousViterbiProb.setSize(cluster.getCount());
for(unsigned long iprob=0;iprob<cluster.getCount();iprob++){
	actualViterbiProb[iprob]=0;
	previousViterbiProb[iprob]=0;
}

if(verbose) cout<<">> First Decoding Pass : "<< endl;
previousSeg=segmentServer;
viterbiDecoding(config,actualHMM,cluster,segmentServer,ss,fs,labelServer,actualViterbiProb);
if(verbose){
	cout << "Segmentation after decoding" << endl;
	displayAllSegments(config, segmentServer);
}
unsigned long indice=0;
if(indice<iterationNb){
	do{
		if(verbose) cout<<">> Loop Decoding Pass - Adaptation : "<<indice<<" <<" << endl;
		segAdaptationClass(config,mapCfg, actualHMM,classModelTab,previousSeg,segmentServer,ss,fs,ms);
		previousSeg=segmentServer;
		viterbiDecoding(config,actualHMM,cluster,segmentServer,ss,fs,labelServer,actualViterbiProb);
		if(verbose){
			cout << "Segmentation after decoding" << endl;
			displayAllSegments(config, segmentServer);
		}
		indice++;
	}while((isStopAcousticSegmentation(config,previousSeg,segmentServer,actualViterbiProb,previousViterbiProb,actualHMM) && (indice<=iterationNb)));// adaptation-decoding loop
}

applyLengthRules(config,segmentServer,actualHMM);

resegOutputServer=segmentServer;
	
}//AcousticSegmentationProcess

/*! \brief Launching of the iterative SAD (Speech Activity Detection) process
 */

void launchAcousticSegmentationProcess(Config & config){

	String outputFilesPath=config.getParam("outputFilesPath");

	String inputListFileName = config.getParam("listFileToSegment");	//file including the list of files to segment
	XLine classToAnalyse;	//Array of labels to analyze
	classToAnalyse.reset();

	if(verbose){
		cout << "*********** Current Configuration ***************" << endl;
		for(unsigned long i=0; i<config.getParamCount(); i++){
			cout << config.getParamName(i) << " => " <<  config.getParamContent(i) << endl;
		}
		cout << "*************************************************" << endl;
	}

	try{
		XList listLabel;
		XList listFileName;
		try{
			listFileName.load(inputListFileName,config);
		}
		catch(FileNotFoundException& e){
			cout<<"There is no files to segment !"<<endl;
		      	exit(-1);
		}
		listFileName.rewind();
		XLine *filep;
		while ((filep=listFileName.getLine()) != NULL){							// For each stream of audio data (in several files in the same line)
			const XLine & listFile=filep->getElements();						// One or several files, as several part of the same stream
		      	MixtureServer ms(config);
	      		StatServer ss(config, ms);
		      	SegServer Resultat;
	      		FeatureServer fs(config,listFile);							// Reading the features (one or more files) 
		      	SegServer segmentsServer;								// Create the segment server for managing the segments/clusters
	      		LabelServer labelServer;								// Create the label server, for indexing the segments/clusters
		      	initializeClusters(listFile,segmentsServer,labelServer,config);				// Reading the segmentation files for each feature input file
	      		verifyClusterFile(segmentsServer,fs,config);						// Verify if the segment terminates before the end of the feature file
			String fileInit=listFile.getElement(0);
			config.setParam("fileSize", String::valueOf(fs.getFeatureCountOfASource(fileInit)));	

			if(config.existsParam("fileRefPath")){
				// assumption: all the segments in the segment server come from the same source file !!!
				displayAllSegmentsFromRef(config, fileInit, fs.getFeatureCountOfASource(fileInit));
			}

			SegServer segOutputServer;
			AcousticSegmentationProcess(config,segOutputServer,segmentsServer,ss,ms,fs,labelServer);
			// saving the segmentation
			for(unsigned long i=0;i<segOutputServer.getSegCount();i++){
				Seg& segment=segOutputServer.getSeg(i);
				Resultat.createSeg(segment.begin(),segment.length(),segment.labelCode(),segment.string(),segment.sourceName());
			}
			saveAcousticSegmentation(config,Resultat,fs,outputFilesPath,0);
		}// while
	} // end try
	catch (Exception& e){ 
		cout << e.toString().c_str() << endl;
	}
}//launchAcousticSegmentationProcess

#endif // !defined(ALIZE_AcousticSegmentation_cpp)

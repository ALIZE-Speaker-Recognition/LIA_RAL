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
 * \file ReSegmentation.cpp
 * \author LIA
 * \version 2.0
 * \date june 2011
 *
 * \brief Definition of speaker re-segmentation behavior
 *
 * The resegmentation step is normally applied in order to refine the segmentation outputs (issued from the segmentation step) and to remove irrelevant speakers.\n
 * A new HMM is generated from the segmentation output and an iterative speaker model training/Viterbi decoding loop is launched.
 * In contrast to the segmentation stage, here the MAP adaptation (coupled with a generic speech/background model) replaces the EM/ML algorithm for speaker model estimation since the segmentation step provides an initial distribution of speech segments among the different speakers detected.\n
 * For the resegmentation process, all the boundaries (except speech/non-speech boundaries) and segment labels are re-examined.
 *
**/

#if !defined(ALIZE_ReSegmentation_cpp)
#define ALIZE_ReSegmentation_cpp

#include "Tools.h"
#include "ReSegmentation.h"

using namespace alize;
using namespace std;

/*! \brief Check :
 *	- whether some global speaker durations (over the segmentation) are under a limit value fixed in the configuration file. In this case, speakers are removed from the segmentation
 *	- whether the stop criteria (for adaptation/decoding step) are reached
 *
 * \return True if Comparison between the last two segmentations (previous and current) are identical
 * otherwise false
 *
 * \remark The comparison between both segmentations is defined by :
 * -# the difference in the number of speakers
 * -# the difference in values of the Viterbi paths according to an epsilon value given in the configuration file
 *
 * \attention The parameter \a decodingEpsilon is required for step 2.\n
 * The value of the parameter \a decodingEpsilon defines the absolute margin of acceptable difference in values of Viterbi paths between both segmentations.
 */

bool isStopReSeg(Config& config,SegServer& previousSeg,SegServer& segmentServer,DoubleVector& actualViterbiProb,DoubleVector& previousViterbiProb,hmm& actualHMM){

unsigned long speakerMinTime=config.getParam("speakerMinTime").toLong();
bool deletedSpeakerFlag=false;
double epsilon=config.getParam("decodingEpsilon").toDouble();

//Checking of the speaker turn length according to the speakerMinTime parameter
unsigned long icluster,indCorr=0, speakerTime;

for(icluster=0;icluster<segmentServer.getClusterCount();icluster++){ 
	SegCluster& cluster=segmentServer.getCluster(icluster);
	speakerTime=totalFrame(cluster);
	
	if (verbose) cout<<" Time for speaker "<<icluster<<" "<<cluster.string()<<" is "<<speakerTime<<endl;
	if(speakerTime<speakerMinTime){//Speaker turn length is under the limit value => the speaker is removed !
		if((verbose) && (verboseLevel == 2)) cout<<" Delete of speaker "<<cluster.string()<<endl;
		actualHMM.deleteState(icluster-indCorr);
		indCorr++;
		// State transition computation to take the speaker deletion into account
		String transitionMethod=config.getParam("transitionMethod");
		double gamma=config.getParam("gammaTransition").toDouble();
		computeTransitions(config,actualHMM,transitionMethod,gamma,segmentServer);		
		
		// Deletion and labeling of cluster (toDelete) associated with the deleted speaker and 
		while(cluster.getCount()!=0){
			segmentServer.remove(cluster.get(0));		
		}
		cluster.setString("ToDelete");
		deletedSpeakerFlag=true;
	}

}//for icluster

// Deletion of segments associated with the deleted speakers 
if(deletedSpeakerFlag){
	// 
	bool flagDelete=false;
	do{
		flagDelete=false;
		icluster=0;
		while((icluster<segmentServer.getClusterCount())&&(!flagDelete)){
			SegCluster& cluster=segmentServer.getCluster(icluster);
			if(cluster.string()=="ToDelete"){
				segmentServer.remove(cluster);
				flagDelete=true;
			}
			icluster++;
		}//while icluster
	}while(flagDelete);
	
	return true;
}//if deletedSpeakerFlag

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
}//isStopReSeg

/*! \brief Build a new cluster and associated segments for the decoding step from an existing segment Server given as input
 *
 * \remark In this function, :
 * - Segments named "nonspeech" are removed
 * - the set of the remaining segments is re-named 'ToDecode0' and grouped together in a same cluster for further decoding
 * In this way, each decoding iteration being based on this cluster, labels and boundaries are re-examined at each of them
 */

void createClusterDecoding(Config& config,SegServer& decoding,SegServer& segmentServer, bool merge=true){

decoding=segmentServer;
decoding.removeAllClusters();//delete all the clusters
decoding.createCluster();
SegCluster& cluster=decoding.getCluster(0);

for(unsigned long i=0;i<decoding.getSegCount();i++)
	if(decoding.getSeg(i).string() == "nonspeech"){
		decoding.remove(decoding.getSeg(i));
		i--;
	}
	else{
		cluster.add(decoding.getSeg(i));
	}

unsigned long i=0;
bool flag=false;

if(merge){
	if((verbose) && (verboseLevel==2)){
		cout<<"Segmentation before Merge"<<endl;
		displayAllSegments(config,decoding);
	}
	//merge the segments
	do{
		i=0;
		flag=false;
		do{
			Seg& seg1=decoding.getSeg(i);
			seg1.setString("ToDecode0");
			for(unsigned long j=i+1;j<decoding.getSegCount();j++){
				Seg& seg2=decoding.getSeg(j);
				if(seg1.sourceName()==seg2.sourceName())
					if((seg1.begin()+seg1.length())==seg2.begin()){
						seg1.merge(seg2);
						flag=true;
						break;
					}
			}//for j
			i++;
		}while((!flag)&&(i<decoding.getSegCount()));
	}while(flag);
	if((verbose) && (verboseLevel==2)){
		cout<<"Segmentation after Merge"<<endl;
		displayAllSegments(config,decoding);
	}
}

}//createClusterDecoding

/*! \brief This procedure initializes the clusters
 */

void initClusterReSeg(Config& config,SegServer& segmentServer,
						  LabelServer& labelServer,StatServer& ss,FeatureServer &fs){
  
String loadSegmentationExtension=config.getParam("loadSegmentationExtension");
String labelPath=config.getParam("labelFilesPath");//path of the label files
segmentServer.removeAllClusters();
segmentServer.removeAllSegs();	
for (unsigned long i=0;i<fs.getSourceCount();i++){                 // For each of the feature file loaded in the fetaure server
	String file=fs.getNameOfASource(i);
   
	config.setParam("labelFilesExtension",loadSegmentationExtension);

	initializeClusters(file,segmentServer,labelServer,config);     // Read the corresponding labels 
}
}//initClusterReSeg

/*! \brief The speaker re-segmentation process
 */

void ReSegmentationProcess(Config& config,SegServer& resegOutputServer,SegServer& segmentServer,
		    StatServer& ss,MixtureServer& ms,FeatureServer &fs,LabelServer& labelServer){

String trainAlgo="MAP";
if(config.existsParam("trainAlgo")) trainAlgo=config.getParam("trainAlgo");
String testAlgo="Viterbi";
if(config.existsParam("testAlgo")) testAlgo=config.getParam("testAlgo");

DoubleVector actualViterbiProb,previousViterbiProb;
SegServer previousSeg;	

resegOutputServer.removeAllClusters();
resegOutputServer.removeAllSegs();

MAPCfg mapCfg(config);
  
// Read the file and the current segmentation
if(verbose) cout<<">> Read the segmentation <<"<<endl;
  
labelServer.clear();
initClusterReSeg(config,segmentServer,labelServer,ss,fs);
SegServer decoding;
createClusterDecoding(config,decoding,segmentServer, true);
SegCluster& cluster=decoding.getCluster(0);

MixtureGD &world=createWorld(config,cluster,ss,fs,ms);
String outputWorldFilename=fs.getNameOfASource(0);
world.save(outputWorldFilename, config);  
  
hmm actualHMM(ms,config); //Create the HMM
InitHMM(config,segmentServer,world,actualHMM); // init the Hmm
previousSeg.removeAllSegs();
previousSeg.removeAllClusters();

if((verbose) && (verboseLevel==2))  cout<<"HMM has "<<actualHMM.getNbState()<<" states"<<endl;
actualViterbiProb.clear();
previousViterbiProb.clear();
actualViterbiProb.setSize(cluster.getCount());
previousViterbiProb.setSize(cluster.getCount());
for(unsigned long iprob=0;iprob<cluster.getCount();iprob++){ 
	actualViterbiProb[iprob]=0;
	previousViterbiProb[iprob]=0;
}

int indice=0;
unsigned long iterationNb=20;
if (config.existsParam("iterationNb")) iterationNb=config.getParam("iterationNb").toLong();


do{
	indice++; 
	if(verbose) cout<<">> Loop Decoding Pass - Adaptation : "<<indice<<" <<" << endl;
	if(trainAlgo == "MAP")
		segAdaptation(config,mapCfg, actualHMM,world,previousSeg,segmentServer,ss,fs,ms);
	else {
//	if(trainAlgo == "FA")
//		segAdaptationFA(config,mapCfg, actualHMM,world,previousSeg,segmentServer,ss,fs,ms);
//	else {

		cout << "Unknown train Algo ! "<< endl;
		exit(-1);
	}
//	}
	previousSeg=segmentServer;
	if(testAlgo == "Viterbi")
		viterbiDecoding(config,actualHMM,cluster,segmentServer,ss,fs,labelServer,actualViterbiProb);
	else {
		cout << "Unknown test Algo ! "<< endl;
		exit(-1);
	}
		
	if(verbose){
		cout << "Segmentation after decoding" << endl;
		displayAllSegments(config, segmentServer);
	}	
}while((isStopReSeg(config,previousSeg,segmentServer,actualViterbiProb,previousViterbiProb,actualHMM)) && (indice < iterationNb));//boucle adaptation-decoding

resegOutputServer=segmentServer;
}//ReSegmentation

/*! \brief Launching of the process speaker re-segmentation
 */

void launchReSegmentationProcess(Config & config){

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
	      		LabelServer labelServer;								// Create the lable server, for indexing the segments/clusters
		      	initializeClusters(listFile,segmentsServer,labelServer,config);				// Reading the segmentation files for each feature input file
	      		verifyClusterFile(segmentsServer,fs,config);						// Verify if the segments ending before the end of the feature files
			String fileInit=listFile.getElement(0);
			config.setParam("fileSize", String::valueOf(fs.getFeatureCountOfASource(fileInit)));	

			if(config.existsParam("fileRefPath")){
				// assumption: all the segments in the segment server come from the same source file !!!
				displayAllSegmentsFromRef(config, fileInit, fs.getFeatureCountOfASource(fileInit));
			}

			SegServer segOutputServer;
			ReSegmentationProcess(config,segOutputServer,segmentsServer,ss,ms,fs,labelServer);
			// saving the segmentation
			for(unsigned long i=0;i<segOutputServer.getSegCount();i++){
				Seg& segment=segOutputServer.getSeg(i);
				Resultat.createSeg(segment.begin(),segment.length(),segment.labelCode(),segment.string(),segment.sourceName());
			}
			saveSegmentation(config,Resultat,fs,outputFilesPath,0);
		}// while
	} // end try
	catch (Exception& e){ 
		cout << e.toString().c_str() << endl;
	}
}//launchReSegmentationProcess

#endif // !defined(ALIZE_ReSegmentation_cpp)

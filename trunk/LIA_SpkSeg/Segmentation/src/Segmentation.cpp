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
 * \file Segmentation.cpp
 * \author LIA
 * \version 2.0
 * \date june 2011
 *
 * \brief Description of speaker segmentation behavior
 *
 * The speaker segmentation and clustering stage works directly on the SAD output.
 * The segmentation and clustering process relies on a one-step segmentation and clustering algorithm in the form of an evolutive hidden Markov model (E-HMM). 
 * Each E-HMM state aims to characterise a single speaker and the transitions represent the speaker turns. 
 * All possible changes between speakers are authorized.
 * The process for each audio show is as follows:
 * -# Initialization: 
 *  - The E-HMM has only one state L0
 *  - An iterative process is then started where a new speaker is added at each iteration.
 * -# A new speaker Lx is added to the E-HMM: 
 *  - According to the selection strategy used (defined through the selectionMethod parameter in the config file), a segment is chosen among all of the segments currently assigned to L0. 
 *  - The selected segment is attributed to Lx and is used to estimate a new GMM with EM training.
 * -# Adaptation/Decoding loop: 
 *  - The objective is to detect all segments belonging to the new speaker Lx. 
 *  - All speaker models are re-estimated through an adaptation process, according to the current segmentation (EM Algorithm) and a new segmentation is obtained via Viterbi decoding. 
 *  - This adaptation/decoding loop is repeated while some significant changes are observed on the speaker segmentation between two successive iterations.
 * -# Speaker model validation and stop criterion: 
 *  - The current segmentation is analyzed in order to decide if the new added speaker Lx is relevant, according to some heuristical rules on the total duration of speech assigned to speaker Lx.
 *  - The stop criterion is reached if there are no more segments available in L0 with which a new speaker can be added, otherwise the process goes back to step 2.
 *
**/

#if !defined(ALIZE_Segmentation_cpp)
#define ALIZE_Segmentation_cpp

#include "Tools.h"
#include "Segmentation.h"

using namespace alize;
using namespace std;

/*! \brief Selection of suitable segment for adding a new speaker
 *
 * \return True if suitable segment is found
 * otherwise false
 *
 * \remark 
 *
 * The method of selection depends on the \a selectionMethod parameter which the following values can be associated with:
 *  - \a 'firstSegmentFound' => Select the first segment which duration is over the \a selectionLength parameter
 *  - \a 'firstLimitedLengthSegmentFound' => Select the first segment which duration is over than the \a selectionLength parameter, but truncated to reflect the \a limitedLength parameter
 *  - \a 'longerSegmentFound' => Select the longest segment which duration is over the \a selectionLength parameter
 *  - \a 'longerLimitedLengthSegmentFound' => Select the longest segment which duration is over than the \a selectionLength parameter, but truncated to reflect the \a limitedLength parameter
 *
 * A segment is considered as suitable if :
 *  - the length of segment is superior or equal to the value of the parameter \a selectionLength
 *  - the pre-selected segment is not already used
 */

bool segmentSelection(Config& config, SegCluster& clusterSeg0,SegCluster& clusterAlreadyUsedSeg, 
			unsigned long &startSeg, unsigned long &stopSeg, String &fileLabel){

String selectionMethod="firstSegmentFound";
if(config.existsParam("selectionMethod")) (selectionMethod=config.getParam("selectionMethod"));
unsigned long selectionLength=(config.getParam("selectionLength")).toLong();// minimum size of the selected segment
bool find=false;
Seg *segment;
clusterSeg0.rewind();
DoubleVector length, starts, ends;                                        // The LLR value and the startindex of each of the potential trial
if(verbose)cout<<">> Segment Selection << "<< selectionMethod << endl;

/********* Select the first segment on the temporal scale which duration is over the selectionLength parameter *********/
if(selectionMethod=="firstSegmentFound"){
	while(!find && (segment=clusterSeg0.getSeg())!=NULL){
		if((segment->length()>=selectionLength) && (!findAlreadyUsed(clusterAlreadyUsedSeg,segment->begin(),segment->length(),segment->sourceName()))){                           // At least the segment is long enough to be selected
			startSeg=segment->begin();
			stopSeg=segment->begin()+segment->length()-1;
			fileLabel=segment->sourceName();
			find=true;
		} // if(segment
	}//while segment
}
/********* Select the first segment on the temporal scale which duration is over the selectionLength parameter, but truncated to reflect the limitedLength parameter *********/
else if(selectionMethod=="firstLimitedLengthSegmentFound"){
	unsigned long limitedLength=(config.getParam("limitedLength")).toLong();// maximum size of the selected segment

	while(!find && (segment=clusterSeg0.getSeg())!=NULL){
		if((segment->length()>=selectionLength) && (!findAlreadyUsed(clusterAlreadyUsedSeg,segment->begin(),segment->length(),segment->sourceName()))){                          // At least the segment is long enough to be selected
			startSeg=segment->begin();
			stopSeg=segment->begin()+segment->length()-1;
			/******** Limited in length *****/
			if((stopSeg-startSeg) > limitedLength) stopSeg=startSeg+limitedLength-1;
			fileLabel=segment->sourceName();
			find=true;
		} // if(segment
	}//while segment
}
/********* Select the longest segment which duration is over the selectionLength parameter *********/
else if(selectionMethod=="longestSegmentFound"){	      
	while((segment=clusterSeg0.getSeg())!=NULL){
		if((segment->length()>=selectionLength) && (!findAlreadyUsed(clusterAlreadyUsedSeg,segment->begin(),segment->length(),segment->sourceName()))){                           // At least the segment is long enough to be selected
			startSeg=segment->begin();
			stopSeg=segment->begin()+segment->length()-1;
			fileLabel=segment->sourceName();
			starts.addValue(startSeg);
			ends.addValue(stopSeg);
			length.addValue(segment->length());
			cout << startSeg*0.01<<" "<<stopSeg*0.01<<" " << endl;		
			find=true;
		} // if(segment
	}//while segment
	if(find){
		unsigned long imax = length.getIndexOfLargestValue();
		startSeg = (unsigned long) starts[imax];
		stopSeg = (unsigned long) ends[imax];	
	}
}
/********* Select the longest segment which duration is over the selectionLength parameter, but truncated to reflect the limitedLength parameter  *********/
else if(selectionMethod=="longestLimitedLengthSegmentFound"){	      
	unsigned long limitedLength=(config.getParam("limitedLength")).toLong();// maximum size of the selected segment

	while((segment=clusterSeg0.getSeg())!=NULL){
		if((segment->length()>=selectionLength) && (!findAlreadyUsed(clusterAlreadyUsedSeg,segment->begin(),segment->length(),segment->sourceName()))){  // A trial is defined{                         // At least the segment is long enough to be selected
			startSeg=segment->begin();
			stopSeg=segment->begin()+segment->length()-1;
			fileLabel=segment->sourceName();
			starts.addValue(startSeg);
			ends.addValue(stopSeg);
			length.addValue(segment->length());
			cout << startSeg*0.01<<" "<<stopSeg*0.01<<" " << endl;		
			find=true;
		} // if(segment
	}//while segment
	if(find){
		unsigned long imax = length.getIndexOfLargestValue();
		startSeg = (unsigned long) starts[imax];
		stopSeg = (unsigned long) ends[imax];
		/******** Limited in length *****/
		if((stopSeg-startSeg) > limitedLength) stopSeg=startSeg+limitedLength-1;
	}
}
else{
	cout << "Selection Method unknown" << endl;
}

return find;  
}//segmentSelection

/*! \brief Add a new speaker after selecting a segment
 *
 * \return 0 if suitable segment is found
 * otherwise -1 if no more segment available
 *
 * \remark This function performs :
 * - the deletion of the speakers for which the length in speech time is zero
 * - the research of an available segment
 * - the creation of a new speaker if a segment is available 
 */

int addSpeaker(Config& config, SegServer& actualSeg,SegServer& alreadyUsedSeg,
		     LabelServer& labelServer, int& numSpeaker,hmm& actualHMM,String clusterName){

unsigned long startSeg, stopSeg;
String fileLabel;
SegServer segTemp=actualSeg;

NoDataSpeakerVerification(config,actualHMM,actualSeg);
SegCluster& clusterSeg0=actualSeg.getCluster(0);
SegCluster& clusterAlreadyUsedSeg=alreadyUsedSeg.getCluster(0);     // The cluster of the not to use segments (already selected)

bool find=segmentSelection(config, clusterSeg0,clusterAlreadyUsedSeg,startSeg,stopSeg, fileLabel);

if(!find){
	if(verbose) cout << "No more segment available"<< endl;
	return -1;
}

numSpeaker++;
unsigned long nbrCluster=numSpeaker;
String nsuiv=clusterName+"_L"+String::valueOf(nbrCluster);
Label lsuiv(nsuiv);
unsigned long codeL=labelServer.addLabel(lsuiv);

// find a segment for initializing a new speaker model
if((verbose) && (verboseLevel==2)){
	cout<<"Cluster add for "<<nsuiv<<endl;
	cout<<"Chosen segment "<<startSeg*0.01<<" :"<<(stopSeg)*0.01<<endl;
}

//add cluster in segTemp
SegCluster& clustersuiv=segTemp.createCluster(codeL,nsuiv," ");
clustersuiv.add(segTemp.createSeg((unsigned long)startSeg,stopSeg-startSeg+1,codeL,nsuiv,fileLabel));
displayOneCluster(config, clustersuiv);

clusterAlreadyUsedSeg.add(alreadyUsedSeg.createSeg((unsigned long)startSeg,stopSeg-startSeg+1,codeL,nsuiv,fileLabel));
SegCluster& clusterSegT0=segTemp.getCluster(0);
removePartOfSeg(clusterSegT0,segTemp,labelServer,fileLabel,(unsigned long)(startSeg),stopSeg-startSeg+1);//Remove the trial in L0 cluster

if(verbose)
	displayAllClusters(config, segTemp);

// Modify the HMM structure
/////////////////////////
actualHMM.addState(nsuiv);
/////////////////////////
if((verbose) && (verboseLevel==2))  cout<<"HMM has "<<actualHMM.getNbState()<<endl;
computeTransitions(config,actualHMM,actualSeg);

actualSeg.removeAllClusters();
actualSeg.removeAllSegs();
actualSeg=segTemp;
 
return (0);
}//addSpeaker

/*! \brief Different tests done to verify the validity of the new speaker and the hmm
 *
 * \return True if the speaker is valid
 * otherwise false
 *
 * \remark A speaker is declared valid if the time of speech allocated to the last speaker is superior to the value of the parameter \a minSpeakerTime.
 */

bool isValidSpeaker(Config& config,hmm& actualHMM,hmm& previousHMM,SegServer& previousSeg,SegServer& actualSeg,SegServer& alreadyUsedSeg, unsigned long initSpeakerTime){

unsigned long minSpeakerTime=config.getParam("speakerMinTime").toLong();
//compute the time allocated to the last speaker N
unsigned long lastSpeakerTime=totalFrame(actualSeg.getCluster(actualSeg.getClusterCount()-1));

if(lastSpeakerTime<minSpeakerTime){
// This second condition may be useful for BN data
	if(verbose) cout<<" Time of Last Speaker is (in frames)"<<lastSpeakerTime<<" < "<<minSpeakerTime<<endl;
	actualSeg=previousSeg;
	actualHMM=previousHMM;
	return false;
}
  
if(verbose)cout<<"Time of Last Speaker is "<<lastSpeakerTime<<endl;
//compute the time allocated of the speaker N-1

return true;
}//isValidSpeaker

/*! \brief Verifies whether a pre-selected segment is not already used
 *
 * \return True if the pre-selected segment is not already used
 * otherwise false
 *
 */

bool isSegAvailable(Config& config,SegServer& alreadyUsedSeg,SegServer& actualSeg){

if(verbose) cout<<">> Segment available <<"<<endl;
unsigned long longSelection=(config.getParam("selectionLength")).toLong();
Seg *segment;
unsigned long ifeature;

SegCluster& clusterSeg0=actualSeg.getCluster(0);
SegCluster& clusterAlreadyUsedSeg=alreadyUsedSeg.getCluster(0);
clusterSeg0.rewind();
while((segment=clusterSeg0.getSeg())!=NULL){ // loop on the files
	if(segment->length()>=longSelection){
		for(ifeature=segment->begin();ifeature<segment->begin()+segment->length()-longSelection;ifeature+=longSelection){
			if(!findAlreadyUsed(clusterAlreadyUsedSeg,ifeature,longSelection,segment->sourceName())){
				return true;
			}//if(!find
		}//for ifeature
	}
}
return false;
}

/*! \brief Verifies that the stop criterion of a segmentation for a speaker is reached or not
 *
 * \return True if the stop criterion is reached
 * otherwise false
 *
 * \remark The stop criterion depend on the validity of the last speaker and the availability of new segments to be selected.
 */

bool isStopCriterionReached(Config& config,hmm& actualHMM,hmm& previousHMM,SegServer& actualSeg,SegServer& previousSeg,SegServer& alreadyUsedSeg,
			    StatServer& ss,FeatureServer &fs,unsigned long initSpeakerTime){
  
if(verbose) cout<<">> Speaker validation <<"<<endl;

if(isValidSpeaker(config,actualHMM,previousHMM,previousSeg,actualSeg,alreadyUsedSeg,initSpeakerTime)){
	return(!isSegAvailable(config,alreadyUsedSeg,actualSeg));
}
return false;
}

/*! \brief Set the parameter \a selectionMethod to the value \a firstSegmentFound if the speaker has only one turn of speech.
 */

void test_segmentation(Config &config, SegCluster& cluster){
	if(cluster.getCount() == 1){
		cout << "Caution: selection method has changed because there is only one segment" << endl;
		config.setParam("selectionMethod", "firstSegmentFound");
	}
}

/*! \brief The speaker segmentation process
 */

void segmentationProcess(Config& config, SegCluster& cluster,SegServer& segOutputServer,
		  StatServer& ss,FeatureServer &fs,MixtureServer&
		  ms,LabelServer& labelServer){

// variable definitions
segOutputServer.removeAllClusters();
segOutputServer.removeAllSegs();
String outputFilesPath=config.getParam("outputFilesPath");
unsigned long limitSpeaker=99999;
if(config.existsParam("limitSpeaker")) limitSpeaker=config.getParam("limitSpeaker").toLong();
hmm previousHMM(ms,config);	
DoubleVector actualViterbiProb,previousViterbiProb;
unsigned long initSpeakerTime=0;

test_segmentation(config, cluster);

// Train variable definition
String trainAlgo="EM";
if(config.existsParam("trainAlgo")) trainAlgo=config.getParam("trainAlgo");
String testAlgo="Viterbi";
if(config.existsParam("testAlgo")) testAlgo=config.getParam("testAlgo");

MAPCfg mapCfg(config);
MixtureGD& world=ms.createMixtureGD();
if(trainAlgo == "MAP"){
	world=createWorld(config,cluster,ss,fs,ms,config.getParam("mixtureDistribCount").toLong());
}

// Creation of L0 and corresponding cluster
if(verbose) cout<<"Create speaker L0 model"<<endl;
int numSpeaker = 0;
DoubleVector transitions;
hmm actualHMM(ms,config);
String et_temp=cluster.string()+"_L"+String::valueOf(numSpeaker);
Label l0(et_temp);
actualHMM.addState(et_temp);
transitions.addValue(1);

SegServer previousSeg,actualSeg;	
SegCluster& clusterSeg0=actualSeg.createCluster(labelServer.addLabel(l0),et_temp," "); //Create the cluster L0
SegServer alreadyUsedSeg;                                     //Create interdits cluster (it is the cluster with the segments to not be used)
alreadyUsedSeg.createCluster();
				  			
Seg *segment;
cluster.rewind();
while((segment=cluster.getSeg())!=NULL)
	clusterSeg0.add(actualSeg.createSeg(segment->begin(),segment->length(),labelServer.addLabel(l0),et_temp,segment->sourceName()));

unsigned long nbDistrib;

unsigned long iterationNb=10;
if (config.existsParam("iterationNb")) iterationNb=config.getParam("iterationNb").toLong();

// Begin of iterative process  
if(isSegAvailable(config,alreadyUsedSeg,actualSeg)){                          // verify if a segment is existing
	do{                                                                     // Segmentation main loop
		previousSeg=actualSeg;
		previousHMM=actualHMM;
      		SegServer segTemp;
		int indice=0;
		actualViterbiProb.clear();
		previousViterbiProb.clear();
		actualViterbiProb.setSize(cluster.getCount());
		previousViterbiProb.setSize(cluster.getCount());
		for(unsigned long iprob=0;iprob<cluster.getCount();iprob++){
			actualViterbiProb[iprob]=0;
			previousViterbiProb[iprob]=0;
		}
		
		//Use of current L0 for segment selection
		int status;
		do{	
			status=addSpeaker(config,actualSeg,alreadyUsedSeg,labelServer, numSpeaker,actualHMM,cluster.string());
			// No more segment available => End of the process
			if(status==-1){
				cout << "End of Processus " << endl;
				segOutputServer=actualSeg;
				return;
			}
			initSpeakerTime=totalFrame(actualSeg.getCluster(actualSeg.getClusterCount()-1));

			nbDistrib=config.getParam("mixtureDistribCount").toLong();
		}while(status!=0);

		
		do{
			indice++;
			if(verbose) cout<<">> Loop Training " << trainAlgo << " - Decoding Pass: "<<indice<< " <<" << endl;
			if(trainAlgo == "EM"){
				segEM(config,actualHMM,segTemp,actualSeg,ss,fs,ms,nbDistrib);
			}
			else if(trainAlgo == "MAP"){
				segAdaptation(config,mapCfg,actualHMM,world,segTemp,actualSeg,ss,fs,ms);
			}

			else {
				cout << "Unknown train Algo ! "<< endl;
				exit(-1);
			}

			segTemp=actualSeg;								
						
			if(testAlgo == "Viterbi")
				viterbiDecoding(config,actualHMM,cluster,actualSeg,ss,fs,labelServer,actualViterbiProb);
			else {
				cout << "Unknown test Algo ! "<< endl;
				exit(-1);
			}
			if((verbose) && (verboseLevel == 2)){
				cout << "Segmentation after adaptation/decoding " << endl;
				displayAllSegments(config,actualSeg);
			}

		}while((!isComparableAndVerifSpeaker(config,segTemp,actualSeg,actualViterbiProb,previousViterbiProb)) && (indice < iterationNb));//Training-decoding loop

				
	}while((!isStopCriterionReached(config,actualHMM,previousHMM,actualSeg,previousSeg,alreadyUsedSeg,ss,fs,initSpeakerTime)) && (actualSeg.getClusterCount()<limitSpeaker));// segmentation loop
}//if isSegAvai


segOutputServer=actualSeg;
cout << "Segmentation process finished for this label" << endl;

}

/*! \brief Launching of the process speaker segmentation
 */

void launchSegmentationProcess(Config & config){

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

	String initialSelectionMethod="firstSegmentFound";
	if(config.existsParam("selectionMethod")) (initialSelectionMethod=config.getParam("selectionMethod"));// size of the selected trials

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

			for(unsigned long icluster=0;icluster<segmentsServer.getClusterCount();icluster++){	// for each cluster
				SegCluster& cluster=segmentsServer.getCluster(icluster);
  				SegServer segOutputServer;
  				//if(verbose) cout<<"Launch Segmentation process "<<endl; 
  				segmentationProcess(config, cluster,segOutputServer,ss,fs,ms,labelServer);
				// Attention: the selection method may change if there is only one segment in the cluster
				config.setParam("selectionMethod", initialSelectionMethod);
  				for(unsigned long i=0;i<segOutputServer.getSegCount();i++){
    					Seg& segment=segOutputServer.getSeg(i);
	    				Resultat.createSeg(segment.begin(),segment.length(),segment.labelCode(),segment.string(),segment.sourceName());
	  			}		
			}//for icluster
			saveSegmentation(config,Resultat,fs,outputFilesPath,0);
		}// while
	} // end try
	catch (Exception& e){ 
		cout << e.toString().c_str() << endl;
	}
}//launchSegmentationProcess

#endif // !defined(ALIZE_Segmentation_cpp)

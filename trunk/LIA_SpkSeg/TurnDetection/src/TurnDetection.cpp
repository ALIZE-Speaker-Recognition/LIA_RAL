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
 * \file TurnDetection.cpp
 * \author LIA
 * \version 2.0
 * \date june 2011
 *
 * \brief Description of turn speaker detection behavior
 *
 * This performs a basic speaker turn detection and a local speaker clustering. It can be used as a preliminary step for the segmentation phase if segments
 * are finally assigned the general speech label (only the frontiers will be used to facilitate the segmentation process, not the labels !)\n
 * The algorithm is applied only to segments with length superior to 2 * \a winSize (parameter, default 50) :
 * - Speaker turn detection
 *   - use of classical criterion as BIC, GLR, DGLR or DELTABIC (parameter \a clusteringCrit, default DGLR)
 *   - criterion applied on 2 consecutive 0.5s long windows (parameter \a winSize) at 0.05s step (parameter \a winStep)
 *   - All the values Xi calculated on the set of segments can be viewed as a curve on which selected peaks are considered as speaker turns
 *   - Xi is a relevant maximum peak if \n
 *     Xi - MinLeft(Xi) > \a alpha * Std-Dev(X) and\n
 *     Xi + MinRight(Xi) > \a alpha * Std-Dev(X)
 *     - MinLeft(Xi) = smaller value localized on the left of Xi
 *     - MinRight(Xi) = smaller value localized on the right of Xi
 *     - Std-Dev(X) = standard deviation calculated with all values X on the curve
 *     - parameter \a alpha (default value 0.7)
 * - Local clustering
 *   - aggregation of consecutive segments is also based on a specified criterion associated with a decision threshold (parameter \a clusteringCritThresh)
 *
**/

#if !defined(ALIZE_TurnDetection_cpp)
#define ALIZE_TurnDetection_cpp

#include "Tools.h"
#include "TurnDetection.h"


using namespace alize;
using namespace std;

double myabs(double value){
double x=value;
	if (value < 0) x*=-1.0;
	return x;
}

/*! \brief The turn speaker detection process
 */

void TurnDetection(Config& config, SegCluster& cluster,SegServer& segOutputServer,
		  StatServer& ss,FeatureServer &fs,MixtureServer&
		  ms,LabelServer& labelServer){

SegServer segTemp;	

segOutputServer.removeAllClusters();
segOutputServer.removeAllSegs();

SegServer actualSeg;	
String et_temp="speech";
Label l(et_temp);
SegCluster& clusterSeg=actualSeg.createCluster(labelServer.addLabel(l),et_temp," "); //Create the cluster L


String crit="DGLR";
if(config.existsParam("clusteringCrit")) 
	crit=config.getParam("clusteringCrit");

double threshold=0.0;
if(config.existsParam("clusteringCritThresh"))
	threshold=config.getParam("clusteringCritThresh").toDouble();

unsigned long winSize=50;
if(config.existsParam("winSize")) winSize=config.getParam("winSize").toLong();
unsigned long winStep=5;
if(config.existsParam("winStep")) winStep=config.getParam("winStep").toLong();
double alpha=0.7;
if(config.existsParam("alpha")) alpha=config.getParam("alpha").toDouble();

unsigned long start1=0, end1=0;
unsigned long start2=0, end2=0;
unsigned long accu=0;


for(unsigned long iseg=0; iseg<cluster.getCount(); iseg++){
	
	
	Seg& segment=(Seg&)cluster.get(iseg);
	if(verbose)
		cout << "Segment" << iseg << ": " << segment.begin() << " " << endSeg(&segment) << endl; 
	if(segment.length() <= 2*winSize){
		clusterSeg.add(actualSeg.createSeg(segment.begin(),endSeg(&segment)-segment.begin()+1,0,segment.string(),segment.sourceName()));
		if(debug) cout << "add: " << segment.begin() << " " << endSeg(&segment) << endl;		
	}
	else{
		ObjectRefVector res;
		start1=segment.begin();
		end1=start1+winSize-1;
		start2=end1+1;
		end2=start2+winSize-1;
		accu = start1;
	
		while(end2 < endSeg(&segment)){
			if(verbose){
				cout << "Computation between: " << start1 << " " << end1; 
				cout << " and " << start2 << " " << end2 << endl; 
			}
			SegCluster& c1=segTemp.createCluster();
			c1.add(segTemp.createSeg(start1,winSize,0,"null",segment.sourceName()));
			SegCluster& c2=segTemp.createCluster();
			c2.add(segTemp.createSeg(start2,winSize,0,"null",segment.sourceName()));
			CritInfo *resCrit=new CritInfo(clusteringCriterionWithoutWorldInitOneGaus(config, c1, c2, ss, fs,crit),false,end1);
			
			res.addObject((Object&)*resCrit);	
			start1+=winStep;
			end1+=winStep;
			start2+=winStep;
			end2+=winStep;	
			
		}	
		
	
		/* smoothing */
	/*	for(unsigned long i=1; i<res.size()-1; i++){
			CritInfo &resCrit=(CritInfo&)(res.getObject(i));
			CritInfo &resCritP=(CritInfo&)(res.getObject(i-1));
			CritInfo &resCritN=(CritInfo&)(res.getObject(i+1));
		
			resCrit.setValue(0.25*resCritP.getValue()+0.25*resCritN.getValue()+0.5*resCrit.getValue());
		}	
	*/

           	DoubleVector score_buffer;
           	score_buffer.setSize(2);
           	score_buffer[0 % 2]=((CritInfo&)(res.getObject(0))).getValue();

           	for(unsigned long i=1; i<res.size()-1; i++)
           	{
               		CritInfo &resCrit=(CritInfo&)(res.getObject(i));
               		CritInfo &resCritN=(CritInfo&)(res.getObject(i+1));//right window

               		score_buffer[i % 2]=resCrit.getValue();
               		resCrit.setValue(0.25*score_buffer[(i-1) % 2]+0.25*resCritN.getValue()+0.5*resCrit.getValue());

           	}

		/* to look for maxima in the criterion value curve */
		/* if difference on left and right of a point with neighboor points is over alpha*standard deviation => maxima is found ! */
	
		double sum=0.0;
		double sum2=0.0;
		for(unsigned long i=0; i<res.size(); i++){
			CritInfo &resCrit=(CritInfo&)res.getObject(i);
			sum += resCrit.getValue();
			sum2+=resCrit.getValue()*resCrit.getValue();
		}
		double mean=sum/(double)res.size();
		double std=sqrt((sum2/(double)(res.size())-(mean*mean)));
	
		if(verbose){
			cout << "Mean and std: " << mean << " " << std << endl;
		}
	
		CritInfo &resCrit=(CritInfo&)res.getObject(0);
		resCrit.setDec(false);
		for(unsigned long i=1, j=0; i<res.size()-1; i++){
			/* for each value */
			/* search left min */
			j=i-1;
			double minL=((CritInfo&)res.getObject(i)).getValue();
			bool ok=true;
			while(ok && (j > 0)){
				if(((CritInfo&)res.getObject(j)).getValue() < minL){
					minL =((CritInfo&)res.getObject(j)).getValue();
					j--;
				}else{	
					ok = false;
				}
			}	
		
			if(myabs(((CritInfo&)res.getObject(i)).getValue()-minL) > alpha*std){
			// search right min 
		
				j=i+1;
				double minR=((CritInfo&)res.getObject(i)).getValue();
				ok=true;
				while(ok && (j < res.size())){
					if(((CritInfo&)res.getObject(j)).getValue() < minR){
						minR = ((CritInfo&)res.getObject(j)).getValue();
						j++;
					}else{
						ok = false;
					}
				}
		
				if(myabs(((CritInfo&)res.getObject(i)).getValue()-minR) > alpha*std){
					((CritInfo&)res.getObject(i)).setDec(true);				
				}else{
					((CritInfo&)res.getObject(i)).setDec(false);				
				}
			}else{
				((CritInfo&)res.getObject(i)).setDec(false);					
			}
			/*double max=((CritInfo&)res.getObject(i)).getValue();
			double minL=((CritInfo&)res.getObject(i-1)).getValue();
			double minR=((CritInfo&)res.getObject(i+1)).getValue();
			if((minL < max) && (minR < max) && (max > mean+std))
				((CritInfo&)res.getObject(i)).setDec(true);
			else
				((CritInfo&)res.getObject(i)).setDec(false);
			*/
		}
	
		start1 = segment.begin();
		for(unsigned long i=0; i<res.size(); i++){
			cout << ((CritInfo&)res.getObject(i)).getFrame() << " " << ((CritInfo&)res.getObject(i)).getValue() << " => " << ((CritInfo&)res.getObject(i)).getDec() << endl;
			if(((CritInfo&)res.getObject(i)).getDec()){
				clusterSeg.add(actualSeg.createSeg(start1,((CritInfo&)res.getObject(i)).getFrame()-start1+1,0,segment.string(),segment.sourceName()));
				if(verbose) cout << "add: " << start1 << " " << ((CritInfo&)res.getObject(i)).getFrame() << endl;
				start1=((CritInfo&)res.getObject(i)).getFrame()+1;
			}
		}
		// last point
		cout << "last point: " << start1 << " fin segment: " <<  endSeg(&segment) << endl;
		if(start1 < endSeg(&segment)){
			clusterSeg.add(actualSeg.createSeg(start1,endSeg(&segment)-start1+1,0,segment.string(),segment.sourceName()));
			if(verbose) cout << "add: " << start1 << " " << endSeg(&segment) << endl;	
		}

	}
}


displayAllClusters(config, actualSeg);
segOutputServer=actualSeg;
}

/*! \brief Launching of the turn speaker detection process
 */

void launchTurnDetectionProcess(Config & config){

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

			for(unsigned long icluster=0;icluster<segmentsServer.getClusterCount();icluster++){	// for each cluster
				SegCluster& cluster=segmentsServer.getCluster(icluster);
  				SegServer segOutputServer;
  				TurnDetection(config,cluster,segOutputServer,ss,fs,ms,labelServer);
				displayAllSegments(config,segOutputServer); 
  				for(unsigned long i=0;i<segOutputServer.getSegCount();i++){
    					Seg& segment=segOutputServer.getSeg(i);
	    				Resultat.createSeg(segment.begin(),segment.length(),segment.labelCode(),segment.string(),segment.sourceName());
	  			}
			}//for icluster
			saveSegmentation(config,Resultat,fs,outputFilesPath,1);
		}// while
	} // end try
	catch (Exception& e){ 
		cout << e.toString().c_str() << endl;
	}
}//launchTurnDetectionProcess

#endif // !defined(ALIZE_TurnDetection_cpp)

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

#if !defined(ALIZE_ShiftedDeltaFeat_cpp)
#define ALIZE_ShiftedDeltaFeat_cpp

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <liatools.h>
#include "ShiftedDeltaFeat.h"

using namespace alize;
using namespace std;

//---------------------------------------------------------------------------------------------------------------
/*!
 * \brief Reads a parametric file (or list of prm files) containing cepstral features and
 *        writes a file containing the shifted delta features.
 * \param config	Config object to fetch parameters from
 * \date 11 april 2008
 * 
 * More detailed description can be found in ShiftedDeltaFeatMain
 * 
 */
int ShiftedDeltaFeat (Config & config) {
  unsigned int numCeps = 0 ;	// N
  unsigned int deltaDelay ;		// d
  unsigned int timeAdvance ;	// P
  unsigned int deltaBlocks ;	// k
  bool keepCepstra = true ;		//default
 
  String inputFeatureFileName =config.getParam("inputFeatureFilename");          // input feature - could be a simple feature file or a list of filename
  if (config.existsParam("vectSize")) numCeps = config.getParam("vectSize").toLong();
  deltaDelay  = config.getParam("SDCdeltaDelay").toLong();
  timeAdvance = config.getParam("SDCtimeAdvance").toLong();
  deltaBlocks = config.getParam("SDCdeltaBlocks").toLong();
  if (config.existsParam("SDCkeepCepstra")) keepCepstra = config.getParam("SDCkeepCepstra").toBool();
  if (verbose){
    cout << "ShiftedDeltaFeat: SDC parameters (N-d-P-k) are: ["<<numCeps<<"-"<<deltaDelay<<"-"<<timeAdvance<<"-"<<deltaBlocks<<"]" << endl;
  }
  if(debug) cout << "SDC parameters: SDCkeepCepstra["<<keepCepstra<<"] vectSize=numCeps["<<numCeps<<"] SDCdeltaDelay["<<deltaDelay<<"] SDCtimeAdvance["<<timeAdvance<<"] SDCdeltaBlocks["<<deltaBlocks<<"]" << endl;
  
  
  // Begin !!
  XLine inputFeatureFileNameList;                                                // The (feature) input filename list
  if ( inputFeatureFileName.endsWith(".lst")){                                   // If the file parameter is the name of a XList file
	XList inputFileNameXList(inputFeatureFileName,config);                       // Read the filename list file
	inputFeatureFileNameList=inputFileNameXList.getAllElements();          		 // And put the filename in a list if the file is a list of feature filenames
  }
  else {                                                                         // It was a simple feature file and not a filename list
    inputFeatureFileNameList.addElement(inputFeatureFileName);                   // add the filename in the list
  }
  try{
	String *file;
	while ((file=inputFeatureFileNameList.getElement())!= NULL){                 // Loop on each feature file
		
		String & featureFileName=(*file);                                        // Current file basename
		FeatureServer fs(config,featureFileName);
		if(!numCeps) numCeps = fs.getVectSize();
		unsigned long numFeats = fs.getFeatureCount() ;

		Feature lowerFeat ;	//feature -d of a delta block
		Feature upperFeat ;	//feature +d of a delta block
		unsigned int cepFeatsInBase = (keepCepstra) ? numCeps : 0 ;		//number of additional feature values (if keeping cepstra)
		unsigned int shiftedDeltaVectSize = cepFeatsInBase + (deltaBlocks*numCeps) ;

		Config configOut( config) ;
		configOut.setParam("vectSize", String::valueOf(shiftedDeltaVectSize) ) ;
		configOut.setParam("featureServerMask", "0-"+String::valueOf(shiftedDeltaVectSize-1)) ;
		cout << "ShiftedDeltaFeat: numCeps["<<numCeps<<"] vectSize["<<fs.getVectSize()<<"] numFeats["<<fs.getFeatureCount()<<"] keepCeps["<<keepCepstra<<"] cepsInBase["<<cepFeatsInBase<<"] OUTnumFeats["<<shiftedDeltaVectSize<<"] OUTcfgVectSize["<<configOut.getParam("vectSize")<<"]" <<endl;
		//FeatureServer fsOut(configOut,featureFileName);
		FeatureServer fsOut(configOut);
		if(debug) cout << "FeatureServer fsOut instantiated" << endl ;

		// Output the SDC features - Take care: only the saveFeatureFileExtension parameter make a difference between the input and output files
		if (!(configOut.existsParam("featureFlags")))
			configOut.setParam("featureFlags",fs.getFeatureFlags().getString());// Put the file flag in the config (could be different for each file   
		cout << "Writing to: " << featureFileName << endl;
		FeatureFileWriter w(featureFileName, configOut);                             // build a featurefile writer to output the features (real features)
		
		// Begin the real work, file by file

		/*
		 * SDC formula is:
		 * 
		 * ShDeltaO_h(t) = O_h( t + i*P + d) - O_h( t + i*P - d)   ;  i = 0, 1, 2 ... k-1  ;  O_h is the input featureVector
		 */
		
		/*for cropping at start and end of file:*/// since we have -d, first t is d and last is Total-d; to avoid looking into nowhere
		/*for cropping at start and end of file:*///for( unsigned long currentTime=deltaDelay; (currentTime+((deltaBlocks-1)*timeAdvance)+deltaDelay)<numFeats; currentTime++)
		// LOOP over t - for all features in the input file - generates the same number of feature in output
		for( unsigned long currentTime=0; currentTime<numFeats; currentTime++) {

			// initialize a (big) featureVector with k*N or k*N+N features
			Feature shiftedDeltaFeature( shiftedDeltaVectSize) ;

			// TODO ev include feature to keep whole input vector in output - not only masked part
			// TODO   this adds the possibility to have ceps, energy, (ev. acceleration) and SDCs (without energy/accel) all together 

			if( cepFeatsInBase > 0) {
				assert( cepFeatsInBase==numCeps) ;
				Feature cepFeats;
				fs.seekFeature( currentTime) ;
				fs.readFeature( cepFeats,0) ;	//TODO chkret
				//TODO replace loop by fct call - new fct
				for( unsigned int i=0; i<numCeps; i++) {
					shiftedDeltaFeature[i] = cepFeats[i] ;	//TODO optimize with memcopy
				}
			}

			// LOOP over i from 0 to (k-1)*P
			for( unsigned int currentDeltaBlock=0; currentDeltaBlock<deltaBlocks; currentDeltaBlock++) {
				// time of lower (t-d) and upper (t+d) feature
				         int lowerTime = currentTime + (currentDeltaBlock*timeAdvance) - deltaDelay ;	// t + i*P - d
				unsigned int upperTime = currentTime + (currentDeltaBlock*timeAdvance) + deltaDelay ; 	// t + i*P + d

				// for begin and end of file, correct times as to simulate repeating edge vectors beyond limits
				if( lowerTime < 0        ) lowerTime = 0 ;
				if( (unsigned long)lowerTime >= numFeats) lowerTime = numFeats-1 ;
				if( upperTime                >= numFeats) upperTime = numFeats-1 ;  //last vector is numFeats-1
				//if (debug || (verboseLevel>2)) cout << "t["<<currentTime<<"] block["<<currentDeltaBlock<<"] lowerT["<<lowerTime<<"] 2d["<<(2*deltaDelay)<<"]" <<endl ;
				if (debug || (verboseLevel>2)) cout << "t["<<currentTime<<"] block["<<currentDeltaBlock<<"] lowerT["<<lowerTime<<"] upperT["<<upperTime<<"]" <<endl ;

				// fetch cepstral features
				fs.seekFeature( lowerTime) ;
				//fs.readFeature( lowerFeat, 2*deltaDelay);	// O_h( t + i*P - d)	//TODO chkret
				fs.readFeature( lowerFeat, 0);				// O_h( t + i*P - d)	//TODO chkret
				if (debug || (verboseLevel>2)) { cout << " lowerFeat[" ; for(unsigned int i=0; i<lowerFeat.getVectSize();i++){ cout<<lowerFeat[i]<<" " ;} ; cout <<"]"<< endl; }
				fs.seekFeature( upperTime) ;
				fs.readFeature( upperFeat);  				// O_h( t + i*P + d)	//TODO chkret
				if (debug || (verboseLevel>2)) { cout << " upperFeat[" ; for(unsigned int i=0; i<upperFeat.getVectSize();i++){ cout<<upperFeat[i]<<" " ;} ; cout <<"]"<< endl; } 

				// calculate delta
				Feature deltaBlockFeature( numCeps) ;
				for( unsigned int i=0; i<numCeps; i++) {
					deltaBlockFeature[i] = upperFeat[i] - lowerFeat[i] ;
				}
				if (debug || (verboseLevel>2)) { cout << "  deltaFeat[" ; for(unsigned int i=0; i<deltaBlockFeature.getVectSize();i++){ cout<<deltaBlockFeature[i]<<" " ;} ; cout <<"]"<< endl; }

				// calculate position to append deltas to output vector: i*N (or N + i*N)
				unsigned int blockBase = cepFeatsInBase + (currentDeltaBlock*numCeps) ;
				if (debug || (verboseLevel>2)) cout << " blkBase["<<blockBase<<"] cepsInBase["<<cepFeatsInBase<<"] curD	Blk["<<currentDeltaBlock<<"]" << endl;

				// append features to new (big) featureVector (at position i*N (or N + i*N) )
				//TODO replace loop by fct call - same new fct as above
				for( unsigned int i=0; i<deltaBlockFeature.getVectSize(); i++) {
					shiftedDeltaFeature[blockBase+i] = deltaBlockFeature[i] ;	//TODO optimize with memcopy
				}
				if (debug || (verboseLevel>2)) { cout << "  ShDltFeat[" ; for(unsigned int i=0; i<shiftedDeltaFeature.getVectSize();i++){ cout<<shiftedDeltaFeature[i]<<" " ;} ; cout <<"]"<< endl; }
			}

			// save (big) featureVector to file
			if (debug || (verboseLevel>2)){ cout << "WRITE: t["<<currentTime<<"] ShDltFeat[" ; for(unsigned int i=0; i<shiftedDeltaFeature.getVectSize();i++){ cout<<shiftedDeltaFeature[i]<<" " ;} ; cout <<"]"<< endl; }
			w.writeFeature( shiftedDeltaFeature) ;
		}
		///not using://outputFeatureFile(config,fsOut,0,fsOut.getFeatureCount(), w);
		if (debug || (verboseLevel>2)){ cout <<"SDC, number of frames= ["<<fsOut.getFeatureCount()<<"]"<<endl; }

	}// End feature file loop

  }// end try
  
  catch (Exception & e)
      {
	  cout << e.toString ().c_str () << endl;
      }
  return 0;
}


#endif // !defined(ALIZE_ShiftedDeltaFeat_cpp)
  

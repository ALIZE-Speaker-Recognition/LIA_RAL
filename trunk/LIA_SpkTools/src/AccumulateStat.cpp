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

#if !defined(ALIZE_GeneralTools_cpp)
#define ALIZE_GeneralTools_cpp

#include <iostream>
#include <fstream>  // pour outFile
#include <cstdio>   // pour printf()
#include <cassert> // pour le debug pratique
#include <cmath>
#include <liatools.h>
#if defined(THREAD)
#include <pthread.h>
#endif

//-------------------------------------------------------------------------
//-- Accumulate the log likelihood for the selected frames and a given model
void accumulateStatLLK(StatServer &ss,FeatureServer &fs,MixtureStat &llkAcc, unsigned long idxBeginFrame,unsigned long nbFrames,
		       Config &config)  
{
  fs.seekFeature(idxBeginFrame);                                           // go to the frame in the buffer (and load it if needed)
  for (unsigned long n=0;n<nbFrames;n++){
    Feature f;
    if(fs.readFeature(f)==false) cout<<"No more features"<<endl;
    
    double toto=llkAcc.computeAndAccumulateLLK(f);
    if (debug) cout << "likelihood Frame["<<idxBeginFrame+n<<"]="<<toto<<endl;
  }
}
// one a Segment
void accumulateStatLLK(StatServer &ss,FeatureServer &fs,MixtureStat &llkAcc,Seg* seg,Config &config)
{
  unsigned long begin=seg->begin()+fs.getFirstFeatureIndexOfASource(seg->sourceName()); // Find the index of the first frame of the file in the buffer
  accumulateStatLLK(ss,fs,llkAcc,begin,seg->length(),config);
}
// One on Cluster
void accumulateStatLLK(StatServer &ss,FeatureServer &fs,MixtureStat &llkAcc,SegCluster &selectedSegments,Config &config)
{
  Seg* seg;                                                     // reset the reader at the begin of the input stream
  selectedSegments.rewind();      
  while((seg=selectedSegments.getSeg())!=NULL)                  // For each of the selected segments
    accumulateStatLLK(ss,fs,llkAcc,seg,config);
}
 
 

//-------------------------------------------------------------------------
//-- accumulate the statistic for EM, using a current accumulator (wordl)
//-- CAUTION: THE ACCUMULATOR SHOULD BE INITIALIZED (resetEM) BEFORE THE CALL
//--          A GET CALL SHOULD BE DONE AFTER THE CALL  
// One a part of the feature stream
double accumulateStatEM(StatServer &ss,FeatureServer &fs,MixtureStat &emAcc,unsigned long idxBeginFrame,unsigned long nbFrames,Config &config){
  double llkAcc=0.0;
  fs.seekFeature(idxBeginFrame);                                           // go to the frame in the buffer (and load it if needed)
  for (unsigned long n=0;n<nbFrames;n++){
    Feature f;
    fs.readFeature(f);
    llkAcc+=log(emAcc.computeAndAccumulateEM(f));
    //for (unsigned long i=0;i<fs.getVectSize();i++)
    //  cout <<"f["<<i<<"]="<<f[i]<<endl;
  } 
  return llkAcc;
}
// one a Segment
double accumulateStatEM(StatServer &ss,FeatureServer &fs,MixtureStat &emAcc,Seg* seg,Config &config){
  unsigned long begin=seg->begin()+fs.getFirstFeatureIndexOfASource(seg->sourceName());              // Find the index of the first frame of the file in the buffer
  return accumulateStatEM(ss,fs,emAcc,begin,seg->length(),config);
}
// One on Cluster
double accumulateStatEMUnThreaded(StatServer &ss,FeatureServer &fs,MixtureStat &emAcc,SegCluster &selectedSegments,Config &config){
  double llkAcc=0.0;
  Seg* seg;                                                     // reset the reader at the begin of the input stream
  selectedSegments.rewind();      
  while((seg=selectedSegments.getSeg())!=NULL)                  // For each of the selected segments
    llkAcc+=accumulateStatEM(ss,fs,emAcc,seg,config);
  return llkAcc;
}

// Choose if its threaded or not here
double accumulateStatEM(StatServer &ss,FeatureServer &fs,MixtureStat &emAcc,SegCluster &selectedSegments,Config &config){  
    double llkAcc=0.0;
    #if defined(THREAD)          
    if (config.existsParam("numThread") && config.getParam("numThread").toLong() >0) llkAcc=accumulateStatEMThreaded(ss,fs,emAcc,selectedSegments,config);                        // Compute EM statistics
    else llkAcc=accumulateStatEMUnThreaded(ss,fs,emAcc,selectedSegments,config);
    #else
   llkAcc= accumulateStatEMUnThreaded(ss,fs,emAcc,selectedSegments,config);            // Compute EM statistics
    #endif
  return llkAcc;
  }

  // With weight on a cluster
double accumulateStatEM(StatServer &ss,FeatureServer &fs,MixtureStat &emAcc,SegCluster &selectedSegments,double & weight, Config &config){  
    double llkAcc=0.0;
    char cW[200];
    sprintf(cW,"%f",weight);
	String sW(cW);
    config.setParam("weightedEM","weight");
    if (verbose) cout << "Weighted EM with weight" << sW << endl;
    llkAcc= accumulateStatEM(ss,fs,emAcc,selectedSegments,config);
    return (weight*llkAcc);
  }
// *****************************************************************************************************
// ****************************** Threaded Version of Accumulate StatEM *********************************
// **************************** Nicolas Scheffer, 16/02/2007 **********************************************
#if defined(THREAD)
pthread_mutex_t mutexsum; // lock variable
bool stop=false; // flag variable once end of SegCluster is reached
//******************** Data strucutre of thread **************************
struct EMthread_data{
        SegCluster *selectedSegments;
        MixtureStat *emAcc;
        FeatureServer *fs;
        Config *config;
        double *llkAcc;
        RefVector <Feature> *featThreadBuff;
        unsigned long nThread;
};
// *********************** Routine **************************************
static void *EMthread(void *threadarg) {
	struct EMthread_data *my_data;
	my_data = (struct EMthread_data *) threadarg;
        SegCluster &selectedSegments=*(my_data->selectedSegments);
        MixtureStat &emAcc=*(my_data->emAcc);
        FeatureServer &fs=*(my_data->fs);
        Config &config=*(my_data->config);
        unsigned long nT=my_data->nThread;
        double weight=1.0;
        if (config.existsParam("weightedEM")) weight=config.getParam("weightedEM").toDouble();
        // **************** Main loop
        Seg* seg;
        unsigned long  cnt=0;
        while(!stop)  {                // For each of the selected segments
             // ***************** Reading features is locked
            RefVector <Feature> featThreadBuff;          
            pthread_mutex_lock (&mutexsum);          
            if (stop) {pthread_mutex_unlock (&mutexsum);break;}                    
            cnt++;
            if ((seg=selectedSegments.getSeg())==NULL) {
              pthread_mutex_unlock (&mutexsum);
              if (verboseLevel >1) cout<<"Thread:"<<nT << " broke" << endl;
              stop=true;
              break;
            }
            unsigned long begin=seg->begin()+fs.getFirstFeatureIndexOfASource(seg->sourceName());
            fs.seekFeature(begin);             
            for (unsigned long j=0;j<seg->length();j++){            
              featThreadBuff.addObject(*(new Feature),j);
              fs.readFeature(featThreadBuff[j]);
            }
            pthread_mutex_unlock (&mutexsum);
             // ***************** End lock
            
            // ***************** Accumulate EM on the frameTab
            for (unsigned long j=0;j<seg->length();j++)
              (*(my_data->llkAcc))+=log(emAcc.computeAndAccumulateEM(featThreadBuff[j],weight));
            featThreadBuff.deleteAllObjects();                    
        }
        if (verboseLevel > 2) cout << "(AccumulateStatEM) Number of segments treated by thread ["<<nT<<"]="<<cnt<<endl;
	pthread_exit((void*) 0);
	return(void*)0;
}

// Split selected cluster in nSplit clusters
void splitSegCluster(SegCluster & selectedSegments,unsigned long nSplit,RefVector <SegCluster> &vSelected) {
  Seg* seg;
  unsigned long t=0;
  unsigned long offset=(totalFrame(selectedSegments))/nSplit;
  if (verbose) cout << "(AccumulateStatEM) Splitting cluster: Total["<<totalFrame(selectedSegments)<<"] "<<offset<<" frames/cluster"<<endl;
  unsigned long limit=offset-1;
  selectedSegments.rewind();   
  unsigned long cnt=0;
  while((seg=selectedSegments.getSeg())!=NULL) { 
    cnt+=seg->length();
    if (cnt>=limit) {
        if (verbose) cout<<"SegCluster["<<t<<"] full with "<<totalFrame(vSelected[t])<<" frames"<<endl;
        if (t!=nSplit-1) t++;cnt=0;
      }    
    vSelected[t].addCopy(*seg);
  }
}

//***************************Threaded version ***************************
double accumulateStatEMThreaded(StatServer &ss,FeatureServer &fs,MixtureStat &emAcc,SegCluster &selectedSegments,Config &config){
  unsigned long NUM_THREADS=1;
  if(config.existsParam("numThread"))	NUM_THREADS = config.getParam("numThread").toLong();

  stop=false;
  if (verbose) cout << "(AccumulateStatEM) Threaded version of EM with "<<NUM_THREADS<<" threads"<<endl;
  double llkAcc=0.0;
  SegServer segServer;  
  //RefVector <SegCluster> vSelected;
  RefVector <MixtureStat> vEmAcc;    
  RealVector <double> vLlkAcc;  
  vLlkAcc.setSize(NUM_THREADS);
  vLlkAcc.setAllValues(0.0);
  selectedSegments.rewind();
  for(unsigned long t=0; t<NUM_THREADS; t++){
    /*vSelected.addObject(segServer.createCluster(selectedSegments.labelCode(),"",""),t);  
    vSelected[t].rewind();*/
    vEmAcc.addObject(ss.createAndStoreMixtureStat(emAcc.getMixture()),t);
    vEmAcc[t].resetEM();    
  }
  //splitSegCluster(selectedSegments,NUM_THREADS,vSelected);  
//  struct EMthread_data thread_data_array[NUM_THREADS];
	struct EMthread_data *thread_data_array = new EMthread_data[NUM_THREADS];
//  pthread_t threads[NUM_THREADS];	
	pthread_t *threads = new pthread_t[NUM_THREADS];


  pthread_attr_t attr;
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  pthread_mutex_init(&mutexsum, NULL);
  int rc,status;
  for(unsigned long t=0; t<NUM_THREADS; t++){
    thread_data_array[t].selectedSegments=&selectedSegments;
    thread_data_array[t].emAcc=&vEmAcc[t];
    thread_data_array[t].fs=&fs;   
    thread_data_array[t].config=&config;
    thread_data_array[t].llkAcc=&vLlkAcc[t];
    thread_data_array[t].nThread=t;	
    if (verbose) cout<<"(AccumulateStatEM) Creating thread n["<< t<< "]"<<endl;    
    rc = pthread_create(&threads[t], &attr, EMthread, (void *)&thread_data_array[t]);		
    if (rc) throw Exception("ERROR; return code from pthread_create() is ",__FILE__,rc);
    }
    if (verbose) cout<<"(AccumulateStatEM) Computing on thread"<<endl;    
    pthread_attr_destroy(&attr);
    for(unsigned long t=0; t<NUM_THREADS; t++) {
      rc = pthread_join(threads[t], (void **)&status);
      if (rc)  throw Exception("ERROR; return code from pthread_join() is ",__FILE__,rc);
      if (verbose) cout <<"(AccumulateStatEM) Completed join with thread ["<<t<<"] status["<<status<<"]"<<endl;
    }
    if (verbose) cout <<"(AccumulateStatEM) Fuse EM Accs "<<endl;
      unsigned long total=0;
    for(unsigned long t=0; t<NUM_THREADS; t++) {
      if (verbose) cout << "Number of frames treated by thread ["<<t<<"]="<<vEmAcc[t].getEMFeatureCount()<<endl;
      total+=vEmAcc[t].getEMFeatureCount();
      emAcc.addAccEM(vEmAcc[t]);
      ss.deleteMixtureStat(vEmAcc[t]);
      llkAcc+=vLlkAcc[t];
    }
    cout << "Total Number of frames in threads: "<<total<<endl;
    pthread_mutex_destroy(&mutexsum);
	free(thread_data_array);
	free(threads);

return llkAcc;
}
#endif
// **************************** End ********************************

/*Alex Preti things
double accumulateStatEM(StatServer & ss, FeatureServer & fs,
  MixtureStat & emAcc, unsigned long idxBeginFrame, unsigned long nbFrames,
  double &weight, Config & config)
{
  double llkAcc = 0.0;
  fs.seekFeature(idxBeginFrame);	// go to the frame in the buffer (and load it if needed)
  for (unsigned long n = 0; n < nbFrames; n++)
    {
      Feature f;
      fs.readFeature(f);
      llkAcc += log(emAcc.computeAndAccumulateEM(f, weight));
    }
  return (llkAcc * weight);
}

// one a Segment
double accumulateStatEM(StatServer & ss, FeatureServer & fs,
  MixtureStat & emAcc, Seg * seg, double &weight, Config & config)
{
  unsigned long begin = seg->begin() + fs.getFirstFeatureIndexOfASource(seg->sourceName());	// Find the index of the first frame of the file in the buffer
  return accumulateStatEM(ss, fs, emAcc, begin, seg->length(), weight,
    config);
}

// One on Cluster
double accumulateStatEM(StatServer & ss, FeatureServer & fs,
  MixtureStat & emAcc, SegCluster & selectedSegments, double &weight,
  Config & config)
{
  double llkAcc = 0.0;
  Seg *seg;			// reset the reader at the begin of the input stream
  selectedSegments.rewind();
  while ((seg = selectedSegments.getSeg()) != NULL)	// For each of the selected segments
    llkAcc += accumulateStatEM(ss, fs, emAcc, seg, weight, config);
  return llkAcc;
}*/


//-------------------------------------------------------------------------
//-- Accumulate the log likelihood for the selected frames and a given model, support weighted Feature Server (A.P)
void accumulateStatLLK(StatServer & ss, FeatureServer & fs,
  MixtureStat & llkAcc, unsigned long idxBeginFrame, unsigned long nbFrames,
  double &weight, Config & config)
{
  fs.seekFeature(idxBeginFrame);	// go to the frame in the buffer (and load it if needed)
  for (unsigned long n = 0; n < nbFrames; n++)
    {
      Feature f;
      if (fs.readFeature(f) == false)
	cout << "No more features" << endl;

      double toto = llkAcc.computeAndAccumulateLLK(f, weight);
      if (debug)
	cout << "likelihood Frame[" << idxBeginFrame +
	  n << "]=" << toto << endl;
    }
}

// one a Segment, support weighted Feature Server (A.P)
void accumulateStatLLK(StatServer & ss, FeatureServer & fs,
  MixtureStat & llkAcc, Seg * seg, double &weight, Config & config)
{
  unsigned long begin = seg->begin() + fs.getFirstFeatureIndexOfASource(seg->sourceName());	// Find the index of the first frame of the file in the buffer
  accumulateStatLLK(ss, fs, llkAcc, begin, seg->length(), weight, config);
}

// One on Cluster, support weighted Feature Server (A.P)
void accumulateStatLLK(StatServer & ss, FeatureServer & fs,
  MixtureStat & llkAcc, SegCluster & selectedSegments, double &weight,
  Config & config)
{
  Seg *seg;			// reset the reader at the begin of the input stream
  selectedSegments.rewind();
  while ((seg = selectedSegments.getSeg()) != NULL)	// For each of the selected segments
    accumulateStatLLK(ss, fs, llkAcc, seg, weight, config);
}


//-------------------------------------------------------------------------
//-- accumulate the statistic on the frames (mean and cov), using a current 
//-- FrameAcc
//--
//-- CAUTION: A COMPUTE_ALL AND A GET CALL SHOULD BE DONE AFTER THE CALLS  
void accumulateStatFrame(FrameAcc & frameAcc,FeatureServer &fs,
		    unsigned long idxBeginFrame,unsigned long nbFrames,Config &config)
{ 	
  fs.seekFeature(idxBeginFrame);                                           // go to the frame in the buffer (and load it if needed)
  Feature f;
  for (unsigned long n=0;n<nbFrames;n++){
    fs.readFeature(f);
    frameAcc.accumulate(f);
  }
}
// one a Segment
void accumulateStatFrame(FrameAcc & frameAcc,FeatureServer &fs,Seg* seg,Config &config)
{
  unsigned long begin=seg->begin()+fs.getFirstFeatureIndexOfASource(seg->sourceName());              // Find the index of the first frame of the file in the buffer
  accumulateStatFrame(frameAcc,fs,begin,seg->length(),config);
}
// One on Cluster
void accumulateStatFrame(FrameAcc &frameAcc,FeatureServer &fs,SegCluster &selectedSegments,Config &config)
{
  Seg* seg;                                                     // reset the reader at the begin of the input stream
  selectedSegments.rewind();      
  while((seg=selectedSegments.getSeg())!=NULL)                  // For each of the selected segments
    accumulateStatFrame(frameAcc,fs,seg,config);
}


//-------------------------------------------------------------------------
//-- accumulate the statistic on the frames raw distribution of each coefficient
//-- CAUTION: 
//         *THE ACCUMULATOR SHOULD BE INITIALIZED BEFORE THE FIRST CALL
//          initHistoTab()
//         *THE HISTO SHOULD BE COMPUTED BEFORE TO USE THE STAT
//          computeHistoTab()
//         *The histoTab should be freezen after use
//          freezeHistoTab();
//      
// Init the Histo array (one by coeff)
double areaHisto(const Histo & histo,unsigned long bin)
{
  return  histo.count(bin)*(histo.higherBound(bin)-histo.lowerBound(bin));
}
double areaHisto(const Histo & histo,unsigned long bin, double nonObserved)
{
  return areaHisto(histo,bin)*(1-nonObserved) ;
}
double linearInterpolation(double val,double lower,double higher){
  const double EPS=0.000000000000000000001;
  double interval=higher-lower;
  if (interval<EPS) return 1;
  return (val-lower)/interval;
} 
void freezeHistoTab(Histo* &histoT)
{
  delete []histoT;
}
void initHistoTab(Histo* &histoT,unsigned long size, unsigned long nbBins)
{
  Histo tmp(nbBins);
  histoT=new Histo[size];
  for (unsigned long i=0;i<size;i++) histoT[i]=tmp;
}
void computeHistoTab(Histo* histoT,unsigned long size)
{
  for (unsigned long i=0;i<size;i++) histoT[i].computeHisto();
}

void accumulateHistoFrame(Histo  *histoT,FeatureServer &fs,
		    unsigned long idxBeginFrame,unsigned long nbFrames,Config &config)
{ 	
  fs.seekFeature(idxBeginFrame);                                           // go to the frame in the buffer (and load it if needed)
  Feature f;
  for (unsigned long n=0;n<nbFrames;n++){
    fs.readFeature(f);
    for (unsigned long c=0;c<fs.getVectSize();c++)
      histoT[c].accumulateValue(f[c]);
  }
}
// one a Segment
void accumulateHistoFrame(Histo *histoT,FeatureServer &fs,Seg* seg,Config &config)
{
  unsigned long begin=seg->begin()+fs.getFirstFeatureIndexOfASource(seg->sourceName());// Find the index of the first frame of the file in the buffer
  accumulateHistoFrame(histoT,fs,begin,seg->length(),config);
}
// One on Cluster
void accumulateHistoFrame(Histo *histoT,FeatureServer &fs,SegCluster &selectedSegments,Config &config)
{
  Seg* seg;                                                     // reset the reader at the begin of the input stream
  selectedSegments.rewind();      
  while((seg=selectedSegments.getSeg())!=NULL)                  // For each of the selected segments
    accumulateHistoFrame(histoT,fs,seg,config);
}


#endif //!defined(ALIZE_AccumulateStat_cpp)



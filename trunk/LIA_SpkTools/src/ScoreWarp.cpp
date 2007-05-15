#if !defined(ALIZE_ScoreWarp_cpp)
#define ALIZE_ScoreWarp_cpp
 
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cassert>
#include <cmath>
#include "TrainTools.h"
#include "ScoreWarp.h"
const double PI2 = 3.14159265358979323846*2;
static double x1;
static double x2;
void boxMullerGeneratorInit(){
 x1= (rand()/ (float)RAND_MAX);
}
double boxMullerGenerator(double mean, double cov){
  x2=x1;
  x1= (rand()/ (float)RAND_MAX);
  if (debug) cout <<"boxMullerGenerator x1["<<x1<<"] x2["<<x2<<"]"<<endl;
  double y= sqrt(-2.0*log(x1))*cos(PI2*x2);   // Box-Muller gaussian generator from 2 sample of [0,1] indep. random numbers
  double ret= (y*cov)+mean;                        // fit to the N(mean,cov) distribution
  if (debug) cout <<"boxMullerGenerator["<<ret<<"]"<<endl;
  return ret;
}


// Build the Gaussian target distribution (mean,cov), by generating nbSample data, specifying the number of bins
Histo makeGausHisto(unsigned long nbSample,double mean, double cov,unsigned long nbBins){
  if (verbose) cout << "makeGaussHisto, nbSample["<<nbSample<<"] mean["<<mean<<"] cov["<<cov<<"] nbBin["<<nbBins<<"]"<<endl;
  Histo histo(nbBins);
  boxMullerGeneratorInit();       // Init The box-muller gaussian number generator
  for (unsigned long i=0;i<nbSample;i++)
    histo.accumulateValue(boxMullerGenerator(mean,cov));
  histo.computeHisto(); 
  if (verbose){
    double tot=0;
    for (unsigned long i=0;i<histo.size();i++) 
      tot+=histo.count(i)*(histo.higherBound(i)-histo.lowerBound(i));
    cout <<"makeGaussHisto tot of the final pdf["<<tot<<"]"<<endl;
  }
  return histo;
}


double centralSpace(const Histo &warpH,double a){
  if (a==0) return 0;
  unsigned long inf,sup;
  double t=(1.0-a)/2.0;
  inf=0;
  for (double tot=0.0;(inf<warpH.size()) && (tot<t);inf++)
    tot+=areaHisto(warpH,inf);
  sup=warpH.size()-1;
  for (double tot=1.0;(sup>=0) && (tot>1-t);sup--)
    tot-=areaHisto(warpH,sup);
  return (warpH.lowerBound(sup)-warpH.higherBound(inf));
}

// Compute the warped score, using the raw score distribution warH and the destination distribution destH
double scoreWarping(double score, const Histo& warpH, const Histo& destH, double nonObserved,double refArea)
{
  if (debug) cout << "scoreWarping score ["<<score<<"] nonObserved["<<nonObserved<<"] refArea["<<refArea<<"] ";
  // Value before the min, after the max
  if (score<warpH.lowerBound(0)){                // the value is less than the minimum of the raw score distribution 
    double infBound=warpH.lowerBound(0)-centralSpace(warpH,refArea);
    return destH.lowerBound(0)-(linearInterpolation(score,infBound,warpH.lowerBound(0))*centralSpace(destH,refArea));
  } 
  if (score>warpH.higherBound(warpH.size()-1)){  // the value is more than the maximum of the raw distribution
    double supBound=warpH.higherBound(warpH.size()-1)+centralSpace(warpH,refArea);
    return destH.higherBound(destH.size()-1)+(linearInterpolation(score,warpH.higherBound(warpH.size()-1),supBound)*centralSpace(destH,refArea)); 
  }
  
  double totalWarp=0.0;
  // Compute the area between -infinite and the score (raw distrib) - totalWarp
  unsigned long idxW=0;                          // Will be the bin number where is the score (in the raw distrib)
  for(;(idxW<warpH.size())&&(warpH.higherBound(idxW)<score); idxW++);
  if (debug) cout << " idxW["<<idxW<<"]";
  totalWarp=nonObserved/2;                       // nonObserved is the estimated amount of data non seen in the set
  for (unsigned long idx=0;idx<idxW;idx++)       // Accumulate the area for the raw distrib (numerical integration)
    totalWarp+=areaHisto(warpH,idx);             // Without the last bin
  double percentW=linearInterpolation(score,warpH.lowerBound(idxW),warpH.higherBound(idxW));
  if (debug) cout << " percentW["<<percentW<<"]";
  totalWarp+=areaHisto(warpH,idxW)*percentW;    // Add a percentage of the last bin (linear interpolation)
  if (debug) cout << " totalW["<<totalWarp<<"]";

  // Find idxD, the index of the corresponding bin (area bin[-infinite] until bin[idxD]=totalWarp
  double totalDest;
  unsigned long idxD; 
  for (idxD=0,totalDest=0;(idxD<destH.size())&&(totalDest<totalWarp);idxD++)
    totalDest+=areaHisto(destH,idxD);
  if (idxD==destH.size()){
    idxD--;
    totalDest-=areaHisto(destH,idxD);
  }
  else
    if (idxD>0){
      totalDest-=areaHisto(destH,idxD);
      idxD--;
    }    
  // Compute the final score
  double ret=destH.lowerBound(idxD);                                                       // Set the ret to the lowerbound of the good bin
  double percentH=linearInterpolation(totalWarp,totalDest,areaHisto(destH,idxD)+totalDest);// Compute the linear interpolation % for the last bin
  if (debug) cout << " idxD["<<idxD<<"] TotalDet["<<totalDest<<"] PercentD["<<percentH<<"]";
  ret+=(destH.higherBound(idxD)-destH.lowerBound(idxD))*percentH;                          // apply the interpolation
  if (debug) cout <<" final["<<ret<<"]"<<endl;
  return ret;
}
 

// Compute the warped value using score, data histo warpH,
// target histo destH
// it assumes that all the values are in warpH bounds

double warping(double score, const Histo& warpH, const Histo& destH)
{
  if (debug) cout << "warping input ["<<score<<"] ";
  
  double totalWarp=0.0;
  // Compute the area between -infinite and the score (raw distrib) - totalWarp
  unsigned long idxW=0;                          // Will be the bin number where is the score (in the raw distrib)
  for(;(idxW<warpH.size())&&(warpH.higherBound(idxW)<score); idxW++);
  if (debug) cout << " idxW["<<idxW<<"] ";

  for (unsigned long idx=0;idx<idxW;idx++)       // Accumulate the area for the raw distrib (numerical integration)
    totalWarp+=areaHisto(warpH,idx);             // Without the last bin
  double percentW=linearInterpolation(score,warpH.lowerBound(idxW),warpH.higherBound(idxW));
  if (debug) cout << " percentW["<<percentW<<"] ";
  totalWarp+=areaHisto(warpH,idxW)*percentW;    // Add a percentage of the last bin (linear interpolation)
  if (debug) cout << " totalW["<<totalWarp<<"] ";

  // Find idxD, the index of the corresponding bin (area bin[-infinite] until bin[idxD]=totalWarp
  double totalDest;
  unsigned long idxD; 
  for (idxD=0,totalDest=0;(idxD<destH.size())&&(totalDest<totalWarp);idxD++)
    totalDest+=areaHisto(destH,idxD);
  if (idxD==destH.size()){
    idxD--;
    totalDest-=areaHisto(destH,idxD);
  }
  else
    if (idxD>0){
      totalDest-=areaHisto(destH,idxD);
      idxD--;
    }    
  // Compute the final value
  double ret=destH.lowerBound(idxD);                                                       // Set the ret to the lowerbound of the good bin
  double percentH=linearInterpolation(totalWarp,totalDest,areaHisto(destH,idxD)+totalDest);// Compute the linear interpolation % for the last bin
  if (debug) cout << " idxD["<<idxD<<"] TotalDet["<<totalDest<<"] Percent last bin["<<percentH<<"] ";
  ret+=(destH.higherBound(idxD)-destH.lowerBound(idxD))*percentH;                          // apply the interpolation
  if (debug) cout <<" warping output["<<ret<<"] "<<endl;
  return ret;
}

#endif // !defined(ALIZE_ScoreWarp_cpp)

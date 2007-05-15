#if !defined(ALIZE_ScoreWarp_h)
#define ALIZE_ScoreWarp_h
#include "alize.h"


using namespace alize;
using namespace std;

// BoxMuller Gaussian random generator
void boxMullerGeneratorInit();
double boxMullerGenerator(double mean, double cov);

// Build the Gaussian target distribution (mean,cov), by generating nbSample data, specifying the number of bins
Histo makeGausHisto(unsigned long nbSample,double mean, double cov,unsigned long nbBins);
// Compute the warped score, using the raw score distribution warH and the destination distribution destH
double scoreWarping(double, const Histo& , const Histo& , double=0.1,double=0.5);
double warping(double score, const Histo& warpH, const Histo& destH);

#endif // !defined(Hist)

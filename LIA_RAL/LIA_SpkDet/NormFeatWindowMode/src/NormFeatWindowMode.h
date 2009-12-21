#if !defined(ALIZE_NormFeatWindowMode_h)
#define ALIZE_NormFeatWindowMode_h

#include "alize.h"
/*!
NormFeatWindowMode is dedicated to the normalisation of features often used in 
Automatic Speech/Speaker Recognition.
It can work in segmental mode to compute mean, cov and normalisation
only on interesting segments.
*/

int normFeatOnlineMode(alize::Config &);
#endif // !defined(NormFeatWindowMode)

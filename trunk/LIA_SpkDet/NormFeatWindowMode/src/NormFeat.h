#if !defined(ALIZE_NormFeat_h)
#define ALIZE_NormFeat_h

#include "alize.h"
/*!
NormFeat is dedicated to the normalisation of features often used in 
Automatic Speech/Speaker Recognition.
It can work in segmental mode to compute mean, cov and normalisation
only on interesting segments.
*/

int normFeatOnlineMode(alize::Config &);
#endif // !defined(NormFeat)

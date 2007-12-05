#if !defined(ALIZE_GMMTokenizer_h)
#define ALIZE_GMMTokenizer_h

#include "alize.h"
extern bool debug;
extern bool verbose;

        int GMMTokenizer(alize::Config&);
        int GaussianConfusionMatrix(alize::Config&);
#endif 

#if !defined(ALIZE_ShiftedDeltaFeat_h)
#define ALIZE_ShiftedDeltaFeat_h

#include "alize.h"

/*!
ShiftedDeltaFeat transforms features (cepstra) into shifted delta features,
a feature type often used in Automatic Speech/Speaker/Language Recognition.
*/
int ShiftedDeltaFeat (alize::Config &);

#endif // !defined(ShiftedDeltaFeat)

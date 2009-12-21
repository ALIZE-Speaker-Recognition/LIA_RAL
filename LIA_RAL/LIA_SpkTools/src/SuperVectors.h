#if !defined(ALIZE_SuperVectors_h)
#define ALIZE_SuperVectors_h

#include <liatools.h>

using namespace alize;
using namespace std;

// Get a mixture and stack its mean parameters into a vector
void modelToSv(const MixtureGD &,RealVector<double> &);

// Get a vector and put values in mean parameters of a mixture
void svToModel(RealVector<double> &,MixtureGD &);

// Project vector X on subspace of U
void projectOnSubSpace(RealVector<double> &,Matrix <double> &,RealVector<double> &);

// Project supervector of M on complementary subspace of U (operator x'=(I-U'U)x)
void computeNap(MixtureGD &,Matrix <double> &);

// Compute supervector channel NAP effect
void computeNAPChannelEffect(MixtureGD &,MixtureGD &,Matrix<double>&);

void getFisherWeightVector(const MixtureGD&,const MixtureGD&, RealVector<double> &,Config&);

void getKLVector(MixtureGD&, RealVector<double> &,Config&); 

void getSuperVector(RealVector<double> &,MixtureGD &,MixtureGD &,Config &);

#endif

#if !defined(ALIZE_ReadModel_cpp)
#define ALIZE_ReadModel_cpp

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cassert>
#include <cmath>
#include "ReadModel.h"

using namespace alize;
using namespace std;

int readModel(Config &config)
{
	String inputFilename = config.getParam("inputModelFilename");
	try{
		MixtureServer ms(config);
		MixtureGD & model=ms.loadMixtureGD(inputFilename);
		cout << model.toString() << endl;
		for (unsigned long i=0;i< model.getDistribCount();i++) {
			cout << "Distribution " << i << endl;
			cout << model.getDistrib(i).toString() << endl;
		}
	}	
	catch (Exception& e){cout << e.toString().c_str() << endl;}
return 0;
}

int outputWeightVector(Config &config)
{
	String inputFilename = config.getParam("inputModelFilename");
	try{
		MixtureServer ms(config);
		MixtureGD & model=ms.loadMixtureGD(inputFilename);
		String format = config.getParam("outputWeightVectorFormat");
		if (format=="LIA") {
			for (unsigned int i=0;i< model.getDistribCount();i++)
				cout<<model.weight(i)<<" ";
		}
		if (format=="SVMLight") {
			cout << "1 "; // target type
			for (unsigned int i=0;i< model.getDistribCount();i++)
				cout<<i+1<<":"<<model.weight(i)<<" ";
		}
	cout << endl;
	}	
	catch (Exception& e){cout << e.toString().c_str() << endl;}
return 0;
}

#endif //!defined(ALIZE_ReadModel_cpp)

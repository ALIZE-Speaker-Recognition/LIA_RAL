#if !defined(ALIZE_ExtractParams_cpp)
#define ALIZE_ExtractParams_cpp

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cassert>
#include <cmath>
#include "ExtractParams.h"

using namespace alize;
using namespace std;

int ExtractParams (Config & config)
{
	try
	{
		String featFileName=config.getParam("inputFile");
		FeatureServer fs(config,featFileName);
		FeatureFileWriter w(featFileName, config);
		Feature f;
		while(fs.readFeature(f,0))
			w.writeFeature (f); 
	}
	catch (Exception & e) {cout << e.toString ().c_str () << endl;}
	return 0;
}

#endif // !defined(ALIZE_ExtractParams_cpp)

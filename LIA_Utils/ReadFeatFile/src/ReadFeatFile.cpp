#if !defined(ALIZE_ReadFeatFile_cpp)
#define ALIZE_ReadFeatFile_cpp

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <liatools.h>
#include "ReadFeatFile.h"

using namespace alize;
using namespace std;

int readFeatFile (Config & config)
{
  try{
	String featFileName=config.getParam("inputFile");
	FeatureServer fs(config,featFileName);
	  
	cout << "***********************************" << endl;
	cout << "File Name: " << featFileName << endl;
	cout << "Features Number: " << fs.getFeatureCount() << endl;  
	cout << "FeatureSize: " << fs.getVectSize() << endl;
	cout << "FeatureFlags: " << fs.getFeatureFlags().getString() << endl;  
	cout << "Format: " << config.getParam("loadFeatureFileFormat") << endl;
	cout << "***********************************" << endl;  
	
	Feature f;
	while (fs.readFeature(f)) {
		for (unsigned long i=0;i<fs.getVectSize();i++)
			cout << f[i] << " ";	
		cout << endl;
	}
	  
  }    
  catch (Exception & e) {
      cout << e.toString ().c_str () << endl;
    }
  return 0;
}
#endif // !defined(ALIZE_ReadFeatFile_cpp)

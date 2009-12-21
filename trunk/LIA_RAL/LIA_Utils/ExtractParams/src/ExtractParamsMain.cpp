#include <iostream>

#include "ExtractParams.h"

int main (int argc, char *argv[])
{
  using namespace std;
  using namespace alize;

  try
  {
    Config config;
    CmdLine cmdLine (argc, argv);
    if (cmdLine.displayHelpRequired ())	// --help
      {
	cout << "ExtractParams" << endl;
	cout <<
	  "ExtractParams.exe --i <inputFile> --featureServerMask <0-m,m+1-N> --saveFeatureFileExtension .red.prm"
	  << endl;
      }
    else if (cmdLine.displayVersionRequired ())	// --version
      cout << "Version 0.0" << endl;
    else
      {
	cmdLine.copyIntoConfig (config);
	if (config.getParamCount() == 0)
		{cerr << "Error: Type ExtractParams.exe --help to get usage" << endl;exit(1);}
	config.setParam("featureServerBufferSize","1");
	config.setParam("saveFeatureFileFormat",config.getParam("loadFeatureFileFormat"));
	ExtractParams (config);
      }
  }
  catch (alize::Exception & e)
  {
    cout << e.toString () << endl;
  }

  return 0;
}

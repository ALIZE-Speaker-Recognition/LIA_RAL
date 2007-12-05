#include <iostream>
#include "FusionScore.h"
#include <liatools.h>


int main(int argc, char* argv[])
{
	using namespace std;
	using namespace alize;

	try
	{
		Config tmp;
		CmdLine cmdLine(argc, argv);
		if (cmdLine.displayHelpRequired()) // --help
		{
			cout << "************************************" << endl;
			cout << "********** FusionScore.exe **********" << endl;
			cout << "************************************" << endl;
			cout << endl;
			cout << "File Fusion" << endl;
			cout << endl;
			cout << "Command Line example: " << endl;
			cout << "FusionScore.exe --inputFileList foo.lst --outputFile foo2.nist --weights <coeffs of each file in a ascii file separate by a white space> --fusionMethod <ArithMean|GeoMean> --format <lia|nist|etf> --Nbest <1...N> (optional)" << endl;
		}
		else
		{
			cmdLine.copyIntoConfig(tmp);
			if (tmp.getParamCount() == 0)
			{cerr << "Error: Type FusionScore.exe --help to get usage" << endl;exit(1);} 
			Config config(tmp.getParam("config"));
			config.setParam(tmp);
			
		        if (config.existsParam("classFusion"))        // if the parameter is present, we work on one NiIST file only by merging class scores 
				FusionClassScore(config);
			else
				FusionScore(config);  
		}
	}
	catch (alize::Exception& e)
	{
		cout << e.toString() << endl;
	}

	return 0;
}

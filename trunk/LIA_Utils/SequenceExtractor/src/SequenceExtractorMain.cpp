#include <iostream>

#include "SequenceExtractor.h"


int main(int argc, char* argv[])
{
    using namespace std;
    using namespace alize;
    

    try
    {

    
      ConfigChecker cc;
      cc.addStringParam("config", false, true, "default config filename");
      cc.addIntegerParam("maxOrder",true,true,"maximum order of the input ngrams (from order 1 to maxOrder)");
      cc.addIntegerParam("maxNgram",true,true,"maximum number of ngram by file");
      cc.addIntegerParam("nbInputSymb",true,true,"number of input symbols in the ngrams. I.E. number of gaussian indexes");
      cc.addIntegerParam("nbOutputSymb",true,true,"number of required output symbols, decided by the program. couls be less");
      cc.addStringParam("ngramFilename",true,true,"Base filename for the ngram. N order filename have the form ngramFilenameN.ext");
      cc.addStringParam("ngramExt",true,true,"Ext for the ngram files, should include the .");
      cc.addStringParam("outputFilename",false,true,"The complete output filename for the output decoder tree");
      cc.addStringParam("outputInfoFilename",false,true,"The complete output filename for the info file");
      cc.addStringParam("equalInputInfo",false,true,"Try to find seequence representing the same amount of input symbols and not only same counts");
     
      CmdLine cmdLine(argc, argv);
      if (cmdLine.displayHelpRequired()){
	cout <<"SequenceExtravtor.exe"<<endl<<"This program uses in input ngram files and outputs the decoder tree"
	     <<" with a set of (variable length) symbol sequences with the *same* probabilty"<<endl<<cc.getParamList()<<endl;
	return 0;
      }
      if (cmdLine.displayVersionRequired())
	{
	  cout <<"Version 1.0"<<endl;
	}
      Config tmp;
      cmdLine.copyIntoConfig(tmp);
      Config config;
      if (tmp.existsParam("config")) config.load(tmp.getParam("config"));
      cmdLine.copyIntoConfig(config);
      cc.check(config);
       
     
      debug=config.getParam_debug();
      if (config.existsParam("verbose")) verbose=config.getParam("verbose").toBool();
      else verbose=false;
      sequenceExtractor(config);
	
     
    }
    catch (alize::Exception& e)
    {
        cout << e.toString() << endl;
    }

    return 0;
}

#include <iostream>
#include "LabelNGram.h"
#include "liatools.h"


int main(int argc, char* argv[])
{
    using namespace std;
    using namespace alize;
    ConfigChecker cc;
    cc.addStringParam("config", false, true, "default config filename");
    cc.addStringParam("inputFilename", false, true, "input symbol file");
    cc.addStringParam("NGramFilename", false, true, "bag of ngram codebook file ");
    cc.addIntegerParam("NGramOrder", false, true, "ngram order (default 3) ");
    cc.addIntegerParam("NGramSelected", false, true, "selected ngram in bag (default 16) ");
    cc.addStringParam("saveLabelFileExtension", false, true, "extension (default .sym.lbl)");
    cc.addStringParam("labelOutputPath", false, true, "path (default ./)");
    cc.addStringParam("symbolFileExtension", false, true, "extension (default .sym)");
    cc.addStringParam("symbolPath", false, true, "path (default ./)");
    cc.addStringParam("symbolFormat", false, true, "format (default ascii)");

  try
  {
    CmdLine cmdLine (argc, argv);
    if (cmdLine.displayHelpRequired ()){	// --help
      cout << "LabelNgram" << endl;
      cout << "Use to transform a stream of tokens in a label file knowing a bag of ngram"<< endl;
      
    
      cout<<cc.getParamList()<<endl;
      }
    else if (cmdLine.displayVersionRequired ())	// --version
      cout << "Version 1.0" << endl;
    else{
      // copy parameters from command line into config
      Config tmp;
      cmdLine.copyIntoConfig (tmp);
      Config config;
      if (tmp.existsParam("config")) config.load(tmp.getParam("config"));
      cmdLine.copyIntoConfig(config);
      cc.check(config);
      debug=config.getParam_debug();
      if (config.existsParam("verbose")) verbose=config.getParam("verbose").toBool();
      else verbose=false;
      if (verbose) verboseLevel=1;else verboseLevel=0;
      if (config.existsParam("verboseLevel"))verboseLevel=config.getParam("verboseLevel").toLong();
      if (verboseLevel>0) verbose=true;
      labelNGram(config);
    }
  }
  catch (alize::Exception & e)
    {
    cout << e.toString () << endl << cc.getParamList() << endl;
  }

  return 0;
}


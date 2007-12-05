#include <iostream>
 
#include "SequenceExtractor.h"


int main(int argc, char* argv[])
{    
    try    
    {    
      // String inputTreeFilename=config.getParam("decoderFilename");  // The complete filename for the decoder tree
      // String outputFilename=config.getParam("outputFilename");      // The complete output filename
      // String inputFilename=config.getParam("inputFilename");        // The complete input filename 
    
      ConfigChecker cc;
      cc.addStringParam("config", false, true, "default config filename");
      cc.addStringParam("decoderFilename",true,true,"filename of the decoder tree");
      cc.addStringParam("inputFilename",true,true,"filename of the input symbol stream");
      cc.addStringParam("outputFilename",true,true,"filename of the output symbol file");
      cc.addStringParam("labelFilename",false,true,"if segment processing is needed, name of the label file");
      cc.addStringParam("labelSelectedFrames",false,true,"if segment processsing is needed, label of the selected segments");        
      cc.addFloatParam("frameLength",false,true,"If segment processing is required, it gives the length in s of an input symbol");
      cc.addBooleanParam("overlapMode",false,true,"if set to true (default=false), find overlap sequences. i.e decode one sequance by input symbol");
      cc.addStringParam("ngramOutputFilename",false,true,"If set, computes the ngram on the output seq and save it in this file");
      cc.addIntegerParam("ngramOutputOrder",false,true,"If output ngram is needed, gives the max order of the ngram");
      CmdLine cmdLine(argc, argv);
      if (cmdLine.displayHelpRequired())
	{ 
	  cout <<"SequenceDecoder.exe"<<endl<<"This program uses in input a decoder tree and an input stream"
	       <<endl<<"It decodes the input stream and outputs the finded sequences"<<endl<<cc.getParamList()<<endl;
	  cout <<"HelpRequired"<<endl; 
	  return 0;  
	}
      if (cmdLine.displayVersionRequired())
	{
	  cout <<"Version Beta "<<endl;
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
      sequenceDecoder(config);
	
     
    }
    catch (alize::Exception& e)
    {
        cout << e.toString() << endl;
    }

    return 0;
}

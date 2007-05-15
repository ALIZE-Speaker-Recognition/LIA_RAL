#include <iostream>

#include "GmmTokenizer.h"
#include "liatools.h"

int main(int argc, char* argv[])
{
    using namespace std;
    using namespace alize;

    try
    {
    ConfigChecker cc;
    cc.addStringParam("config", false, true, "default config filename");

    CmdLine cmdLine(argc, argv);
        if (cmdLine.displayHelpRequired()) // --help
        {cout << "Symbol Generation by outputing best gaussians" << endl;
	cout << "usage: AnalyseModel.exe --config <foo.cfg> <--symbolsExtract|--confusionMatrix> opt:<duration> (to get repeated symbols or not)" << endl;}
        else if (cmdLine.displayVersionRequired()) // --version
            cout << "Version 0.0" << endl;
        else
        {
        Config tmp;
        cmdLine.copyIntoConfig(tmp);
    Config config;
    if (tmp.existsParam("config")) config.load(tmp.getParam("config"));
    cmdLine.copyIntoConfig(config);
    cc.check(config);
    debug=config.getParam_debug();
    if (config.existsParam("verbose"))verbose=config.getParam("verbose").toBool();else verbose=false;
    if (verbose) verboseLevel=1;else verboseLevel=0;
    if (config.existsParam("verboseLevel"))verboseLevel=config.getParam("verboseLevel").toLong();
    if (verboseLevel>0) verbose=true;
        // Two modes available Symboles output or confusion Matrix
        bool confusionMatrixMode=false;
    bool symbolsExtractMode=true;
    if (config.existsParam("confusionMatrix")) {confusionMatrixMode=true;symbolsExtractMode=false;}
        if (verbose){
                if (confusionMatrixMode) { cout << " Compute the confusion Matrix"<<endl; GaussianConfusionMatrix(config);}
                if (symbolsExtractMode) {cout << " GMM Tokenizer mode"<<endl; GMMTokenizer(config);}
        }
        if ((confusionMatrixMode || symbolsExtractMode ) == 0) {cerr << "No Mode found <confusionMatrix|symbolsExtract>" << endl ; exit(1);}
        }  
   }
    catch (alize::Exception& e)
    {
        cout << e.toString() << endl;
    }

    return 0;
}

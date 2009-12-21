#include <iostream>
#include "TrainTools.h"
#include "ScoreWarp.h"

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
	cout << "ScoreWarp.exe" << endl;
	cout <<  " "  << endl;   
      }
    else if (cmdLine.displayVersionRequired ())	// --version
      cout << "Version 0.0" << endl;
    else
      {
	cmdLine.copyIntoConfig (config);
	if (config.getParamCount() == 0)
		{cerr << "Error: Type ScoreWarp --help to get usage" << endl;exit(1);} 
	if (config.existsParam("debug")) debug=true;
	if (config.existsParam("verbose")) verbose=true;
	if (config.existsParam("generateGauss")){
	  unsigned long nbBin=50;
	  if (config.existsParam("nbBinGauss")) nbBin=config.getParam("nbBinGauss").toLong();
	  unsigned long nbSample=1000000;
	  if (config.existsParam("nbSample")) nbSample=config.getParam("nbSample").toLong();
	  double mean=0;
	  if (config.existsParam("mean")) mean=config.getParam("mean").toDouble();
	  double cov=1;
	  if (config.existsParam("cov")) cov=config.getParam("cov").toDouble();
	  String outputFile="foo";
	  if (config.existsParam("outputFilename")) outputFile=config.getParam("outputFilename");
	  Histo gaussH=makeGausHisto(nbSample,mean, cov,nbBin);
	  gaussH.saveGnuplot(outputFile+".gnuplot");
	  gaussH.save(outputFile+".histo");
	  }
	if (config.existsParam("testWarp")){
	  String outputFile="foo.out";
	  if (config.existsParam("outputFilename")) outputFile=config.getParam("outputFilename");
	  ofstream ofile(outputFile.c_str(),ios::out);
	  String fileDestH="destH";
	  if (config.existsParam("filenameDestH")) fileDestH=config.getParam("filenameDestH");
	  String fileRawH="rawH";
	  if (config.existsParam("filenameRawH")) fileRawH=config.getParam("filenameRawH");
	  unsigned long nbSample=1000000;
	  if (config.existsParam("nbSample")) nbSample=config.getParam("nbSample").toLong();
	  double mean=0;
	  if (config.existsParam("mean")) mean=config.getParam("mean").toDouble();
	  double cov=1;
	  if (config.existsParam("cov")) cov=config.getParam("cov").toDouble();
	  double nonObserved=0.1;
	  if (config.existsParam("nonObserved")) nonObserved=config.getParam("nonObserved").toDouble();
	  double area=0.5;
	  if (config.existsParam("refArea")) area=config.getParam("refArea").toDouble();	  
	  Histo gaussH;
	  Histo gaussR;
	  gaussH.load(fileDestH+".histo");
	  gaussR.load(fileRawH+".histo");
	  boxMullerGeneratorInit();
	  for (unsigned long i=0;i<nbSample;i++){
	    double val=boxMullerGenerator(mean, cov);
	    ofile << val <<" "<<scoreWarping(val,gaussR ,gaussH,nonObserved,area)<<endl;
	  }
	}
	if (config.existsParam("normalizeHisto")){
	  String outputFile="foo";
	  if (config.existsParam("outputFilename")) outputFile=config.getParam("outputFilename");
	  ofstream ofile(outputFile.c_str(),ios::out);
	  String fileInputH="";
	  if (config.existsParam("filenameInputH")) fileInputH=config.getParam("filenameInputH");  
	  Histo inputH;
	  inputH.load(fileInputH+".histo");   
	  double area=0;
	  for (unsigned long i=0;i<inputH.size();i++) area+=areaHisto(inputH,i);
	  if (debug) cout <<"normalizeHisto, area["<<area<<"]"<<endl;
	  inputH.div(area); 
	  inputH.save(outputFile+".histo"); 
	  inputH.saveGnuplot(outputFile+".gnuplot");
	}
	if (config.existsParam("testWarpLinear")){
	  String outputFile="foo.out";
	  if (config.existsParam("outputFilename")) outputFile=config.getParam("outputFilename");
	  ofstream ofile(outputFile.c_str(),ios::out);
	  String fileDestH="destH";
	  if (config.existsParam("filenameDestH")) fileDestH=config.getParam("filenameDestH");
	  String fileRawH="rawH";
	  if (config.existsParam("filenameRawH")) fileRawH=config.getParam("filenameRawH");
	  unsigned long nbSample=1000000;
	  if (config.existsParam("nbSample")) nbSample=config.getParam("nbSample").toLong();
	  double min=-10; 
	  if (config.existsParam("min")) min=config.getParam("min").toDouble();
	  double max=10;
	  if (config.existsParam("max")) max=config.getParam("max").toDouble();
	  double nonObserved=0.1;
	  if (config.existsParam("nonObserved")) nonObserved=config.getParam("nonObserved").toDouble();
	  double area=0.5;
	  if (config.existsParam("refArea")) area=config.getParam("refArea").toDouble();
	  Histo gaussH; 
	  Histo gaussR; 
	  gaussH.load(fileDestH+".histo");
	  gaussR.load(fileRawH+".histo");
	  for (unsigned long i=0;i<nbSample;i++){
	    double val=min+(i*((max-min)/(double)nbSample)); 
	    ofile << val <<" "<<scoreWarping(val,gaussR ,gaussH,nonObserved,area)<<endl;
	  }
	}
      }
  }
  catch (alize::Exception & e)
  {
    cout << e.toString () << endl;
  }

  return 0;
}

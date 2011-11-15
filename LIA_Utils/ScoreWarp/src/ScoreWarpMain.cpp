/*
This file is part of LIA_RAL which is a set of software based on ALIZE
toolkit for speaker recognition. ALIZE toolkit is required to use LIA_RAL.

LIA_RAL project is a development project was initiated by the computer
science laboratory of Avignon / France (Laboratoire Informatique d'Avignon -
LIA) [http://lia.univ-avignon.fr <http://lia.univ-avignon.fr/>]. Then it
was supported by two national projects of the French Research Ministry:
	- TECHNOLANGUE program [http://www.technolangue.net]
	- MISTRAL program [http://mistral.univ-avignon.fr]

LIA_RAL is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation, either version 3 of
the License, or any later version.

LIA_RAL is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with LIA_RAL.
If not, see [http://www.gnu.org/licenses/].

The LIA team as well as the LIA_RAL project team wants to highlight the
limits of voice authentication in a forensic context.
The "Person Authentification by Voice: A Need of Caution" paper
proposes a good overview of this point (cf. "Person
Authentification by Voice: A Need of Caution", Bonastre J.F.,
Bimbot F., Boe L.J., Campbell J.P., Douglas D.A., Magrin-
chagnolleau I., Eurospeech 2003, Genova].
The conclusion of the paper of the paper is proposed bellow:
[Currently, it is not possible to completely determine whether the
similarity between two recordings is due to the speaker or to other
factors, especially when: (a) the speaker does not cooperate, (b) there
is no control over recording equipment, (c) recording conditions are not
known, (d) one does not know whether the voice was disguised and, to a
lesser extent, (e) the linguistic content of the message is not
controlled. Caution and judgment must be exercised when applying speaker
recognition techniques, whether human or automatic, to account for these
uncontrolled factors. Under more constrained or calibrated situations,
or as an aid for investigative purposes, judicious application of these
techniques may be suitable, provided they are not considered as infallible.
At the present time, there is no scientific process that enables one to
uniquely characterize a persones voice or to identify with absolute
certainty an individual from his or her voice.]

Copyright (C) 2004-2010
Laboratoire d'informatique d'Avignon [http://lia.univ-avignon.fr]
LIA_RAL admin [alize@univ-avignon.fr]
Jean-Francois Bonastre [jean-francois.bonastre@univ-avignon.fr]
*/

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

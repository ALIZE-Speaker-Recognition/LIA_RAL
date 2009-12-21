#if !defined(ALIZE_Hist_cpp)
#define ALIZE_Hist_cpp

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cassert>
#include <cmath>
#include "Hist.h"

using namespace alize;
using namespace std;

int Hist (Config & config)
{
	try {
	String datalist=config.getParam("dataFile");
	unsigned long nbBins=config.getParam("nbBins").toLong();
	String format="gnuplot";
	if (config.existsParam("format")) format=config.getParam("format");
	XList data;
	data.load(datalist,config);
	XLine * linep;
	Histo histDATA(nbBins);

	while((linep=data.getLine())!=NULL) {
		histDATA.accumulateValue(linep->getElement(0).toDouble());
	}
	histDATA.computeHisto();
	if (format=="txt")
	  histDATA.save(config.getParam("outFile"));
	else
	  histDATA.saveGnuplot(config.getParam("outFile"));
	}			// fin try

	catch (Exception & e) {cout << e.toString ().c_str () << endl;}
return 0;
}

#endif // !defined(ALIZE_Hist_cpp)

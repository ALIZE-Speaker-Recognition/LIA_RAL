// ReadFeatFile.cpp
//
// This file is a part of Mistral Software MISTRAL_SpkDet, based on ALIZE library 
// MISTRAL_SpkDet  is a free, open tool for speaker recognition
// MISTRAL_SpkDet is a development project initiated and funded by the ANR (agence nationale 
// de la recherche).
// See http://mistral.univ-avignon.fr
// 
// ALIZE is needed for MISTRAL_SpkDet
//
// Copyright (C) 2004, 2005, 2006, 2007, 2008
//  Laboratoire d'informatique d'Avignon [www.lia.univ-avignon.fr]
//  Jean-Francois Bonastre [jean-francois.bonastre@lia.univ-avignon.fr]
//      
// MISTRAL_SpkDet is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
// You should have received a copy of the GNU General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// Contact Jean-Francois Bonastre (jean-francois.bonastre@lia.univ-avignon.fr) for
// more information about the licence or the use of LIA_SpkDet
// First version 15/07/2004
// This version 28/01/2008 - (edited by e.c / eric.charton@univ-avignon.fr)


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

/*!
	This metho read a feature file and display informations about it. You will get the Features numbers,
	size, flags and format. 

	@param config	give a config file for evaluation by config checker
	@return 		send back a 0 value in all cases
*/
int readFeatFile (Config & config)
{
  try{
	
	/** Get as param a feature file name 
	 *  with --inputFeatureFileName cfg option 
      */
	String featFileName=config.getParam("inputFeatureFileName"); 

	FeatureServer fs(config,featFileName);
	  
	//! Display informations about browsed file
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

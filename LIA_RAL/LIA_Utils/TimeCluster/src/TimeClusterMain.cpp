// TimeClusterMain.cpp
// This file is a part of LIA Softwares LIA_SpkDet and LIA_SpkSeg, based on ALIZE toolkit 
// LIA_SpkDet is a free, open tool for speaker recognition
// LIA_SpkSeg is a free, open tool for speaker segmentation
// This project is a development project initiated and funded by the LIA lab.
// See www.lia.univ-avignon.fr
// 
// ALIZE is needed for LIA_SpkDet and LIA_SpkSeg
// for more information about ALIZE, see http://www.lia.univ-avignon.fr/heberges/ALIZE/
// First version March 26 2005
//
// Copyright (C) 2005
//  Laboratoire d'informatique d'Avignon [www.lia.univ-avignon.fr]
// Main authors
//  Jean-Francois Bonastre [jean-francois.bonastre@lia.univ-avignon.fr]
//      
// LIA_SpkDet and LIA_SpkSeg are free software; you can redistribute it and/or
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
// The LIA team as well as the ALIZE project want to highlight the limits of voice authentication
// in a forensic context. 
// The following paper proposes a good overview of this point:
// [Bonastre J.F., Bimbot F., Boe L.J., Campbell J.P., Douglas D.A., Magrin-chagnolleau I.,
//  Person  Authentification by Voice: A Need of Caution,
//  Eurospeech 2003, Genova]
// The conclusion of the paper of the paper is proposed bellow:
// [Currently, it is not possible to completely determine whether the
//  similarity between two recordings is due to the speaker or to other
//  factors, especially when: (a) the speaker does not cooperate, (b) there
//  is no control over recording equipment, (c) recording conditions are not 
//  known, (d) one does not know whether the voice was disguised and, to a
//  lesser extent, (e) the linguistic content of the message is not
//  controlled. Caution and judgment must be exercised when applying speaker
//  recognition techniques, whether human or automatic, to account for these
//  uncontrolled factors. Under more constrained or calibrated situations,
//  or as an aid for investigative purposes, judicious application of these
//  techniques may be suitable, provided they are not considered as infallible.
//  At the present time, there is no scientific process that enables one to
//  uniquely characterize a person=92s voice or to identify with absolute
//  certainty an individual from his or her voice.]
//
// Contact Jean-Francois Bonastre (jean-francois.bonastre@lia.univ-avignon.fr) for
// more information about the licence or the use of LIA_SpkDet or LIA_SpkSeg

#include <iostream>
#include <iostream>
#include <fstream>  // pour outFile
#include <cstdio>   // pour printf()
#include <cassert> // pour le debug pratique
#include <cmath>
#include <liatools.h>

using namespace alize;
using namespace std;


int main(int argc, char* argv[])
{

    try
    {
      Config tmp;
      CmdLine cmdLine(argc, argv);
      cmdLine.copyIntoConfig(tmp);
      if (tmp.getParamCount()==0)
	{ cerr << "Error: type TimeCluster.exe --help to get usage" << endl;}
      Config config(tmp.getParam("config"));
      config.setParam(tmp);
     
      if (config.existsParam("debug"))debug=true; else debug=false;  
      if (config.existsParam("verbose"))verbose=true; else verbose=false;
      double frameLength = 0.01;
      if (config.existsParam("frameLength")) 
	frameLength=config.getParam("frameLength").toDouble();                                // length in s of a frame
      bool timeMode=false;
      bool frmMode=true;
      if (config.existsParam("timeMode")) timeMode=!(config.getParam("timeMode")=="false");
      if (config.existsParam("frmMode")) frmMode=!(config.getParam("frmMode")=="false");
      String file=config.getParam("filename");
      String labelSelectedFrames=config.getParam("labelSelectedFrames");
      SegServer segServer;                
      LabelServer labelServer;
      loadClusterFile(file,segServer,labelServer,config);
      long codeSelectedFrame=labelServer.getLabelIndexByString(labelSelectedFrames);          // Get the index of the selected cluster
      if (codeSelectedFrame==-1){                                                             // No data for this model !!!!!!!!!!!!!!
	cout << " WARNING - NO DATA with the label["<<labelSelectedFrames<<"] in file ["<<file<<"]"<<endl;
	exit(0);
      } 
      SegCluster& cluster=segServer.getCluster(codeSelectedFrame);                             // Gives the cluster of the selected segs     
      unsigned long time=totalFrame(cluster);
      if(verbose) cout << file;
      if (timeMode){
      	if(verbose) cout << " time:";
	cout << frameIdxToTime(time,frameLength);
      }
      if (frmMode){
       if(verbose) cout << " nb frames:";
       cout << time;
      }
      cout <<endl;
    }
    catch (alize::Exception& e)
    {
        cout << e.toString() << endl;
    }

    return 0;
}

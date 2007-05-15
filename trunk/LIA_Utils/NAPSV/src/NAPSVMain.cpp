// NAPSVMain.cpp
// This file is a part of LIA Software LIA_SpkDet, based on ALIZE toolkit 
// LIA_SpkDet  is a free, open tool for speaker recognition
// LIA_SpkDet is a development project initiated and funded by the LIA lab.
// See www.lia.univ-avignon.fr
// 
// ALIZE is needed for LIA_SpkDet
// for more information about ALIZE, see http://www.lia.univ-avignon.fr/heberges/ALIZE/
//
// Copyright (C) 2004
//  Laboratoire d'informatique d'Avignon [www.lia.univ-avignon.fr]
//  Jean-Francois Bonastre [jean-francois.bonastre@lia.univ-avignon.fr]
//      
// LIA_SpkDet is free software; you can redistribute it and/or
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
// more information about the licence or the use of LIA_SpkDet
// Author: Nicolas Scheffer

#include <iostream>
#include "liatools.h"

void nap(Matrix <double>& sv,Matrix <double>& NAP) {
    RealVector <double> in;
    in.setSize(sv.cols());
    RealVector <double> proj;
    proj.setSize(sv.cols());    
    for (unsigned long i=0;i<in.size();i++)
        in[i]=sv(0,i);
    projectOnSubSpace(in,NAP,proj);
    if (debug) cout <<"* sv:"<<sv(0,50)<<" proj:"<<proj[50]<<endl;    
    in-=proj;
    sv=(Matrix <double>)in;
    if (debug) cout <<"** sv:"<<sv(0,50)<<" proj:"<<proj[50]<<endl;        
}

void project(Matrix <double>& sv,Matrix <double>& NAP) {
    RealVector <double> in;
    in.setSize(sv.cols());
    RealVector <double> proj;
    proj.setSize(sv.cols());    
    for (unsigned long i=0;i<in.size();i++)
        in[i]=sv(0,i);
    projectOnSubSpace(in,NAP,proj);
    sv=(Matrix <double>)proj;
    if (debug) cout <<"** sv:"<<sv(0,50)<<" proj:"<<proj[50]<<endl;        
}

int main(int argc, char* argv[]){

    ConfigChecker cc;
    cc.addStringParam("config", false, true, "default config filename");
    try{
    CmdLine cmdLine(argc, argv);
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

    // Begins here
    XList inputList(config.getParam("inputFilename"));
    String vPath=config.getParam("vectorFilesPath");
    String vExt=config.getParam("vectorFilesExt");        
    String vOExt=config.getParam("napvectorFilesExt");  
    Matrix <double> NAP;
    NAP.load(config.getParam("NAP"),config);
    if (verbose) cout << "NAP matrix dimensions: ("<<NAP.rows()<<","<<NAP.cols()<<")"<<endl;
     for (unsigned long r=0;r<inputList.getLineCount();r++) {
            String in=vPath+inputList.getLine(r).getElement(0)+vExt;
            String out=vPath+inputList.getLine(r).getElement(0)+vOExt;         
            if (verbose) cout << "Processing ["<<in<<"] to ["<<out<<"]"<<endl;         
            Matrix <double> sv;        
            sv.load(in,config);
            if (verboseLevel > 1) cout << "Vector dimensions: ("<<sv.rows()<<","<<sv.cols()<<")"<<endl;    
            if (config.existsParam("project")) project(sv,NAP);
            nap(sv,NAP);
            sv.save(out,config);
    }
    return 0;
  }
  catch(alize::Exception & e){ cout <<"TrainTarget "<< e.toString() << endl << cc.getParamList()<<endl;}
}

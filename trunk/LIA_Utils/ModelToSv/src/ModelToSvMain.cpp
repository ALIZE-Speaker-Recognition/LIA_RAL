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
#include "liatools.h"

void getMeanNorm(RealVector <double>& norm,MixtureGD & UBM) {
    unsigned long mSize=UBM.getDistribCount();
    unsigned long vSize=UBM.getVectSize();            
    norm.setSize(mSize*vSize);    
    for (unsigned long i=0;i<mSize;i++) {
        DistribGD & d=UBM.getDistrib(i);
        double w=UBM.weight(i);
        for (unsigned long j=0;j<vSize;j++)
            norm[i*vSize+j]=sqrt(w*d.getCovInv(j)); 
    }
}

void getWeightNorm(RealVector <double>& norm,MixtureGD & UBM) {
    unsigned long mSize=UBM.getDistribCount();    
    norm.setSize(mSize);       
    for (unsigned long i=0;i<mSize;i++) 
        norm[i]=1/sqrt(UBM.weight(i));
}

int main(int argc, char* argv[]){

    ConfigChecker cc;
    cc.addStringParam("config", false, true, "default config filename");    
    cc.addBooleanParam("meanSv", true, true, "output Mean vector");
    cc.addBooleanParam("weightSv", true, true, "output Weight vector");    
    cc.addBooleanParam("normSv", true, true, "normalize vector (need UBM)");
    cc.addBooleanParam("vectors", false, true, "read from vectors");       
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
    bool mean=config.getParam("meanSv").toBool();
    bool weight=config.getParam("weightSv").toBool();
    if (!(weight^mean)) throw Exception("weightSv and meanSv have same value !",__FILE__,__LINE__);
    bool norm=config.getParam("normSv").toBool();    
    bool vectors=false;
    if (config.existsParam("vectors")) vectors=true; // read either from models or vectors
    bool gmm=!vectors;
    if (verbose) {
        if (mean) cout << "(modelToSv) MeanSv"<<endl;
        if (weight) cout << "(modelToSv) WeightSv"<<endl;
        if (norm) cout << "(modelToSv) NormSv: Normalizing vector with UBM model "<<endl;
        if (norm && !(weight||mean)) cout << "(modelToSv) NormSv: only normalize vectors" << endl;
        }
    XList inputList(config.getParam("inputFilename"));
    String vPath=config.getParam("vectorFilesPath");
    String vExt=config.getParam("vectorFilesExtension");        
    MixtureServer ms(config);
    MixtureGD UBM=ms.loadMixtureGD(config.getParam("inputWorldFilename"));
    RealVector <double> normVector;
    if (norm) {
        if (verbose)  cout << "(modelToSv) Compute normalisation vector"<<endl;
        if(mean) getMeanNorm(normVector,UBM); // do seometing better, if only want to normalize assume it's mean
        if (weight) getWeightNorm(normVector,UBM);     
        String normFile=vPath+"norm"+vExt; 
        if (debug ) {
            ((Matrix<double>)normVector).save(normFile,config);
            cout<<"(modelToSv) Saving normalisation vector to ["<<normFile<<"]"<<endl;
        }
    }
         
    for (unsigned long r=0;r<inputList.getLineCount();r++) {
            String out=vPath+inputList.getLine(r).getElement(0)+vExt;         
            if (verbose) cout << "(modelToSv) Processing ["<<inputList.getLine(r).getElement(0)<<"] to ["<<out<<"]"<<endl;         
            RealVector <double> sv;        
            if (gmm) { // read from models
                if (verbose) cout << "(modelToSv) Read from GMM models"<<endl;
                MixtureGD &curr=ms.loadMixtureGD(inputList.getLine(r).getElement(0));
                if (mean) modelToSv(curr,sv);
                if (weight) {
                    unsigned long mSize=curr.getDistribCount();
                    sv.setSize(mSize);        
                    for (unsigned long i=0;i<mSize;i++)
                        sv[i]=curr.weight(i);
                }
                ms.deleteMixture(curr);
                ms.deleteUnusedDistribs();  
            }
            else if (vectors) {//read from vectors (mean and weight at false)
                if (verbose) cout << "(modelToSv) Read from SV from vectors"<<endl;     
                String vIExt=config.getParam("inputVectorFilesExtension"); 
                String in=vPath+inputList.getLine(r).getElement(0)+vIExt;   
                Matrix <double> a;
                a.load(in,config);   
                sv.setSize(a.cols());
                for (unsigned long i=0;i<sv.size();i++)
                    sv[i]=a(0,i);
            }
            if(norm) {
                for (unsigned long i=0;i<sv.size();i++) 
                    sv[i]*=normVector[i];
            }
            ((Matrix<double>)sv).save(out,config);        
    }
    return 0;
  }
  catch(alize::Exception & e){ cout <<"modelToSv"<< e.toString() << endl << cc.getParamList()<<endl;}
}

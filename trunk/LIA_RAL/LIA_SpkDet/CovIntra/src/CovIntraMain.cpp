#if !defined(ALIZE_CovIntraMain_cpp)
#define ALIZE_CovIntraMain_cpp

#include <iostream>
#include "liatools.h"
#include "CovIntra.h"
extern "C" {
#include "svdlib.h" 
}

/* ******************************
Cov Intra: Compute the NAP matrix with an eigenvalue decomposition of CC'.
Find eigenvalues of C with svdlibc, thanks authors: C=Ut.S.V, where S*S are eigenvalues of CC'
Ut left eigenvectors are saved.
C is found by computing channel covariance (remove supervector mean for each spealer)
*/

int main(int argc, char* argv[]){

	ConfigChecker cc;
try {
	cc.addStringParam("ndx",true,true,"NDX of multiple GMM speaker recordings");
	cc.addStringParam("channelMatrix",true,true,"Output file for matrix");	
	cc.addIntegerParam("nbEigenVectors",true,true,"nb factors");	
	cc.addStringParam("svd",false,false,"do the svd");	
	cc.addStringParam("gmm",false,false,"get mean super vector from gmms");
	cc.addIntegerParam("vsize",false,false,"size of vectors");	
	
	CmdLine cmdLine(argc, argv);
	if (cmdLine.displayHelpRequired()){
		cout <<"CovIntra.exe"<<endl<<"This program is used to compute covariance channel matrix (and optionally do the svd)"
	   <<endl<<cc.getParamList()<<endl;
	return 0;  
	}
	if (cmdLine.displayVersionRequired()) cout <<"Version 2-beta"<<endl;
	Config tmp;
	cmdLine.copyIntoConfig(tmp);
	Config config;
	if (tmp.existsParam("config")) config.load(tmp.getParam("config"));
	cmdLine.copyIntoConfig(config);
	cc.check(config);
	debug=config.getParam_debug();
        if (config.existsParam("verbose"))verbose=config.getParam("verbose").toBool();
        else verbose=false;
        if (verbose) verboseLevel=1;else verboseLevel=0;
        if (config.existsParam("verboseLevel"))verboseLevel=config.getParam("verboseLevel").toLong();
        if (verboseLevel>0) verbose=true;        
        //if (config.existsParam("binary")) CovIntraBin(config);
	CovIntra(config);
	}
catch (Exception& e) {cout << e.toString() << cc.getParamList() << endl;}
} 
#endif

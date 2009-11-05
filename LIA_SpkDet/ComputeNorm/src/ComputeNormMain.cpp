#include <iostream>

#include "ComputeNorm.h"
#include <liatools.h>
using namespace std;
using namespace alize;

int main(int argc, char* argv[])
{
  ConfigChecker cc;
  cc.addStringParam("testNistFile",true, true,"target_seg score file");
  cc.addStringParam("normType",true, true,"tnorm|znorm|ztnorm|tznorm, the normalization method (ztnorm is outputing both tnorm and ztnorm scores, tznorm both tznorm and znorm scores)");
  cc.addStringParam("tnormNistFile",false,true,"imp_seg score file, impostor scores used for tnorm and ztnorm");
  cc.addStringParam("znormNistFile",false,true,"target_imp score file, impostor scores used for znorm and ztnorm");
  cc.addStringParam("ztnormNistFile",false,true,"imp_imp score file, impostor scores used for ztnorm and tznorm");
  cc.addStringParam("outputFileBaseName",true, true,"the output file(s) basename");
  cc.addStringParam("znormFilesExtension",false, true,"znorm output file extension (default='.znorm'");
  cc.addStringParam("ztnormFilesExtension",false, true,"ztnorm output file extension (default='.znorm'");
  cc.addStringParam("tnormFilesExtension",false, true,"tnorm output file extension (default='.tnorm'");
  cc.addStringParam("tznormFilesExtension",false, true,"tznorm output file extension (default='.tznorm'");
  cc.addStringParam("cohortFilePath",false, true,"cohort files path, for selectTargetDependentCohortInFile");
  cc.addStringParam("cohortFileExt",false, true,"cohort files extension, for selectTargetDependentCohortInFile");
  cc.addIntegerParam("maxIdNb",false, true,"Max target speakers - use to fix the max number of znorm score distributions (default=1000)");
  cc.addIntegerParam("maxSegNb",false, true,"Max test segments - use to fix the max number of tnorm score distributions (default=1000)");
  cc.addIntegerParam("maxScoreDistribNb",false, true,"Max scores per distribution - use to fix the max number of score in a distribution (default=1000)");
  cc.addStringParam("selectType",false, true,"Define the score selection method,'noSelect|selectNBestByTestSegment|selectTargetDependentCohortInFile'  (default='noSelect')");
  cc.addIntegerParam("fieldGender",false, true,"The field for gender in the nist file format (default=0)");
  cc.addIntegerParam("fieldName",false, true,"The field for gender in the nist file format (default=1)");
  cc.addIntegerParam("fieldDecision",false, true,"The field for gender in the nist file format (default=2)");
  cc.addIntegerParam("fieldSeg",false, true,"The field for gender in the nist file format (default=3)");
  cc.addIntegerParam("fieldLLR",false, true,"The field for gender in the nist file format (default=4)");
  try {
      CmdLine cmdLine(argc, argv);
      if (cmdLine.displayHelpRequired()){
	cout << "************************************" << endl;
	cout << "********** ComputeNorm.exe *************" << endl;
	cout << "************************************" << endl;
	cout << endl;
	cout << "Apply Z-T-ZT or ZT Norm to a Score File" << endl<< "ztnorm includes tnorm"<<endl;
        cout <<endl<<cc.getParamList()<<endl;
        return 0;  
      }
      if (cmdLine.displayVersionRequired()){
        cout <<"Version 3"<<endl;
      } 
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
      // Initialize the default values for the parameters TODO : add this option in the addParam functions...
 	  if(!config.existsParam("znormFilesExtension")) config.setParam("znormFilesExtension",".znorm");
 	  if(!config.existsParam("tnormFilesExtension")) config.setParam("tnormFilesExtension",".tnorm");
	  if(!config.existsParam("ztnormFilesExtension")) config.setParam("ztnormFilesExtension",".ztnorm");
	  if(!config.existsParam("tznormFilesExtension")) config.setParam("tznormFilesExtension",".tznorm");
	  if(!config.existsParam("maxIdNb")) config.setParam("maxIdNb","1000");
	  if(!config.existsParam("maxSegNb")) config.setParam("maxSegNb","1000");
	  if(!config.existsParam("maxScoreDistribNb")) config.setParam("maxScoreDistribNb","1000");
	  if(!config.existsParam("selectType")) config.setParam("selectType","noSelect");
	  if(!config.existsParam("fieldGender")) config.setParam("fieldGender","0");
	  if(!config.existsParam("fieldName")) config.setParam("fieldName","1");
	  if(!config.existsParam("fieldDecision")) config.setParam("fieldDecision","2");	  	  	  	  	  	  	   	  
	  if(!config.existsParam("fieldSeg")) config.setParam("fieldSeg","3");	  	  	  	  	  	  	   	  
	  if(!config.existsParam("fieldLLR")) config.setParam("fieldLLR","4");	  	  	  	  	  	  	   	  

      // start the prog!
      ComputeNorm(config);
    }
catch (alize::Exception& e) {cout << e.toString() << endl << cc.getParamList()<< endl;}
return 0;     
}



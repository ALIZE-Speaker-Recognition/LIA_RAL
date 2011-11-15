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
#include "ShiftedDeltaFeat.h"
#include <alize.h>
#include <liatools.h>

//---------------------------------------------------------------------------------------------------------------
/*!
 * \brief Reads a parametric file (or list of prm files) containing cepstral features and
 *        writes a file containing the shifted delta features.
 * \param argc	number or arguments
 * \param argv	argument array - in usual Mistral format [see detailed description below for config params]
 * \date 11 april 2008
 * 
 * Config-Params are:
 *  inputFeatureFilename	could be a simple feature file or a list of filename
 *  featureServerMode		should really be FEATURE_WRITABLE to write features to outpu file
 *  vectSize				SDC param N (below)
 *  SDCdeltaDelay			SDC param d (below)
 *  SDCtimeAdvance			SDC param P (below)
 *  SDCdeltaBlocks			SDC param k (below)
 *  [SDCkeepCepstra]		bool, if we want to keep the cepstral features in the output
 * 
 * SDC parameters are:
 *  N	number of cepstral coefficients
 *  d	advance and delay for the delta computation (in each delta block)
 *  P	time shift between consecutive delta blocks
 *  k	number of delta blocks whose delta coefficients are concatenated to form the final feature vector
 *
 * For each input vector at "time" t:
 *   ShDeltaO_h(t,i) = O_h( t + i*P + d) - O_h( t + i*P - d)   ;  i = 0, 1, 2 ... k-1  ;  O_h is the input featureVector
 *   ShiftedDeltaCepstraVectorO_h(t) = concatenation of all ShDeltaO_h(t,i) ; i = 0, 1, 2 ... k-1
 *
 *      t        t+P       t+2P    ...     t+(k-1)P
 *      |         |         |                 |
 * t----------------------------------------------------->
 *    |   |     |   |     |   |    ...      |   |
 *   -d  +d    -d  +d    -d  +d            -d  +d
 *     '-'       '-'       '-'               '-'
 *      |         |         |                 |
 *    block0    block1    block2   ..    block(k-1)
 * 
 * SDC-vector: [ block0' block1' block2' ... block(k-1)' ]'  ; Thus, its dimension is k*N
 * 
 * There is the possibility to keep the cepstral features in the output vector:
 *   SDC-vector: [ cepstra' block0' block1' block2' ... block(k-1)' ]'  ; with dimension N+(k*N)
 * 
 * Often used SDC parameters are:  7-1-3-7  or 7-1-3-4 or 10-1-3-3
 * 
 * Refs (param description and formula):
 * - Matějka Pavel, Burget Lukáš, Schwarz Petr, Černocký Jan: Brno University of Technology System 
 *    for NIST 2005 Language Recognition Evaluation, In: Proceedings of Odyssey 2006: The Speaker 
 *    and Language Recognition Workshop, San Juan, PR, 2006, p. 57-64, ISBN 1-4244-0472-X
 *    http://www.fit.vutbr.cz/research/view_pub.php?id=8135
 *    http://www.fit.vutbr.cz/~matejkap/pubs.php
 * 
 */
int main (int argc, char *argv[])
{
  using namespace std; 
  using namespace alize;

    ConfigChecker cc;
    cc.addStringParam("config", false, true, "default config filename");
    cc.addStringParam("inputFeatureFilename",true,true,"input feature - could be a simple feature file or a list of filename");
    // SDCnumCeps is param vectSize
    cc.addIntegerParam("SDCdeltaDelay",true,true,"Advance and delay for delta computation: F(t+d)-F(t-d); parameter d in N-d-P-k");  
    cc.addIntegerParam("SDCtimeAdvance",true,true,"Time shift between consecutive delta blocks; parameter P in N-d-P-k");  
    cc.addIntegerParam("SDCdeltaBlocks",true,true,"Number of delta blocks that will be concatenated to form the SDC vector; parameter k in N-d-P-k");  
    cc.addBooleanParam("SDCkeepCepstra",false,true,"Keep orignial cepstral features and append SDC features; feature vector will be that bigger (default: true)");  
    cc.addStringParam("featureServerMode",true,true,"should really be FEATURE_WRITABLE to write features to outpu file"); // TODO code it implicitely  
  try
  {
    CmdLine cmdLine (argc, argv);
    if (cmdLine.displayHelpRequired ()){	// --help
      cout << "ShiftedDeltaFeat" << endl;
      cout << "ShiftedDeltaFeat.exe --config <foo.cfg> --inputFeatureFileName <foo.prm> --SDCdeltaDelay <d> --SDCtimeAdvance <P> --SDCdeltaBlocks <k>"<< endl;
      cout<<cc.getParamList()<<endl;
      }
    else if (cmdLine.displayVersionRequired ())	// --version
      cout << "Version 1.0" << endl;
    else{
      Config tmp;
      cmdLine.copyIntoConfig (tmp);
      Config config;
      if (tmp.existsParam("config")) config.load(tmp.getParam("config"));
      cmdLine.copyIntoConfig(config);
      cc.check(config);
      debug=config.getParam_debug();
      if (config.existsParam("verbose")) verbose=config.getParam("verbose").toBool();
      else verbose=false;
      if (verbose) verboseLevel=1;else verboseLevel=0;
      if (config.existsParam("verboseLevel"))verboseLevel=config.getParam("verboseLevel").toLong();
      if (verboseLevel>0) verbose=true;

      // DO IT
      ShiftedDeltaFeat(config);

    }
  }
  catch (alize::Exception & e) {cout << e.toString () << endl << cc.getParamList() << endl;}
  #ifdef NDEBUG 
  cout<<"*** Objects created and destroyed **"<<Object::getCreationCounter()<<"-"<<Object::getDestructionCounter()<<endl;    
  #endif
  return 0;
}

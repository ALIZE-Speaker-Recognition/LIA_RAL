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

#if !defined(ALIZE_Hmm_cpp)
#define ALIZE_Hmm_cpp

#include <alize.h>
#include "Hmm.h"

using namespace alize;
using namespace std;
	
// Class hmm
	
	// --------------------------------------------------------------------------------
	
unsigned long hmm::getNbState(){
	 return (tabState.size());
}

unsigned long hmm::getStateIndex(String name){
	return tabStateName.getIndex(name);
}

	// --------------------------------------------------------------------------------
void hmm::setDensity(String fileName, unsigned long nModel)
{
	 tabState.setObject(ms->loadMixture(fileName),nModel);
}

	// --------------------------------------------------------------------------------
Mixture& hmm::getDensity(unsigned long nModel)
{
	 return (static_cast<Mixture&>(tabState.getObject(nModel)));
}

	// --------------------------------------------------------------------------------
void hmm::setDensity(Mixture& m, unsigned long nModel)
{
	// A MODIFIER POUR DETRUIRE SI NECESSAIRE LA MIXTURE CREE AUPARAVANT !
	unsigned long temp;

	temp=ms->getMixtureIndex(static_cast<Mixture&>(tabState.getObject(nModel)).getId());
	ms->deleteMixtures(temp,temp);
	ms->deleteUnusedDistribs();
 	tabState.setObject(ms->duplicateMixture(m,DUPL_DISTRIB),nModel);
}

	// --------------------------------------------------------------------------------
void hmm::setDensity(MixtureGD& m, unsigned long nModel)
{
	unsigned long temp;

	temp=ms->getMixtureIndex(static_cast<MixtureGD&>(tabState.getObject(nModel)).getId());
	ms->deleteMixtures(temp,temp);
	ms->deleteUnusedDistribs();
 	tabState.setObject(ms->duplicateMixture(m,DUPL_DISTRIB),nModel);
}

//-------------------------------------------------------------------------------
const String &hmm::getStateName(unsigned long nModel){
	 return tabStateName.getElement(nModel,false);
}

//-------------------------------------------------------------------------------
void hmm::setStateName(unsigned long nModel, String s){
	 tabStateName.getElement(nModel)=s;
}

//--------------------------------------------------------------------------
unsigned long hmm::LoadState(String fileName)
{
tabState.addObject(ms->loadMixture(fileName));
tabStateName.addElement(fileName);		
transitions.setSize(tabState.size()*tabState.size());
return (tabState.size());
}

//--------------------------------------------------------------------------
unsigned long hmm::LoadState(String fileName,String name)
{
tabState.addObject(ms->loadMixture(fileName));
tabStateName.addElement(name);		
transitions.setSize(tabState.size()*tabState.size());

return (tabState.size());
}


//--------------------------------------------------------------------------
unsigned long hmm::LoadState(Mixture& m,String name)
{
tabState.addObject(ms->duplicateMixture(m,DUPL_DISTRIB));
tabStateName.addElement(name);		
transitions.setSize(tabState.size()*tabState.size());

return (tabState.size());
}

//--------------------------------------------------------------------------
unsigned long hmm::LoadState(MixtureGD& m,String name)
{
tabState.addObject(ms->duplicateMixture(m,DUPL_DISTRIB));
tabStateName.addElement(name);		
transitions.setSize(tabState.size()*tabState.size());

return (tabState.size());
}

//--------------------------------------------------------------------------
unsigned long hmm::LoadState(Mixture& m,String name, unsigned long l)
{
tabState.addObject(ms->duplicateMixture(m,DUPL_DISTRIB));
tabStateName.addElement(name);		
transitions.setSize(tabState.size()*tabState.size());
length.addValue(l);

return (tabState.size());
}

//--------------------------------------------------------------------------
unsigned long hmm::LoadState(MixtureGD& m,String name, unsigned long l)
{
tabState.addObject(ms->duplicateMixture(m,DUPL_DISTRIB));
tabStateName.addElement(name);		
transitions.setSize(tabState.size()*tabState.size());
length.addValue(l);

return (tabState.size());
}


//--------------------------------------------------------------------------
unsigned long hmm::LoadState(Mixture& m,String name, unsigned long l, unsigned long r)
{
tabState.addObject(ms->duplicateMixture(m,DUPL_DISTRIB));
tabStateName.addElement(name);		
transitions.setSize(tabState.size()*tabState.size());
length.addValue(l);
replacement.addValue(r);

return (tabState.size());
}

//--------------------------------------------------------------------------
unsigned long hmm::LoadState(MixtureGD& m,String name, unsigned long l, unsigned long r)
{
tabState.addObject(ms->duplicateMixture(m,DUPL_DISTRIB));
tabStateName.addElement(name);		
transitions.setSize(tabState.size()*tabState.size());
length.addValue(l);
replacement.addValue(r);

return (tabState.size());
}

//--------------------------------------------------------------------------
unsigned long hmm::addState(String nameState)
{
tabState.addObject(ms->createMixture());
tabStateName.addElement(nameState);		
transitions.setSize(tabState.size()*tabState.size());

return (tabState.size());
}
//-----------------------------------------------------------------------------------
double hmm::getTransition(int i,int j){
	return transitions[i+j*tabState.size()];

}
//-----------------------------------------------------------------------------------
void hmm::getTransition(DoubleVector &tmpTrans){

tmpTrans.setSize(tabState.size()*tabState.size());
for(unsigned long i=0;i<tabState.size();i++)
	for(unsigned long j=0;j<tabState.size();j++){
		tmpTrans[i+j*tabState.size()]=transitions[i+j*tabState.size()];
	}
}

//-------------------------------------------------------------------------------------
void hmm::setTransition(double a,int i,int j){
  transitions[i+j*tabState.size()]=a;
}

//-----------------------------------------------------------------------------------
unsigned long hmm::getLength(unsigned long i){

return length[i];
}

//-----------------------------------------------------------------------------------
unsigned long hmm::getLength(String name){

return length[getStateIndex(name)];
}

//-----------------------------------------------------------------------------------
unsigned long hmm::getReplacement(String name){

return replacement[getStateIndex(name)];
}

//-----------------------------------------------------------------------------------
unsigned long hmm::getReplacement(unsigned long i){

return replacement[i];
}



//-----------------------------------------------------------------------------------
void hmm::setLength(unsigned long i, unsigned long l){

length[i]=l;
}

//-----------------------------------------------------------------------------------
void hmm::setReplacement(unsigned long i, unsigned long r){

replacement[i]=r;
}

//--------------------------------------------------------------------------------------
void hmm::reset(void){
  unsigned long temp;
  for(unsigned long i=0;i<tabState.size();i++){
    temp=ms->getMixtureIndex(static_cast<Mixture&>(tabState.getObject(i)).getId());
    ms->deleteMixtures(temp,temp);
  }
  ms->deleteUnusedDistribs();
  tabState.clear();
  tabStateName.reset();
  transitions.clear();
}
//----------------------------------------------------------------------------------
unsigned long hmm::deleteState(unsigned long indice){
  unsigned long temp;
  temp=ms->getMixtureIndex(static_cast<Mixture&>(tabState.getObject(indice)).getId()); // TODO, change for the mixture delete function
  ms->deleteMixtures(temp,temp);
  ms->deleteUnusedDistribs();
  tabState.removeObjects(indice,indice);
  tabStateName.deleteElement(tabStateName.getElement(indice));
  transitions.setSize(tabState.size()*tabState.size());
  if(length.size() != 0)
  	length.removeValues(indice, indice);
  return temp; // TODO VERIFY THE RETURN VALUE - JFB
}


	// --------------------------------------------------------------------------------
hmm::hmm(MixtureServer & m, Config &config)
:ms(&m),conf(&config)
{
}
//---------------------------------------------------------------------
const hmm& hmm::operator=(const hmm& hmmc)
{
	assign(hmmc);
	return *this;
}
//-------------------------------------------------------------------------
void hmm::assign(const hmm& hmmc){
  this->reset();
  ms=hmmc.ms;
  conf=hmmc.conf;
  tabStateName=hmmc.tabStateName;
  //copie des données
  transitions=hmmc.transitions;
  for(unsigned long i=0;i<hmmc.tabState.size();i++) {
    tabState.addObject(ms->duplicateMixture(static_cast<Mixture&>( hmmc.tabState.getObject(i)),DUPL_DISTRIB));
  }
}

//-----------------------------------------------------------------------	
hmm::~hmm() { }

//-------------------------------------------------------------------------
hmm::hmm(const hmm& h){
  assign(h);
}
//------------------------------------------------------------------------
void hmm::load()
{
//TO DO
}

//-------------------------------------------------------------------------
void hmm::save()
{
//TO DO
}


// Return a copy of the transition matrix - should be included into hmm class - TODO
/*DoubleVector &copyTransition(hmm& cHmm)
{
  DoubleVector transitions;
  transitions.setSize(cHmm.getNbState()*cHmm.getNbState());
  for(int i=0;i<cHmm.getNbState();i++)
    for(int j=0;j<cHmm.getNbState();j++){
      transitions[i+j*cHmm.getNbState()]=cHmm.getTransition(i,j);
    }
  return transitions;
}
*/


#endif


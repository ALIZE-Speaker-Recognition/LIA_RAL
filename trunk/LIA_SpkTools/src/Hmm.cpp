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

	// --------------------------------------------------------------------------------
void hmm::setDensity(String fileName, unsigned long nModel)
{
	 tabState.setObject(ms->loadMixtureGD(fileName),nModel);
}

	// --------------------------------------------------------------------------------
MixtureGD& hmm::getDensity(unsigned long nModel)
{
	 return (static_cast<MixtureGD&>(tabState.getObject(nModel)));
}

	// --------------------------------------------------------------------------------
void hmm::setDensity(MixtureGD& m, unsigned long nModel)
{
	// A MODIFIER POUR DETRUIRE SI NECESSAIRE LA MIXTURE CREE AUPARAVANT !
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
tabState.addObject(ms->loadMixtureGD(fileName));
tabStateName.addElement(fileName);		
transitions.setSize(tabState.size()*tabState.size());

return (tabState.size());
}

//--------------------------------------------------------------------------
unsigned long hmm::LoadState(String fileName,String name)
{
tabState.addObject(ms->loadMixtureGD(fileName));
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
unsigned long hmm::addState(String nameState)
{
tabState.addObject(ms->createMixtureGD());
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

//--------------------------------------------------------------------------------------
void hmm::reset(void){
  unsigned long temp;
  for(unsigned long i=0;i<tabState.size();i++){
    temp=ms->getMixtureIndex(static_cast<MixtureGD&>(tabState.getObject(i)).getId());
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
  temp=ms->getMixtureIndex(static_cast<MixtureGD&>(tabState.getObject(indice)).getId()); // TODO, change for the mixture delete function
  ms->deleteMixtures(temp,temp);
  ms->deleteUnusedDistribs();
  tabState.removeObjects(indice,indice);
  tabStateName.deleteElement(tabStateName.getElement(indice));
  transitions.setSize(tabState.size()*tabState.size());
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
    tabState.addObject(ms->duplicateMixtureGD(static_cast<MixtureGD&>( hmmc.tabState.getObject(i)),DUPL_DISTRIB));
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


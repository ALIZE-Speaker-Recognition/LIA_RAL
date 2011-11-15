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

#if !defined(ALIZE_Hmm_h)
#define ALIZE_Hmm_h

#if defined(_WIN32)
#if defined(LIA_SPKTOOLS_EXPORTS)
#define LIA_SPKTOOLS_API __declspec(dllexport)
#else
#define LIA_SPKTOOLS_API __declspec(dllimport)
#endif
#else
#define LIA_SPKTOOLS_API
#endif

#include <alize.h>


using namespace alize;
using namespace std;

class LIA_SPKTOOLS_API hmm{
	MixtureServer *ms;
	Config *conf;
	ObjectRefVector tabState;
	XLine tabStateName;
	DoubleVector transitions;
	ULongVector	length;	
	ULongVector	replacement;	
	
	public:
	unsigned long getNbState();
	unsigned long getStateIndex(String);
	unsigned long LoadState(String);
	unsigned long LoadState(String,String);
	unsigned long LoadState(Mixture&,String);	
	unsigned long LoadState(MixtureGD&,String);	
	unsigned long LoadState(Mixture&,String,unsigned long);	
	unsigned long LoadState(Mixture&,String,unsigned long, unsigned long);	
	unsigned long LoadState(MixtureGD&,String,unsigned long);	
	unsigned long LoadState(MixtureGD&,String,unsigned long, unsigned long);	
	unsigned long addState(String);
	unsigned long deleteState(unsigned long indice);
	double getTransition(int i,int j);
	void getTransition(DoubleVector &);
	void setTransition(double,int , int );
	void reset();
	const hmm& operator=(const hmm& hmmc);
	Mixture& getDensity(unsigned long);
	void setDensity(Mixture& m, unsigned long nModel);
	void setDensity(MixtureGD& m, unsigned long nModel);
	void setDensity(String, unsigned long);
	const String &getStateName(unsigned long);
	void setStateName(unsigned long, String);
	void setLength(unsigned long, unsigned long);
	void setReplacement(unsigned long, unsigned long);
	unsigned long getLength(unsigned long);
	unsigned long getLength(String);
	unsigned long getReplacement(unsigned long);
	unsigned long getReplacement(String);
	hmm(MixtureServer &, Config &);
	hmm(const hmm&);
	~hmm();
	
	private:
	void assign(const hmm& hmmc);
	void save();
	void load();
};

// Return a copy of the transition matrix - should be included into hmm class - TODO
//DoubleVector &copyTransition(hmm& cHmm);
#endif

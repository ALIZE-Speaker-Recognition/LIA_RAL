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

#if !defined(ALIZE_TopGauss_h)
#define ALIZE_TopGauss_h

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
#include "liatools.h"

using namespace alize;
using namespace std;

class LIA_SPKTOOLS_API TopGauss {
	private:	
		ULongVector _idx; // gaussian idxs
		ULongVector _nbg; // nb gauss per frame
		DoubleVector _snsw; // resting weight
		DoubleVector _snsl; // resting lk
	
		unsigned long _nt; //nb frames
		unsigned long _nbgcnt; //nb gauss per file
	
	public:
		TopGauss() {
		};
		~TopGauss() {
		};
		double compute(MixtureGD &,FeatureServer &,String &,Config &);
		double get(MixtureGD &,FeatureServer &,String &,Config &);
		
		double compute(MixtureGD &,String &,Config &);
		RealVector <double> compute(MixtureGD &,FeatureServer &,Config &);
		
		double get(MixtureGD &,String &,Config &);
		RealVector <double> get(MixtureGD &,FeatureServer &,Config &);		
		
		void read(String &,Config &);
		void write(String &,Config &);
		
		ULongVector & nbg() {
			return _nbg;}
		ULongVector & idx() {
			return _idx;}			
		unsigned long frameToIdx(unsigned long &);
		unsigned long nbgcnt() {
			return _nbgcnt;}			
		unsigned long nt() {
			return _nt;}			
};
	
#endif
	

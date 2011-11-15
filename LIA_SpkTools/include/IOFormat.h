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

#if !defined(ALIZE_IOFormat_h)
#define ALIZE_IOFormat_h

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

//******************************************************
// These functions are used to ouput Objects into a vector
//******************************************************
//------------------------------------------------------------------------
// Output Weight of a  model in a vector -- Nicolas SCHEFFER (depracated)
LIA_SPKTOOLS_API void outputWeightVector(MixtureGD& world,MixtureGD& clientMixture, double frameCount, Config& config); 
//------------------------------------------------------------------------
// Output vector in SVMLight Format -- Nicolas SCHEFFER
LIA_SPKTOOLS_API void outputSVMLightVector(RealVector<double	>&,const String &,Config&); 
//------------------------------------------------------------------------
// Output Matrix in SVMLight Format
LIA_SPKTOOLS_API void outputDSM(DoubleSquareMatrix &,MixtureGD&,Config&);

//******************************************************
// These functions are related to file format
//******************************************************

//-------------------------------------------------------------------------
// LIA Internal Format
//Output a result line in a file in LIA format: F client 1 test 0 20 -0.02
LIA_SPKTOOLS_API void outputResultLine(double LLR, const String & clientName, const String & seg, const String &gender, int decision, ostream& strm);

LIA_SPKTOOLS_API void outputResultLine(double LLR, const String & clientName, const String & seg,real_t begin,real_t end, const String &gender, int decision, ostream& strm);

LIA_SPKTOOLS_API void outputResultLIARALLine(const String &gender, const String & clientName, const String & channel, const String & seg, const String & start, const String & duration, double LLR, ostream& strm);

//-------------------------------------------------------------------------
//Output a result line in NIST 04 format
LIA_SPKTOOLS_API void outputResultNIST04Line(const String &trainTypeTest, const String &adaptationMode, const String &segTypeTest, const String &gender, const String & clientName, const String & seg, const String &decision, double LLR, ostream& strm);

//-------------------------------------------------------------------------
//Output a result line in ETF format
LIA_SPKTOOLS_API void outputResultETFLine(const String &source, const String &channel, const String &start, double duration, const String &type, const String &sub, const String &event, double score, const String &decision, ostream& strm);

//-------------------------------------------------------------------------
//Output a result line in MDTM format
LIA_SPKTOOLS_API void outputResultMDTMLine(const String &source, const String &channel, const String &start, double duration, const String &type, double conf, const String &sub, ostream& strm);
//-------------------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------------------
// Affect the fields values depending on the format
//-------------------------------------------------------------------------------------------------------

// Return LIA info fields
LIA_SPKTOOLS_API void setLIAInfoFields(unsigned long & genderField, unsigned long & nameField, unsigned long & decisionField, unsigned long & segField, unsigned long & scoreField);

// Return NIST info fields
LIA_SPKTOOLS_API void setNIST04InfoFields(unsigned long & genderField, unsigned long & nameField, unsigned long & decisionField, unsigned long & segField, unsigned long & scoreField);

#endif //!defined(ALIZE_IOFormat_h)

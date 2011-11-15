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

#if !defined(ALIZE_FileInfo_cpp)
#define ALIZE_FileInfo_cpp

#if defined(_WIN32)
#pragma warning( disable : 4127 4702)
#endif

#include "FileInfo.h"
#include "Exception.h"

using namespace alize;

//-------------------------------------------------------------------------
FileInfo::FileInfo(const FileName & f):Object(), _pFileStruct(NULL),
_fileName(f), _swap(false)
{
}

//-------------------------------------------------------------------------
bool FileInfo::isClosed() const 
{
  return _pFileStruct == NULL;
}

//-------------------------------------------------------------------------
bool FileInfo::isOpen() const 
{
  return _pFileStruct != NULL;
}

//-------------------------------------------------------------------------
void FileInfo::open(int mode)
{
  if (isOpen())
    close();
  if (_fileName.isEmpty())
    throw Exception("empty file name", __FILE__, __LINE__);
  if (mode == 1)
    _pFileStruct =::fopen(_fileName.c_str(), "wb");
  else if (mode == 2)
    _pFileStruct =::fopen(_fileName.c_str(), "rb");
  if (_pFileStruct == NULL)
    throw IOException("Cannot create new file", __FILE__, __LINE__,
    _fileName);
}

//-------------------------------------------------------------------------
void FileInfo::close()
{
  if (isOpen())
    if (::fclose(_pFileStruct) == EOF)
      throw IOException("Cannot close file", __FILE__, __LINE__, _fileName);
  _pFileStruct = NULL;
}

//-------------------------------------------------------------------------
void FileInfo::writeTopInfo(const LKVector & v, Config & c)
{
  if (isClosed())
    open(1);			// can throw Exception if file name = ""
  assert(_pFileStruct != NULL);
  for (int i = 0; i < c.getParam("topDistribsCount").toLong(); i++)
    if (::fwrite(&(v.getArray()[i].idx), sizeof(v.getArray()[i].idx), 1,
	_pFileStruct) != 1)
      throw IOException("Cannot write in file", __FILE__, __LINE__,
    _fileName);

  if (::fwrite(&v.sumNonTopDistribLK, sizeof(v.sumNonTopDistribLK), 1,
      _pFileStruct) != 1)
    throw IOException("Cannot write in file", __FILE__, __LINE__, _fileName);

  if (::fwrite(&v.sumNonTopDistribWeights, sizeof(v.sumNonTopDistribWeights),
      1, _pFileStruct) != 1)
    throw IOException("Cannot write in file", __FILE__, __LINE__, _fileName);

}

//-------------------------------------------------------------------------
void FileInfo::writeTopInfo(RealVector <real_t> & v, Config & c)
{
  if (isClosed())
    open(1);			// can throw Exception if file name = ""
  assert(_pFileStruct != NULL);

  for (unsigned long i = 0; i < v.size(); i++){
		if (i<c.getParam("topDistribsCount").toULong()){			
			unsigned long tmp=(unsigned long)v[i];
			if (::fwrite(&tmp, sizeof(tmp), 1,_pFileStruct) != 1)
			throw IOException("Cannot write in file", __FILE__, __LINE__, _fileName);
		}
		else {
			 if (::fwrite(&v[i], sizeof(v[i]), 1,_pFileStruct) != 1)
			 throw IOException("Cannot write in file", __FILE__, __LINE__, _fileName);
		}
		   
}

}

//-------------------------------------------------------------------------
void FileInfo::loadTopInfo(StatServer & ss, unsigned long &numLigne,
  Config & c)
{
  ULongVector index;
  double sumNonSelectedWeights;
  double sumNonSelectedLLK;
  unsigned long s;
  if (isClosed())
    {
      open(2);			// can throw Exception if file name = ""
    }

  assert(_pFileStruct != NULL);
  unsigned long pos =    numLigne * (c.getParam("topDistribsCount").toLong() *    sizeof(unsigned long) + 2 * sizeof(real_t));


  if (::fseek(_pFileStruct, pos, SEEK_SET) != 0)
    throw IOException("seek out of bounds", __FILE__, __LINE__, _fileName);

  for (int i = 0; i < c.getParam("topDistribsCount").toLong(); i++)
    {
      ::fread(&s, 1, sizeof(unsigned long), _pFileStruct);
      index.addValue(s);
    }

  ::fread(&sumNonSelectedLLK, 1, sizeof(double), _pFileStruct);

  ::fread(&sumNonSelectedWeights, 1, sizeof(double), _pFileStruct);
  ss.setTopDistribIndexVector(index, sumNonSelectedWeights, sumNonSelectedLLK);	//SET the top Component



}

//-------------------------------------------------------------------------
String FileInfo::toString() const 
{
  return Object::toString() + "\n  file name = '" + _fileName;
}

//-------------------------------------------------------------------------
String FileInfo::getClassName() const 
{
  return "FileInfo";
}

//-------------------------------------------------------------------------
FileInfo::~FileInfo()
{
  close();
}

//-------------------------------------------------------------------------

#endif // !defined(ALIZE_FileInfo_cpp)

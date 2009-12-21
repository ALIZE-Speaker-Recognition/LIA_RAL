/*
Alize is a free, open tool for speaker recognition

Alize is a development project initiated by the ELISA consortium
  [www.lia.univ-avignon.fr/heberges/ALIZE/ELISA] and funded by the
  French Research Ministry in the framework of the
  TECHNOLANGUE program [www.technolangue.net]
  [www.technolangue.net]

The Alize project team wants to highlight the limits of voice 
  authentication in a forensic context.
  The following paper proposes a good overview of this point:
  [Bonastre J.F., Bimbot F., Boe L.J., Campbell J.P., Douglas D.A., 
  Magrin-chagnolleau I., Person  Authentification by Voice: A Need of 
  Caution, Eurospeech 2003, Genova]
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
  uniquely characterize a person=92s voice or to identify with absolute 
  certainty an individual from his or her voice.]
  Contact Jean-Francois Bonastre for more information about the licence or
  the use of Alize

Copyright (C) 2003-2005
  Laboratoire d'informatique d'Avignon [www.lia.univ-avignon.fr]
  Frederic Wils [frederic.wils@lia.univ-avignon.fr]
  Jean-Francois Bonastre [jean-francois.bonastre@lia.univ-avignon.fr]
      
This file is part of Alize.

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
//Author : Alexandre PRETI.
#if !defined(ALIZE_FileInfo_cpp)
#define ALIZE_FileInfo_cpp

#ifdef WIN32
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

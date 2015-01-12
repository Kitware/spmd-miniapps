#ifndef __LoadImage3D_h
#define __LoadImage3D_h

#include "ConvertBuffer.h"
#include "Image3D.h"
#include "TypeInfo.h"

#include <algorithm>
#include <cassert>
#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

class bad_format : public std::exception
{
public:
  bad_format(const char* const message) : message(message) {}

private:
  const char* const message;
};

class LineStream
{
public:
  LineStream(std::istream &in);
  std::istream& stream();

  void readline();

private:
  std::istream &in;
  std::stringstream sstream;
  std::string line;
};

inline LineStream::LineStream(std::istream &in)
  : in(in)
{
}

inline std::istream& LineStream::stream()
{
  return this->sstream;
}

inline void LineStream::readline()
{
  std::getline(this->in, this->line);
  this->sstream.str(this->line);
  this->sstream.seekg(0);
}


template <typename T>
void loadImage3D(const char *vtkFileName, Image3D<T> *image)
{
  std::ifstream stream;
  stream.open(vtkFileName);

  LineStream reader(stream);

  // read and discard header
  reader.readline();

  reader.readline();
  std::string name;
  reader.stream() >> name;

  reader.readline();
  std::string format;
  reader.stream() >> format;
  if (format != "BINARY")
    {
    throw bad_format("Only 'BINARY' format supported");
    }

  std::string tag;

  reader.readline();
  std::string dataset;
  reader.stream() >> tag >> dataset;
  if (tag != "DATASET" || dataset != "STRUCTURED_POINTS")
    {
    throw bad_format("Expecting STRUCTURED_POINTS dataset");
    }

  int count = 3;
  int xdim, ydim, zdim;
  T spacing[3];
  T origin[3];
  while (count)
    {
    reader.readline();
    reader.stream() >> tag;

    if (tag == "DIMENSIONS")
      {
      reader.stream() >> xdim >> ydim >> zdim;
      if (reader.stream().bad())
        {
        throw bad_format("Expecting DIMENSIONS [3]");
        }
      --count;
      }
    else if (tag == "SPACING")
      {
      reader.stream() >> spacing[0] >> spacing[1] >> spacing[2];
      if (reader.stream().bad())
        {
        throw bad_format("Expecting SPACING [3]");
        }
      --count;
      }
    else if (tag == "ORIGIN")
      {
      reader.stream() >> origin[0] >> origin[1] >> origin[2];
      if (reader.stream().bad())
        {
        throw bad_format("Expecting ORIGIN [3]");
        }
      --count;
      }
    else
      {
      throw bad_format("Expecting DIMENSIONS, SPACING and ORIGIN");
      }
    }

  reader.readline();
  int npoints = 0;
  reader.stream() >> tag >> npoints;
  if (tag != "POINT_DATA" || reader.stream().bad())
    {
    throw bad_format("Expecting POINT_DATA <npoints>");
    }

  reader.readline();
  std::string scalName, typeName;
  reader.stream() >> tag >> scalName >> typeName;
  if (tag != "SCALARS" || reader.stream().bad())
    {
    throw bad_format("Expecting SCALARS <name> <type>");
    }
  TypeInfo ti = createTypeInfo(typeName.c_str());
  if (ti.getId() == TypeInfo::ID_UNKNOWN)
    {
    throw bad_format("Unsupported datatype");
    }

  reader.readline();
  reader.stream() >> tag >> name;
  if (tag != "LOOKUP_TABLE" || reader.stream().bad())
    {
    throw bad_format("Expecting LOOKUP_TABLE name");
    }
  if (name != "default")
    {
    int size;
    reader.stream() >> size;
    }

  int bufsize = npoints * ti.size();
  std::vector<char> rbuf(bufsize);
  stream.read(&rbuf[0], bufsize);

  image->setDimension(xdim, ydim, zdim);
  image->setSpacing(spacing[0], spacing[1], spacing[2]);
  image->setOrigin(origin[0], origin[1], origin[2]);
  image->allocate();

  convertBufferWithTypeInfo(&rbuf[0], ti, npoints, image->getData());

  stream.close();
}

#endif


#ifndef __LoadTetrahedronMesh_h
#define __LoadTetrahedronMesh_h

#include "ConvertBuffer.h"
#include "TetrahedronMesh.h"
#include "TypeInfo.h"

#include <fstream>
#include <iostream>
#include <sstream>

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
  this->sstream.clear();
  this->sstream.str(this->line);
  this->sstream.seekg(0);
}


template <typename T>
void loadTetrahedronMesh(const char *vtkFileName, TetrahedronMesh<T> *tetmesh)
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
  if (tag != "DATASET" || dataset != "UNSTRUCTURED_GRID")
    {
    throw bad_format("Expecting UNSTRUCTURED_GRID dataset");
    }

  std::vector<char> rbuf;
  size_t buffsize;

  std::vector<T> points, values;
  std::vector<unsigned> indexes;


  reader.readline();
  std::string type;
  unsigned npoints;
  reader.stream() >> tag >> npoints >> type;
  if (tag != "POINTS" || reader.stream().bad())
    {
    throw bad_format("Expecting POINTS <npoints> <type>");
    }
  TypeInfo ti = createTypeInfo(type.c_str());
  if (ti.getId() == TypeInfo::ID_UNKNOWN)
    {
    throw bad_format("Unsupported datatype");
    }

  buffsize = npoints * 3 * ti.size();
  rbuf.resize(buffsize);
  stream.read(&rbuf[0], rbuf.size());
  stream >> std::ws;

  points.resize(npoints * 3);
  convertBufferWithTypeInfo(&rbuf[0], ti, npoints * 3, &points[0]);


  reader.readline();
  unsigned ncells, nints;
  reader.stream() >> tag >> ncells >> nints;
  if (tag != "CELLS" || reader.stream().bad())
    {
    throw bad_format("Expecting CELLS <ncells> <ninds>");
    }
  
  buffsize = nints * sizeof(int);
  rbuf.resize(buffsize);
  stream.read(&rbuf[0], buffsize);
  stream >> std::ws;

  indexes.resize(nints);
  convertBufferWithTypeInfo(&rbuf[0], createTypeInfo(TypeInfo::ID_INT), nints,
                            &indexes[0]);
  unsigned ninds = 0;
  for (unsigned i = 0; i < nints; i += 5)
    {
    if (indexes[i] != 4)
      {
      throw bad_format("Expecting all cells to be tetrahedrons");
      }
    for (int ii = 1; ii < 5; ++ii)
      {
      indexes[ninds++] = indexes[i + ii];
      }
    }
  indexes.resize(ninds);


  reader.readline();
  reader.stream() >> tag >> ncells;
  if (tag != "CELL_TYPES" || reader.stream().bad())
    {
    throw bad_format("Expecting CELL_TYPES <ncells>");
    }
  buffsize = ncells * sizeof(unsigned);
  stream.seekg(buffsize, std::istream::cur); // assumes all are tetrahedra
  stream >> std::ws;


  reader.readline();
  reader.stream() >> tag >> npoints;
  if (tag != "POINT_DATA" || reader.stream().bad())
    {
    throw bad_format("Expecting POINT_DATA <npoints>");
    }
  reader.readline();
  reader.stream() >> tag >> name >> type;
  if (tag != "SCALARS" || reader.stream().bad())
    {
    throw bad_format("Expecting SCALARS <name> <type>");
    }
  ti = createTypeInfo(type.c_str());
  if (ti.getId() == TypeInfo::ID_UNKNOWN)
    {
    throw bad_format("Unsupported datatype");
    }
  reader.readline();
  reader.stream() >> tag >> name;
  if (tag != "LOOKUP_TABLE" || reader.stream().bad())
    {
    throw bad_format("Expecting LOOKUP_TABLE <name>");
    }
  
  buffsize = npoints * ti.size();
  rbuf.resize(buffsize);
  stream.read(&rbuf[0], rbuf.size());
  stream >> std::ws;

  values.resize(npoints);
  convertBufferWithTypeInfo(&rbuf[0], ti, npoints, &values[0]);

  stream.close();

  tetmesh->points.swap(points);
  tetmesh->values.swap(values);
  tetmesh->indexes.swap(indexes);
}

#endif


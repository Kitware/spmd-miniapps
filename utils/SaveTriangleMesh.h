#ifndef __SaveTriangleMesh_h
#define __SaveTriangleMesh_h

#include "ConvertBuffer.h"
#include "TriangleMesh.h"
#include "TypeInfo.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>


template <typename T>
void saveTriangleMesh(const TriangleMesh<T> &mesh, const char *vtkFileName)
{
  std::ofstream stream;
  stream.open(vtkFileName);

  unsigned nverts = mesh.points.size()/3;
  unsigned ntriangles = mesh.indexes.size()/3;
  TypeInfo ti = createTypeInfo<T>();

  size_t bufsize = 0;
  std::vector<char> wbuff;

  stream << "# vtk DataFile Version 3.0" << std::endl;
  stream << "Isosurface Mesh" << std::endl;
  stream << "BINARY" << std::endl;
  stream << "DATASET POLYDATA" << std::endl;

  bufsize = mesh.points.size() * sizeof(T);
  wbuff.resize(bufsize);
  convertBuffer(&mesh.points[0], mesh.points.size(),
                reinterpret_cast<T*>(&wbuff[0]));

  stream << "POINTS " << nverts << " " << ti.name() << std::endl;
  stream.write(&wbuff[0], wbuff.size());
  stream << std::endl;

  bufsize = ntriangles * 4 * sizeof(unsigned);
  wbuff.resize(bufsize);
  unsigned *ind = reinterpret_cast<unsigned*>(&wbuff[0]);
  const unsigned *src = &mesh.indexes[0];
  for (unsigned i = 0; i < ntriangles; ++i)
    {
    *ind = 3;
    flipEndianness(*ind++);
    for (int j = 0; j < 3; ++j)
      {
      *ind = *src++;
      flipEndianness(*ind++);
      }
    }

  stream << "POLYGONS " << ntriangles << " " << ntriangles * 4 << std::endl;
  stream.write(&wbuff[0], wbuff.size());
  stream << std::endl;

  bufsize = mesh.normals.size() * sizeof(T);
  wbuff.resize(bufsize);
  convertBuffer(&mesh.normals[0], mesh.normals.size(), 
                reinterpret_cast<T*>(&wbuff[0]));

  stream << "POINT_DATA " << nverts << std::endl;
  stream << "NORMALS Normals " << ti.name() << std::endl;
  stream.write(&wbuff[0], wbuff.size());
  stream << std::endl;

  stream.close();
}

#endif


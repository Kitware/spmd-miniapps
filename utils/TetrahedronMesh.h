#ifndef __TetrahedronMesh_h
#define __TetrahedronMesh_h

#include <vector>

template <typename T>
struct TetrahedronMesh
{
public:
  std::vector<T> points;
  std::vector<T> values;
  std::vector<unsigned> indexes;

  unsigned numberOfPoints() const
  {
    return this->values.size();
  }

  unsigned numberOfTetrahedra() const
  {
    return this->indexes.size()/4;
  }
};

#endif


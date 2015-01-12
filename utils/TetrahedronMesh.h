#ifndef __TetrahedronMesh_h
#define __TetrahedronMesh_h

#include <vector>

template <typename T>
struct TetrahedronMesh
{
public:
  std::vector<T> points;
  std::vector<T> values;
  std::vector<int> indexes;

  int numberOfPoints() const
  {
    return this->values.size();
  }

  int numberOfTetrahedra() const
  {
    return this->indexes.size()/4;
  }
};

#endif


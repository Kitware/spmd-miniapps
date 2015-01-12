#ifndef __TriangleMesh_h
#define __TriangleMesh_h

#include <vector>

template <typename T>
struct TriangleMesh
{
public:
  std::vector<T> points;
  std::vector<T> normals;
  std::vector<int> indexes;

  int numberOfVertices() const;
  int numberOfTriangles() const;
};

template <typename T>
inline int TriangleMesh<T>::numberOfVertices() const
{
  return points.size()/3;
}

template <typename T>
inline int TriangleMesh<T>::numberOfTriangles() const
{
  return indexes.size()/3;
}

#endif


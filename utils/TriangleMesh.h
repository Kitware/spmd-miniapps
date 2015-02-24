#ifndef __TriangleMesh_h
#define __TriangleMesh_h

#include <vector>

using std::size_t;

template <typename T>
struct TriangleMesh
{
public:
  std::vector<T> points;
  std::vector<T> normals;
  std::vector<unsigned> indexes;

  unsigned numberOfVertices() const;
  unsigned numberOfTriangles() const;
};

template <typename T>
inline unsigned TriangleMesh<T>::numberOfVertices() const
{
  return points.size()/3;
}

template <typename T>
inline unsigned TriangleMesh<T>::numberOfTriangles() const
{
  return indexes.size()/3;
}


template <typename T, typename MeshListIterator>
void mergeTriangleMeshes(MeshListIterator from, MeshListIterator to,
                         TriangleMesh<T> *outMesh)
{
  size_t ptsSize = 0, indsSize = 0;
  for (MeshListIterator i = from; i != to; ++i)
    {
    ptsSize += i->points.size();
    indsSize += i->indexes.size();
    }

  outMesh->points.reserve(ptsSize);
  outMesh->normals.reserve(ptsSize);
  outMesh->indexes.reserve(indsSize);

  for (MeshListIterator i = from; i != to; ++i)
    {
    unsigned lastInd = outMesh->points.size()/3;

    outMesh->points.insert(outMesh->points.end(), i->points.begin(),
                           i->points.end());
    outMesh->normals.insert(outMesh->normals.end(), i->normals.begin(),
                            i->normals.end());

    for (size_t j = 0; j < i->indexes.size(); ++j)
      {
      outMesh->indexes.push_back(i->indexes[j] + lastInd);
      }
    }
}

#endif


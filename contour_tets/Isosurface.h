#ifndef __Isosurface_h
#define __Isosurface_h

#include "type.h"

namespace scalar
{
void extractIsosurface(const TetrahedronMesh_t &tetmesh, Float_t isoval,
                       TriangleMesh_t *trimesh);
};

namespace simd
{
void extractIsosurface(const TetrahedronMesh_t &tetmesh, Float_t isoval,
                       TriangleMesh_t *trimesh);
};

namespace simd_2
{
void extractIsosurface(const TetrahedronMesh_t &tetmesh, Float_t isoval,
                       TriangleMesh_t *trimesh);
};


inline unsigned estimatedNumberOfActiveTets(unsigned numTets)
{
  // assume 20% tests produce triangles;
  return numTets/5;
}

inline unsigned estimatedNumberOfTriangles(unsigned numTets)
{
  // assume active tests produce an average of 3/2 triangles
  return (estimatedNumberOfActiveTets(numTets) * 3)/2;
}

inline unsigned estimatedNumberOfPoints(unsigned numTets)
{
  // assume each edge is shared by around 3 tets
  return (estimatedNumberOfTriangles(numTets) * 3)/3;
}

#endif


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

#endif


#ifndef __MarchingCubes_h
#define __MarchingCubes_h

#include "type.h"

namespace scalar {
void extractIsosurface(const Image3D_t &vol, Float_t isoval,
                       TriangleMesh_t *mesh);
};

namespace shortvec {
void extractIsosurface(const Image3D_t &vol, Float_t isoval,
                       TriangleMesh_t *mesh);
};

namespace simd {
void extractIsosurface(const Image3D_t &vol, Float_t isoval,
                       TriangleMesh_t *mesh);
};

namespace scalar_2 {
void extractIsosurface(const Image3D_t &vol, Float_t isoval,
                       TriangleMesh_t *mesh);
};

namespace simd_2 {
void extractIsosurface(const Image3D_t &vol, Float_t isoval,
                       TriangleMesh_t *mesh);
};

namespace scalar_2_1 {
void extractIsosurface(const Image3D_t &vol, Float_t isoval,
                       TriangleMesh_t *mesh);
};

namespace simd_2_1 {
void extractIsosurface(const Image3D_t &vol, Float_t isoval,
                       TriangleMesh_t *mesh);
};


inline unsigned estimatedNumberOfActiveCells(const unsigned indim[3])
{
  // assume about 30% of cells produce triangles.
  return (indim[0] * indim[1] * indim[2])/3;
}

inline unsigned estimatedNumberOfTriangles(const unsigned indim[3])
{
  // each active cell has about 5/2 triangles
  return (estimatedNumberOfActiveCells(indim) * 5)/2;
}

inline unsigned estimatedNumberOfPoints(const unsigned indim[3])
{
  // 3 vertices per triangle, each edge is shared by upto 4 cells
  return (estimatedNumberOfTriangles(indim) * 3)/4;
}

#endif


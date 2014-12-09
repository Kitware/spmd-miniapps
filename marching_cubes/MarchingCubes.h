#ifndef __MarchingCubes_h
#define __MarchingCubes_h

struct TriangleMesh3D;

namespace serial {
void extractIsosurface(const Image3D &vol, double isoval,
                       TriangleMesh3D *mesh);
};

namespace simd {
void extractIsosurface(const Image3D &vol, double isoval,
                       TriangleMesh3D *mesh);
};

namespace ispc {
void extractIsosurface(const Image3D &vol, double isoval,
                       TriangleMesh3D *mesh);
};

#endif


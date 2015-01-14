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

#endif


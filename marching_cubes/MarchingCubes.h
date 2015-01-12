#ifndef __MarchingCubes_h
#define __MarchingCubes_h

#include "type.h"

namespace scalar {
void extractIsosurface(const Image3D_t &vol, Float_t isoval,
                       TriangleMesh_t *mesh);
};

namespace mixed {
void extractIsosurface(const Image3D_t &vol, Float_t isoval,
                       TriangleMesh_t *mesh);
};

namespace vectorized {
void extractIsosurface(const Image3D_t &vol, Float_t isoval,
                       TriangleMesh_t *mesh);
};

namespace scalar_2 {
void extractIsosurface(const Image3D_t &vol, Float_t isoval,
                       TriangleMesh_t *mesh);
};

namespace scalar_3 {
void extractIsosurface(const Image3D_t &vol, Float_t isoval,
                       TriangleMesh_t *mesh);
};

namespace mixed_3 {
void extractIsosurface(const Image3D_t &vol, Float_t isoval,
                       TriangleMesh_t *mesh);
};

namespace scalar_3_1 {
void extractIsosurface(const Image3D_t &vol, Float_t isoval,
                       TriangleMesh_t *mesh);
};

namespace mixed_3_1 {
void extractIsosurface(const Image3D_t &vol, Float_t isoval,
                       TriangleMesh_t *mesh);
};

#endif


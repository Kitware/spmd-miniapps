#ifndef __type_h
#define __type_h

#if !ISPC
# include <Image3D.h>
# include <TriangleMesh.h>
#endif

typedef double Float_t;

#if !ISPC
typedef Image3D<Float_t> Image3D_t;
typedef TriangleMesh<Float_t> TriangleMesh_t;
#endif

#endif


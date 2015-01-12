#ifndef __type_h
#define __type_h

#if !ISPC
# include <TetrahedronMesh.h>
# include <TriangleMesh.h>
#endif

typedef double Float_t;

#if !ISPC
typedef TetrahedronMesh<Float_t> TetrahedronMesh_t;
typedef TriangleMesh<Float_t> TriangleMesh_t;
#endif

#endif


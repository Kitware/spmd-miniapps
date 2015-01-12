
namespace vectorized_2
{

void extractIsosurface(const TetrahedronMesh_t &tetmesh, Float_t isoval,
                       TriangleMesh_t *trimesh)
{
  const Float_t *tetPoints = &tetmesh.points[0];
  const Float_t *tetScalars = &tetmesh.values[0];
  const Float_t *tetIndexes = &tetmesh.indexes[0];

  
}

}


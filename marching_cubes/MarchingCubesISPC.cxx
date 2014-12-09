#include "Image3D.h"
#include "MarchingCubes.h"
#include "MarchingCubesISPC.h"
#include "MarchingCubesTables.h"
#include "PointLocator3D.h"
#include "TriangleMesh3D.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

namespace ispc {

struct VertexId
{
  int idx, id;
};

class GetPosition
{
public:
  GetPosition(const std::vector<Vertex> &vertices)
    : vertices(vertices) {}

  const double* operator()(const VertexId &v) const
  {
    return vertices[v.idx].pos;
  }

private:
  const std::vector<Vertex> &vertices;
};

void extractIsosurface(const Image3D &vol, double isoval, TriangleMesh3D *mesh)
{
  const int *dims = vol.getDimension();
  const double *origin = vol.getOrigin();
  const double *spacing = vol.getSpacing();
  const double *buffer = vol.getData();

  int totalCells = (dims[0] - 1) * (dims[1] - 1) * (dims[2] - 1);
  std::vector<int> cases(totalCells);
  int activeCells = ispc::getCellCases(buffer, dims, isoval, &cases[0]);

  int ntris = 0;
  std::vector<int> cellIndex(activeCells), vertInds(activeCells);
  for (int s = 0, d = 0; d < totalCells; ++d)
    {
    if (cases[d] != 0 && cases[d] != 255)
      {
      cellIndex[s] = d;
      vertInds[s] = ntris * 3;
      cases[s++] = cases[d];
      ntris += MarchingCubesTables::getTrianglesPerCase()[cases[d]];
      }
    }
  cases.resize(activeCells);

  int nverts = ntris * 3;
  std::vector<Vertex> verts(nverts);

  ispc::extractIsosurface_impl(buffer, dims, origin, spacing, isoval,
                               MarchingCubesTables::getTriangleCases()[0],
                               activeCells, &(cellIndex[0]), &(cases[0]),
                               &(vertInds[0]), &(verts[0]));

  // remove duplicates
  double range[6];
  for (int i = 0; i < 3; ++i)
    {
    range[i * 2] = origin[i];
    range[(i * 2) + 1] = origin[i] +
                         (static_cast<double>(dims[i] - 1) * spacing[i]);
    }
  PointLocator3D<VertexId, GetPosition> pointLocator(range, 32, 32, 16);
  GetPosition getPosition(verts);

  mesh->indexes.reserve(nverts);
  int tri[3];
  for (int i = 0, j = 0; i < nverts; ++i, ++j)
    {
    bool exists = false;

    VertexId put;
    put.idx = i;
    put.id = mesh->vertices.size();

    VertexId *get = pointLocator.insert(put, &exists, getPosition);
    if (!exists)
      {
        TriangleMesh3D::Vertex v;
        std::copy(verts[i].pos, verts[i].pos + 3, v.pos);
        std::copy(verts[i].norm, verts[i].norm + 3, v.norm);
        mesh->vertices.push_back(v);
      }
    tri[j] = get->id;

    if (j == 2)
      {
        if (tri[0] != tri[1] && tri[1] != tri[2] && tri[2] != tri[0])
          {
          mesh->indexes.push_back(tri[0]);
          mesh->indexes.push_back(tri[1]);
          mesh->indexes.push_back(tri[2]);
          }
        j = -1;
      }
    }
}

}; // namespace ispc


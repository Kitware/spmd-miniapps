#include "MarchingCubes.h"
#include "MarchingCubes.ispc.h"
#include "MarchingCubesShortVec.ispc.h"
#include "MarchingCubesTables.h"

#include <PointLocator3D.h>

#include <iostream>

struct Position
{
  Float_t xyz[3];
  int idx;

  const Float_t* getPosition() const
  {
    return xyz;
  }
};


namespace scalar {

static void computeGradient(int xidx, int yidx, int zidx,
                            const Float_t *buffer, const int dims[3],
                            const Float_t spacing[3], Float_t grad[3])
{
  int xysize = dims[0] * dims[1];
  int ptidx = xidx + yidx * dims[0] + zidx * xysize;

  if (xidx == 0)
    {
    Float_t x1 = buffer[ptidx + 1];
    Float_t x2 = buffer[ptidx];
    grad[0] = (x2 - x1)/spacing[0];
    }
  else if (xidx == (dims[0] - 1))
    {
    Float_t x1 = buffer[ptidx];
    Float_t x2 = buffer[ptidx - 1];
    grad[0] = (x2 - x1)/spacing[0];
    }
  else
    {
    Float_t x1 = buffer[ptidx + 1];
    Float_t x2 = buffer[ptidx - 1];
    grad[0] = (0.5 * (x2 - x1))/spacing[0];
    }

  if (yidx == 0)
    {
    Float_t y1 = buffer[ptidx + dims[0]];
    Float_t y2 = buffer[ptidx];
    grad[1] = (y2 - y1)/spacing[1];
    }
  else if (yidx == (dims[1] - 1))
    {
    Float_t y1 = buffer[ptidx];
    Float_t y2 = buffer[ptidx - dims[0]];
    grad[1] = (y2 - y1)/spacing[1];
    }
  else
    {
    Float_t y1 = buffer[ptidx + dims[0]];
    Float_t y2 = buffer[ptidx - dims[0]];
    grad[1] = (0.5 * (y2 - y1))/spacing[1];
    }

  if (zidx == 0)
    {
    Float_t z1 = buffer[ptidx + xysize];
    Float_t z2 = buffer[ptidx];
    grad[2] = (z2 - z1)/spacing[2];
    }
  else if (zidx == (dims[2] - 1))
    {
    Float_t z1 = buffer[ptidx];
    Float_t z2 = buffer[ptidx - xysize];
    grad[2] = (z2 - z1)/spacing[2];
    }
  else
    {
    Float_t z1 = buffer[ptidx + xysize];
    Float_t z2 = buffer[ptidx - xysize];
    grad[2] = (0.5 * (z2 - z1))/spacing[2];
    }
}


static inline Float_t lerp(Float_t a, Float_t b, Float_t w)
{
  //return ((1.0 - w) * a) + (w * b);
  return a + (w * (b - a));
}

void extractIsosurface(const Image3D_t &vol, Float_t isoval,
                       TriangleMesh_t *mesh)
{
  static const int caseMask[] = { 1, 2, 4, 8, 16, 32, 64, 128 };

  const int *dims = vol.getDimension();
  const Float_t *origin = vol.getOrigin();
  const Float_t *spacing = vol.getSpacing();
  const Float_t *buffer = vol.getData();

  int sliceSize = dims[0] * dims[1];
  Float_t range[6];

  // initialize point locator data-structure
  for (int i = 0; i < 3; ++i)
    {
    range[i * 2] = origin[i];
    range[(i * 2) + 1] = origin[i] +
                         (static_cast<Float_t>(dims[i] - 1) * spacing[i]);
    }
  PointLocator3D<Float_t, Position> pointLocator(range, 32, 32, 16);

  int ptIdx = 0;

  // march through each cell
  Float_t zpos = origin[2];
  for (int zidx = 0; zidx < (dims[2] - 1); ++zidx, zpos += spacing[2])
    {
    Float_t ypos = origin[1];
    for (int yidx = 0; yidx < (dims[1] - 1); ++yidx, ypos += spacing[1])
      {
      Float_t xpos = origin[0];
      for (int xidx = 0; xidx < (dims[0] - 1); ++xidx, xpos += spacing[0])
        {
        Float_t pos[8][3], grad[8][3];
        Float_t val[8];

        // get cell-points values
        int idx = xidx + (yidx * dims[0]) + (zidx * sliceSize);
        val[0] = buffer[idx];
        val[1] = buffer[idx + 1];
        val[2] = buffer[idx + dims[0] + 1];
        val[3] = buffer[idx + dims[0]];
        val[4] = buffer[idx + sliceSize];
        val[5] = buffer[idx + 1 + sliceSize];
        val[6] = buffer[idx + dims[0] + 1 + sliceSize];
        val[7] = buffer[idx + dims[0] + sliceSize];

        int cellId = 0;
        for (int i = 0; i < 8; ++i)
          {
          cellId |= (val[i] >= isoval) ? caseMask[i] : 0;
          }

        // no intersections
        if (cellId == 0 || cellId == 255)
          {
          continue;
          }

        // get physical position and gradient of the points
        pos[0][0] = xpos;
        pos[0][1] = ypos;
        pos[0][2] = zpos;

        pos[1][0] = xpos + spacing[0];
        pos[1][1] = ypos;
        pos[1][2] = zpos;

        pos[2][0] = xpos + spacing[0];
        pos[2][1] = ypos + spacing[1];
        pos[2][2] = zpos;

        pos[3][0] = xpos;
        pos[3][1] = ypos + spacing[1];
        pos[3][2] = zpos;

        pos[4][0] = xpos;
        pos[4][1] = ypos;
        pos[4][2] = zpos + spacing[2];

        pos[5][0] = xpos + spacing[0];
        pos[5][1] = ypos;
        pos[5][2] = zpos + spacing[2];

        pos[6][0] = xpos + spacing[0];
        pos[6][1] = ypos + spacing[1];
        pos[6][2] = zpos + spacing[2];

        pos[7][0] = xpos;
        pos[7][1] = ypos + spacing[1];
        pos[7][2] = zpos + spacing[2];

        computeGradient(xidx, yidx, zidx, buffer, dims, spacing, grad[0]);
        computeGradient(xidx + 1, yidx, zidx, buffer, dims, spacing,
                        grad[1]);
        computeGradient(xidx + 1, yidx + 1, zidx, buffer, dims, spacing,
                        grad[2]);
        computeGradient(xidx, yidx + 1, zidx, buffer, dims, spacing,
                        grad[3]);
        computeGradient(xidx, yidx, zidx + 1, buffer, dims, spacing,
                        grad[4]);
        computeGradient(xidx + 1, yidx, zidx + 1, buffer, dims, spacing,
                        grad[5]);
        computeGradient(xidx + 1, yidx + 1, zidx + 1, buffer, dims, spacing,
                        grad[6]);
        computeGradient(xidx, yidx + 1, zidx + 1, buffer, dims, spacing,
                        grad[7]);

        // get the triangles to generate
        const int *edges = MarchingCubesTables::getCaseTrianglesEdges(cellId);
        for (; *edges != -1; edges += 3)
          {
          int tri[3];
          for (int i = 0; i < 3; ++i)
            {
            int v1 = MarchingCubesTables::getEdgeVertices(edges[i])[0];
            int v2 = MarchingCubesTables::getEdgeVertices(edges[i])[1];
            Float_t w = (isoval - val[v1])/(val[v2] - val[v1]);

            // interpolate vertex position
            Position pt;
            pt.xyz[0] = lerp(pos[v1][0], pos[v2][0], w);
            pt.xyz[1] = lerp(pos[v1][1], pos[v2][1], w);
            pt.xyz[2] = lerp(pos[v1][2], pos[v2][2], w);
            pt.idx = ptIdx;

            bool exists = false;
            Position *pt_ = pointLocator.insert(pt, &exists);
            if (!exists)
              {
              Float_t norm[3];

              // interpolate vertex normal
              norm[0] = lerp(grad[v1][0], grad[v2][0], w);
              norm[1] = lerp(grad[v1][1], grad[v2][1], w);
              norm[2] = lerp(grad[v1][2], grad[v2][2], w);

              for (int ii = 0; ii < 3; ++ii)
                {
                mesh->points.push_back(pt.xyz[ii]);
                mesh->normals.push_back(norm[ii]);
                }
              tri[i] = ptIdx++;
              }
            else
              {
              tri[i] = pt_->idx;
              }
            }

          if (tri[0] == tri[1] || tri[1] == tri[2] || tri[2] == tri[0])
            {
            //std::cout << "dgenerate triangle" << std::endl;
            }
          else
            {
            mesh->indexes.push_back(tri[0]);
            mesh->indexes.push_back(tri[1]);
            mesh->indexes.push_back(tri[2]);
            }
          }
        }
      }
    }
}

}; //namespace scalar

namespace shortvec {

void extractIsosurface(const Image3D_t &vol, Float_t isoval,
                       TriangleMesh_t *mesh)
{
  static const int caseMask[] = { 1, 2, 4, 8, 16, 32, 64, 128 };

  const int *dims = vol.getDimension();
  const Float_t *origin = vol.getOrigin();
  const Float_t *spacing = vol.getSpacing();
  const Float_t *buffer = vol.getData();

  int sliceSize = dims[0] * dims[1];
  Float_t range[6];

  // initialize point locator data-structure
  for (int i = 0; i < 3; ++i)
    {
    range[i * 2] = origin[i];
    range[(i * 2) + 1] = origin[i] +
                         (static_cast<Float_t>(dims[i] - 1) * spacing[i]);
    }
  PointLocator3D<Float_t, Position> pointLocator(range, 32, 32, 16);

  int ptIdx = 0;

  // march through each cell
  Float_t zpos = origin[2];
  for (int zidx = 0; zidx < (dims[2] - 1); ++zidx, zpos += spacing[2])
    {
    Float_t ypos = origin[1];
    for (int yidx = 0; yidx < (dims[1] - 1); ++yidx, ypos += spacing[1])
      {
      Float_t xpos = origin[0];
      for (int xidx = 0; xidx < (dims[0] - 1); ++xidx, xpos += spacing[0])
        {
        Float_t pos[8][3], grad[8][3];
        Float_t val[8];

        // get cell-points values
        int idx = xidx + (yidx * dims[0]) + (zidx * sliceSize);
        val[0] = buffer[idx];
        val[1] = buffer[idx + 1];
        val[2] = buffer[idx + dims[0] + 1];
        val[3] = buffer[idx + dims[0]];
        val[4] = buffer[idx + sliceSize];
        val[5] = buffer[idx + 1 + sliceSize];
        val[6] = buffer[idx + dims[0] + 1 + sliceSize];
        val[7] = buffer[idx + dims[0] + sliceSize];

        int cellId = 0;
        for (int i = 0; i < 8; ++i)
          {
          cellId |= (val[i] >= isoval) ? caseMask[i] : 0;
          }

        // no intersections
        if (cellId == 0 || cellId == 255)
          {
          continue;
          }

        // get physical position and gradient of the points
        pos[0][0] = xpos;
        pos[0][1] = ypos;
        pos[0][2] = zpos;

        pos[1][0] = xpos + spacing[0];
        pos[1][1] = ypos;
        pos[1][2] = zpos;

        pos[2][0] = xpos + spacing[0];
        pos[2][1] = ypos + spacing[1];
        pos[2][2] = zpos;

        pos[3][0] = xpos;
        pos[3][1] = ypos + spacing[1];
        pos[3][2] = zpos;

        pos[4][0] = xpos;
        pos[4][1] = ypos;
        pos[4][2] = zpos + spacing[2];

        pos[5][0] = xpos + spacing[0];
        pos[5][1] = ypos;
        pos[5][2] = zpos + spacing[2];

        pos[6][0] = xpos + spacing[0];
        pos[6][1] = ypos + spacing[1];
        pos[6][2] = zpos + spacing[2];

        pos[7][0] = xpos;
        pos[7][1] = ypos + spacing[1];
        pos[7][2] = zpos + spacing[2];

        int xinds[8] = { xidx, xidx + 1, xidx + 1, xidx,
                         xidx, xidx + 1, xidx + 1, xidx };
        int yinds[8] = { yidx, yidx, yidx + 1, yidx + 1,
                         yidx, yidx, yidx + 1, yidx + 1 };
        int zinds[8] = { zidx, zidx, zidx, zidx,
                         zidx + 1, zidx + 1, zidx + 1, zidx + 1 };

        ispc::computeGradient(xinds, yinds, zinds, buffer, dims, spacing,
                              grad);

        // get the triangles to generate
        const int *edges = MarchingCubesTables::getCaseTrianglesEdges(cellId);
        for (; *edges != -1; edges += 3)
          {
          int tri[3];
          for (int i = 0; i < 3; ++i)
            {
            int v1 = MarchingCubesTables::getEdgeVertices(edges[i])[0];
            int v2 = MarchingCubesTables::getEdgeVertices(edges[i])[1];
            Float_t w = (isoval - val[v1])/(val[v2] - val[v1]);

            // interpolate vertex position
            Position pt;
            ispc::lerp(3, pos[v1], pos[v2], w, pt.xyz);
            pt.idx = ptIdx;

            bool exists = false;
            Position *pt_ = pointLocator.insert(pt, &exists);
            if (!exists)
              {
              Float_t norm[3];

              // interpolate vertex normal
              ispc::lerp(3, grad[v1], grad[v2], w, norm);
              for (int ii = 0; ii < 3; ++ii)
                {
                mesh->points.push_back(pt.xyz[ii]);
                mesh->normals.push_back(norm[ii]);
                }
              tri[i] = ptIdx++;
              }
            else
              {
              tri[i] = pt_->idx;
              }
            }

            if (tri[0] == tri[1] || tri[1] == tri[2] || tri[2] == tri[0])
              {
                //std::cout << "degenerate triangle" << std::endl;
              }
            else
              {
                mesh->indexes.push_back(tri[0]);
                mesh->indexes.push_back(tri[1]);
                mesh->indexes.push_back(tri[2]);
              }
          }
        }
      }
    }
}

}; // namespace shortvec

namespace simd {

struct VertexId
{
  int idx, id;
};

class GetPosition
{
public:
  GetPosition(const std::vector<Float_t> &points)
    : points(points) {}

  const Float_t* operator()(const VertexId &v) const
  {
    return &points[v.idx * 3];
  }

private:
  const std::vector<Float_t> &points;
};

void extractIsosurface(const Image3D_t &vol, Float_t isoval,
                       TriangleMesh_t *mesh)
{
  const int *dims = vol.getDimension();
  const Float_t *origin = vol.getOrigin();
  const Float_t *spacing = vol.getSpacing();
  const Float_t *buffer = vol.getData();

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
      vertInds[s] = ntris * 3 * 3;
      cases[s++] = cases[d];
      ntris += MarchingCubesTables::getNumberOfTriangles(cases[d]);
      }
    }
  cases.resize(activeCells);

  int nverts = ntris * 3;
  std::vector<Float_t> points(nverts * 3), normals(nverts * 3);

  ispc::extractIsosurface_impl(buffer, dims, origin, spacing, isoval,
                               MarchingCubesTables::getCaseTrianglesEdges(0),
                               MarchingCubesTables::getEdgeVertices(0),
                               activeCells, &(cellIndex[0]), &(cases[0]),
                               &(vertInds[0]), &(points[0]), &normals[0]);

  // remove duplicates
  Float_t range[6];
  for (int i = 0; i < 3; ++i)
    {
    range[i * 2] = origin[i];
    range[(i * 2) + 1] = origin[i] +
                         (static_cast<Float_t>(dims[i] - 1) * spacing[i]);
    }
  PointLocator3D<Float_t, VertexId, GetPosition>
    pointLocator(range, 32, 32, 16);
  GetPosition getPosition(points);

  mesh->indexes.reserve(nverts);
  int tri[3];
  for (int i = 0, j = 0, ptid = 0; i < nverts; ++i, ++j)
    {
    bool exists = false;

    VertexId put;
    put.idx = i;
    put.id = ptid;

    VertexId *get = pointLocator.insert(put, &exists, getPosition);
    if (!exists)
      {
      for (int ii = 0; ii < 3; ++ii)
        {
        mesh->points.push_back(points[(i * 3) + ii]);
        mesh->normals.push_back(normals[(i * 3) + ii]);
        ++ptid;
        }
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

}; // namespace simd


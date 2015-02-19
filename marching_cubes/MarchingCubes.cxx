#include "MarchingCubes.h"
#include "MarchingCubes.ispc.h"
#include "MarchingCubesShortVec.ispc.h"
#include "MarchingCubesTables.h"

#include <PointLocator3D.h>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range3d.h>

#include <omp.h>

#include <iostream>

class GetPointPosition
{
public:
  GetPointPosition(const std::vector<Float_t> &points) : points(points)
  {
  }

  const Float_t* operator()(int idx) const
  {
    return &points[idx * 3];
  }

private:
  const std::vector<Float_t> &points;
};

typedef PointLocator3D<Float_t, int, GetPointPosition> PointLocator_t;


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


void extractIsosurfaceFromBlock(const Image3D_t &vol, const int ext[6],
  Float_t isoval, PointLocator_t &pointLocator, TriangleMesh_t *mesh)
{
  static const int caseMask[] = { 1, 2, 4, 8, 16, 32, 64, 128 };

  const int *dims = vol.getDimension();
  const Float_t *origin = vol.getOrigin();
  const Float_t *spacing = vol.getSpacing();
  const Float_t *buffer = vol.getData();

  int sliceSize = dims[0] * dims[1];

  int ptIdx = pointLocator.numberOfPoints();
  GetPointPosition getPosition(mesh->points);

  // march through each cell
  Float_t zpos = origin[2] + (Float_t(ext[4]) * spacing[2]);
  for (int zidx = ext[4]; zidx <= ext[5]; ++zidx, zpos += spacing[2])
    {
    Float_t ypos = origin[1] + (Float_t(ext[2]) * spacing[1]);
    for (int yidx = ext[2]; yidx <= ext[3]; ++yidx, ypos += spacing[1])
      {
      Float_t xpos = origin[0] + (Float_t(ext[0]) * spacing[0]);
      for (int xidx = ext[0]; xidx <= ext[1]; ++xidx, xpos += spacing[0])
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

            // interpolate vertex position and temporarily add
            mesh->points.push_back(lerp(pos[v1][0], pos[v2][0], w));
            mesh->points.push_back(lerp(pos[v1][1], pos[v2][1], w));
            mesh->points.push_back(lerp(pos[v1][2], pos[v2][2], w));

            bool exists = false;
            int eid = *pointLocator.insert(ptIdx, &exists, getPosition);
            if (!exists)
              {
              Float_t norm[3];

              // interpolate vertex normal
              mesh->normals.push_back(lerp(grad[v1][0], grad[v2][0], w));
              mesh->normals.push_back(lerp(grad[v1][1], grad[v2][1], w));
              mesh->normals.push_back(lerp(grad[v1][2], grad[v2][2], w));

              tri[i] = ptIdx++;
              }
            else
              {
              // duplicate, remove it
              mesh->points.resize(mesh->points.size() - 3);
              tri[i] = eid;
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

class IsosurfaceFunctor
{
public:
  typedef tbb::enumerable_thread_specific<PointLocator_t> TLS_pl;
  typedef tbb::enumerable_thread_specific<TriangleMesh_t> TLS_tm;
  typedef tbb::blocked_range3d<int, int, int> Range_t;

  IsosurfaceFunctor(const Image3D_t &vol, Float_t isoval,
                    TLS_pl &pointLocators, TLS_tm &meshPieces)
    : input(vol), isoval(isoval), pointLocators(pointLocators),
      meshPieces(meshPieces)
  {
  }

  void operator()(const Range_t &range) const
  {
    int extent[6] = { range.cols().begin(), range.cols().end() - 1,
                      range.rows().begin(), range.rows().end() - 1,
                      range.pages().begin(), range.pages().end() - 1 };
    //std::cout << "Processing block: [" << extent[0] << "," << extent[1] << ","
    //          << extent[2] << "," << extent[3] << "," << extent[4] << ","
    //          << extent[5] << "] dim: (" << extent[1] - extent[0] + 1 << ","
    //          << extent[3] - extent[2] + 1 << "," << extent[5] - extent[4] + 1
    //          << ")" << std::endl;
    PointLocator_t &ptLocator = this->pointLocators.local();
    TriangleMesh_t &meshPiece = this->meshPieces.local();
    extractIsosurfaceFromBlock(this->input, extent, this->isoval, ptLocator,
                               &meshPiece);
  }

private:
  const Image3D_t &input;
  Float_t isoval;
  TLS_pl &pointLocators;
  TLS_tm &meshPieces;
};


void extractIsosurface(const Image3D_t &vol, Float_t isoval,
                       TriangleMesh_t *mesh)
{
  const int *dims = vol.getDimension();
  const Float_t *origin = vol.getOrigin();
  const Float_t *spacing = vol.getSpacing();

  Float_t range[6];
  for (int i = 0; i < 3; ++i)
    {
    range[i * 2] = origin[i];
    range[(i * 2) + 1] = origin[i] +
                         (static_cast<Float_t>(dims[i] - 1) * spacing[i]);
    }

  PointLocator_t pl(range, dims[0]/8, dims[1]/8, dims[2]/8);
  IsosurfaceFunctor::TLS_pl pointLocators(pl);
  IsosurfaceFunctor::TLS_tm meshPieces;
  IsosurfaceFunctor::Range_t cellRange(0, dims[2] - 1, 64, 0, dims[1] - 1, 64,
                                       0, dims[0] - 1, 64);

  IsosurfaceFunctor func(vol, isoval, pointLocators, meshPieces);
  tbb::parallel_for(cellRange, func);
  mergeTriangleMeshes(meshPieces.begin(), meshPieces.end(), mesh);

  //std::cout << "Computed using " << meshPieces.size() << " threads"
  //          << std::endl;
}

}; //namespace scalar

namespace ompimp
{

void extractIsosurface(const Image3D_t &vol, Float_t isoval,
                       TriangleMesh_t *mesh)
{
  const int *dims = vol.getDimension();
  const Float_t *origin = vol.getOrigin();
  const Float_t *spacing = vol.getSpacing();

  Float_t range[6];
  for (int i = 0; i < 3; ++i)
    {
    range[i * 2] = origin[i];
    range[(i * 2) + 1] = origin[i] +
                         (static_cast<Float_t>(dims[i] - 1) * spacing[i]);
    }

  PointLocator_t pl(range, dims[0]/8, dims[1]/8, dims[2]/8);
  std::vector<TriangleMesh_t> meshPieces;

  int blocksDim[] = {(dims[0] - 1 + 63)/64, (dims[1] - 1 + 63)/64,
                     (dims[2] - 1 + 63)/64};
  int nblocks = blocksDim[0] * blocksDim[1] * blocksDim[2];

# pragma omp parallel for firstprivate(pl) shared(meshPieces)
  for (int i = 0; i < nblocks; ++i)
    {
#     pragma omp critical
      if (meshPieces.empty())
        {
        meshPieces.resize(omp_get_num_threads());
        }

      int blockIdx[3];
      blockIdx[2] = i/(blocksDim[0] * blocksDim[1]);
      blockIdx[1] = (i%(blocksDim[0] * blocksDim[1]))/blocksDim[0];
      blockIdx[0] = (i%(blocksDim[0] * blocksDim[1]))%blocksDim[0];

      int extent[6];
      extent[0] = blockIdx[0] * 64;
      extent[1] = std::min(extent[0] + 63, dims[0] - 2);
      extent[2] = blockIdx[1] * 64;
      extent[3] = std::min(extent[2] + 63, dims[1] - 2);
      extent[4] = blockIdx[2] * 64;
      extent[5] = std::min(extent[4] + 63, dims[2] - 2);

      scalar::extractIsosurfaceFromBlock(vol, extent, isoval, pl,
        &meshPieces[omp_get_thread_num()]);
    }

  mergeTriangleMeshes(meshPieces.begin(), meshPieces.end(), mesh);
}

}

namespace shortvec {

void extractIsosurfaceFromBlock(const Image3D_t &vol, const int ext[6],
  Float_t isoval, PointLocator_t &pointLocator, TriangleMesh_t *mesh)
{
  static const int caseMask[] = { 1, 2, 4, 8, 16, 32, 64, 128 };

  const int *dims = vol.getDimension();
  const Float_t *origin = vol.getOrigin();
  const Float_t *spacing = vol.getSpacing();
  const Float_t *buffer = vol.getData();

  int sliceSize = dims[0] * dims[1];

  int ptIdx = pointLocator.numberOfPoints();
  GetPointPosition getPosition(mesh->points);

  // march through each cell
  Float_t zpos = origin[2] + (static_cast<Float_t>(ext[4]) * spacing[2]);
  for (int zidx = ext[4]; zidx <= ext[5]; ++zidx, zpos += spacing[2])
    {
    Float_t ypos = origin[1] + (static_cast<Float_t>(ext[2]) * spacing[1]);
    for (int yidx = ext[2]; yidx <= ext[3]; ++yidx, ypos += spacing[1])
      {
      Float_t xpos = origin[0] + (static_cast<Float_t>(ext[0]) * spacing[0]);
      for (int xidx = ext[0]; xidx <= ext[1]; ++xidx, xpos += spacing[0])
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
            Float_t ptpos[3];
            ispc::lerp(3, pos[v1], pos[v2], w, ptpos);

            // temporarily add the point
            mesh->points.push_back(ptpos[0]);
            mesh->points.push_back(ptpos[1]);
            mesh->points.push_back(ptpos[2]);

            bool exists = false;
            int eid = *pointLocator.insert(ptIdx, &exists, getPosition);
            if (!exists)
              {
              // interpolate vertex normal
              Float_t norm[3];
              ispc::lerp(3, grad[v1], grad[v2], w, norm);

              mesh->normals.push_back(norm[0]);
              mesh->normals.push_back(norm[1]);
              mesh->normals.push_back(norm[2]);

              tri[i] = ptIdx++;
              }
            else
              {
              // remove the point
              mesh->points.resize(mesh->points.size() - 3);
              tri[i] = eid;
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

class IsosurfaceFunctor
{
public:
  typedef tbb::enumerable_thread_specific<PointLocator_t> TLS_pl;
  typedef tbb::enumerable_thread_specific<TriangleMesh_t> TLS_tm;
  typedef tbb::blocked_range3d<int, int, int> Range_t;

  IsosurfaceFunctor(const Image3D_t &vol, Float_t isoval,
                    TLS_pl &pointLocators, TLS_tm &meshPieces)
    : input(vol), isoval(isoval), pointLocators(pointLocators),
      meshPieces(meshPieces)
  {
  }

  void operator()(const Range_t &range) const
  {
    int extent[6] = { range.cols().begin(), range.cols().end() - 1,
                      range.rows().begin(), range.rows().end() - 1,
                      range.pages().begin(), range.pages().end() - 1 };
    PointLocator_t &ptLocator = this->pointLocators.local();
    TriangleMesh_t &meshPiece = this->meshPieces.local();
    extractIsosurfaceFromBlock(this->input, extent, this->isoval, ptLocator,
                               &meshPiece);
  }

private:
  const Image3D_t &input;
  Float_t isoval;
  TLS_pl &pointLocators;
  TLS_tm &meshPieces;
};


void extractIsosurface(const Image3D_t &vol, Float_t isoval,
                       TriangleMesh_t *mesh)
{
  const int *dims = vol.getDimension();
  const Float_t *origin = vol.getOrigin();
  const Float_t *spacing = vol.getSpacing();

  Float_t range[6];
  for (int i = 0; i < 3; ++i)
    {
    range[i * 2] = origin[i];
    range[(i * 2) + 1] = origin[i] +
                         (static_cast<Float_t>(dims[i] - 1) * spacing[i]);
    }

  PointLocator_t pl(range, dims[0]/8, dims[1]/8, dims[2]/8);
  IsosurfaceFunctor::TLS_pl pointLocators(pl);
  IsosurfaceFunctor::TLS_tm meshPieces;
  IsosurfaceFunctor::Range_t cellRange(0, dims[2] - 1, 64, 0, dims[1] - 1, 64,
                                       0, dims[0] - 1, 64);

  IsosurfaceFunctor func(vol, isoval, pointLocators, meshPieces);
  tbb::parallel_for(cellRange, func);
  mergeTriangleMeshes(meshPieces.begin(), meshPieces.end(), mesh);
}

}; // namespace shortvec

namespace simd {

void extractIsosurfaceFromBlock(const Image3D_t &vol, const int ext[6],
  Float_t isoval, PointLocator_t &pointLocator, TriangleMesh_t *mesh)
{
  const int *dims = vol.getDimension();
  const Float_t *origin = vol.getOrigin();
  const Float_t *spacing = vol.getSpacing();
  const Float_t *buffer = vol.getData();

  int blockDim[3] = { ext[1] - ext[0] + 1, ext[3] - ext[2] + 1,
                      ext[5] - ext[4] + 1 };

  int totalCells = blockDim[0] * blockDim[1] * blockDim[2];
  std::vector<int> cases(totalCells);
  int activeCells = ispc::getCellCases(buffer, dims, ext, blockDim, isoval,
                                       &cases[0]);

  int ntris = 0;
  std::vector<int> cellIndex(activeCells), vertInds(activeCells);
  for (int z = ext[4], s = 0, d = 0; z <= ext[5]; ++z)
    {
    for (int y = ext[2]; y <= ext[3]; ++y)
      {
      for (int x = ext[0]; x <= ext[1]; ++x, ++s)
        {
        if (cases[s] != 0 && cases[s] != 255)
          {
          cellIndex[d] = z * (dims[0] - 1) * (dims[1] - 1) +
                         y * (dims[0] - 1) + x;
          vertInds[d] = ntris * 3 * 3;
          cases[d++] = cases[s];
          ntris += MarchingCubesTables::getNumberOfTriangles(cases[s]);
          }
        }
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
  GetPointPosition getPosition(mesh->points);

  int tri[3];
  for (int i = 0, j = 0, ptid = pointLocator.numberOfPoints(); i < nverts;
       ++i, ++j)
    {
    bool exists = false;

    // temporarily add the point
    mesh->points.push_back(points[(i * 3) + 0]);
    mesh->points.push_back(points[(i * 3) + 1]);
    mesh->points.push_back(points[(i * 3) + 2]);

    int eid = *pointLocator.insert(ptid, &exists, getPosition);
    if (!exists)
      {
      mesh->normals.push_back(normals[(i * 3) + 0]);
      mesh->normals.push_back(normals[(i * 3) + 1]);
      mesh->normals.push_back(normals[(i * 3) + 2]);
      ++ptid;
      }
    else
      {
      // remove the point
      mesh->points.resize(mesh->points.size() - 3);
      }
    tri[j] = eid;

    if (j == 2)
      {
      // remove degenerate triangles
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

class IsosurfaceFunctor
{
public:
  typedef tbb::enumerable_thread_specific<PointLocator_t> TLS_pl;
  typedef tbb::enumerable_thread_specific<TriangleMesh_t> TLS_tm;
  typedef tbb::blocked_range3d<int, int, int> Range_t;

  IsosurfaceFunctor(const Image3D_t &vol, Float_t isoval,
                    TLS_pl &pointLocators, TLS_tm &meshPieces)
    : input(vol), isoval(isoval), pointLocators(pointLocators),
      meshPieces(meshPieces)
  {
  }

  void operator()(const Range_t &range) const
  {
    int extent[6] = { range.cols().begin(), range.cols().end() - 1,
                      range.rows().begin(), range.rows().end() - 1,
                      range.pages().begin(), range.pages().end() - 1 };
    PointLocator_t &ptLocator = this->pointLocators.local();
    TriangleMesh_t &meshPiece = this->meshPieces.local();
    extractIsosurfaceFromBlock(this->input, extent, this->isoval, ptLocator,
                               &meshPiece);
  }

private:
  const Image3D_t &input;
  Float_t isoval;
  TLS_pl &pointLocators;
  TLS_tm &meshPieces;
};


void extractIsosurface(const Image3D_t &vol, Float_t isoval,
                       TriangleMesh_t *mesh)
{
  const int *dims = vol.getDimension();
  const Float_t *origin = vol.getOrigin();
  const Float_t *spacing = vol.getSpacing();

  Float_t range[6];
  for (int i = 0; i < 3; ++i)
    {
    range[i * 2] = origin[i];
    range[(i * 2) + 1] = origin[i] +
                         (static_cast<Float_t>(dims[i] - 1) * spacing[i]);
    }

  PointLocator_t pl(range, dims[0]/8, dims[1]/8, dims[2]/8);
  IsosurfaceFunctor::TLS_pl pointLocators(pl);
  IsosurfaceFunctor::TLS_tm meshPieces;
  IsosurfaceFunctor::Range_t cellRange(0, dims[2] - 1, 64, 0, dims[1] - 1, 64,
                                       0, dims[0] - 1, 64);

  IsosurfaceFunctor func(vol, isoval, pointLocators, meshPieces);
  tbb::parallel_for(cellRange, func);
  mergeTriangleMeshes(meshPieces.begin(), meshPieces.end(), mesh);
}

}; // namespace simd


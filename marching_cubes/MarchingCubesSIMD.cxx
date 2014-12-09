#include "Image3D.h"
#include "MarchingCubesTables.h"
#include "PointLocator3D.h"
#include "TriangleMesh3D.h"
#include "Gradient.h"

#include <iostream>

struct Position
{
  double xyz[3];
  int idx;

  const double* getPosition() const
  {
    return xyz;
  }
};


namespace simd {

void extractIsosurface(const Image3D &vol, double isoval, TriangleMesh3D *mesh)
{
  static const int caseMask[] = { 1, 2, 4, 8, 16, 32, 64, 128 };
  static const int everts[12][2] = {{0,1}, {1,2}, {3,2}, {0,3}, {4,5}, {5,6},
                                    {7,6}, {4,7}, {0,4}, {1,5}, {3,7}, {2,6}};
  const int *dims = vol.getDimension();
  const double *origin = vol.getOrigin();
  const double *spacing = vol.getSpacing();
  const double *buffer = vol.getData();

  int sliceSize = dims[0] * dims[1];
  double range[6];

  // initialize point locator data-structure
  for (int i = 0; i < 3; ++i)
    {
    range[i * 2] = origin[i];
    range[(i * 2) + 1] = origin[i] +
                         (static_cast<double>(dims[i] - 1) * spacing[i]);
    }
  PointLocator3D<Position> pointLocator(range, 32, 32, 16);

  int ptIdx = 0;

  // march through each cell
  double zpos = origin[2];
  for (int zidx = 0; zidx < (dims[2] - 1); ++zidx, zpos += spacing[2])
    {
    double ypos = origin[1];
    for (int yidx = 0; yidx < (dims[1] - 1); ++yidx, ypos += spacing[1])
      {
      double xpos = origin[0];
      for (int xidx = 0; xidx < (dims[0] - 1); ++xidx, xpos += spacing[0])
        {
        double pos[8][3], grad[8][3];
        double val[8];

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

        // compute cell id based on isosurface intersection with the edges
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

        // get physical position and gadient at the cell points
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
        const int *edges = MarchingCubesTables::getTriangleCases()[cellId];
        for (; *edges != -1; edges += 3)
          {
          int tri[3];
          for (int i = 0; i < 3; ++i)
            {
            int v1 = everts[edges[i]][0], v2 = everts[edges[i]][1];
            double w = (isoval - val[v1])/(val[v2] - val[v1]);

            // interpolate vertex position
            TriangleMesh3D::Vertex v;
            ispc::lerp(3, pos[v1], pos[v2], w, v.pos);

            Position pt;
            pt.xyz[0] = v.pos[0];
            pt.xyz[1] = v.pos[1];
            pt.xyz[2] = v.pos[2];
            pt.idx = ptIdx;

            bool exists = false;
            Position *pt_ = pointLocator.insert(pt, &exists);
            if (!exists)
              {
              // interpolate vertex normal
              ispc::lerp(3, grad[v1], grad[v2], w, v.norm);

              mesh->vertices.push_back(v);
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

}; // namespace simd


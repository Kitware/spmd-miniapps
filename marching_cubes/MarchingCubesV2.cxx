#include "MarchingCubes.h"
#include "MarchingCubesTables.h"

#include <boost/chrono.hpp>
#include <boost/unordered_map.hpp>

#include <iostream>
#include <vector>

typedef boost::chrono::steady_clock Clock;

inline void generateEdgeIdOffsets(int xdim, int xydim, int offsets[12])
{
  offsets[0] = 0;
  offsets[1] = 4;
  offsets[2] = 3 * xdim;
  offsets[3] = 1;
  offsets[4] = 3 * xydim;
  offsets[5] = (3 * (xydim + 1)) + 1;
  offsets[6] = 3 * (xdim + xydim);
  offsets[7] = (3 * xydim) + 1;
  offsets[8] = 2;
  offsets[9] = 5;
  offsets[10] = (3 * xdim) + 2;
  offsets[11] = (3 * (xdim + 1)) + 2;
}

inline size_t estimateNumberOfBuckets(const int dims[3], float loadFactor)
{
  size_t estimatedNumOfVertexes = (dims[0] * dims[1] * dims[2])/32;
  return estimatedNumOfVertexes/loadFactor;
}

struct PointInfo
{
  Float_t pos[3];
  Float_t val;
  Float_t grad[3];
};

struct EdgeInfo
{
  PointInfo pt1, pt2;
};

typedef boost::unordered_map<int, int> EdgeUnorderedMap;

inline void createPointInfo(Float_t xpos, Float_t ypos, Float_t zpos,
                            Float_t val, PointInfo *p)
{
  p->pos[0] = xpos;
  p->pos[1] = ypos;
  p->pos[2] = zpos;
  p->val = val;
}

static void computeGradient(int xidx, int yidx, int zidx, int idx,
                            const Float_t cachedVals[8],
                            const Float_t *buffer, const int dims[3],
                            const Float_t spacing[3],
                            Float_t grad[3])
{
  Float_t v1, v2, fac;

  v1 = (xidx > 0) ? buffer[idx - 1] : cachedVals[0];
  v2 = (xidx < (dims[0] - 1)) ? cachedVals[1] : cachedVals[0];
  fac = (xidx > 0 && xidx < (dims[0] - 1)) ? 2.0 : 1.0;
  grad[0] = (v1 - v2)/(fac * spacing[0]);

  v1 = (yidx > 0) ? buffer[idx - dims[0]] : cachedVals[0];
  v2 = (yidx < (dims[1] - 1)) ? cachedVals[3] : cachedVals[0];
  fac = (yidx > 0 && yidx < (dims[1] - 1)) ? 2.0 : 1.0;
  grad[1] = (v1 - v2)/(fac * spacing[1]);

  v1 = (zidx > 0) ? buffer[idx - (dims[0] * dims[1])] : cachedVals[0];
  v2 = (zidx < (dims[2] - 1)) ? cachedVals[4] : cachedVals[0];
  fac = (zidx > 0 && zidx < (dims[2] - 1)) ? 2.0 : 1.0;
  grad[2] = (v1 - v2)/(fac * spacing[2]);
}

inline bool intersects(Float_t p1val, Float_t p2val, Float_t isoval)
{
  return ((p1val >= isoval && p2val < isoval) ||
          (p2val >= isoval && p1val < isoval));
}

inline Float_t lerp(Float_t a, Float_t b, Float_t w)
{
  return a + (w * (b - a));
}

namespace scalar_2 {

void extractIsosurface(const Image3D_t &vol, Float_t isoval, TriangleMesh_t *mesh)
{
  static const int caseMask[] = { 1, 2, 4, 8, 16, 32, 64, 128 };
  static const int everts[12][2] = {{0,1}, {1,2}, {3,2}, {0,3}, {4,5}, {5,6},
                                    {7,6}, {4,7}, {0,4}, {1,5}, {3,7}, {2,6}};

  const int *dims = vol.getDimension();
  const Float_t *origin = vol.getOrigin();
  const Float_t *spacing = vol.getSpacing();
  const Float_t *buffer = vol.getData();
  int sliceSize = dims[0] * dims[1];

  int edgeIdOffsets[12];
  generateEdgeIdOffsets(dims[0], sliceSize, edgeIdOffsets);

  Clock::time_point start, middle, finish;
  start = Clock::now();

  std::vector<int> indexes;
  std::vector<Float_t> points, normals;

  std::vector<EdgeInfo> edges;

  const Float_t loadFactor = 1.0;
  EdgeUnorderedMap edgeMap(estimateNumberOfBuckets(dims, loadFactor));
  edgeMap.max_load_factor(loadFactor);

#if 0
  std::cout << "nbuckets: " << edgeMap.bucket_count()
            << ", estimate: " << estimateNumberOfBuckets(dims, loadFactor)
            << ", maxbuckets: " << edgeMap.max_bucket_count()
            << ", loadFactor: " << edgeMap.load_factor()
            << ", maxloadfactor: " << edgeMap.max_load_factor()
            << ", maxsize: " << edgeMap.max_size()
            << std::endl;
#endif

  Float_t zpos = origin[2];
  for (int zidx = 0, idx = 0; zidx < dims[2]; ++zidx, zpos += spacing[2])
    {
    Float_t ypos = origin[1];
    for (int yidx = 0; yidx < dims[1]; ++yidx, ypos += spacing[1])
      {
      Float_t xpos = origin[0];
      for (int xidx = 0; xidx < dims[0]; ++xidx, ++idx, xpos += spacing[0])
        {
        Float_t val[8];
        val[0] = buffer[idx];

        int edgeId = idx * 3;
        int edgexsects[12] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

        bool isCell = false, skipOutEdges = false;
        if (xidx < (dims[0] - 1) && yidx < (dims[1] - 1) &&
            zidx < (dims[2] - 1))
          {
          isCell = true;
          val[1] = buffer[idx + 1];
          val[2] = buffer[idx + dims[0] + 1];
          val[3] = buffer[idx + dims[0]];
          val[4] = buffer[idx + sliceSize];
          val[5] = buffer[idx + 1 + sliceSize];
          val[6] = buffer[idx + dims[0] + 1 + sliceSize];
          val[7] = buffer[idx + dims[0] + sliceSize];

          // compute cell case based on isosurface intersection with the edges
          int caseId = 0;
          for (int i = 0; i < 8; ++i)
            {
            caseId |= (val[i] >= isoval) ? caseMask[i] : 0;
            }

           if (caseId != 0 && caseId != 255)
            {
            const int *eids = MarchingCubesTables::getTriangleCases()[caseId];
            for (; *eids != -1; ++eids)
              {
              indexes.push_back(edgeId + edgeIdOffsets[*eids]);
              edgexsects[*eids] = 1;
              }
            }
          else
            {
            skipOutEdges = true;
            }
          }

        bool needPointGradient = false;
        int outEdges[3] = { -1, -1, -1 };
        if (!skipOutEdges)
          {
          bool xsects = false;
          if (xidx < (dims[0] - 1))
            {
            if (isCell)
              {
              xsects = edgexsects[0];
              }
            else
              {
              val[1] = buffer[idx + 1];
              xsects = intersects(val[0], val[1], isoval);
              }
            if (xsects)
              {
              EdgeInfo e;
              createPointInfo(xpos, ypos, zpos, val[0], &e.pt1);
              createPointInfo(xpos + spacing[0], ypos, zpos, val[1], &e.pt2);

              outEdges[0] = edgeMap[edgeId] = edges.size();
              edges.push_back(e);
              needPointGradient = true;
              }
            }
          if (yidx < (dims[1] - 1))
            {
            if (isCell)
              {
              xsects = edgexsects[3];
              }
            else
              {
              val[3] = buffer[idx + dims[0]];
              xsects = intersects(val[0], val[3], isoval);
              }
            if (xsects)
              {
              EdgeInfo e;
              createPointInfo(xpos, ypos, zpos, val[0], &e.pt1);
              createPointInfo(xpos, ypos + spacing[1], zpos, val[3], &e.pt2);

              outEdges[1] = edgeMap[edgeId + 1] = edges.size();
              edges.push_back(e);
              needPointGradient = true;
              }
            }
          if (zidx < (dims[2] - 1))
            {
            if (isCell)
              {
              xsects = edgexsects[8];
              }
            else
              {
              val[4] = buffer[idx + sliceSize];
              xsects = intersects(val[0], val[4], isoval);
              }
            if (xsects)
              {
              EdgeInfo e;
              createPointInfo(xpos, ypos, zpos, val[0], &e.pt1);
              createPointInfo(xpos, ypos, zpos + spacing[2], val[4], &e.pt2);

              outEdges[2] = edgeMap[edgeId + 2] = edges.size();
              edges.push_back(e);
              needPointGradient = true;
              }
            }
          }

        int inEdges[3] = { -1, -1, -1 };
        int inEdgeIds[3] = { (idx - 1) * 3,
                             ((idx - dims[0]) * 3) + 1,
                             ((idx - sliceSize) * 3) + 2 };
        for (int i = 0; i < 3; ++i)
          {
          EdgeUnorderedMap::iterator loc = edgeMap.find(inEdgeIds[i]);
          if (loc != edgeMap.end())
            {
            inEdges[i] = loc->second;
            needPointGradient = true;
            }
          }

        if (needPointGradient)
          {
          Float_t grad[3];
          computeGradient(xidx, yidx, zidx, idx, val, buffer, dims, spacing,
                          grad);
          for (int i = 0; i < 3; ++i)
            {
            if (inEdges[i] != -1)
              {
              edges[inEdges[i]].pt2.grad[0] = grad[0];
              edges[inEdges[i]].pt2.grad[1] = grad[1];
              edges[inEdges[i]].pt2.grad[2] = grad[2];
              }
            if (outEdges[i] != -1)
              {
              edges[outEdges[i]].pt1.grad[0] = grad[0];
              edges[outEdges[i]].pt1.grad[1] = grad[1];
              edges[outEdges[i]].pt1.grad[2] = grad[2];
              }
            }
          }
        }
      }
    }

  for (unsigned i = 0; i < indexes.size(); ++i)
    {
    indexes[i] = edgeMap.at(indexes[i]);
    }

  middle = Clock::now();

  points.reserve(edges.size() * 3);
  normals.reserve(edges.size() * 3);
  for (unsigned i = 0; i < edges.size(); ++i)
    {
    const EdgeInfo &edge = edges[i];

    Float_t w = (isoval - edge.pt1.val)/(edge.pt2.val - edge.pt1.val);
    for (int ii = 0; ii < 3; ++ii)
      {
      points.push_back(lerp(edge.pt1.pos[ii], edge.pt2.pos[ii], w));
      normals.push_back(lerp(edge.pt1.grad[ii], edge.pt2.grad[ii], w));
      }
    }

  mesh->points.swap(points);
  mesh->normals.swap(normals);
  mesh->indexes.swap(indexes);

  finish = Clock::now();

  boost::chrono::duration<double> part1 = middle - start;
  boost::chrono::duration<double> part2 = finish - middle;
  std::cout << "part1 cost: " << part1 << ", part2 cost: " << part2
            << std::endl;

#if 0
  std::cout << "nbuckets: " << edgeMap.bucket_count()
            << ", loadFactor: " << edgeMap.load_factor()
            << ", maxloadfactor: " << edgeMap.max_load_factor()
            << ", size: " << edgeMap.size()
            << std::endl;
#endif
}

}; //namespace scalar_2


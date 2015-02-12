#include "Isosurface.h"
#include "Isosurface.ispc.h"
#include "MarchingTetsTables.h"

#include <PointLocator3D.h>

#include <boost/unordered_map.hpp>
#include <algorithm>
#include <cmath>
#include <iostream>

class GetPositionFromIndex
{
public:
  GetPositionFromIndex(const std::vector<Float_t> &points,
                       const Float_t *thispt)
    : points(points), thispt(thispt)
  {
  }

  const Float_t * operator()(int ptind) const
  {
  return (this->points.size() == static_cast<size_t>(ptind * 3)) ?
          this->thispt : &(this->points[ptind * 3]);
  }

private:
  const std::vector<Float_t> &points;
  const Float_t *thispt;
};


struct EdgeKey
{
  int v1, v2;

  EdgeKey(int v1, int v2) : v1(v1), v2(v2) {}
  bool operator==(const EdgeKey &e) const
  {
    return (v1 == e.v1) && (v2 == e.v2);
  }
};

inline size_t hash_value(const EdgeKey &k)
{
  size_t seed = 0;
  boost::hash_combine(seed, k.v1);
  boost::hash_combine(seed, k.v2);
  return seed;
}

typedef boost::unordered_map<EdgeKey, int> EdgeUnorderedMap;


namespace scalar
{

inline Float_t lerp(Float_t a, Float_t b, Float_t w)
{
  return a + (w * (b - a));
}

inline void computeTriangleNormal(const Float_t verts[3][3], Float_t normal[3])
{
  Float_t v2v1[3] = { verts[1][0] - verts[0][0], verts[1][1] - verts[0][1],
                      verts[1][2] - verts[0][2] };
  Float_t v3v1[3] = { verts[2][0] - verts[0][0], verts[2][1] - verts[0][1],
                      verts[2][2] - verts[0][2] };

  normal[0] = (v2v1[1] * v3v1[2]) - (v2v1[2] * v3v1[1]);
  normal[1] = (v2v1[2] * v3v1[0]) - (v2v1[0] * v3v1[2]);
  normal[2] = (v2v1[0] * v3v1[1]) - (v2v1[1] * v3v1[0]);
}

inline void normalize(Float_t normal[3])
{
  Float_t mag = sqrt((normal[0] * normal[0]) + (normal[1] * normal[1]) +
                     (normal[2] * normal[2]));
  if (mag != 0)
    {
    normal[0] /= mag;
    normal[1] /= mag;
    normal[2] /= mag;
    }
}

inline void computeTetrahedronMeshBoundingBox(const TetrahedronMesh_t &tetmesh,
                                              Float_t *range)
{
  const Float_t *pt = &tetmesh.points.front();
  const Float_t *end = pt + tetmesh.points.size();
  for (; pt != end; pt += 3)
    {
    range[0] = std::min(range[0], pt[0]);
    range[1] = std::max(range[1], pt[0]);
    range[2] = std::min(range[2], pt[1]);
    range[3] = std::max(range[3], pt[1]);
    range[4] = std::min(range[4], pt[2]);
    range[5] = std::max(range[5], pt[2]);
    }
}


#define USE_HASHING_LOCATOR 1

void extractIsosurface(const TetrahedronMesh_t &tetmesh, Float_t isoval,
                       TriangleMesh_t *trimesh)
{
  static const int caseMask[4] = { 1, 2, 4, 8 };

#if !USE_HASHING_LOCATOR
  Float_t range[6] = { 0, 0, 0, 0, 0, 0 };
  computeTetrahedronMeshBoundingBox(tetmesh, range);

  PointLocator3D<Float_t, int, GetPositionFromIndex> ptlocator(range,
                                                               32, 32, 32);
#else
  EdgeUnorderedMap emap;
#endif

  for (int i = 0, npts = 0; i < tetmesh.numberOfTetrahedra(); ++i)
    {
    int ptinds[4];
    std::copy(tetmesh.indexes.begin() + (i * 4),
              tetmesh.indexes.begin() + (i * 4) + 4, ptinds);

    Float_t val[4];
    val[0] = tetmesh.values[ptinds[0]];
    val[1] = tetmesh.values[ptinds[1]];
    val[2] = tetmesh.values[ptinds[2]];
    val[3] = tetmesh.values[ptinds[3]];

    int caseId = 0;
    for (int ii = 0; ii < 4; ++ii)
      {
      if (val[ii] >= isoval)
        {
        caseId |= caseMask[ii];
        }
      }

    const int *triEdges = MarchingTetsTables::getCaseTrianglesEdges(caseId);
    for (; *triEdges != -1; triEdges += 3)
      {
      Float_t verts[3][3];
      int vinds[3];
      for (int e = 0; e < 3; ++e)
        {
        int v1 = MarchingTetsTables::getEdgeVertices(triEdges[e])[0];
        int v2 = MarchingTetsTables::getEdgeVertices(triEdges[e])[1];

        if (val[v1] > val[v2])
          {
          std::swap(v1, v2);
          }

        Float_t dv = val[v2] - val[v1];
        Float_t w = (dv == 0.0) ? 0.0 : (isoval - val[v1])/dv;

        for (int ii = 0; ii < 3; ++ii)
          {
          Float_t x1 = tetmesh.points[(ptinds[v1] * 3) + ii];
          Float_t x2 = tetmesh.points[(ptinds[v2] * 3) + ii];
          verts[e][ii] = lerp(x1, x2, w);
          }

#if !USE_HASHING_LOCATOR
        GetPositionFromIndex gp(trimesh->points, verts[e]);
        bool exists = false;
        vinds[e] = *ptlocator.insert(npts, &exists, gp);
#else
        EdgeKey ekey(ptinds[v1], ptinds[v2]);
        EdgeUnorderedMap::iterator loc = emap.find(ekey);
        bool exists = (loc != emap.end());
        vinds[e] = exists ? loc->second : npts;
        if (!exists)
          {
          emap[ekey] = npts;
          }
#endif
        if (!exists)
          {
          for (int ii = 0; ii < 3; ++ii)
            {
            trimesh->points.push_back(verts[e][ii]);
            trimesh->normals.push_back(0);
            }
          ++npts;
          }
        trimesh->indexes.push_back(vinds[e]);
        }

      Float_t normal[3];
      computeTriangleNormal(verts, normal);
      for (int e = 0; e < 3; ++e)
        {
        Float_t *vp = &(trimesh->normals[vinds[e] * 3]);
        vp[0] += normal[0];
        vp[1] += normal[1];
        vp[2] += normal[2];
        }
      }
    }

  for (size_t i = 0; i < trimesh->normals.size(); i += 3)
    {
    normalize(&(trimesh->normals[i]));
    }
}

}; // namespace scalar


namespace simd
{

const int GANG_SIZE = ISPC_GANG_SIZE;

inline size_t getSOASize(size_t aos_size, size_t gangSize)
{
  return (aos_size + gangSize - (aos_size % gangSize));
}

inline size_t getNumberOfGangs(size_t aos_size, size_t gangSize)
{
  return (aos_size + gangSize - 1)/gangSize;
}

struct TetInfoA
{
  int tetIdx;
  int caseId;
};

class CaseIdIsLess
{
public:
  bool operator()(const TetInfoA &t1, const TetInfoA &t2) const
  {
    return (t1.caseId < t2.caseId);
  }
};

struct TetInfoB
{
  int ptidx[4];
  Float_t val[4];
};

#define USE_VECTORIZED_NORMALIZE 1

using ispc::TetInfo_soa;

void extractIsosurface(const TetrahedronMesh_t &tetmesh, Float_t isoval,
                       TriangleMesh_t *trimesh)
{
  static const int caseMask[4] = { 1, 2, 4, 8 };

  std::vector<TetInfoA> tetInfoA;
  std::vector<TetInfoB> tetInfoB;
  for (int i = 0, idx = 0; i < tetmesh.numberOfTetrahedra(); ++i)
    {
    TetInfoA tiA;
    TetInfoB tiB;

    tiA.tetIdx = idx;
    tiA.caseId = 0;
    for (int ii = 0; ii < 4; ++ii)
      {
      tiB.ptidx[ii] = tetmesh.indexes[(i * 4) + ii];
      tiB.val[ii] = tetmesh.values[tiB.ptidx[ii]];

      if (tiB.val[ii] >= isoval)
        {
        tiA.caseId |= caseMask[ii];
        }
      }

    if (tiA.caseId != 0 && tiA.caseId != 15)
      {
      tetInfoA.push_back(tiA);
      tetInfoB.push_back(tiB);
      ++idx;
      }
    }

  std::sort(tetInfoA.begin(), tetInfoA.end(), CaseIdIsLess());

  int numContributingTets = tetInfoA.size();
  std::vector<TetInfo_soa> tetinfo_soa(getNumberOfGangs(numContributingTets,
                                                        GANG_SIZE));
  int numTriangles = 0;
  for (int i = 0, g = 0, gm = 0; i < numContributingTets; ++i, ++gm)
    {
    if (gm == GANG_SIZE)
      {
      ++g;
      gm = 0;
      }

    int tetIdx = tetinfo_soa[g].tetIdx[gm] = tetInfoA[i].tetIdx;
    int caseId = tetinfo_soa[g].caseId[gm] = tetInfoA[i].caseId;
    for (int ii = 0; ii < 4; ++ii)
      {
      int ptidx = tetinfo_soa[g].ptidx[ii][gm] = tetInfoB[tetIdx].ptidx[ii];
      tetinfo_soa[g].val[ii][gm] = tetInfoB[tetIdx].val[ii];

      ptidx *= 3;
      tetinfo_soa[g].pts[ii][0][gm] = tetmesh.points[ptidx];
      tetinfo_soa[g].pts[ii][1][gm] = tetmesh.points[ptidx + 1];
      tetinfo_soa[g].pts[ii][2][gm] = tetmesh.points[ptidx + 2];
      }

    numTriangles += MarchingTetsTables::getNumberOfTriangles(caseId);
    }

  int numVerts = numTriangles * 3;
  std::vector<Float_t> triPoints(numVerts * 3);
  std::vector<int> triPointKeys(numVerts * 2);
  std::vector<Float_t>
    triCellNorms(getSOASize(numTriangles * 3, GANG_SIZE * 3));

  ispc::extractIsosurface_impl(&tetinfo_soa[0], numContributingTets, isoval,
                               MarchingTetsTables::getCaseTrianglesEdges(0),
                               MarchingTetsTables::getEdgeVertices(0),
                               &triPoints[0], &triPointKeys[0],
                               &triCellNorms[0]);

  EdgeUnorderedMap emap;

  std::vector<int> triIndexes(numTriangles * 3);
  std::vector<Float_t>
    triPtNorms(getSOASize(triPoints.size(), GANG_SIZE * 3), 0);

  int numPts = 0;
  for (int t = 0, i = 0; t < numTriangles; ++t)
    {
    for (int v = 0; v < 3; ++v, ++i)
      {
      int s = i * 3, d = numPts * 3;
      int ind = 0;
      EdgeKey k(triPointKeys[i * 2], triPointKeys[i * 2 + 1]);
      EdgeUnorderedMap::iterator loc = emap.find(k);
      if (loc == emap.end())
        {
        if (i > numPts)
          {
          triPoints[d++] = triPoints[s++];
          triPoints[d++] = triPoints[s++];
          triPoints[d++] = triPoints[s++];
          }
        ind = numPts++;
        emap[k] = ind;
        }
      else
        {
        ind = loc->second;
        }
      triIndexes[i] = ind;

      int didx = ((ind/GANG_SIZE) * GANG_SIZE * 3) + (ind % GANG_SIZE);
      int sidx = ((t/GANG_SIZE) * GANG_SIZE * 3) + (t % GANG_SIZE);
      for (int ii = 0; ii < 3; ++ii)
        {
#if USE_VECTORIZED_NORMALIZE
        triPtNorms[didx + (ii * GANG_SIZE)] +=
          triCellNorms[sidx + (ii * GANG_SIZE)];
#else
        triPtNorms[ind * 3 + ii] +=
          triCellNorms[sidx + (ii * GANG_SIZE)];
#endif
        }
      }
    }
  triPoints.resize(numPts * 3);

#if USE_VECTORIZED_NORMALIZE
  ispc::normalizeNormals(&triPtNorms[0], numPts);
  triPtNorms.resize(numPts * 3);
#else
  triPtNorms.resize(numPts * 3);
  for (size_t i = 0; i < triPtNorms.size(); i += 3)
    {
    scalar::normalize(&triPtNorms[i]);
    }
#endif

  trimesh->points.swap(triPoints);
  trimesh->normals.swap(triPtNorms);
  trimesh->indexes.swap(triIndexes);
}

}; // namespace simd



namespace simd_2
{

struct TriMeshHandle
{
  TriangleMesh_t *trimesh;
  EdgeUnorderedMap emap;
  int npts;
};


extern "C" void addTriangleToOutput(void* trimeshHandle, Float_t verts[3][3],
                                    int keys[3][2], Float_t normal[3])
{
  TriMeshHandle &th = *reinterpret_cast<TriMeshHandle*>(trimeshHandle);

  for (int v = 0; v < 3; ++v)
    {
    int vind = 0;
    EdgeKey ekey(keys[v][0], keys[v][1]);
    EdgeUnorderedMap::iterator loc = th.emap.find(ekey);
    if (loc != th.emap.end())
      {
      vind = loc->second;
      }
    else
      {
      for (int ii = 0; ii < 3; ++ii)
        {
        th.trimesh->points.push_back(verts[v][ii]);
        th.trimesh->normals.push_back(0);
        }
      vind = th.npts++;
      th.emap[ekey] = vind;
      }

    th.trimesh->indexes.push_back(vind);
    th.trimesh->normals[vind * 3 + 0] += normal[0];
    th.trimesh->normals[vind * 3 + 1] += normal[1];
    th.trimesh->normals[vind * 3 + 2] += normal[2];
    }
}


void extractIsosurface(const TetrahedronMesh_t &tetmesh, Float_t isoval,
                       TriangleMesh_t *trimesh)
{
  TriMeshHandle th;
  th.trimesh = trimesh;
  th.npts = 0;

  ispc::extractIosurface_impl_v2(&tetmesh.points[0], &tetmesh.values[0],
                                 &tetmesh.indexes[0], tetmesh.numberOfPoints(),
                                 tetmesh.numberOfTetrahedra(), isoval,
                                 MarchingTetsTables::getCaseTrianglesEdges(0),
                                 MarchingTetsTables::getEdgeVertices(0),
                                 &th);

  for (size_t i = 0; i < trimesh->normals.size(); i += 3)
    {
    scalar::normalize(&(trimesh->normals[i]));
    }
}

}; // namespace simd_2



#include "Isosurface.h"
#include "Isosurface.ispc.h"
#include "MarchingTetsTables.h"

#include <smp.h>

#include <boost/unordered_map.hpp>
#include <algorithm>
#include <cmath>
#include <iostream>

static const unsigned grainSize = 64 * 1024;

struct EdgeKey
{
  unsigned v1, v2;

  EdgeKey(unsigned v1, unsigned v2) : v1(v1), v2(v2) {}
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

typedef boost::unordered_map<EdgeKey, unsigned> EdgeUnorderedMap;


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


void extractIsosurfaceFromTetsRange(const TetrahedronMesh_t &tetmesh,
  unsigned from, unsigned to, Float_t isoval, EdgeUnorderedMap &emap,
  TriangleMesh_t *trimesh)
{
  static const int caseMask[4] = { 1, 2, 4, 8 };

  unsigned npts = trimesh->points.size()/3;
  for (unsigned i = from; i < to; ++i)
    {
    unsigned ptinds[4];
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
      unsigned vinds[3];
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

        EdgeKey ekey(ptinds[v1], ptinds[v2]);
        EdgeUnorderedMap::iterator loc = emap.find(ekey);
        bool exists = (loc != emap.end());
        vinds[e] = exists ? loc->second : npts;
        if (!exists)
          {
          emap[ekey] = npts;
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
}

class NormalizeFunctor
{
public:
  NormalizeFunctor(std::vector<Float_t> *normals) : normals(normals)
  {
  }

  void operator()(const smp::Range &range) const
  {
    for (unsigned i = range.begin(); i < range.end(); ++i)
      {
      normalize(&(this->normals->at(i * 3)));
      }
  }

private:
  std::vector<Float_t> *normals;
};

class IsosurfaceFunctor
{
public:
  typedef smp::ThreadLocalStorage<TriangleMesh_t> TLS_tm;
  typedef smp::ThreadLocalStorage<EdgeUnorderedMap> TLS_em;

  IsosurfaceFunctor(const TetrahedronMesh_t &tetmesh, Float_t isoval,
    TLS_em &edgeMaps, TLS_tm &meshPieces)
    : tetmesh(tetmesh), isoval(isoval), edgeMaps(edgeMaps),
      meshPieces(meshPieces)
  {
  }

  void operator()(const smp::Range &range) const
  {
    //std::cout << "Processing range: [" << range.begin() << "," << range.end()
    //          << "] ntets=" << range.size() << std::endl;
    TriangleMesh_t &meshPiece = this->meshPieces.local();
    EdgeUnorderedMap &edgeMap = this->edgeMaps.local();
    extractIsosurfaceFromTetsRange(this->tetmesh, range.begin(), range.end(),
      this->isoval, edgeMap, &meshPiece);
  }

private:
  const TetrahedronMesh_t &tetmesh;
  Float_t isoval;
  TLS_em &edgeMaps;
  TLS_tm &meshPieces;
};

void extractIsosurface(const TetrahedronMesh_t &tetmesh, Float_t isoval,
                       TriangleMesh_t *trimesh)
{
  IsosurfaceFunctor::TLS_em edgeMaps;
  IsosurfaceFunctor::TLS_tm meshPieces;
  smp::Range tetsRange(0, tetmesh.numberOfTetrahedra(), grainSize);
  IsosurfaceFunctor func(tetmesh, isoval, edgeMaps, meshPieces);
  smp::parallel_for(tetsRange, func);

  mergeTriangleMeshes(meshPieces.begin(), meshPieces.end(), trimesh);

  smp::Range normRange(0, trimesh->numberOfVertices(), grainSize);
  NormalizeFunctor nfunc(&trimesh->normals);
  smp::parallel_for(normRange, nfunc);

  //std::cout << "Computed using " << meshPieces.size() << " threads"
  //          << std::endl;
}

}; // namespace scalar


namespace simd
{

const unsigned GANG_SIZE = ISPC_GANG_SIZE;

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
  unsigned tetIdx;
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
  unsigned ptidx[4];
  Float_t val[4];
};

#define USE_VECTORIZED_NORMALIZE 1

using ispc::TetInfo_soa;

void extractIsosurfaceFromTetsRange(const TetrahedronMesh_t &tetmesh,
  unsigned from, unsigned to, Float_t isoval, EdgeUnorderedMap &emap,
  TriangleMesh_t *trimesh)
{
  static const int caseMask[4] = { 1, 2, 4, 8 };

  std::vector<TetInfoA> tetInfoA;
  std::vector<TetInfoB> tetInfoB;
  for (unsigned i = from, idx = 0; i < to; ++i)
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

  unsigned numContributingTets = tetInfoA.size();
  std::vector<TetInfo_soa> tetinfo_soa(getNumberOfGangs(numContributingTets,
                                                        GANG_SIZE));
  unsigned numTriangles = 0;
  for (unsigned i = 0, g = 0, gm = 0; i < numContributingTets; ++i, ++gm)
    {
    if (gm == GANG_SIZE)
      {
      ++g;
      gm = 0;
      }

    unsigned tetIdx = tetinfo_soa[g].tetIdx[gm] = tetInfoA[i].tetIdx;
    int caseId = tetinfo_soa[g].caseId[gm] = tetInfoA[i].caseId;
    for (int ii = 0; ii < 4; ++ii)
      {
      unsigned ptidx = tetinfo_soa[g].ptidx[ii][gm]
                     = tetInfoB[tetIdx].ptidx[ii];
      tetinfo_soa[g].val[ii][gm] = tetInfoB[tetIdx].val[ii];

      ptidx *= 3;
      tetinfo_soa[g].pts[ii][0][gm] = tetmesh.points[ptidx];
      tetinfo_soa[g].pts[ii][1][gm] = tetmesh.points[ptidx + 1];
      tetinfo_soa[g].pts[ii][2][gm] = tetmesh.points[ptidx + 2];
      }

    numTriangles += MarchingTetsTables::getNumberOfTriangles(caseId);
    }

  unsigned numVerts = numTriangles * 3;
  std::vector<Float_t> triPoints(numVerts * 3);
  std::vector<unsigned> triPointKeys(numVerts * 2);
  std::vector<Float_t>
    triCellNorms(getSOASize(numTriangles * 3, GANG_SIZE * 3));

  ispc::extractIsosurface_impl(&tetinfo_soa[0], numContributingTets, isoval,
                               MarchingTetsTables::getCaseTrianglesEdges(0),
                               MarchingTetsTables::getEdgeVertices(0),
                               &triPoints[0], &triPointKeys[0],
                               &triCellNorms[0]);


  std::vector<unsigned> triIndexes(numTriangles * 3);

  unsigned curNumPts = trimesh->points.size()/3;
#if USE_VECTORIZED_NORMALIZE
  // resize for potential new normals
  trimesh->normals.resize(getSOASize((curNumPts + numVerts) * 3,
    GANG_SIZE * 3), 0);
#endif

  unsigned numPts = 0, ptIdx = curNumPts;
  for (unsigned t = 0, i = 0; t < numTriangles; ++t)
    {
    for (int v = 0; v < 3; ++v, ++i)
      {
      unsigned s = i * 3, d = numPts * 3;
      unsigned ind = 0;
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
#if !USE_VECTORIZED_NORMALIZE
        trimesh->normals.push_back(0.0);
        trimesh->normals.push_back(0.0);
        trimesh->normals.push_back(0.0);
#endif
        ind = ptIdx++;
        emap[k] = ind;
        ++numPts;
        }
      else
        {
        ind = loc->second;
        }
      triIndexes[i] = ind;

      unsigned didx = ((ind/GANG_SIZE) * GANG_SIZE * 3) + (ind % GANG_SIZE);
      unsigned sidx = ((t/GANG_SIZE) * GANG_SIZE * 3) + (t % GANG_SIZE);
      for (int ii = 0; ii < 3; ++ii)
        {
#if USE_VECTORIZED_NORMALIZE
        trimesh->normals[didx + (ii * GANG_SIZE)] +=
          triCellNorms[sidx + (ii * GANG_SIZE)];
#else
        trimesh->normals[ind * 3 + ii] +=
          triCellNorms[sidx + (ii * GANG_SIZE)];
#endif
        }
      }
    }
  triPoints.resize(numPts * 3);

  trimesh->points.insert(trimesh->points.end(), triPoints.begin(),
                         triPoints.end());
  trimesh->indexes.insert(trimesh->indexes.end(), triIndexes.begin(),
                          triIndexes.end());
#if USE_VECTORIZED_NORMALIZE
  // resize to fit
  trimesh->normals.resize(getSOASize((curNumPts + numPts) * 3, GANG_SIZE * 3));
#endif
}


class IsosurfaceFunctor
{
public:
  typedef smp::ThreadLocalStorage<TriangleMesh_t> TLS_tm;
  typedef smp::ThreadLocalStorage<EdgeUnorderedMap> TLS_em;

  IsosurfaceFunctor(const TetrahedronMesh_t &tetmesh, Float_t isoval,
    TLS_em &edgeMaps, TLS_tm &meshPieces)
    : tetmesh(tetmesh), isoval(isoval), edgeMaps(edgeMaps),
      meshPieces(meshPieces)
  {
  }

  void operator()(const smp::Range &range) const
  {
    TriangleMesh_t &meshPiece = this->meshPieces.local();
    EdgeUnorderedMap &edgeMap = this->edgeMaps.local();
    extractIsosurfaceFromTetsRange(this->tetmesh, range.begin(), range.end(),
      this->isoval, edgeMap, &meshPiece);
  }

private:
  const TetrahedronMesh_t &tetmesh;
  Float_t isoval;
  TLS_em &edgeMaps;
  TLS_tm &meshPieces;
};

class NormalizeFunctor
{
public:
  NormalizeFunctor(const std::vector<TriangleMesh_t*> &meshptrs)
    : meshptrs(meshptrs)
  {
  }

  void operator()(const smp::Range &range) const
  {
    for (unsigned i = range.begin(); i != range.end(); ++i)
      {
      unsigned count = meshptrs[i]->numberOfVertices();
      ispc::normalizeNormals(&(meshptrs[i]->normals[0]), count);
      meshptrs[i]->normals.resize(count * 3);
      }
  }

private:
  const std::vector<TriangleMesh_t*> &meshptrs;
};

void extractIsosurface(const TetrahedronMesh_t &tetmesh, Float_t isoval,
                       TriangleMesh_t *trimesh)
{
  IsosurfaceFunctor::TLS_em edgeMaps;
  IsosurfaceFunctor::TLS_tm meshPieces;
  smp::Range tetsRange(0, tetmesh.numberOfTetrahedra());
  IsosurfaceFunctor func(tetmesh, isoval, edgeMaps, meshPieces);
  smp::parallel_for(tetsRange, func);

#if USE_VECTORIZED_NORMALIZE
  std::vector<TriangleMesh_t*> meshptrs;
  unsigned count = 0;
  for (IsosurfaceFunctor::TLS_tm::iterator i = meshPieces.begin();
       i != meshPieces.end(); ++i)
    {
    meshptrs.push_back(&(*i));
    ++count;
    }
  simd::NormalizeFunctor nfunc(meshptrs);
  smp::parallel_for(smp::Range(0, count), nfunc);
#endif

  mergeTriangleMeshes(meshPieces.begin(), meshPieces.end(), trimesh);

#if !USE_VECTORIZED_NORMALIZE
  smp::Range normRange(0, trimesh->numberOfVertices());
  scalar::NormalizeFunctor nfunc(&trimesh->normals);
  smp::parallel_for(normRange, nfunc);
#endif
}

}; // namespace simd



namespace simd_2
{

struct TriMeshHandle
{
  TriangleMesh_t trimesh;
  EdgeUnorderedMap emap;
  unsigned npts;

  TriMeshHandle() : npts(0)
  {
  }
};


extern "C" void addTriangleToOutput(void* trimeshHandle, Float_t verts[3][3],
                                    unsigned keys[3][2], Float_t normal[3])
{
  TriMeshHandle &th = *reinterpret_cast<TriMeshHandle*>(trimeshHandle);

  for (int v = 0; v < 3; ++v)
    {
    unsigned vind = 0;
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
        th.trimesh.points.push_back(verts[v][ii]);
        th.trimesh.normals.push_back(0);
        }
      vind = th.npts++;
      th.emap[ekey] = vind;
      }

    th.trimesh.indexes.push_back(vind);
    th.trimesh.normals[vind * 3 + 0] += normal[0];
    th.trimesh.normals[vind * 3 + 1] += normal[1];
    th.trimesh.normals[vind * 3 + 2] += normal[2];
    }
}


void extractIsosurfaceFromTetsRange(const TetrahedronMesh_t &tetmesh,
  unsigned from, unsigned to, Float_t isoval, TriMeshHandle &tmh)
{
  const unsigned *tetIdxPtr = &tetmesh.indexes[from * 4];
  unsigned numTets = to - from;
  ispc::extractIosurface_impl_v2(&tetmesh.points[0], &tetmesh.values[0],
                                 tetIdxPtr, tetmesh.numberOfPoints(),
                                 numTets, isoval,
                                 MarchingTetsTables::getCaseTrianglesEdges(0),
                                 MarchingTetsTables::getEdgeVertices(0),
                                 &tmh);
}

class IsosurfaceFunctor
{
public:
  typedef smp::ThreadLocalStorage<TriMeshHandle> TLS_tmh;

  IsosurfaceFunctor(const TetrahedronMesh_t &tetmesh, Float_t isoval,
    TLS_tmh &trimeshHandles)
    : tetmesh(tetmesh), isoval(isoval), trimeshHandles(trimeshHandles)
  {
  }

  void operator()(const smp::Range &range) const
  {
    TriMeshHandle &myMesh = this->trimeshHandles.local();
    extractIsosurfaceFromTetsRange(this->tetmesh, range.begin(), range.end(),
      this->isoval, myMesh);
  }

private:
  const TetrahedronMesh_t &tetmesh;
  Float_t isoval;
  TLS_tmh &trimeshHandles;
};


inline void swapTrianlgeMeshes(TriangleMesh_t &m1, TriangleMesh_t &m2)
{
  m1.points.swap(m2.points);
  m1.normals.swap(m2.normals);
  m1.indexes.swap(m2.indexes);
}

void extractIsosurface(const TetrahedronMesh_t &tetmesh, Float_t isoval,
                       TriangleMesh_t *trimesh)
{
  IsosurfaceFunctor::TLS_tmh trimeshHandles;
  smp::Range tetsRange(0, tetmesh.numberOfTetrahedra());
  IsosurfaceFunctor func(tetmesh, isoval, trimeshHandles);
  smp::parallel_for(tetsRange, func);

  std::vector<TriangleMesh_t> meshes(trimeshHandles.size());
  int i = 0;
  for (IsosurfaceFunctor::TLS_tmh::iterator it = trimeshHandles.begin();
    it != trimeshHandles.end(); ++it)
    {
    swapTrianlgeMeshes(meshes[i++], it->trimesh);
    }
  mergeTriangleMeshes(meshes.begin(), meshes.end(), trimesh);

  smp::Range normRange(0, trimesh->numberOfVertices());
  scalar::NormalizeFunctor nfunc(&trimesh->normals);
  smp::parallel_for(normRange, nfunc);
}

}; // namespace simd_2


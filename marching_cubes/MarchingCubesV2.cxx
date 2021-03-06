#include "MarchingCubes.h"
#include "MarchingCubesTables.h"
#include "MarchingCubesV2.ispc.h"

#include <smp.h>

#include <iostream>
#include <vector>

static const unsigned grainDim = 32;
static const unsigned INVALID_EDGE = ~0u;

inline void generateEdgeIdOffsets(unsigned xdim, unsigned xydim,
                                  unsigned offsets[12])
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

static void computePosition(const unsigned idx[3], const Float_t origin[3],
                            const Float_t spacing[3], Float_t pos[3])
{
  pos[0] = origin[0] + (static_cast<Float_t>(idx[0]) * spacing[0]);
  pos[1] = origin[1] + (static_cast<Float_t>(idx[1]) * spacing[1]);
  pos[2] = origin[2] + (static_cast<Float_t>(idx[2]) * spacing[2]);
}

static void computeGradient(const unsigned idx[3], const Float_t *buffer,
                            const unsigned dims[3], const Float_t spacing[3],
                            Float_t grad[3])
{
  unsigned flatIdx = idx[0] + (idx[1] * dims[0]) + (idx[2] * dims[0] * dims[1]);
  Float_t v1, v2, fac;

  v1 = (idx[0] > 0) ? buffer[flatIdx - 1] : buffer[flatIdx];
  v2 = (idx[0] < (dims[0] - 1)) ? buffer[flatIdx + 1] : buffer[flatIdx];
  fac = (idx[0] > 0 && idx[0] < (dims[0] - 1)) ? 2.0 : 1.0;
  grad[0] = (v1 - v2)/(fac * spacing[0]);

  v1 = (idx[1] > 0) ? buffer[flatIdx - dims[0]] : buffer[flatIdx];
  v2 = (idx[1] < (dims[1] - 1)) ? buffer[flatIdx + dims[0]] : buffer[flatIdx];
  fac = (idx[1] > 0 && idx[1] < (dims[1] - 1)) ? 2.0 : 1.0;
  grad[1] = (v1 - v2)/(fac * spacing[1]);

  v1 = (idx[2] > 0) ? buffer[flatIdx - (dims[0] * dims[1])] : buffer[flatIdx];
  v2 = (idx[2] < (dims[2] - 1)) ? buffer[flatIdx + (dims[0] * dims[1])]
                                : buffer[flatIdx];
  fac = (idx[2] > 0 && idx[2] < (dims[2] - 1)) ? 2.0 : 1.0;
  grad[2] = (v1 - v2)/(fac * spacing[2]);
}

inline Float_t lerp(Float_t a, Float_t b, Float_t w)
{
  return a + (w * (b - a));
}



namespace scalar_2 {

struct CellInfo
{
  unsigned idx[3];
  int caseId;
};

class CaseIdIsLess
{
public:
  bool operator()(const CellInfo &c1, const CellInfo &c2) const
  {
    return (c1.caseId < c2.caseId);
  }
};

class CellIsEmpty
{
public:
  bool operator()(const CellInfo &c) const
  {
    return (c.caseId == 0 || c.caseId == 255);
  }
};

struct EdgeInfo
{
  unsigned edgeId, pos;
  unsigned idx[4];
};

class EdgeIsInvalid
{
public:
  bool operator()(const EdgeInfo &e) const
  {
    return (e.edgeId == INVALID_EDGE);
  }
};

class EdgeIdIsLess
{
public:
  bool operator()(const EdgeInfo &e1, const EdgeInfo &e2) const
  {
    return (e1.edgeId < e2.edgeId);
  }
};

void extractIsosurfaceFromBlock(const Image3D_t &vol, const unsigned ext[6],
  Float_t isoval, TriangleMesh_t *mesh)
{
  static const int caseMask[] = { 1, 2, 4, 8, 16, 32, 64, 128 };
  static const unsigned idxOffset[8][3] = { { 0, 0, 0}, { 1, 0, 0},
                                            { 1, 1, 0}, { 0, 1, 0},
                                            { 0, 0, 1}, { 1, 0, 1},
                                            { 1, 1, 1}, { 0, 1, 1} };

  const unsigned *dims = vol.getDimension();
  const Float_t *origin = vol.getOrigin();
  const Float_t *spacing = vol.getSpacing();
  const Float_t *buffer = vol.getData();
  unsigned xydim = dims[0] * dims[1];

  unsigned edgeIdOffsets[12];
  generateEdgeIdOffsets(dims[0], xydim, edgeIdOffsets);

  // part 1: compute caseId of each cell
  unsigned ncells = (ext[1] - ext[0] + 1) * (ext[3] - ext[2] + 1) *
                    (ext[5] - ext[4] + 1);
  std::vector<CellInfo> cellInfo;
  cellInfo.reserve(ncells);

  for (unsigned zidx = ext[4]; zidx <= ext[5]; ++zidx)
    {
    for (unsigned yidx = ext[2]; yidx <= ext[3]; ++yidx)
      {
      for (unsigned xidx = ext[0]; xidx <= ext[1]; ++xidx)
        {
        Float_t val[8];
        unsigned idx = xidx + (yidx * dims[0]) + (zidx * xydim);

        val[0] = buffer[idx];
        val[1] = buffer[idx + 1];
        val[2] = buffer[idx + 1 + dims[0]];
        val[3] = buffer[idx + dims[0]];
        val[4] = buffer[idx + xydim];
        val[5] = buffer[idx + 1 + xydim];
        val[6] = buffer[idx + 1 + dims[0] + xydim];
        val[7] = buffer[idx + dims[0] + xydim];

        int caseId = 0;
        for (int i = 0; i < 8; ++i)
          {
          caseId |= (val[i] >= isoval) ? caseMask[i] : 0;
          }

        CellInfo ci = { { xidx, yidx, zidx}, caseId };
        cellInfo.push_back(ci);
        }
      }
    }

  // part 2: remove empty cells and sort by caseId
  std::vector<CellInfo>::iterator newEnd
    = std::remove_if(cellInfo.begin(), cellInfo.end(), CellIsEmpty());
  cellInfo.resize(newEnd - cellInfo.begin());
  std::sort(cellInfo.begin(), cellInfo.end(), CaseIdIsLess());

  // part 3: generate edges
  ncells = cellInfo.size();
  std::vector<EdgeInfo> edgeInfo(ncells * 15);
  for (unsigned i = 0, pos = 0; i < ncells; ++i)
    {
    unsigned idx = i * 15;

    const CellInfo &ci = cellInfo[i];
    unsigned baseEdge = (ci.idx[0] + (ci.idx[1] * dims[0]) +
                        (ci.idx[2] * xydim)) * 3;

    const int *edgeIds = MarchingCubesTables::getCaseTrianglesEdges(ci.caseId);
    for (int j = 0; j < 15; ++j, ++idx, ++edgeIds)
      {
      EdgeInfo &ei = edgeInfo[idx];
      if (*edgeIds == -1)
        {
        ei.edgeId = INVALID_EDGE;
        }
      else
        {
        ei.edgeId = baseEdge + edgeIdOffsets[*edgeIds];

        ei.idx[0] = ci.idx[0];
        ei.idx[1] = ci.idx[1];
        ei.idx[2] = ci.idx[2];
        ei.idx[3] = *edgeIds;

        ei.pos = pos++;
        }
      }
    }

  // part 4: remove invalid edges and sort by unique edgeId
  std::vector<EdgeInfo>::iterator newEnd1
    = std::remove_if(edgeInfo.begin(), edgeInfo.end(), EdgeIsInvalid());
  edgeInfo.resize(newEnd1 - edgeInfo.begin());
  std::sort(edgeInfo.begin(), edgeInfo.end(), EdgeIdIsLess());

  // part 5: merge duplicate edges, generate indexes
  std::vector<unsigned> indexes(edgeInfo.size());
  unsigned last = 0, ptsIdxOffset = mesh->points.size()/3;
  for (unsigned i = 0; i < edgeInfo.size(); ++i)
    {
    if (edgeInfo[i].edgeId != edgeInfo[last].edgeId)
      {
      if (++last != i)
        {
        edgeInfo[last] = edgeInfo[i];
        }
      }
    indexes[edgeInfo[i].pos] = last + ptsIdxOffset;
    }
  edgeInfo.resize(last + 1);

  mesh->indexes.insert(mesh->indexes.end(), indexes.begin(), indexes.end());

  // part 6: compute gradients and generate vertices
  for (unsigned i = 0;  i < edgeInfo.size(); ++i)
    {
    const EdgeInfo &ei = edgeInfo[i];

    unsigned p1idx[3], p2idx[3];
    int p1 = MarchingCubesTables::getEdgeVertices(ei.idx[3])[0];
    int p2 = MarchingCubesTables::getEdgeVertices(ei.idx[3])[1];
    p1idx[0] = ei.idx[0] + idxOffset[p1][0];
    p1idx[1] = ei.idx[1] + idxOffset[p1][1];
    p1idx[2] = ei.idx[2] + idxOffset[p1][2];
    p2idx[0] = ei.idx[0] + idxOffset[p2][0];
    p2idx[1] = ei.idx[1] + idxOffset[p2][1];
    p2idx[2] = ei.idx[2] + idxOffset[p2][2];

    unsigned p1FlatIdx = p1idx[0] + (p1idx[1] * dims[0]) + (p1idx[2] * xydim);
    unsigned p2FlatIdx = p2idx[0] + (p2idx[1] * dims[0]) + (p2idx[2] * xydim);
    Float_t w = (isoval - buffer[p1FlatIdx]) /
               (buffer[p2FlatIdx] - buffer[p1FlatIdx]);

    Float_t pos1[3], pos2[3], grad1[3], grad2[3];
    computePosition(p1idx, origin, spacing, pos1);
    computeGradient(p1idx, buffer, dims, spacing, grad1);
    computePosition(p2idx, origin, spacing, pos2);
    computeGradient(p2idx, buffer, dims, spacing, grad2);

    for (int ii = 0; ii < 3; ++ii)
      {
      mesh->points.push_back(lerp(pos1[ii], pos2[ii], w));
      mesh->normals.push_back(lerp(grad1[ii], grad2[ii], w));
      }
    }
}

class IsosurfaceFunctor
{
public:
  typedef smp::ThreadLocalStorage<TriangleMesh_t> TLS_tm;

  IsosurfaceFunctor(const Image3D_t &vol, Float_t isoval, TLS_tm &meshPieces)
    : input(vol), isoval(isoval), meshPieces(meshPieces)
  {
  }

  void operator()(const smp::Range3D &range) const
  {
    unsigned extent[6] = { range.cols().begin(), range.cols().end() - 1,
                           range.rows().begin(), range.rows().end() - 1,
                           range.pages().begin(), range.pages().end() - 1 };
    TriangleMesh_t &meshPiece = this->meshPieces.local();
    extractIsosurfaceFromBlock(this->input, extent, this->isoval, &meshPiece);
  }

private:
  const Image3D_t &input;
  Float_t isoval;
  TLS_tm &meshPieces;
};

void extractIsosurface(const Image3D_t &vol, Float_t isoval,
                       TriangleMesh_t *mesh)
{
  const unsigned *dims = vol.getDimension();
  const Float_t *origin = vol.getOrigin();
  const Float_t *spacing = vol.getSpacing();

  IsosurfaceFunctor::TLS_tm meshPieces;
  smp::Range3D cellRange(0, dims[2] - 1, grainDim,
                         0, dims[1] - 1, grainDim,
                         0, dims[0] - 1, grainDim);

  IsosurfaceFunctor func(vol, isoval, meshPieces);
  smp::parallel_for(cellRange, func);
  mergeTriangleMeshes(meshPieces.begin(), meshPieces.end(), mesh);
}

}; //namespace scalar_2

namespace simd_2 {

struct CellInfo
{
  unsigned idx[3];
  int caseId;
};

class CaseIdIsLess
{
public:
  bool operator()(const CellInfo &c1, const CellInfo &c2) const
  {
    return (c1.caseId < c2.caseId);
  }
};

class CellIsEmpty
{
public:
  bool operator()(const CellInfo &c) const
  {
    return (c.caseId == 0 || c.caseId == 255);
  }
};

struct EdgeInfo
{
  unsigned edgeId, pos;
  unsigned idx[4];
};

class EdgeIsInvalid
{
public:
  bool operator()(const EdgeInfo &e) const
  {
    return (e.edgeId == INVALID_EDGE);
  }
};

class EdgeIdIsLess
{
public:
  bool operator()(const EdgeInfo &e1, const EdgeInfo &e2) const
  {
    return (e1.edgeId < e2.edgeId);
  }
};

void extractIsosurfaceFromBlock(const Image3D_t &vol, const unsigned ext[6],
  Float_t isoval, TriangleMesh_t *mesh)
{
  static const int caseMask[] = { 1, 2, 4, 8, 16, 32, 64, 128 };
  static const unsigned idxOffset[8][3] = { { 0, 0, 0}, { 1, 0, 0},
                                            { 1, 1, 0}, { 0, 1, 0},
                                            { 0, 0, 1}, { 1, 0, 1},
                                            { 1, 1, 1}, { 0, 1, 1} };

  const unsigned *dims = vol.getDimension();
  const Float_t *origin = vol.getOrigin();
  const Float_t *spacing = vol.getSpacing();
  const Float_t *buffer = vol.getData();
  unsigned xydim = dims[0] * dims[1];

  unsigned edgeIdOffsets[12];
  generateEdgeIdOffsets(dims[0], xydim, edgeIdOffsets);

  // part 1: compute caseId of each cell
  unsigned ncells = (ext[1] - ext[0] + 1) * (ext[3] - ext[2] + 1) *
                    (ext[5] - ext[4] + 1);
  std::vector<int> caseIds(ncells);
  ispc::getCellCaseIds(buffer, dims, ext, isoval, &caseIds[0]);

  std::vector<CellInfo> cellInfo(ncells);
  for (unsigned z = ext[4], i = 0; z <= ext[5]; ++z)
    {
    for (unsigned y = ext[2]; y <= ext[3]; ++y)
      {
      for (unsigned x = ext[0]; x <= ext[1]; ++x, ++i)
        {
        cellInfo[i].idx[0] = x;
        cellInfo[i].idx[1] = y;
        cellInfo[i].idx[2] = z;
        cellInfo[i].caseId = caseIds[i];
        }
      }
    }

  // part 2: remove empty cells and sort by caseId
  std::vector<CellInfo>::iterator newEnd
    = std::remove_if(cellInfo.begin(), cellInfo.end(), CellIsEmpty());
  cellInfo.resize(newEnd - cellInfo.begin());
  std::sort(cellInfo.begin(), cellInfo.end(), CaseIdIsLess());

  // part 3: generate edges
  ncells = cellInfo.size();
  std::vector<EdgeInfo> edgeInfo(ncells * 15);
  for (unsigned i = 0, pos = 0; i < ncells; ++i)
    {
    unsigned idx = i * 15;

    const CellInfo &ci = cellInfo[i];
    unsigned baseEdge = (ci.idx[0] + (ci.idx[1] * dims[0]) +
                        (ci.idx[2] * xydim)) * 3;

    const int *edgeIds = MarchingCubesTables::getCaseTrianglesEdges(ci.caseId);
    for (int j = 0; j < 15; ++j, ++idx, ++edgeIds)
      {
      EdgeInfo &ei = edgeInfo[idx];
      if (*edgeIds == -1)
        {
        ei.edgeId = INVALID_EDGE;
        }
      else
        {
        ei.edgeId = baseEdge + edgeIdOffsets[*edgeIds];

        ei.idx[0] = ci.idx[0];
        ei.idx[1] = ci.idx[1];
        ei.idx[2] = ci.idx[2];
        ei.idx[3] = *edgeIds;

        ei.pos = pos++;
        }
      }
    }

  // part 4: remove invalid edges and sort by unique edgeId
  std::vector<EdgeInfo>::iterator newEnd1
    = std::remove_if(edgeInfo.begin(), edgeInfo.end(), EdgeIsInvalid());
  edgeInfo.resize(newEnd1 - edgeInfo.begin());
  std::sort(edgeInfo.begin(), edgeInfo.end(), EdgeIdIsLess());

  // part 5: merge duplicate edges, generate indexes
  std::vector<unsigned> indexes(edgeInfo.size());
  unsigned last = 0, ptsIdxOffset = mesh->points.size()/3;
  for (unsigned i = 0; i < edgeInfo.size(); ++i)
    {
    if (edgeInfo[i].edgeId != edgeInfo[last].edgeId)
      {
      if (++last != i)
        {
        edgeInfo[last] = edgeInfo[i];
        }
      }
    indexes[edgeInfo[i].pos] = last + ptsIdxOffset;
    }
  edgeInfo.resize(last + 1);

  mesh->indexes.insert(mesh->indexes.end(), indexes.begin(), indexes.end());

  // part 6: compute gradients and generate vertices
  for (unsigned i = 0;  i < edgeInfo.size(); ++i)
    {
    const EdgeInfo &ei = edgeInfo[i];

    unsigned p1idx[3], p2idx[3];
    int p1 = MarchingCubesTables::getEdgeVertices(ei.idx[3])[0];
    int p2 = MarchingCubesTables::getEdgeVertices(ei.idx[3])[1];
    p1idx[0] = ei.idx[0] + idxOffset[p1][0];
    p1idx[1] = ei.idx[1] + idxOffset[p1][1];
    p1idx[2] = ei.idx[2] + idxOffset[p1][2];
    p2idx[0] = ei.idx[0] + idxOffset[p2][0];
    p2idx[1] = ei.idx[1] + idxOffset[p2][1];
    p2idx[2] = ei.idx[2] + idxOffset[p2][2];

    unsigned p1FlatIdx = p1idx[0] + (p1idx[1] * dims[0]) + (p1idx[2] * xydim);
    unsigned p2FlatIdx = p2idx[0] + (p2idx[1] * dims[0]) + (p2idx[2] * xydim);
    Float_t w = (isoval - buffer[p1FlatIdx]) /
               (buffer[p2FlatIdx] - buffer[p1FlatIdx]);

    Float_t pos1[3], pos2[3], grad1[3], grad2[3];
    computePosition(p1idx, origin, spacing, pos1);
    computeGradient(p1idx, buffer, dims, spacing, grad1);
    computePosition(p2idx, origin, spacing, pos2);
    computeGradient(p2idx, buffer, dims, spacing, grad2);

    for (int ii = 0; ii < 3; ++ii)
      {
      mesh->points.push_back(lerp(pos1[ii], pos2[ii], w));
      mesh->normals.push_back(lerp(grad1[ii], grad2[ii], w));
      }
    }
}

class IsosurfaceFunctor
{
public:
  typedef smp::ThreadLocalStorage<TriangleMesh_t> TLS_tm;

  IsosurfaceFunctor(const Image3D_t &vol, Float_t isoval, TLS_tm &meshPieces)
    : input(vol), isoval(isoval), meshPieces(meshPieces)
  {
  }

  void operator()(const smp::Range3D &range) const
  {
    unsigned extent[6] = { range.cols().begin(), range.cols().end() - 1,
                           range.rows().begin(), range.rows().end() - 1,
                           range.pages().begin(), range.pages().end() - 1 };
    TriangleMesh_t &meshPiece = this->meshPieces.local();
    extractIsosurfaceFromBlock(this->input, extent, this->isoval, &meshPiece);
  }

private:
  const Image3D_t &input;
  Float_t isoval;
  TLS_tm &meshPieces;
};

void extractIsosurface(const Image3D_t &vol, Float_t isoval,
                       TriangleMesh_t *mesh)
{
  const unsigned *dims = vol.getDimension();
  const Float_t *origin = vol.getOrigin();
  const Float_t *spacing = vol.getSpacing();

  IsosurfaceFunctor::TLS_tm meshPieces;
  smp::Range3D cellRange(0, dims[2] - 1, grainDim,
                         0, dims[1] - 1, grainDim,
                         0, dims[0] - 1, grainDim);

  IsosurfaceFunctor func(vol, isoval, meshPieces);
  smp::parallel_for(cellRange, func);
  mergeTriangleMeshes(meshPieces.begin(), meshPieces.end(), mesh);
}

}; //namespace simd_2


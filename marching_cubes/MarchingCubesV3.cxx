#include "MarchingCubes.h"
#include "MarchingCubesTables.h"
#include "MarchingCubesV3.ispc.h"

#include <boost/chrono.hpp>

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

static void computePosition(const int idx[3], const Float_t origin[3],
                            const Float_t spacing[3], Float_t pos[3])
{
  pos[0] = origin[0] + (static_cast<Float_t>(idx[0]) * spacing[0]);
  pos[1] = origin[1] + (static_cast<Float_t>(idx[1]) * spacing[1]);
  pos[2] = origin[2] + (static_cast<Float_t>(idx[2]) * spacing[2]);
}

static void computeGradient(const int idx[3], const Float_t *buffer,
                            const int dims[3], const Float_t spacing[3],
                            Float_t grad[3])
{
  int flatIdx = idx[0] + (idx[1] * dims[0]) + (idx[2] * dims[0] * dims[1]);
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

static const char * partDesc[] = { "Begin", "Compute case ids", 
                                   "Remove empty cells and sort by caseId",
                                   "Generate edges", 
                                   "Remove invalid edges and sort by edge id",
                                   "Merge duplicate edges, generate indexes",
                                   "Compute gradients and generate vertices" };

namespace scalar_3 {

struct CellInfo
{
  int idx[3];
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
  int edgeId, pos;
  int idx[4];
};

class EdgeIsInvalid
{
public:
  bool operator()(const EdgeInfo &e) const
  {
    return (e.edgeId == -1);
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

void extractIsosurface(const Image3D_t &vol, Float_t isoval,
                       TriangleMesh_t *mesh)
{
  static const int caseMask[] = { 1, 2, 4, 8, 16, 32, 64, 128 };
  static const int everts[12][2] = {{0,1}, {1,2}, {3,2}, {0,3}, {4,5}, {5,6},
                                    {7,6}, {4,7}, {0,4}, {1,5}, {3,7}, {2,6}};
  static const int idxOffset[8][3] = { { 0, 0, 0}, { 1, 0, 0},
                                       { 1, 1, 0}, { 0, 1, 0},
                                       { 0, 0, 1}, { 1, 0, 1},
                                       { 1, 1, 1}, { 0, 1, 1} };

  const int *dims = vol.getDimension();
  const Float_t *origin = vol.getOrigin();
  const Float_t *spacing = vol.getSpacing();
  const Float_t *buffer = vol.getData();
  int xydim = dims[0] * dims[1];

  int edgeIdOffsets[12];
  generateEdgeIdOffsets(dims[0], xydim, edgeIdOffsets);

  Clock::time_point checkpoints[7];

  // part 1: compute caseId of each cell
  checkpoints[0] = Clock::now();

  int ncells = (dims[0] - 1) * (dims[1] - 1) * (dims[2] - 1);
  std::vector<CellInfo> cellInfo;
  cellInfo.reserve(ncells);

  for (int zidx = 0; zidx < (dims[2] - 1); ++zidx)
    {
    for (int yidx = 0; yidx < (dims[1] - 1); ++yidx)
      {
      for (int xidx = 0; xidx < (dims[0] - 1); ++xidx)
        {
        Float_t val[8];
        int idx = xidx + (yidx * dims[0]) + (zidx * xydim);

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
  checkpoints[1] = Clock::now();

  std::vector<CellInfo>::iterator newEnd
    = std::remove_if(cellInfo.begin(), cellInfo.end(), CellIsEmpty());
  cellInfo.resize(newEnd - cellInfo.begin());
  std::sort(cellInfo.begin(), cellInfo.end(), CaseIdIsLess());

  // part 3: generate edges
  checkpoints[2] = Clock::now();

  ncells = cellInfo.size();
  std::vector<EdgeInfo> edgeInfo(ncells * 15);
  for (int i = 0, pos = 0; i < ncells; ++i)
    {
    int idx = i * 15;

    const CellInfo &ci = cellInfo[i];
    int baseEdge = (ci.idx[0] + (ci.idx[1] * dims[0]) + (ci.idx[2] * xydim)) *
                   3;

    const int *edgeIds =
      MarchingCubesTables::getTriangleCases()[ci.caseId];
    for (int j = 0; j < 15; ++j, ++idx, ++edgeIds)
      {
      EdgeInfo &ei = edgeInfo[idx];
      if (*edgeIds == -1)
        {
        ei.edgeId = -1;
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
  checkpoints[3] = Clock::now();

  std::vector<EdgeInfo>::iterator newEnd1
    = std::remove_if(edgeInfo.begin(), edgeInfo.end(), EdgeIsInvalid());
  edgeInfo.resize(newEnd1 - edgeInfo.begin());
  std::sort(edgeInfo.begin(), edgeInfo.end(), EdgeIdIsLess());

  // part 5: merge duplicate edges, generate indexes
  checkpoints[4] = Clock::now();

  std::vector<int> indexes(edgeInfo.size());
  size_t last = 0;
  for (size_t i = 0; i < edgeInfo.size(); ++i)
    {
    if (edgeInfo[i].edgeId != edgeInfo[last].edgeId)
      {
      if (++last != i)
        {
        edgeInfo[last] = edgeInfo[i];
        }
      }
    indexes[edgeInfo[i].pos] = last;
    }
  edgeInfo.resize(last + 1);

  // part 6: compute gradients and generate vertices
  checkpoints[5] = Clock::now();

  std::vector<Float_t> points, normals;
  points.reserve(edgeInfo.size() * 3);
  normals.reserve(edgeInfo.size() * 3);

  for (size_t i = 0;  i < edgeInfo.size(); ++i)
    {
    const EdgeInfo &ei = edgeInfo[i];

    int p1idx[3], p2idx[3];
    int p1 = everts[ei.idx[3]][0], p2 = everts[ei.idx[3]][1];
    p1idx[0] = ei.idx[0] + idxOffset[p1][0];
    p1idx[1] = ei.idx[1] + idxOffset[p1][1];
    p1idx[2] = ei.idx[2] + idxOffset[p1][2];
    p2idx[0] = ei.idx[0] + idxOffset[p2][0];
    p2idx[1] = ei.idx[1] + idxOffset[p2][1];
    p2idx[2] = ei.idx[2] + idxOffset[p2][2];

    int p1FlatIdx = p1idx[0] + (p1idx[1] * dims[0]) + (p1idx[2] * xydim);
    int p2FlatIdx = p2idx[0] + (p2idx[1] * dims[0]) + (p2idx[2] * xydim);
    Float_t w = (isoval - buffer[p1FlatIdx]) /
               (buffer[p2FlatIdx] - buffer[p1FlatIdx]);

    Float_t pos1[3], pos2[3], grad1[3], grad2[3];
    computePosition(p1idx, origin, spacing, pos1);
    computeGradient(p1idx, buffer, dims, spacing, grad1);
    computePosition(p2idx, origin, spacing, pos2);
    computeGradient(p2idx, buffer, dims, spacing, grad2);

    for (int ii = 0; ii < 3; ++ii)
      {
      points.push_back(lerp(pos1[ii], pos2[ii], w));
      normals.push_back(lerp(grad1[ii], grad2[ii], w));
      }
    }

  checkpoints[6] = Clock::now();

  mesh->points.swap(points);
  mesh->normals.swap(normals);
  mesh->indexes.swap(indexes);

  boost::chrono::duration<double> total = checkpoints[6] - checkpoints[0];
  std::cout << "Timings: (total = " << total.count() << ")" << std::endl;
  for (int i = 1; i < 7; ++i)
    {
    boost::chrono::duration<double> dur = checkpoints[i] - checkpoints[i - 1];
    double pc = (100.0 * dur.count())/total.count();
    std::cout << "Part " << i << " - " << partDesc[i] << ": "
              << dur.count() << " seconds, " << pc << "%" << std::endl;
    }
}

}; //namespace scalar_3

namespace mixed_3 {

using ispc::CellInfo;

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
  int edgeId, pos;
  int idx[4];
};

class EdgeIsInvalid
{
public:
  bool operator()(const EdgeInfo &e) const
  {
    return (e.edgeId == -1);
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

void extractIsosurface(const Image3D_t &vol, Float_t isoval,
                       TriangleMesh_t *mesh)
{
  static const int caseMask[] = { 1, 2, 4, 8, 16, 32, 64, 128 };
  static const int everts[12][2] = {{0,1}, {1,2}, {3,2}, {0,3}, {4,5}, {5,6},
                                    {7,6}, {4,7}, {0,4}, {1,5}, {3,7}, {2,6}};
  static const int idxOffset[8][3] = { { 0, 0, 0}, { 1, 0, 0},
                                       { 1, 1, 0}, { 0, 1, 0},
                                       { 0, 0, 1}, { 1, 0, 1},
                                       { 1, 1, 1}, { 0, 1, 1} };

  const int *dims = vol.getDimension();
  const Float_t *origin = vol.getOrigin();
  const Float_t *spacing = vol.getSpacing();
  const Float_t *buffer = vol.getData();
  int xydim = dims[0] * dims[1];

  int edgeIdOffsets[12];
  generateEdgeIdOffsets(dims[0], xydim, edgeIdOffsets);

  Clock::time_point checkpoints[7];

  // part 1: compute caseId of each cell
  checkpoints[0] = Clock::now();

  int ncells = (dims[0] - 1) * (dims[1] - 1) * (dims[2] - 1);
  std::vector<CellInfo> cellInfo(ncells);

  ispc::getCellCaseIds(buffer, dims, isoval, &cellInfo[0]);


  // part 2: remove empty cells and sort by caseId
  checkpoints[1] = Clock::now();

  std::vector<CellInfo>::iterator newEnd
    = std::remove_if(cellInfo.begin(), cellInfo.end(), CellIsEmpty());
  cellInfo.resize(newEnd - cellInfo.begin());
  std::sort(cellInfo.begin(), cellInfo.end(), CaseIdIsLess());

  // part 3: generate edges
  checkpoints[2] = Clock::now();

  ncells = cellInfo.size();
  std::vector<EdgeInfo> edgeInfo(ncells * 15);
  for (int i = 0, pos = 0; i < ncells; ++i)
    {
    int idx = i * 15;

    const CellInfo &ci = cellInfo[i];
    int baseEdge = (ci.idx[0] + (ci.idx[1] * dims[0]) + (ci.idx[2] * xydim)) *
                   3;

    const int *edgeIds =
      MarchingCubesTables::getTriangleCases()[ci.caseId];
    for (int j = 0; j < 15; ++j, ++idx, ++edgeIds)
      {
      EdgeInfo &ei = edgeInfo[idx];
      if (*edgeIds == -1)
        {
        ei.edgeId = -1;
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
  checkpoints[3] = Clock::now();

  std::vector<EdgeInfo>::iterator newEnd1
    = std::remove_if(edgeInfo.begin(), edgeInfo.end(), EdgeIsInvalid());
  edgeInfo.resize(newEnd1 - edgeInfo.begin());
  std::sort(edgeInfo.begin(), edgeInfo.end(), EdgeIdIsLess());

  // part 5: merge duplicate edges, generate indexes
  checkpoints[4] = Clock::now();

  std::vector<int> indexes(edgeInfo.size());
  size_t last = 0;
  for (size_t i = 0; i < edgeInfo.size(); ++i)
    {
    if (edgeInfo[i].edgeId != edgeInfo[last].edgeId)
      {
      if (++last != i)
        {
        edgeInfo[last] = edgeInfo[i];
        }
      }
    indexes[edgeInfo[i].pos] = last;
    }
  edgeInfo.resize(last + 1);

  // part 6: compute gradients and generate vertices
  checkpoints[5] = Clock::now();

  std::vector<Float_t> points, normals;
  points.reserve(edgeInfo.size() * 3);
  normals.reserve(edgeInfo.size() * 3);

  for (size_t i = 0;  i < edgeInfo.size(); ++i)
    {
    const EdgeInfo &ei = edgeInfo[i];

    int p1idx[3], p2idx[3];
    int p1 = everts[ei.idx[3]][0], p2 = everts[ei.idx[3]][1];
    p1idx[0] = ei.idx[0] + idxOffset[p1][0];
    p1idx[1] = ei.idx[1] + idxOffset[p1][1];
    p1idx[2] = ei.idx[2] + idxOffset[p1][2];
    p2idx[0] = ei.idx[0] + idxOffset[p2][0];
    p2idx[1] = ei.idx[1] + idxOffset[p2][1];
    p2idx[2] = ei.idx[2] + idxOffset[p2][2];

    int p1FlatIdx = p1idx[0] + (p1idx[1] * dims[0]) + (p1idx[2] * xydim);
    int p2FlatIdx = p2idx[0] + (p2idx[1] * dims[0]) + (p2idx[2] * xydim);
    Float_t w = (isoval - buffer[p1FlatIdx]) /
               (buffer[p2FlatIdx] - buffer[p1FlatIdx]);

    Float_t pos1[3], pos2[3], grad1[3], grad2[3];
    computePosition(p1idx, origin, spacing, pos1);
    computeGradient(p1idx, buffer, dims, spacing, grad1);
    computePosition(p2idx, origin, spacing, pos2);
    computeGradient(p2idx, buffer, dims, spacing, grad2);

    for (int ii = 0; ii < 3; ++ii)
      {
      points.push_back(lerp(pos1[ii], pos2[ii], w));
      normals.push_back(lerp(grad1[ii], grad2[ii], w));
      }
    }

  checkpoints[6] = Clock::now();

  mesh->points.swap(points);
  mesh->normals.swap(normals);
  mesh->indexes.swap(indexes);

  boost::chrono::duration<double> total = checkpoints[6] - checkpoints[0];
  std::cout << "Timings: (total = " << total.count() << ")" << std::endl;
  for (int i = 1; i < 7; ++i)
    {
    boost::chrono::duration<double> dur = checkpoints[i] - checkpoints[i - 1];
    double pc = (100.0 * dur.count())/total.count();
    std::cout << "Part " << i << " - " << partDesc[i] << ": "
              << dur.count() << " seconds, " << pc << "%" << std::endl;
    }
}

}; //namespace mixed_3


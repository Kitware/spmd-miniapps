#include "type.h"

typedef unsigned int uint;

static const int caseMask[] = { 1, 2, 4, 8, 16, 32, 64, 128 };


export uniform uint getCellCases(uniform const Float_t buffer[],
                                 uniform const uint dims[3],
                                 uniform const uint cellExtent[6],
                                 uniform const uint ncells[3],
                                 uniform Float_t isoval,
                                 uniform int cellCase[])
{
  uniform uint ptxyCount = dims[0] * dims[1];
  uniform uint cellxyCount = ncells[0] * ncells[1];

  uint numberOfActiveCells = 0;
  foreach (zc = cellExtent[4] ... cellExtent[5] + 1,
           yc = cellExtent[2] ... cellExtent[3] + 1,
           xc = cellExtent[0] ... cellExtent[1] + 1)
    {
    uint ptIdx = zc * ptxyCount + yc * dims[0] + xc;

    Float_t val[8];
    val[0] = buffer[ptIdx];
    val[1] = buffer[ptIdx + 1];
    val[2] = buffer[ptIdx + 1 + dims[0]];
    val[3] = buffer[ptIdx + dims[0]];
    val[4] = buffer[ptIdx + ptxyCount];
    val[5] = buffer[ptIdx + 1 + ptxyCount];
    val[6] = buffer[ptIdx + 1 + dims[0] + ptxyCount];
    val[7] = buffer[ptIdx + dims[0] + ptxyCount];

    int caseId = 0;
    for (int i = 0; i < 8; ++i)
      {
      caseId |= (val[i] >= isoval) ? caseMask[i] : 0;
      }

    uint cellIdx = ((zc - cellExtent[4]) * cellxyCount) +
                   ((yc - cellExtent[2]) * ncells[0]) +
                   (xc - cellExtent[0]);
    cellCase[cellIdx] = caseId;

    if (caseId != 0 && caseId != 255)
      {
      ++numberOfActiveCells;
      }
    }

  return reduce_add(numberOfActiveCells);
}

inline void computeGradient(varying uint xidx, varying uint yidx,
                            varying uint zidx,
                            uniform const Float_t *buffer,
                            uniform const uint dims[3],
                            uniform const Float_t spacing[3],
                            varying Float_t grad[3])
{
  uniform uint ptxyCount = dims[0] * dims[1];
  uint ptIdx = xidx + (yidx * dims[0]) + (zidx * ptxyCount);

  if (xidx == 0)
    {
    Float_t x1 = buffer[ptIdx + 1];
    Float_t x2 = buffer[ptIdx];
    grad[0] = (x2 - x1)/spacing[0];
    }
  else if (xidx == (dims[0] - 1))
    {
    Float_t x1 = buffer[ptIdx];
    Float_t x2 = buffer[ptIdx - 1];
    grad[0] = (x2 - x1)/spacing[0];
    }
  else
    {
    Float_t x1 = buffer[ptIdx + 1];
    Float_t x2 = buffer[ptIdx - 1];
    grad[0] = (0.5 * (x2 - x1))/spacing[0];
    }

  if (yidx == 0)
    {
    Float_t y1 = buffer[ptIdx + dims[0]];
    Float_t y2 = buffer[ptIdx];
    grad[1] = (y2 - y1)/spacing[1];
    }
  else if (yidx == (dims[1] - 1))
    {
    Float_t y1 = buffer[ptIdx];
    Float_t y2 = buffer[ptIdx - dims[0]];
    grad[1] = (y2 - y1)/spacing[1];
    }
  else
    {
    Float_t y1 = buffer[ptIdx + dims[0]];
    Float_t y2 = buffer[ptIdx - dims[0]];
    grad[1] = (0.5 * (y2 - y1))/spacing[1];
    }

  if (zidx == 0)
    {
    Float_t z1 = buffer[ptIdx + ptxyCount];
    Float_t z2 = buffer[ptIdx];
    grad[2] = (z2 - z1)/spacing[2];
    }
  else if (zidx == (dims[2] - 1))
    {
    Float_t z1 = buffer[ptIdx];
    Float_t z2 = buffer[ptIdx - ptxyCount];
    grad[2] = (z2 - z1)/spacing[2];
    }
  else
    {
    Float_t z1 = buffer[ptIdx + ptxyCount];
    Float_t z2 = buffer[ptIdx - ptxyCount];
    grad[2] = (0.5 * (z2 - z1))/spacing[2];
    }
}

inline Float_t lerp(varying Float_t a, varying Float_t b, varying Float_t w)
{
  return a + (w * (b - a));
}

export void extractIsosurface_impl(uniform const Float_t buffer[],
                                   uniform const uint dims[3],
                                   uniform const Float_t origin[3],
                                   uniform const Float_t spacing[3],
                                   uniform Float_t isoval,
                                   uniform const int triangleCases[],
                                   uniform const int edgeVerts[],
                                   uniform uint numberOfCells,
                                   uniform const uint cellIndex[],
                                   uniform const int caseIds[],
                                   uniform const uint vertIndex[],
                                   uniform Float_t points[],
                                   uniform Float_t normals[])
{
  uniform uint ptxyCount = dims[0] * dims[1];
  uniform const uint idxOff[8] = { 0, 1, dims[0] + 1, dims[0], ptxyCount,
                                   1 + ptxyCount, 1 + dims[0] + ptxyCount,
                                   dims[0] + ptxyCount };

  uniform uint ncells[3] = { dims[0] - 1, dims[1] - 1, dims[2] - 1 };
  uniform uint cellxyCount = ncells[0] * ncells[1];

  foreach (i = 0 ... numberOfCells)
    {
    uint cell = cellIndex[i];

    uint zc = cell/cellxyCount;
    uint yc = (cell%cellxyCount)/ncells[0];
    uint xc = (cell%cellxyCount)%ncells[0];

    uint ptIdx = zc * ptxyCount + yc * dims[0] + xc;

    Float_t pos[8][3], grad[8][3];

    pos[0][0] = origin[0] + (((Float_t)xc) * spacing[0]);
    pos[0][1] = origin[1] + (((Float_t)yc) * spacing[1]);
    pos[0][2] = origin[2] + (((Float_t)zc) * spacing[2]);

    pos[1][0] = pos[0][0] + spacing[0];
    pos[1][1] = pos[0][1];
    pos[1][2] = pos[0][2];

    pos[2][0] = pos[0][0] + spacing[0];
    pos[2][1] = pos[0][1] + spacing[1];
    pos[2][2] = pos[0][2];

    pos[3][0] = pos[0][0];
    pos[3][1] = pos[0][1] + spacing[1];
    pos[3][2] = pos[0][2];

    pos[4][0] = pos[0][0];
    pos[4][1] = pos[0][1];
    pos[4][2] = pos[0][2] + spacing[2];

    pos[5][0] = pos[0][0] + spacing[0];
    pos[5][1] = pos[0][1];
    pos[5][2] = pos[0][2] + spacing[2];

    pos[6][0] = pos[0][0] + spacing[0];
    pos[6][1] = pos[0][1] + spacing[1];
    pos[6][2] = pos[0][2] + spacing[2];

    pos[7][0] = pos[0][0];
    pos[7][1] = pos[0][1] + spacing[1];
    pos[7][2] = pos[0][2] + spacing[2];

    computeGradient(xc, yc, zc, buffer, dims, spacing, grad[0]);
    computeGradient(xc + 1, yc, zc, buffer, dims, spacing, grad[1]);
    computeGradient(xc + 1, yc + 1, zc, buffer, dims, spacing, grad[2]);
    computeGradient(xc, yc + 1, zc, buffer, dims, spacing, grad[3]);
    computeGradient(xc, yc, zc + 1, buffer, dims, spacing, grad[4]);
    computeGradient(xc + 1, yc, zc + 1, buffer, dims, spacing, grad[5]);
    computeGradient(xc + 1, yc + 1, zc + 1, buffer, dims, spacing, grad[6]);
    computeGradient(xc, yc + 1, zc + 1, buffer, dims, spacing, grad[7]);

    const int *edges = triangleCases + (caseIds[i] * 16);
    uint vidx = vertIndex[i];
    for (; *edges != -1; ++edges, vidx += 3)
      {
      int v1 = edgeVerts[*edges * 2 + 0], v2 = edgeVerts[*edges * 2 + 1];
      Float_t val1 = buffer[ptIdx + idxOff[v1]];
      Float_t val2 = buffer[ptIdx + idxOff[v2]];
      Float_t w = (isoval - val1)/(val2 - val1);

      for (int ii = 0; ii < 3; ++ii)
        {
        points[vidx + ii] = lerp(pos[v1][ii], pos[v2][ii], w);
        normals[vidx + ii] = lerp(grad[v1][ii], grad[v2][ii], w);
        }
      }
    }
}


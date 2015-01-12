#include "type.h"

static const int caseMask[] = { 1, 2, 4, 8, 16, 32, 64, 128 };
static const int edgeVerts[12][2] = {{0,1}, {1,2}, {3,2}, {0,3}, {4,5}, {5,6},
                                     {7,6}, {4,7}, {0,4}, {1,5}, {3,7}, {2,6}};


export uniform int getCellCases(uniform const Float_t buffer[],
                                uniform const int dims[3],
                                uniform Float_t isoval, uniform int cellCase[])
{
  uniform int ptxyCount = dims[0] * dims[1];

  uniform int ncells[3] = { dims[0] - 1, dims[1] - 1, dims[2] - 1 };
  uniform int cellxyCount = ncells[0] * ncells[1];

  int numberOfActiveCells = 0;
  foreach (zc = 0 ... ncells[2], yc = 0 ... ncells[1], xc = 0 ... ncells[0])
    {
    int cellIdx = zc * cellxyCount + yc * ncells[0] + xc;
    int ptIdx = zc * ptxyCount + yc * dims[0] + xc;

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

    cellCase[cellIdx] = caseId;

    if (caseId != 0 && caseId != 255)
      {
      ++numberOfActiveCells;
      }
    }

  return reduce_add(numberOfActiveCells);
}

inline void computeGradient(varying int xidx, varying int yidx,
                            varying int zidx,
                            uniform const Float_t *buffer,
                            uniform const int dims[3],
                            uniform const Float_t spacing[3],
                            varying Float_t grad[3])
{
  uniform int ptxyCount = dims[0] * dims[1];
  int ptIdx = xidx + (yidx * dims[0]) + (zidx * ptxyCount);

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
                                   uniform const int dims[3],
                                   uniform const Float_t origin[3],
                                   uniform const Float_t spacing[3],
                                   uniform Float_t isoval,
                                   uniform const int triangleCases[],
                                   uniform int numberOfCells,
                                   uniform const int cellIndex[],
                                   uniform const int caseIds[],
                                   uniform const int vertIndex[],
                                   uniform Float_t points[],
                                   uniform Float_t normals[])
{
  uniform int ptxyCount = dims[0] * dims[1];

  uniform int ncells[3] = { dims[0] - 1, dims[1] - 1, dims[2] - 1 };
  uniform int cellxyCount = ncells[0] * ncells[1];

  foreach (i = 0 ... numberOfCells)
    {
    int cell = cellIndex[i];

    int zc = cell/cellxyCount;
    int yc = (cell%cellxyCount)/ncells[0];
    int xc = (cell%cellxyCount)%ncells[0];

    int ptIdx = zc * ptxyCount + yc * dims[0] + xc;

    Float_t val[8], pos[8][3], grad[8][3];

    val[0] = buffer[ptIdx];
    val[1] = buffer[ptIdx + 1];
    val[2] = buffer[ptIdx + 1 + dims[0]];
    val[3] = buffer[ptIdx + dims[0]];
    val[4] = buffer[ptIdx + ptxyCount];
    val[5] = buffer[ptIdx + 1 + ptxyCount];
    val[6] = buffer[ptIdx + 1 + dims[0] + ptxyCount];
    val[7] = buffer[ptIdx + dims[0] + ptxyCount];

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
    int vidx = vertIndex[i];
    for (; *edges != -1; ++edges, vidx += 3)
      {
      int v1 = edgeVerts[*edges][0], v2 = edgeVerts[*edges][1];
      Float_t w = (isoval - val[v1])/(val[v2] - val[v1]);

      for (int ii = 0; ii < 3; ++ii)
        {
        points[vidx + ii] = lerp(pos[v1][ii], pos[v2][ii], w);
        normals[vidx + ii] = lerp(grad[v1][ii], grad[v2][ii], w);
        }
      }
    }
}


#include "type.h"

static const int caseMask[] = { 1, 2, 4, 8, 16, 32, 64, 128 };

export void getCellCaseIds(uniform const Float_t buffer[],
                           uniform const int dims[3],
                           uniform const int cellExt[6],
                           uniform Float_t isoval,
                           uniform int caseIds[])
{
  uniform int xydim = dims[0] * dims[1];
  uniform int ncellsxyz[3] = { cellExt[1] - cellExt[0] + 1,
                               cellExt[3] - cellExt[2] + 1,
                               cellExt[5] - cellExt[4] + 1 };
  uniform int ncellsxy = ncellsxyz[0] * ncellsxyz[1];

  foreach (z = 0 ... ncellsxyz[2], y = 0 ... ncellsxyz[1],
           x = 0 ... ncellsxyz[0])
    {
    int xidx = x + cellExt[0];
    int yidx = y + cellExt[2];
    int zidx = z + cellExt[4];
    int idx = xidx + (yidx * dims[0]) + (zidx * xydim);

    Float_t val[8];
    val[0] = buffer[idx];
    val[1] = buffer[idx + 1];
    val[2] = buffer[idx + 1 + dims[0]];
    val[3] = buffer[idx + dims[0]];
    val[4] = buffer[idx + xydim];
    val[5] = buffer[idx + 1 + xydim];
    val[6] = buffer[idx + 1 + dims[0] + xydim];
    val[7] = buffer[idx + dims[0] + xydim];

    int caseId = 0;
    for (uniform int i = 0; i < 8; ++i)
      {
      caseId |= (val[i] >= isoval) ? caseMask[i] : 0;
      }

    caseIds[x + (y * ncellsxyz[0]) + (z * ncellsxy)] = caseId;
    }
}

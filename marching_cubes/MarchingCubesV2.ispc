#include "type.h"

#define GANG_SIZE (programCount)

struct CellInfo
{
  int idx[3];
  int caseId;
};

struct CellInfo_soa
{
  int idx[3][GANG_SIZE];
  int caseId[GANG_SIZE];
};

static const int caseMask[] = { 1, 2, 4, 8, 16, 32, 64, 128 };

export void getCellCaseIds(uniform const Float_t buffer[],
                           uniform const int dims[3], uniform Float_t isoval,
                           uniform CellInfo cellInfo[])
{
  uniform int xydim = dims[0] * dims[1];
  uniform int cellDims[3] = { (dims[0] - 1), (dims[1] - 1), (dims[2] - 1) };
  uniform int xycells = cellDims[0] * cellDims[1];
  uniform int ncells = cellDims[0] * cellDims[1] * cellDims[2];

  uniform int nallocs = (ncells + GANG_SIZE - 1)/GANG_SIZE;
  uniform CellInfo_soa * uniform cis = uniform new uniform CellInfo_soa[nallocs];

  for (uniform int zidx = 0; zidx < cellDims[2]; ++zidx)
    {
    for (uniform int yidx = 0; yidx < cellDims[1]; ++yidx)
      {
      uniform int cellIdx = (yidx * cellDims[0]) + (zidx * xycells);
      foreach (xidx = 0 ... cellDims[0])
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

        uniform int index1 = (cellIdx + extract(xidx, 0))/GANG_SIZE;
        int index2 = programIndex;

        cis[index1].idx[0][index2] = xidx;
        cis[index1].idx[1][index2] = yidx;
        cis[index1].idx[2][index2] = zidx;
        cis[index1].caseId[index2] = caseId;
        }
      }
    }

  foreach (i = 0 ... ncells)
    {
    uniform unsigned int i0 = extract(i, 0);
    uniform unsigned int index1 = i0/GANG_SIZE;
    unsigned int index2 = programIndex;

#if 0
    cellInfo[i].idx[0] = cis[index1].idx[0][index2];
    cellInfo[i].idx[1] = cis[index1].idx[1][index2];
    cellInfo[i].idx[2] = cis[index1].idx[2][index2];
    cellInfo[i].caseId = cis[index1].caseId[index2];
#else
    uniform int * uniform cip = (uniform int*)(cellInfo + i0);
    soa_to_aos4(cis[index1].idx[0][index2], cis[index1].idx[1][index2],
                cis[index1].idx[2][index2], cis[index1].caseId[index2],
                cip);
#endif
    }

  delete[] cis;
}


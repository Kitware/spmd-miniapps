#include "type.h"

typedef unsigned int uint;

export inline void computeGradient(uniform uint xidx[],
                                   uniform uint yidx[],
                                   uniform uint zidx[],
                                   uniform const Float_t buffer[],
                                   uniform const uint dims[3],
                                   uniform const Float_t spacing[3],
                                   uniform Float_t grad[][3])
{
  uniform uint xysize = dims[0] * dims[1];

  foreach (i = 0 ... 8) {
    if (xidx[i] == 0)
      {
      Float_t x1 = buffer[(xidx[i] + 1) + (yidx[i] * dims[0]) +
                         (zidx[i] * xysize)];
      Float_t x2 = buffer[xidx[i] + (yidx[i] * dims[0]) + (zidx[i] * xysize)];
      grad[i][0] = (x2 - x1)/spacing[0];
      }
    else if (xidx[i] == (dims[0] - 1))
      {
      Float_t x1 = buffer[xidx[i] + (yidx[i] * dims[0]) + (zidx[i] * xysize)];
      Float_t x2 = buffer[(xidx[i] - 1) + (yidx[i] * dims[0]) +
                         (zidx[i] * xysize)];
      grad[i][0] = (x2 - x1)/spacing[0];
      }
    else
      {
      Float_t x1 = buffer[(xidx[i] + 1) + (yidx[i] * dims[0]) +
                         (zidx[i] * xysize)];
      Float_t x2 = buffer[(xidx[i] - 1) + (yidx[i] * dims[0]) +
                         (zidx[i] * xysize)];
      grad[i][0] = (0.5 * (x2 - x1))/spacing[0];
      }

    if (yidx[i] == 0)
      {
      Float_t y1 = buffer[xidx[i] + ((yidx[i] + 1) * dims[0]) +
                         (zidx[i] * xysize)];
      Float_t y2 = buffer[xidx[i] + (yidx[i] * dims[0]) + (zidx[i] * xysize)];
      grad[i][1] = (y2 - y1)/spacing[1];
      }
    else if (yidx[i] == (dims[1] - 1))
      {
      Float_t y1 = buffer[xidx[i] + (yidx[i] * dims[0]) + (zidx[i] * xysize)];
      Float_t y2 = buffer[xidx[i] + ((yidx[i] - 1) * dims[0]) +
                         (zidx[i] * xysize)];
      grad[i][1] = (y2 - y1)/spacing[1];
      }
    else
      {
      Float_t y1 = buffer[xidx[i] + ((yidx[i] + 1) * dims[0]) +
                         (zidx[i] * xysize)];
      Float_t y2 = buffer[xidx[i] + ((yidx[i] - 1) * dims[0]) +
                         (zidx[i] * xysize)];
      grad[i][1] = (0.5 * (y2 - y1))/spacing[1];
      }

    if (zidx[i] == 0)
      {
      Float_t z1 = buffer[xidx[i] + (yidx[i] * dims[0]) +
                         ((zidx[i] + 1) * xysize)];
      Float_t z2 = buffer[xidx[i] + (yidx[i] * dims[0]) + (zidx[i] * xysize)];
      grad[i][2] = (z2 - z1)/spacing[2];
      }
    else if (zidx[i] == (dims[2] - 1))
      {
      Float_t z1 = buffer[xidx[i] + (yidx[i] * dims[0]) + (zidx[i] * xysize)];
      Float_t z2 = buffer[xidx[i] + (yidx[i] * dims[0]) +
                         ((zidx[i] - 1) * xysize)];
      grad[i][2] = (z2 - z1)/spacing[2];
      }
    else
      {
      Float_t z1 = buffer[xidx[i] + (yidx[i] * dims[0]) +
                         ((zidx[i] + 1) * xysize)];
      Float_t z2 = buffer[xidx[i] + (yidx[i] * dims[0]) +
                         ((zidx[i] - 1) * xysize)];
      grad[i][2] = (0.5 * (z2 - z1))/spacing[2];
      }
  }
}

export inline void lerp(uniform int len,
                        uniform const Float_t as[], uniform const Float_t bs[],
                        uniform Float_t w, uniform Float_t result[])
{
  foreach (i = 0 ... len) {
    result[i] = as[i] + (w * (bs[i] - as[i]));
  }
}


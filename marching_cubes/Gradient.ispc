
export inline void computeGradient(uniform int xidx[], uniform int yidx[],
                                   uniform int zidx[],
                                   uniform const double buffer[],
                                   uniform const int dims[3],
                                   uniform const double spacing[3],
                                   uniform double grad[][3])
{
  uniform int xysize = dims[0] * dims[1];

  foreach (i = 0 ... 8) {
    if (xidx[i] == 0)
      {
      double x1 = buffer[(xidx[i] + 1) + (yidx[i] * dims[0]) +
                         (zidx[i] * xysize)];
      double x2 = buffer[xidx[i] + (yidx[i] * dims[0]) + (zidx[i] * xysize)];
      grad[i][0] = (x2 - x1)/spacing[0];
      }
    else if (xidx[i] == (dims[0] - 1))
      {
      double x1 = buffer[xidx[i] + (yidx[i] * dims[0]) + (zidx[i] * xysize)];
      double x2 = buffer[(xidx[i] - 1) + (yidx[i] * dims[0]) +
                         (zidx[i] * xysize)];
      grad[i][0] = (x2 - x1)/spacing[0];
      }
    else
      {
      double x1 = buffer[(xidx[i] + 1) + (yidx[i] * dims[0]) +
                         (zidx[i] * xysize)];
      double x2 = buffer[(xidx[i] - 1) + (yidx[i] * dims[0]) +
                         (zidx[i] * xysize)];
      grad[i][0] = (0.5 * (x2 - x1))/spacing[0];
      }

    if (yidx[i] == 0)
      {
      double y1 = buffer[xidx[i] + ((yidx[i] + 1) * dims[0]) +
                         (zidx[i] * xysize)];
      double y2 = buffer[xidx[i] + (yidx[i] * dims[0]) + (zidx[i] * xysize)];
      grad[i][1] = (y2 - y1)/spacing[1];
      }
    else if (yidx[i] == (dims[1] - 1))
      {
      double y1 = buffer[xidx[i] + (yidx[i] * dims[0]) + (zidx[i] * xysize)];
      double y2 = buffer[xidx[i] + ((yidx[i] - 1) * dims[0]) +
                         (zidx[i] * xysize)];
      grad[i][1] = (y2 - y1)/spacing[1];
      }
    else
      {
      double y1 = buffer[xidx[i] + ((yidx[i] + 1) * dims[0]) +
                         (zidx[i] * xysize)];
      double y2 = buffer[xidx[i] + ((yidx[i] - 1) * dims[0]) +
                         (zidx[i] * xysize)];
      grad[i][1] = (0.5 * (y2 - y1))/spacing[1];
      }

    if (zidx[i] == 0)
      {
      double z1 = buffer[xidx[i] + (yidx[i] * dims[0]) +
                         ((zidx[i] + 1) * xysize)];
      double z2 = buffer[xidx[i] + (yidx[i] * dims[0]) + (zidx[i] * xysize)];
      grad[i][2] = (z2 - z1)/spacing[2];
      }
    else if (zidx[i] == (dims[2] - 1))
      {
      double z1 = buffer[xidx[i] + (yidx[i] * dims[0]) + (zidx[i] * xysize)];
      double z2 = buffer[xidx[i] + (yidx[i] * dims[0]) +
                         ((zidx[i] - 1) * xysize)];
      grad[i][2] = (z2 - z1)/spacing[2];
      }
    else
      {
      double z1 = buffer[xidx[i] + (yidx[i] * dims[0]) +
                         ((zidx[i] + 1) * xysize)];
      double z2 = buffer[xidx[i] + (yidx[i] * dims[0]) +
                         ((zidx[i] - 1) * xysize)];
      grad[i][2] = (0.5 * (z2 - z1))/spacing[2];
      }
  }
}

export inline void lerp(uniform int len,
                        uniform const double as[], uniform const double bs[],
                        uniform double w, uniform double result[])
{
  foreach (i = 0 ... len) {
    result[i] = as[i] + (w * (bs[i] - as[i]));
  }
}


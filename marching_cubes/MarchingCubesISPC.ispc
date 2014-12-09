static const int caseMask[] = { 1, 2, 4, 8, 16, 32, 64, 128 };
static const int edgeVerts[12][2] = {{0,1}, {1,2}, {3,2}, {0,3}, {4,5}, {5,6},
                                     {7,6}, {4,7}, {0,4}, {1,5}, {3,7}, {2,6}};


export uniform int getCellCases(uniform const double buffer[],
                                uniform const int dims[3],
                                uniform double isoval, uniform int cellCase[])
{
  uniform int ptxyCount = dims[0] * dims[1];

  uniform int ncells[3] = { dims[0] - 1, dims[1] - 1, dims[2] - 1 };
  uniform int cellxyCount = ncells[0] * ncells[1];

  int numberOfActiveCells = 0;
  foreach (zc = 0 ... ncells[2], yc = 0 ... ncells[1], xc = 0 ... ncells[0])
    {
    int cellIdx = zc * cellxyCount + yc * ncells[0] + xc;
    int ptIdx = zc * ptxyCount + yc * dims[0] + xc;

    double val[8];
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
                            uniform const double *buffer,
                            uniform const int dims[3],
                            uniform const double spacing[3],
                            varying double grad[3])
{
  uniform int ptxyCount = dims[0] * dims[1];
  int ptIdx = xidx + (yidx * dims[0]) + (zidx * ptxyCount);

  if (xidx == 0)
    {
    double x1 = buffer[ptIdx + 1];
    double x2 = buffer[ptIdx];
    grad[0] = (x2 - x1)/spacing[0];
    }
  else if (xidx == (dims[0] - 1))
    {
    double x1 = buffer[ptIdx];
    double x2 = buffer[ptIdx - 1];
    grad[0] = (x2 - x1)/spacing[0];
    }
  else
    {
    double x1 = buffer[ptIdx + 1];
    double x2 = buffer[ptIdx - 1];
    grad[0] = (0.5 * (x2 - x1))/spacing[0];
    }

  if (yidx == 0)
    {
    double y1 = buffer[ptIdx + dims[0]];
    double y2 = buffer[ptIdx];
    grad[1] = (y2 - y1)/spacing[1];
    }
  else if (yidx == (dims[1] - 1))
    {
    double y1 = buffer[ptIdx];
    double y2 = buffer[ptIdx - dims[0]];
    grad[1] = (y2 - y1)/spacing[1];
    }
  else
    {
    double y1 = buffer[ptIdx + dims[0]];
    double y2 = buffer[ptIdx - dims[0]];
    grad[1] = (0.5 * (y2 - y1))/spacing[1];
    }

  if (zidx == 0)
    {
    double z1 = buffer[ptIdx + ptxyCount];
    double z2 = buffer[ptIdx];
    grad[2] = (z2 - z1)/spacing[2];
    }
  else if (zidx == (dims[2] - 1))
    {
    double z1 = buffer[ptIdx];
    double z2 = buffer[ptIdx - ptxyCount];
    grad[2] = (z2 - z1)/spacing[2];
    }
  else
    {
    double z1 = buffer[ptIdx + ptxyCount];
    double z2 = buffer[ptIdx - ptxyCount];
    grad[2] = (0.5 * (z2 - z1))/spacing[2];
    }
}

inline double lerp(varying double a, varying double b, varying double w)
{
  return a + (w * (b - a));
}

struct Vertex
{
  double pos[3];
  double norm[3];
};

export void extractIsosurface_impl(uniform const double buffer[],
                                   uniform const int dims[3],
                                   uniform const double origin[3],
                                   uniform const double spacing[3],
                                   uniform double isoval,
                                   uniform const int triangleCases[],
                                   uniform int numberOfCells,
                                   uniform const int cellIndex[],
                                   uniform const int caseIds[],
                                   uniform const int vertIndex[],
                                   uniform Vertex verts[])
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

    double val[8], pos[8][3], grad[8][3];

    val[0] = buffer[ptIdx];
    val[1] = buffer[ptIdx + 1];
    val[2] = buffer[ptIdx + 1 + dims[0]];
    val[3] = buffer[ptIdx + dims[0]];
    val[4] = buffer[ptIdx + ptxyCount];
    val[5] = buffer[ptIdx + 1 + ptxyCount];
    val[6] = buffer[ptIdx + 1 + dims[0] + ptxyCount];
    val[7] = buffer[ptIdx + dims[0] + ptxyCount];

    pos[0][0] = origin[0] + (((double)xc) * spacing[0]);
    pos[0][1] = origin[1] + (((double)yc) * spacing[1]);
    pos[0][2] = origin[2] + (((double)zc) * spacing[2]);

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
    for (; *edges != -1; ++edges, ++vidx)
      {
      int v1 = edgeVerts[*edges][0], v2 = edgeVerts[*edges][1];
      double w = (isoval - val[v1])/(val[v2] - val[v1]);

      verts[vidx].pos[0] = lerp(pos[v1][0], pos[v2][0], w);
      verts[vidx].pos[1] = lerp(pos[v1][1], pos[v2][1], w);
      verts[vidx].pos[2] = lerp(pos[v1][2], pos[v2][2], w);

      verts[vidx].norm[0] = lerp(grad[v1][0], grad[v2][0], w);
      verts[vidx].norm[1] = lerp(grad[v1][1], grad[v2][1], w);
      verts[vidx].norm[2] = lerp(grad[v1][2], grad[v2][2], w);
      }
    }
}


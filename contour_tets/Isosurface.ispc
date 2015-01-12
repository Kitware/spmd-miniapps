#include "type.h"

enum CONSTANTS
{
  GANG_SIZE = programCount,
  EDGES_PER_CASE = 7
};

struct TetInfo_soa
{
  int tetIdx[GANG_SIZE];
  int caseId[GANG_SIZE];
  int ptidx[4][GANG_SIZE];
  Float_t val[4][GANG_SIZE];
  Float_t pts[4][3][GANG_SIZE];
};

inline Float_t lerp(varying Float_t a, varying Float_t b, varying Float_t w)
{
  return a + (w * (b - a));
}

inline void computeTriangleNormal(varying const Float_t verts[3][3],
                                  varying Float_t normal[3])
{
  Float_t v2v1[3] = { verts[1][0] - verts[0][0],
                      verts[1][1] - verts[0][1],
                      verts[1][2] - verts[0][2] };
  Float_t v3v1[3] = { verts[2][0] - verts[0][0],
                      verts[2][1] - verts[0][1],
                      verts[2][2] - verts[0][2] };

  normal[0] = (v2v1[1] * v3v1[2]) - (v2v1[2] * v3v1[1]);
  normal[1] = (v2v1[2] * v3v1[0]) - (v2v1[0] * v3v1[2]);
  normal[2] = (v2v1[0] * v3v1[1]) - (v2v1[1] * v3v1[0]);
}

export void extractIsosurface_impl(uniform const TetInfo_soa tetinfo[],
                                   uniform int ntets,
                                   uniform Float_t isoval,
                                   uniform const int triangleEdges[],
                                   uniform const int edgeVertices[],
                                   uniform Float_t vertices[],
                                   uniform int vertKeys[],
                                   uniform Float_t normals[])
{
  uniform int triInd = 0;
  foreach (i = 0 ... ntets)
    {
    uniform int gang = extract(i, 0)/GANG_SIZE;

    int caseId = tetinfo[gang].caseId[programIndex];
    foreach_unique(u_caseId in caseId)
      {
      uniform const int * uniform triEdges = &triangleEdges[u_caseId *
                                                            EDGES_PER_CASE];
      for (; *triEdges != -1; triEdges += 3)
        {
        Float_t verts[3][3];
        int keys[3][2];
        for (uniform int e = 0; e < 3; ++e)
          {
          uniform int v1 = edgeVertices[triEdges[e] * 2 + 0];
          uniform int v2 = edgeVertices[triEdges[e] * 2 + 1];

          Float_t val1 = tetinfo[gang].val[v1][programIndex];
          Float_t val2 = tetinfo[gang].val[v2][programIndex];

          if (val1 < val2)
            {
            Float_t dv = val2 - val1;
            Float_t w = (isoval - val1)/dv;
            for (uniform int ii = 0; ii < 3; ++ii)
              {
              verts[e][ii] = lerp(tetinfo[gang].pts[v1][ii][programIndex],
                                  tetinfo[gang].pts[v2][ii][programIndex], w);
              }
            keys[e][0] = tetinfo[gang].ptidx[v1][programIndex];
            keys[e][1] = tetinfo[gang].ptidx[v2][programIndex];
            }
          else
            {
            Float_t dv = val1 - val2;
            Float_t w = (dv == 0.0) ? 0.0 : (isoval - val2)/dv;
            for (uniform int ii = 0; ii < 3; ++ii)
              {
              verts[e][ii] = lerp(tetinfo[gang].pts[v2][ii][programIndex],
                                  tetinfo[gang].pts[v1][ii][programIndex], w);
              }
            keys[e][0] = tetinfo[gang].ptidx[v2][programIndex];
            keys[e][1] = tetinfo[gang].ptidx[v1][programIndex];
            }
          }

        Float_t normal[3];
        computeTriangleNormal(verts, normal);

        foreach_active (p)
          {
          uniform int vind = triInd * 3 * 3;
          uniform int kind = triInd * 3 * 2;
          uniform int nind = ((triInd/GANG_SIZE) * GANG_SIZE * 3) +
                             (triInd % GANG_SIZE);
          for (uniform int ii = 0; ii < 3; ++ii)
            {
            vertices[vind++] = extract(verts[ii][0], p);
            vertices[vind++] = extract(verts[ii][1], p);
            vertices[vind++] = extract(verts[ii][2], p);
            vertKeys[kind++] = extract(keys[ii][0], p);
            vertKeys[kind++] = extract(keys[ii][1], p);

            normals[nind + (ii * GANG_SIZE)] = extract(normal[ii], p);
            }
          ++triInd;
          }
        }
      }
    }
}

export void normalizeNormals(uniform Float_t normals[], uniform int count)
{
  foreach (i = 0 ... count)
    {
    uniform int gang = extract(i, 0)/GANG_SIZE;
    uniform int sind = gang * GANG_SIZE * 3;

    Float_t norm[3];
    for (uniform int ii = 0; ii < 3; ++ii)
      {
      norm[ii] = normals[sind + (ii * GANG_SIZE) + programIndex];
      }

    Float_t mag = sqrt((norm[0] * norm[0]) + (norm[1] * norm[1]) +
                       (norm[2] * norm[2]));
    if (mag != 0)
      {
      norm[0] /= mag;
      norm[1] /= mag;
      norm[2] /= mag;
      }

    int dind = i * 3;
    foreach_active (p)
      {
      uniform int u_dind = extract(dind, p);
      normals[u_dind + 0] = extract(norm[0], p);
      normals[u_dind + 1] = extract(norm[1], p);
      normals[u_dind + 2] = extract(norm[2], p);
      }
    }
}


extern "C" void addTriangleToOutput(void* uniform trimeshHandle, 
                                    uniform Float_t verts[3][3],
                                    uniform int keys[3][2],
                                    uniform Float_t normal[3]);

export void extractIosurface_impl_v2(uniform const Float_t tetPoints[],
                                     uniform const Float_t tetScalars[],
                                     uniform const int tetIndexes[],
                                     uniform int numTetPoints,
                                     uniform int numTetCells,
                                     uniform Float_t isoval,
                                     uniform const int triangleEdges[],
                                     uniform const int edgeVertices[],
                                     void* uniform trimeshHandle)
{
  static uniform const int caseMask[4] = { 1, 2, 4, 8 };

  foreach (i = 0 ... numTetCells)
    {
    int ptinds[4];
    Float_t val[4];
    int caseId = 0;
    for (uniform int ii = 0; ii < 4; ++ii)
      {
      ptinds[ii] = tetIndexes[i * 4 + ii];
      val[ii] = tetScalars[ptinds[ii]];
      if (val[ii] >= isoval)
        {
        caseId |= caseMask[ii];
        }
      }

    const uniform int *triEdges = &triangleEdges[caseId * EDGES_PER_CASE];
    for (; *triEdges != -1; triEdges += 3)
      {
      Float_t verts[3][3];
      Float_t normal[3];
      int keys[3][2];

      for (uniform int e = 0; e < 3; ++e)
        {
        int v1 = edgeVertices[triEdges[e] * 2];
        int v2 = edgeVertices[triEdges[e] * 2 + 1];

        Float_t val1 = val[v1];
        Float_t val2 = val[v2];

        Float_t w = 0.0;
        if (val1 < val2)
          {
          Float_t dv = val2 - val1;
          w = (isoval - val1)/dv;
          keys[e][0] = ptinds[v1];
          keys[e][1] = ptinds[v2];
          }
        else
          {
          Float_t dv = val1 - val2;
          w = (isoval - val2)/dv;
          keys[e][0] = ptinds[v2];
          keys[e][1] = ptinds[v1];
          }
        for (uniform int ii = 0; ii < 3; ++ii)
          {
          verts[e][ii] = lerp(tetPoints[keys[e][0] * 3 + ii],
                              tetPoints[keys[e][1] * 3 + ii], w);
          }
        }
      computeTriangleNormal(verts, normal);

      foreach_active (p)
        {
        uniform Float_t v[3][3];
        uniform int k[3][2];
        uniform Float_t n[3];
        for (uniform int ii = 0; ii < 3; ++ii)
          {
          v[ii][0] = extract(verts[ii][0], p);
          v[ii][1] = extract(verts[ii][1], p);
          v[ii][2] = extract(verts[ii][2], p);
          k[ii][0] = extract(keys[ii][0], p);
          k[ii][1] = extract(keys[ii][1], p);
          n[ii] = extract(normal[ii], p);
          }
        addTriangleToOutput(trimeshHandle, v, k, n);
        }
      }
    }
}



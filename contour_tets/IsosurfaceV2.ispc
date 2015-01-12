
export void extractIosurface_impl(uniform const Float_t tetPoints[],
                                  uniform const Float_t tetScalars[],
                                  uniform const int tetIndexes[],
                                  uniform int numTetPoints,
                                  uniform int numTetCells,
                                  uniform Float_t isoval,
                                  uniform const int triangleEdges[],
                                  uniform const int edgeVertices[],
                                  uniform void * uniform trimeshHandle)
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

    const uniform int *triEdges = &triangleEdges[caseId * 7];
    for (; *triEdges != -1; triEdges += 3)
      {
      Float_t verts[3][3];
      Float_t normal[3];
      Float_t keys[3][2];

      for (uniform int e = 0; e < 3; ++e)
        {
        int v1 = edgeVertices[triEdges[e] * 2];
        int v2 = edgeVertices[triEdges[e] * 2 + 1];

        Float_t val1 = val[v1];
        Float_t val2 = val[v2];

        if (val1 < val2)
          {
          Float_t dv = val2 - val1;
          Float_t w = (isoval - val1)/dv;
          keys[e][0] = ptinds[v1];
          keys[e][1] = ptinds[v2];
          }
        else
          {
          Float_t dv = val1 - val2;
          Float_t w = (isoval - val2)/dv;
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
        addTriangleToOutput(trimeshHandle, extract(verts[ii], p),
                            extract(keys[ii], p), extract(normal, p));
        }
      }
    }
}


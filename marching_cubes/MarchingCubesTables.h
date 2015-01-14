#ifndef __MarchingCubesTables_h
#define __MarchingCubesTables_h

class MarchingCubesTables
{
public:
  static const int* getEdgeVertices(int edgeId);
  static const int* getCaseTrianglesEdges(int caseId);
  static int getNumberOfTriangles(int caseId);

private:
  static const int edgeVertices[12][2];
  static const int caseTrianglesEdges[256][16];
  static const int numberOfTriangles[256];
};

inline const int* MarchingCubesTables::getEdgeVertices(int edgeId)
{
  return edgeVertices[edgeId];
}

inline const int* MarchingCubesTables::getCaseTrianglesEdges(int caseId)
{
  return caseTrianglesEdges[caseId];
}

inline int MarchingCubesTables::getNumberOfTriangles(int caseId)
{
  return numberOfTriangles[caseId];
}

#endif


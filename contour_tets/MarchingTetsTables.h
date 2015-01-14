#ifndef __MarchingTetsTables_h
#define __MarchingTetsTables_h

class MarchingTetsTables
{
public:
  static const int* getEdgeVertices(int edgeId);
  static const int* getCaseTrianglesEdges(int caseId);
  static int getNumberOfTriangles(int caseId);

private:
  static const int edgeVertices[6][2];
  static const int caseTrianglesEdges[16][7];
  static const int numberOfTriangles[16];
};

inline const int* MarchingTetsTables::getEdgeVertices(int edgeId)
{
  return edgeVertices[edgeId];
}

inline const int* MarchingTetsTables::getCaseTrianglesEdges(int caseId)
{
  return caseTrianglesEdges[caseId];
}

inline int MarchingTetsTables::getNumberOfTriangles(int caseId)
{
  return numberOfTriangles[caseId];
}

#endif


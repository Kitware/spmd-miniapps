#ifndef __MarchingCubesTables_h
#define __MarchingCubesTables_h

class MarchingCubesTables
{
public:
  typedef int EdgeId;
  typedef EdgeId EdgeList[16];

  static const EdgeList* getTriangleCases();
  static const int* getTrianglesPerCase();
};

#endif


#ifndef __TriangleMesh3D_h
#define __TriangleMesh3D_h

#include <vector>

struct TriangleMesh3D
{
public:
  struct Vertex
  {
    double pos[3];
    double norm[3];
  };

  std::vector<Vertex> vertices;
  std::vector<int> indexes;
};

#endif


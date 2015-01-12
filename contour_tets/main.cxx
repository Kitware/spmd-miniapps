#include "Isosurface.h"

#include <LoadTetrahedronMesh.h>
#include <SaveTriangleMesh.h>
#include <TetrahedronMesh.h>
#include <TriangleMesh.h>

#include <boost/chrono.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <string>

typedef boost::chrono::steady_clock Clock;

static const char* impstr[] = { "scalar", "vectorized", "vectorized_v2" };

enum Implementations
{
  IMP_SCALAR = 0,
  IMP_VECTORIZED,
  IMP_VECTORIZED_2,
  NUM_IMPLEMENTATIONS
};

int main(int argc, char* argv[])
{
  if (argc != 5)
    {
    std::cout << "Usage:" << std::endl;
    std::cout << argv[0]
              << " in_tetmesh.vtk out_poly.vtk isoval <implementation>"
              << std::endl;
    std::cout << "Implementations:";
    for (unsigned i = 0; i < NUM_IMPLEMENTATIONS; ++i)
      {
      std::cout << " " << impstr[i];
      }
    std::cout << std::endl;
    return 1;
    }

  std::string opt(argv[4]);
  Implementations imp = NUM_IMPLEMENTATIONS;
  for (int i = 0; i < NUM_IMPLEMENTATIONS; ++i)
    {
    if (opt == impstr[i])
      {
      imp = static_cast<Implementations>(i);
      }
    }
  if (imp == NUM_IMPLEMENTATIONS)
    {
    std::cout << "Invalid implementation" << std::endl;
    return 1;
    }

  TetrahedronMesh_t tetmesh;
  loadTetrahedronMesh(argv[1], &tetmesh);

  TriangleMesh_t trimesh;

  Float_t isoval = boost::lexical_cast<Float_t>(argv[3]);

  Clock::time_point start, finish;
  start = Clock::now();
  switch (imp)
    {
    case IMP_SCALAR:
      scalar::extractIsosurface(tetmesh, isoval, &trimesh);
      break;
    case IMP_VECTORIZED:
      vectorized::extractIsosurface(tetmesh, isoval, &trimesh);
      break;
    case IMP_VECTORIZED_2:
      vectorized_2::extractIsosurface(tetmesh, isoval, &trimesh);
      break;
    default:
      break;
    }
  finish = Clock::now();
  boost::chrono::duration<double> sec = finish - start;

  std::cout << "nverts: " << trimesh.numberOfVertices() << std::endl;
  std::cout << "ntris: " << trimesh.numberOfTriangles() << std::endl;
  std::cout << "done in " << sec.count() << " seconds\n";

  saveTriangleMesh(trimesh, argv[2]);

  return 0;
}


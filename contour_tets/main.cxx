#include "Isosurface.h"

#include <LoadTetrahedronMesh.h>
#include <SaveTriangleMesh.h>
#include <TetrahedronMesh.h>
#include <TriangleMesh.h>

#include <tbb/task_scheduler_init.h>

#include <boost/chrono.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <string>

typedef boost::chrono::steady_clock Clock;

static const char* impstr[] = { "scalar", "simd", "simd_2" };

enum ImplementationId
{
  IMP_SCALAR = 0,
  IMP_SIMD,
  IMP_SIMD_2,
  NUM_IMPLEMENTATIONS
};

ImplementationId getImplementationId(const char *str)
{
  std::string name(str);
  ImplementationId imp = NUM_IMPLEMENTATIONS;
  for (int i = 0; i < NUM_IMPLEMENTATIONS; ++i)
    {
    if (name == impstr[i])
      {
      imp = static_cast<ImplementationId>(i);
      }
    }

  return imp;
}

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

  ImplementationId choice = getImplementationId(argv[4]);
  if (choice == NUM_IMPLEMENTATIONS)
    {
    std::cout << "Invalid implementation" << std::endl;
    return 1;
    }

  tbb::task_scheduler_init tbbInit(tbb::task_scheduler_init::deferred);
  const char *omp_num_threads_str = getenv("OMP_NUM_THREADS");
  if (omp_num_threads_str)
    {
    int nthreads = boost::lexical_cast<int>(omp_num_threads_str);
    tbbInit.initialize(nthreads);
    }

  TetrahedronMesh_t tetmesh;
  loadTetrahedronMesh(argv[1], &tetmesh);

  TriangleMesh_t trimesh;

  Float_t isoval = boost::lexical_cast<Float_t>(argv[3]);

  Clock::time_point start, finish;
  start = Clock::now();
  switch (choice)
    {
    case IMP_SCALAR:
      scalar::extractIsosurface(tetmesh, isoval, &trimesh);
      break;
    case IMP_SIMD:
      simd::extractIsosurface(tetmesh, isoval, &trimesh);
      break;
    case IMP_SIMD_2:
      simd_2::extractIsosurface(tetmesh, isoval, &trimesh);
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


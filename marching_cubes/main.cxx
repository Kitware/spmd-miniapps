#include "MarchingCubes.h"

#include <Image3D.h>
#include <LoadImage3D.h>
#include <SaveTriangleMesh.h>
#include <TriangleMesh.h>

#if USE_TBB_BACKEND
#include <tbb/task_scheduler_init.h>
#include <boost/thread.hpp>
#endif
#if USE_OMP_BACKEND
#include <omp.h>
#endif

#include <boost/chrono.hpp>
#include <boost/lexical_cast.hpp>

#include <iostream>
#include <string>

typedef boost::chrono::steady_clock Clock;

static const char * impstr[] = { "scalar", "shortvec", "simd", "scalar_2",
                                 "simd_2", "scalar_2.1", "simd_2.1" };

enum ImplementationId
{
  IMP_SCALAR = 0,
  IMP_SHORTVEC,
  IMP_SIMD,
  IMP_SCALAR_2,
  IMP_SIMD_2,
  IMP_SCALAR_2_1,
  IMP_SIMD_2_1,
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

int NumberOfThreads = 1;

int main(int argc, char* argv[])
{
  if (argc != 5)
    {
    std::cout << "Usage:" << std::endl;
    std::cout << argv[0]
              << " in_image.vtk out_poly.vtk isoval <implementation>"
              << std::endl;
    std::cout << "Implementations:";
    for (int i = 0; i < NUM_IMPLEMENTATIONS; ++i)
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
    }

#if USE_TBB_BACKEND
  tbb::task_scheduler_init tbbInit(tbb::task_scheduler_init::deferred);
  const char *omp_num_threads_str = getenv("OMP_NUM_THREADS");
  if (omp_num_threads_str)
    {
    NumberOfThreads = boost::lexical_cast<int>(omp_num_threads_str);
    tbbInit.initialize(NumberOfThreads);
    }
  else
    {
    NumberOfThreads = std::max(1u, boost::thread::hardware_concurrency());
    }
  std::cout << "Using tbb backend" << std::endl;
#else
  NumberOfThreads = omp_get_max_threads();
  std::cout << "Using omp backend" << std::endl;
#endif

  Image3D_t volume;
  loadImage3D(argv[1], &volume);

  Float_t isoval = boost::lexical_cast<Float_t>(argv[3]);

  TriangleMesh_t mesh;

  Clock::time_point start, finish;
  start = Clock::now();
  switch (choice)
    {
    case IMP_SCALAR:
      scalar::extractIsosurface(volume, isoval, &mesh);
      break;
    case IMP_SHORTVEC:
      shortvec::extractIsosurface(volume, isoval, &mesh);
      break;
    case IMP_SIMD:
      simd::extractIsosurface(volume, isoval, &mesh);
      break;
    case IMP_SCALAR_2:
      scalar_2::extractIsosurface(volume, isoval, &mesh);
      break;
    case IMP_SIMD_2:
      simd_2::extractIsosurface(volume, isoval, &mesh);
      break;
    case IMP_SCALAR_2_1:
      scalar_2_1::extractIsosurface(volume, isoval, &mesh);
      break;
    case IMP_SIMD_2_1:
      simd_2_1::extractIsosurface(volume, isoval, &mesh);
      break;
    default:
      break;
    }
  finish = Clock::now();

  boost::chrono::duration<double> sec = finish - start;

  std::cout << "nverts: " << mesh.numberOfVertices() << std::endl;
  std::cout << "ntris: " << mesh.numberOfTriangles() << std::endl;
  std::cout << "done in " << sec.count() << " seconds\n";

#if !SCALING_TEST_BUILD
  saveTriangleMesh(mesh, argv[2]);
#endif

  return 0;
}


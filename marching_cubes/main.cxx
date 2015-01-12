#include "MarchingCubes.h"

#include <Image3D.h>
#include <LoadImage3D.h>
#include <SaveTriangleMesh.h>
#include <TriangleMesh.h>

#include <boost/chrono.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>

typedef boost::chrono::steady_clock Clock;

int main(int argc, char* argv[])
{
  if (argc != 5)
    {
    std::cout << "Usage:" << std::endl;
    std::cout << argv[0]
              << " in_image.vtk out_poly.vtk isoval <implementation>"
              << std::endl;
    std::cout << "Implementations: scalar, mixed, vectorized, scalar_2, "
              << "scalar_3, mixed_3, scalar_3.1, mixed_3.1" << std::endl;
    return 1;
    }

  Image3D_t volume;
  loadImage3D(argv[1], &volume);

  TriangleMesh_t mesh;

  std::string opt(argv[4]);
  Clock::time_point start, finish;

  Float_t isoval = boost::lexical_cast<Float_t>(argv[3]);

  if (opt == "scalar")
    {
    start = Clock::now();
    scalar::extractIsosurface(volume, isoval, &mesh);
    finish = Clock::now();
    }
  else if (opt == "mixed")
    {
    start = Clock::now();
    mixed::extractIsosurface(volume, isoval, &mesh);
    finish = Clock::now();
    }
  else if (opt == "vectorized")
    {
    start = Clock::now();
    vectorized::extractIsosurface(volume, isoval, &mesh);
    finish = Clock::now();
    }
  else if (opt == "scalar_2")
    {
    start = Clock::now();
    scalar_2::extractIsosurface(volume, isoval, &mesh);
    finish = Clock::now();
    }
  else if (opt == "scalar_3")
    {
    start = Clock::now();
    scalar_3::extractIsosurface(volume, isoval, &mesh);
    finish = Clock::now();
    }
  else if (opt == "mixed_3")
    {
    start = Clock::now();
    mixed_3::extractIsosurface(volume, isoval, &mesh);
    finish = Clock::now();
    }
  else if (opt == "scalar_3.1")
    {
    start = Clock::now();
    scalar_3_1::extractIsosurface(volume, isoval, &mesh);
    finish = Clock::now();
    }
  else if (opt == "mixed_3.1")
    {
    start = Clock::now();
    mixed_3_1::extractIsosurface(volume, isoval, &mesh);
    finish = Clock::now();
    }
  else
    {
    std::cout << "Invalid option\n";
    return 1;
    }

  boost::chrono::duration<double> sec = finish - start;

  std::cout << "nverts: " << mesh.numberOfVertices() << std::endl;
  std::cout << "ntris: " << mesh.numberOfTriangles() << std::endl;
  std::cout << "done in " << sec.count() << " seconds\n";

  saveTriangleMesh(mesh, argv[2]);

  return 0;
}


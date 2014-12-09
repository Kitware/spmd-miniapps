#include "Image3D.h"
#include "LoadImage3D.h"
#include "MarchingCubes.h"
#include "SaveTriangleMesh3D.h"
#include "TriangleMesh3D.h"

#include <boost/chrono.hpp>

#include <iostream>

typedef boost::chrono::steady_clock Clock;

int main(int argc, char* argv[])
{
  if (argc != 2)
    {
    std::cout << "Usage: " << argv[0] << " <serial|simd|ispc>" << std::endl;
    return 1;
    }

  Image3D volume;
  loadImage3D("/home/sujinphilip/Temp/PlasticSkull.vti", &volume);

  TriangleMesh3D mesh;

  std::string opt(argv[1]);
  Clock::time_point start, finish;

  if (opt == "serial")
    {
    start = Clock::now();
    serial::extractIsosurface(volume, 979, &mesh);
    finish = Clock::now();
    }
  else if (opt == "simd")
    {
    start = Clock::now();
    simd::extractIsosurface(volume, 979, &mesh);
    finish = Clock::now();
    }
  else if (opt == "ispc")
    {
    start = Clock::now();
    ispc::extractIsosurface(volume, 979, &mesh);
    finish = Clock::now();
    }
  else
    {
    std::cout << "Invalid option\n";
    return 1;
    }

  boost::chrono::duration<double> sec = finish - start;

  std::cout << "nverts: " << mesh.vertices.size() << std::endl;
  std::cout << "ntris: " << mesh.indexes.size() << "/3 = "
            << mesh.indexes.size()/3 << std::endl;
  std::cout << "done in " << sec.count() << " seconds\n";

  saveTriangleMesh3D(mesh, "/home/sujinphilip/Temp/PlasticSkull.vtp");

  return 0;
}


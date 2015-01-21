#include "poisson.h"
#include "type.h"

#include <Image2D.h>
#include <pgm.h>

#include <boost/chrono.hpp>

#include <iostream>
#include <string>

typedef boost::chrono::steady_clock Clock;

int main(int argc, char *argv[])
{
  enum {
    IMPL_SERIAL, IMPL_SIMD
  };

  if (argc != 4)
    {
    std::cout << "usage:\n"
              << argv[0] << "in.pgm out.pgm <implementations>\n"
              << "Implementations: serial, simd\n";
    return 1;
    }

  std::string optstr(argv[3]);
  int opt = IMPL_SERIAL;
  if (optstr == "serial")
    {
    opt = IMPL_SERIAL;
    }
  else if (optstr == "simd")
    {
    opt = IMPL_SIMD;
    }
  else
    {
    std::cout << "Invalid options\n";
    return 1;
    }

  Image2D<Float_t> image;
  load_pgm(argv[1], &image);

  int xdim = image.getXDim();
  int ydim = image.getYDim();
  Image2D<Float_t> divImg(xdim, ydim);
  Image2D<Float_t> solImg(xdim, ydim);

  scalar::computeDivergence(image.getBuffer(), xdim, ydim, divImg.getBuffer());

  Clock::time_point start, finish;

  start = Clock::now();
  switch (opt)
    {
    case IMPL_SERIAL:
      scalar::solvePoisson(divImg.getBuffer(), xdim, ydim, solImg.getBuffer());
      break;
    case IMPL_SIMD:
      simd::solvePoisson(divImg.getBuffer(), xdim, ydim, solImg.getBuffer());
      break;
    }
  finish = Clock::now();

  boost::chrono::duration<double> sec = finish - start;
  std::cout << "done in " << sec.count() << " seconds\n";

  write_pgm(argv[2], solImg);

  return 0;
}


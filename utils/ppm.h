#ifndef __ppm_h
#define __ppm_h

#include "Image2D.h"

#include <fstream>
#include <iostream>
#include <string>

template <typename PType>
void load_ppm(const char *fname, Image2D<PType> *img)
{
  std::ifstream ppm;
  ppm.open(fname);

  std::string magic;
  int xdim, ydim, maxColor;

  ppm >> magic >> xdim >> ydim >> maxColor >> std::ws;

  int npixels = xdim * ydim;
  char *readBuffer = new char[npixels];
  ppm.read(readBuffer, npixels);
  ppm.close();

  img->allocate(xdim, ydim);
  PType *buffer = img->getBuffer();
  for (int j = 0, idx = 0; j < ydim; ++j, buffer += img->getPitch())
    {
    for (int i = 0; i < xdim; ++i, ++idx)
      {
      buffer[i] = static_cast<PType>(readBuffer[idx]);
      }
    }

  delete [] readBuffer;
}

template <typename PType>
void write_ppm(const char *fname, const Image2D<PType> &img)
{
  int npixels = img.getXDim() * img.getYDim();
  char *writeBuffer = new char[npixels];

  const PType *buffer = img.getBuffer();
  for (int j = 0, idx = 0; j < img.getYDim(); ++j, buffer += img.getPitch())
    {
    for (int i = 0; i < img.getXDim(); ++i, ++idx)
      {
      writeBuffer[idx] = static_cast<char>(buffer[i]);
      }
    }

  std::ofstream ppm;
  ppm.open(fname);

  ppm << "P5\n" << img.getXDim() << " " << img.getYDim() << "\n255\n";
  ppm.write(writeBuffer, npixels);
  ppm.close();

  delete [] writeBuffer;
}

#endif


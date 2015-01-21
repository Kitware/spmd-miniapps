#ifndef __pgm_h
#define __pgm_h

#include "Image2D.h"

#include <fstream>
#include <iostream>
#include <string>

template <typename PType>
void load_pgm(const char *fname, Image2D<PType> *img)
{
  std::ifstream pgm;
  pgm.open(fname);

  std::string magic;
  int xdim, ydim, maxColor;

  pgm >> magic >> xdim >> ydim >> maxColor >> std::ws;

  int npixels = xdim * ydim;
  char *readBuffer = new char[npixels];
  pgm.read(readBuffer, npixels);
  pgm.close();

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
void write_pgm(const char *fname, const Image2D<PType> &img)
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

  std::ofstream pgm;
  pgm.open(fname);

  pgm << "P5\n" << img.getXDim() << " " << img.getYDim() << "\n255\n";
  pgm.write(writeBuffer, npixels);
  pgm.close();

  delete [] writeBuffer;
}

#endif


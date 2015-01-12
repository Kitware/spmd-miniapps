#ifndef __Image2D_h
#define __Image2D_h

template <typename PType>
class Image2D
{
public:
  Image2D()
    : buffer(0), xdim(0), ydim(0), pitch(0)
  {}

  Image2D(int xdim, int ydim)
    : xdim(xdim), ydim(ydim), pitch(xdim)
  {
    this->buffer = new PType[xdim * ydim];
  }

  ~Image2D()
  {
    delete [] this->buffer;
  }

  void allocate(int xdim, int ydim)
  {
    delete [] this->buffer;
    this->xdim = xdim;
    this->ydim = ydim;
    this->pitch = xdim;
    this->buffer = new  PType[xdim * ydim];
  }

  PType* getBuffer()
  {
    return this->buffer;
  }

  int getXDim() const
  {
    return this->xdim;
  }

  int getYDim() const
  {
    return this->ydim;
  }

  int getPitch() const
  {
    return this->pitch;
  }

  const PType* getBuffer() const
  {
    return this->buffer;
  }

private:
  PType *buffer;
  int xdim, ydim, pitch;
};

#endif


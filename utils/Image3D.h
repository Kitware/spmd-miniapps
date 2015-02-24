#ifndef __Image3D_h
#define __Image3D_h

template <typename T>
class Image3D
{
public:
  Image3D();
  ~Image3D();

  void setDimension(unsigned xdim, unsigned ydim, unsigned zdim);
  void setSpacing(T xspc, T yspc, T zspc);
  void setOrigin(T x, T y, T z);

  void allocate();
  T *getData();

  const unsigned* getDimension() const;
  unsigned getNumberOfPoints() const;
  const T* getSpacing() const;
  const T* getOrigin() const;
  const T* getData() const;

private:
  unsigned dim[3];
  unsigned npoints;
  T spacing[3];
  T origin[3];
  T *data;
};

template <typename T>
inline Image3D<T>::Image3D()
    : npoints(0), data(0)
{
  dim[0] = dim[1] = dim[2] = 0;
  spacing[0] = spacing[1] = spacing[2] = 0.0;
  origin[0] = origin[1] = origin[2] = 0.0;
}

template <typename T>
inline Image3D<T>::~Image3D()
{
  delete [] this->data;
}

template <typename T>
inline void Image3D<T>::setDimension(unsigned xdim, unsigned ydim,
                                     unsigned zdim)
{
  this->dim[0] = xdim;
  this->dim[1] = ydim;
  this->dim[2] = zdim;
}

template <typename T>
inline void Image3D<T>::setSpacing(T xspc, T yspc, T zspc)
{
  this->spacing[0] = xspc;
  this->spacing[1] = yspc;
  this->spacing[2] = zspc;
}

template <typename T>
inline void Image3D<T>::setOrigin(T x, T y, T z)
{
  this->origin[0] = x;
  this->origin[1] = y;
  this->origin[2] = z;
}

template <typename T>
inline void Image3D<T>::allocate()
{
  this->npoints = dim[0] * dim[1] * dim[2];
  delete [] this->data;
  this->data = new T[npoints];
}

template <typename T>
inline T* Image3D<T>::getData()
{
  return this->data;
}

template <typename T>
inline const unsigned* Image3D<T>::getDimension() const
{
  return this->dim;
}

template <typename T>
inline unsigned Image3D<T>::getNumberOfPoints() const
{
  return this->npoints;
}

template <typename T>
inline const T* Image3D<T>::getSpacing() const
{
  return this->spacing;
}

template <typename T>
inline const T* Image3D<T>::getOrigin() const
{
  return this->origin;
}

template <typename T>
inline const T* Image3D<T>::getData() const
{
  return this->data;
}

#endif


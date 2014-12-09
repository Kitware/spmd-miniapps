#ifndef __Image3D_h
#define __Image3D_h

class Image3D
{
public:
  Image3D();
  ~Image3D();

  void setDimension(int xdim, int ydim, int zdim);
  void setSpacing(double xspc, double yspc, double zspc);
  void setOrigin(double x, double y, double z);

  void allocate();
  double *getData();

  const int* getDimension() const;
  int getNumberOfPoints() const;
  const double* getSpacing() const;
  const double* getOrigin() const;
  const double* getData() const;

private:
  int dim[3];
  int npoints;
  double spacing[3];
  double origin[3];
  double *data;
};

inline Image3D::Image3D()
    : npoints(0), data(0)
{
  dim[0] = dim[1] = dim[2] = 0;
  spacing[0] = spacing[1] = spacing[2] = 0.0;
  origin[0] = origin[1] = origin[2] = 0.0;
}

inline Image3D::~Image3D()
{
  delete [] this->data;
}

inline void Image3D::setDimension(int xdim, int ydim, int zdim)
{
  this->dim[0] = xdim;
  this->dim[1] = ydim;
  this->dim[2] = zdim;
}

inline void Image3D::setSpacing(double xspc, double yspc, double zspc)
{
  this->spacing[0] = xspc;
  this->spacing[1] = yspc;
  this->spacing[2] = zspc;
}

inline void Image3D::setOrigin(double x, double y, double z)
{
  this->origin[0] = x;
  this->origin[1] = y;
  this->origin[2] = z;
}

inline void Image3D::allocate()
{
  this->npoints = dim[0] * dim[1] * dim[2];
  delete [] this->data;
  this->data = new double[npoints];
}

inline double* Image3D::getData()
{
  return this->data;
}

inline const int* Image3D::getDimension() const
{
  return this->dim;
}

inline int Image3D::getNumberOfPoints() const
{
  return this->npoints;
}

inline const double* Image3D::getSpacing() const
{
  return this->spacing;
}

inline const double* Image3D::getOrigin() const
{
  return this->origin;
}

inline const double* Image3D::getData() const
{
  return this->data;
}

#endif


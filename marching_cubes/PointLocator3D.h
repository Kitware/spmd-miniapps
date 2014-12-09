#ifndef __PointLocator3D_h
#define __PointLocator3D_h

#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

template<typename PointType>
class GetPointPosition
{
public:
  const double* operator()(const PointType &p) const
  {
    return p.getPosition();
  }
};

template <typename PointType,
          typename GetPointPositionType = GetPointPosition<PointType> >
class PointLocator3D
{
public:
  PointLocator3D(double range[6], int xdivs, int ydivs, int zdivs);

  PointType* insert(const PointType &point, bool *exists,
                    GetPointPositionType getPosition = GetPointPositionType());

  const int numberOfPoints() const;
  void dumpStats() const;

  template<typename Iterator>
  void getPoints(Iterator destination) const;

private:
  int xdivs, ydivs, zdivs;
  int xydivs, npoints;
  double range[6];

  typedef std::vector<PointType> Bin;
  std::vector<Bin> bins;
};

template <typename PointType, typename GetPointPositionType>
inline PointLocator3D<PointType, GetPointPositionType>::PointLocator3D(
  double range[6], int xdivs, int ydivs, int zdivs)
  : xdivs(xdivs), ydivs(ydivs), zdivs(zdivs), xydivs(xdivs * ydivs),
    npoints(0)
{
  std::copy(range, range + 6, this->range);

  int nbins = xydivs * zdivs;
  this->bins.resize(nbins);
}

template <typename PointType, typename GetPointPositionType>
PointType* PointLocator3D<PointType, GetPointPositionType>::insert(
  const PointType &point, bool *exists, GetPointPositionType getPosition)
{
  const double *pval = getPosition(point);

  assert(pval[0] >= this->range[0] && pval[0] <= this->range[1] &&
         pval[1] >= this->range[2] && pval[1] <= this->range[3] &&
         pval[2] >= this->range[4] && pval[2] <= this->range[5]);

  int binIdx[3];
  binIdx[0] = ((pval[0] - this->range[0])/(this->range[1] - this->range[0])) *
              static_cast<double>(this->xdivs);
  if (binIdx[0] == this->xdivs)
    {
    --binIdx[0];
    }

  binIdx[1] = ((pval[1] - this->range[2])/(this->range[3] - this->range[2])) *
              static_cast<double>(this->ydivs);
  if (binIdx[1] == this->ydivs)
    {
    --binIdx[1];
    }

  binIdx[2] = ((pval[2] - this->range[4])/(this->range[5] - this->range[4])) *
              static_cast<double>(this->zdivs);
  if (binIdx[2] == this->zdivs)
    {
    --binIdx[2];
    }

  int idx = binIdx[0] + (binIdx[1] * this->xdivs) + (binIdx[2] * this->xydivs);
  Bin &thisBin = this->bins[idx];

  PointType *ret = 0;
  for (typename Bin::iterator p = thisBin.begin(); p != thisBin.end(); ++p)
    {
    const double *pv = getPosition(*p);
    if (pv[0] == pval[0] && pv[1] == pval[1] && pv[2] == pval[2])
      {
      ret = &(*p);
      *exists = true;
      break;
      }
    }

  if (!ret)
    {
    thisBin.push_back(point);
    ++npoints;

    ret = &thisBin.back();
    *exists = false;
    }

  return ret;
}

template <typename PointType, typename GetPointPositionType>
inline const int PointLocator3D<PointType, GetPointPositionType>::
  numberOfPoints() const
{
  return this->npoints;
}

template <typename PointType, typename GetPointPositionType>
template <typename Iterator>
void PointLocator3D<PointType, GetPointPositionType>::getPoints(
  Iterator destination) const
{
  for (typename std::vector<Bin>::const_iterator bi = this->bins.begin();
       bi != this->bins.end(); ++bi)
    {
      for (typename Bin::const_iterator pi = bi->begin();
           pi != bi->end(); ++pi)
      {
        *destination = *pi;
        ++destination;
      }
    }
}

const int NBUCKETS = 50;

template <typename PointType, typename GetPointPositionType>
void PointLocator3D<PointType, GetPointPositionType>::dumpStats() const
{
  int maxPoints = 0, nempty = 0;
  for (int i = 0; i <  bins.size(); ++i)
    {
      maxPoints = std::max(maxPoints, static_cast<int>(bins[i].size()));
    }

  int bucketSize = ((maxPoints + NBUCKETS - 1)/NBUCKETS);
  int buckets[NBUCKETS];

  std::fill(buckets, buckets + NBUCKETS, 0);

  for (int i = 0; i < bins.size(); ++i)
    {
    int binsize = bins[i].size();
    if (!binsize)
      {
      ++nempty;
      }
    else
      {
      ++buckets[(binsize - 1)/bucketSize];
      }
    }

  std::cout << "Maximum bin size: " << maxPoints << "\n";
  std::cout << "Number of empty bins: " << nempty << "\n";
  for (int i = 0, j = 1; i < NBUCKETS; ++i, j += bucketSize)
    {
    std::cout << "(" << j << "-" << j + bucketSize - 1 << "): "
              << buckets[i] << "\n";
    }
}

#endif


#ifndef __PointLocator3D_h
#define __PointLocator3D_h

#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

namespace pl_internal
{
template<typename T, typename PointHandle>
class GetPointPosition
{
public:
  const T* operator()(const PointHandle &p) const
  {
    return p.getPosition();
  }
};
}

template <typename T, typename PointHandle, typename GetPointPositionType =
  pl_internal::GetPointPosition<T, PointHandle> >
class PointLocator3D
{
public:
  PointLocator3D(T range[6], unsigned xdivs, unsigned ydivs, unsigned zdivs);

  PointHandle* insert(const PointHandle &point, bool *exists,
    GetPointPositionType getPosition = GetPointPositionType());

  unsigned numberOfPoints() const;
  void dumpStats() const;

private:
  unsigned xdivs, ydivs, zdivs;
  unsigned xydivs, npoints;
  T range[6];

  typedef std::vector<PointHandle> Bin;
  std::vector<Bin> bins;
};

template <typename T, typename PointHandle, typename GetPointPositionType>
inline PointLocator3D<T, PointHandle, GetPointPositionType>::
PointLocator3D(T range[6], unsigned xdivs, unsigned ydivs, unsigned zdivs)
  : npoints(0)
{
  this->xdivs = std::max(1u, xdivs);
  this->ydivs = std::max(1u, ydivs);
  this->zdivs = std::max(1u, zdivs);
  this->xydivs = this->xdivs * this->ydivs;

  std::copy(range, range + 6, this->range);

  unsigned nbins = this->xydivs * this->zdivs;
  this->bins.resize(nbins);
}

template <typename T>
inline T clampValue(const T &minVal, const T &maxVal, const T& val)
{
  return std::min(maxVal, std::max(minVal, val));
}

template <typename T, typename PointHandle, typename GetPointPositionType>
PointHandle* PointLocator3D<T, PointHandle, GetPointPositionType>::
insert(const PointHandle &point, bool *exists,
       GetPointPositionType getPosition)
{
  const T e = 1e-9;
  const T *pval = getPosition(point);

  assert(pval[0] >= this->range[0] - e && pval[0] <= this->range[1] + e &&
         pval[1] >= this->range[2] - e && pval[1] <= this->range[3] + e &&
         pval[2] >= this->range[4] - e && pval[2] <= this->range[5] + e);

  unsigned binIdx[3];
  binIdx[0] = ((pval[0] - this->range[0])/(this->range[1] - this->range[0])) *
              static_cast<T>(this->xdivs);
  binIdx[0] = clampValue(0u, xdivs - 1, binIdx[0]);

  binIdx[1] = ((pval[1] - this->range[2])/(this->range[3] - this->range[2])) *
              static_cast<T>(this->ydivs);
  binIdx[1] = clampValue(0u, ydivs - 1, binIdx[1]);

  binIdx[2] = ((pval[2] - this->range[4])/(this->range[5] - this->range[4])) *
              static_cast<T>(this->zdivs);
  binIdx[2] = clampValue(0u, zdivs - 1, binIdx[2]);

  unsigned idx = binIdx[0] + (binIdx[1] * this->xdivs) +
                 (binIdx[2] * this->xydivs);
  Bin &thisBin = this->bins[idx];

  PointHandle *ret = 0;
  for (typename Bin::iterator p = thisBin.begin(); p != thisBin.end(); ++p)
    {
    const T *pv = getPosition(*p);
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

template <typename T, typename PointHandle, typename GetPointPositionType>
inline unsigned PointLocator3D<T, PointHandle, GetPointPositionType>::
numberOfPoints() const
{
  return this->npoints;
}

const unsigned NBUCKETS = 50;

template <typename T, typename PointHandle, typename GetPointPositionType>
void PointLocator3D<T, PointHandle, GetPointPositionType>::dumpStats() const
{
  unsigned maxPoints = 0, nempty = 0;
  for (unsigned i = 0; i <  bins.size(); ++i)
    {
      maxPoints = std::max(maxPoints, static_cast<unsigned>(bins[i].size()));
    }

  unsigned bucketSize = ((maxPoints + NBUCKETS - 1)/NBUCKETS);
  unsigned buckets[NBUCKETS];

  std::fill(buckets, buckets + NBUCKETS, 0);

  for (unsigned i = 0; i < bins.size(); ++i)
    {
    unsigned binsize = bins[i].size();
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
  for (unsigned i = 0, j = 1; i < NBUCKETS; ++i, j += bucketSize)
    {
    std::cout << "(" << j << "-" << j + bucketSize - 1 << "): "
              << buckets[i] << "\n";
    }
}

#endif


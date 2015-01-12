#include "poisson.h"
#include "poisson.ispc.h"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace scalar
{

void computeDivergence(const Float_t *in, int xdim, int ydim, Float_t *out)
{
  for (int j = 0; j < ydim; ++j)
    {
    for (int i = 0; i < xdim; ++i)
      {
      int idx = (j * xdim) + i;

      Float_t count = 0;
      Float_t div = 0;

      if (i > 0)
        {
        count += 1.0;
        div += in[idx - 1];
        }
      if (i < (xdim - 1))
        {
        count += 1.0;
        div += in[idx + 1];
        }
      if (j > 0)
        {
        count += 1.0;
        div += in[idx - xdim];
        }
      if (j < (ydim - 1))
        {
        count += 1.0;
        div += in[idx + xdim];
        }

      out[idx] = div - (count * in[idx]);
      }
    }
}

inline Float_t clamp(Float_t in, Float_t minVal, Float_t maxVal)
{
  return std::max(minVal, std::min(maxVal, in));
}

bool converged(const Float_t *curr, const Float_t *prev, int xdim, int ydim)
{
  const double esqr = 1.0e-6;

  double norm = 0.0;
  for (int j = 0; j < ydim; ++j)
    {
    for (int i = 0; i < xdim; ++i)
      {
        int idx = (j * xdim) + i;
        double diff = (curr[idx] - prev[idx]);
        norm += diff * diff;
      }
    }

  norm /= static_cast<double>(xdim * ydim);
  if (norm <= esqr)
    {
    return true;
    }

  return false;
}

void solvePoisson(const Float_t *in, int xdim, int ydim, Float_t *out)
{
  int npixels = xdim * ydim;
  std::fill(out, out + npixels, 127.0);

  Float_t *med = new Float_t[npixels];

  Float_t w =
    2.0/(1.0 + sin(M_PI/(static_cast<Float_t>(std::min(xdim, ydim)) + 1.0)));

  for (int k = 0; k < 1000; ++k)
    {
    bool checkForConvergence = ((k % 100) == 99);
    if (checkForConvergence)
      {
        std::copy(out, out + npixels, med);
      }

    for (int j = 0; j < ydim; ++j)
      {
      for (int i = 0; i < xdim; ++i)
        {
        int idx = (j * xdim) + i;

        Float_t count = 0;
        Float_t col = 0;

        if (i > 0)
          {
          count += 1;
          col += out[idx - 1];
          }
        if (i < (xdim - 1))
          {
          count += 1;
          col += out[idx + 1];
          }
        if (j > 0)
          {
          count += 1;
          col += out[idx - xdim];
          }
        if (j < (ydim - 1))
          {
          count += 1;
          col += out[idx + xdim];
          }

        //out[idx] = ((1.0 - w) * out[idx]) + (w * (col - in[idx])/count);
        out[idx] = (col - in[idx])/count;
        }
      }

    if (checkForConvergence && converged(out, med, xdim, ydim))
      {
      std::cout << "converged in " << k << " iterations\n";
      break;
      }
    }

  for (int j = 0; j < ydim; ++j)
    {
    for (int i = 0; i < xdim; ++i)
      {
        int idx = (j * xdim) + i;
        out[idx] = clamp(out[idx], 0.0, 255.0);
      }
    }
}

}; // namespace scalar

namespace simd
{

void solvePoisson(const Float_t *in, int xdim, int ydim, Float_t *out)
{
  Float_t *med = new Float_t[xdim * ydim];
  ispc::solvePoisson_impl(in, xdim, ydim, out, med);
  delete [] med;
}

};


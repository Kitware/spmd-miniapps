#include "type.h"

#if 0
inline void fill(uniform Float_t *from, uniform Float_t * uniform to,
                 uniform Float_t val)
{
  from += programIndex;
  while (from < to)
    {
      *from = val;
      from += programCount;
    }
}
#endif

inline void fill(uniform Float_t buf[], uniform int count, uniform Float_t val)
{
  foreach (i = 0 ... count)
    {
    buf[i] = val;
    }
}

inline void copy(uniform const Float_t src[], uniform int count,
                 uniform Float_t dst[])
{
  foreach (i = 0 ... count)
    {
    dst[i] = src[i];
    }
}

inline uniform bool converged(uniform const Float_t curr[],
                              uniform const Float_t prev[],
                              uniform int xdim, uniform int ydim)
{
  uniform const double esqr = 1.0e-6;

  double lnorm = 0.0;
  foreach (j = 0 ... ydim, i = 0 ... xdim)
    {
    int idx = (j * xdim) + i;
    double diff = (curr[idx] - prev[idx]);
    lnorm += diff * diff;
    }

  uniform double norm = reduce_add(lnorm)/((uniform double)(xdim * ydim));
  if (norm <= esqr)
    {
    return true;
    }

  return false;
}

export void solvePoisson_impl(uniform const Float_t div[], uniform int xdim,
                              uniform int ydim, uniform Float_t sol[],
                              uniform Float_t med[])
{
  uniform int npixels = xdim * ydim;
  fill(sol, npixels, 127.0);

  //print("count = %\n", programCount);
  //print("sizeof(float) = %\n", sizeof(float));
  //print("sizeof(uniform float) = %\n", sizeof(uniform float));

  for (uniform int k = 0; k < 1000; ++k)
    {
    uniform bool checkForConvergence = ((k % 100) == 99);
    if (checkForConvergence)
      {
      copy(sol, npixels, med);
      }

    foreach (j = 0 ... ydim, i = 0 ... xdim)
      {
      int idx = (j * xdim) + i;

      Float_t count = 0;
      Float_t col = 0;

      if (i > 0)
        {
        count += 1;
        col += sol[idx - 1];
        }
      if (i < (xdim - 1))
        {
        count += 1;
        col += sol[idx + 1];
        }
      if (j > 0)
        {
        count += 1;
        col += sol[idx - xdim];
        }
      if (j < (ydim - 1))
        {
        count += 1;
        col += sol[idx + xdim];
        }

      //sol[idx] = ((1.0 - w) * sol[idx]) + (w * (col - div[idx])/count);
      sol[idx] = (col - div[idx])/count;
      }

    if (checkForConvergence && converged(sol, med, xdim, ydim))
      {
      print ("converged in % iterations", k);
      break;
      }
    }

  foreach (j = 0 ... ydim, i = 0 ... xdim)
    {
    int idx = (j * xdim) + i;
    sol[idx] = clamp(sol[idx], (Float_t)0.0, (Float_t)255.0);
    }
}


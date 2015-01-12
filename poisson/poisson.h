#ifndef __poisson_h
#define __poisson_h

#include "type.h"

namespace scalar
{

void computeDivergence(const Float_t *in, int xdim, int ydim, Float_t *out);
void solvePoisson(const Float_t *div, int xdim, int ydim, Float_t *sol);

};

namespace simd
{
void solvePoisson(const Float_t *div, int xdim, int ydim, Float_t *sol);
};

#endif


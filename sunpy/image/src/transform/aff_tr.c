/*
Copyright (c) 2011 The SunPy developers
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.
Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/* Rotate, scale and shift in plain old unoptimised C */
/*
## well OK - 32 bit ints on this box


 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <assert.h>

#include "aff_tr.h"


// A couple of typedefs to avoid horrendous declarations
typedef double Kernfun(double, double);
typedef double Intfun(int *, INTYPE *, double, double, int, Kernfun *, double, char *, double);

// Some decl just to catch type errors
static Kernfun k_bicub;
static Kernfun k_bilin;
static Intfun interpol_nearest;
static Intfun interpol_kernel;

/* If compile with sinc() code, needs a couple of arrays to be bigger
   plus the extra code is added.  The extra array space is only used
   if needed (i.e. by sinc(), but does exist. */
#ifndef HAVESINC
#define PATCHLEN 4  // max size for kernel - 4 should be enough for bicub...
#else
static Kernfun k_sinc;
#define PATCHLEN 8  //  max size for kernel - sinc is quite big
#endif


/* Nearest "interp" */
static double interpol_nearest(int *dims, INTYPE *img, double row, double col, int dummy1, Kernfun dummy2, double dummy3, char *mode, double missing)
{
  int ri, ci;  // must be signed as calc can give negative

  ri = (int)round(row);
  ci = (int)round(col);
  if (ri < 0 || ci < 0 || ri >= dims[0] || ci >= dims[1])
    return missing;
  else
    return img[ri*dims[1]+ci];
}


/* Basic bicubic as convolution with separable kernel. 

While you can get exactly five points in a 1D convolution, then the
end ones are zero, so just do 4 - start arbitrarily at one end.

Also arb., we'll convolve along rows first.
 */

// 1D bicubic kernel (for 2D separable)
static double k_bicub(double x, double a)
{

  x = fabs(x);
  if (x <= 1.0)
    return (a+2.0)*x*x*x - (a+3.0)*x*x + 1.0;
  else if (x < 2.0)
    return a*x*x*x - 5.0*a*x*x + 8.0*a*x - 4.0*a;
  else
    return 0.0;
}

/* 1D bilinear kernel (for 2D separable).  Second arg only placeholder */
static double k_bilin(double x, double a)
{

  x = fabs(x);
  if (x <= 1.0)
    return 1 - x;
  else
    return 0;
}

#ifdef HAVESINC
/* A sinc fn which can be used for test but which is slow and takes up
   space.  As for others second arg is placeholder. */
static double k_sinc(double x, double a)
{

  if (fabs(x) <= 0.00001)
    return 1.0;
  else
    x = x*M_PI;
  return sin(x)/x;
}
#endif


/* For a point [row,col] between pixel values in an array img
calculate from a kernel an interpolated value.  That is, this uses
fractional indices on the input array.

k_size is the size actually used for the kernel, so some parts of
declared arrays might not in fact be used.

See affine_transform_kc() for more on args.

  */
static double interpol_kernel(int *dims, INTYPE *img, double row, double col, int k_size, Kernfun k_fun, double intparam, char *mode, double missing)
{
  int cols[PATCHLEN];
  double colw[PATCHLEN];
  int rows[PATCHLEN];
  double roww[PATCHLEN];
  double csum;
  double rsum;
  
  int c0;    // must be signed as calc can give negative
  int r0;
  int i,j;

  assert(k_size <= PATCHLEN);


// $$ could prob do N-dim
/* Tabulate the values to process for this point */
  c0 = (int)floor(col) - k_size/2 + 1;
  r0 = (int)floor(row) - k_size/2 + 1;
  for (j=0; j < k_size; j++)
  {
    cols[j] = c0 + j;  // $$ do we want to save this or recalc...
    colw[j] = k_fun((double)(c0 + j) - col, intparam); 
  }
  for (i=0; i < k_size; i++)
  {
    rows[i] = r0 + i;
    roww[i] = k_fun((double)(r0 + i) - row, intparam); 
  }
// convolve by cols - can be fn $$
/* Each step for separable reduces the dims by one, so NxN->N->scalar */
  rsum = 0.0;
  for (i=0; i<k_size; i++)
  {
    csum = 0.0;
    for (j=0; j<k_size; j++)
    {
      if (rows[i] < 0 || rows[i] >= dims[0] || cols[j] < 0 || cols[j] >= dims[1])
        csum = csum + colw[j]*missing;
      else
        csum = csum + colw[j]*img[rows[i]*dims[1] + cols[j]];
    }

    rsum = rsum + roww[i]*csum;
  }
  return rsum;
}


/* Run over the *output* image and get the input values.

Should behave like Python affine_transform, only requiring the top row
of the rotation matrix divided by mag.  Since by convention
r11,r12,r21,r22<->r[0],r[1],r[2],r[3] and we pass &r, the matrix
dimensions do not matter.

The int_type arg is an integer enum to NEAREST, BILINEAR or BICUBIC,
with int_param only used if need be.

## mode ignored for now
## not sure what speed implications of function pointers is
## likewise, implications of fixed index array refs

*/
int affine_transform_kc(int *dims, OUTTYPE *out_arr, INTYPE *in_arr, double *mat, double *tr, int int_type, double int_param, char *mode, double miss_val)
{

  int out1, out2;       // counters   
  double in1, in2;      // fractional indexes
  double o1, o2; 
  Intfun *i_fun;        // the interpolation function
  Kernfun *k_fun;       // the evt. kernel for the function
  int k_size;           // the evt. kernel size

  switch (int_type) {
  case NEAREST:
    i_fun = interpol_nearest;
    k_fun = NULL;
    k_size = 0;
    break;
  case BILINEAR:
    i_fun = interpol_kernel;
    k_fun = k_bilin;
    k_size = 2;
    break;
  case BICUBIC:
    i_fun = interpol_kernel;
    k_fun = k_bicub;
    k_size = 4;
    break;
#ifdef HAVESINC
  case SINC:
    i_fun = interpol_kernel;
    k_fun = k_sinc;
    k_size = 8;
    break;
#endif
  default:
    return -1;
  }

  for (out1=0; out1<dims[0]; out1++)    // rows
    for (out2=0; out2<dims[1]; out2++)  // cols
    {
      o1 = (double)out1;
      o2 = (double)out2;

      in1 = mat[0] * o1 + mat[1]  * o2 + tr[0];
      in2 = -mat[1] * o1 + mat[0]  * o2 + tr[1];
      out_arr[out1*dims[1]+out2] = i_fun(dims, in_arr, in1, in2, k_size, k_fun, int_param, mode, miss_val);
    }

return 0;
}




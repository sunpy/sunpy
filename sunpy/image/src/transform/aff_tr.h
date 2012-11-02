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

/* Header providing affine transform functions. */

#ifndef Z_AFFTR_H
#define Z_AFFTR_H 1

/* Image data types - to get these from include by caller, so can it
   cast if required types do not match. */
#define INTYPE float
#define OUTTYPE float

/* define if want sinc() interp option */
//#define HAVESINC

// Allowed values for int_param arg
#define NEAREST 0
#define BILINEAR 1
#define BICUBIC 2
#ifdef HAVESINC
#define SINC 3
#endif

#ifndef M_PI   // prob the case for c99
#define M_PI 3.14159265358979323846
#endif

// Quoter
#define QU(A) xQU(A)
#define xQU(A) #A


/* Rotate/scale/shift as seen from output image (i.e. like Scipy fn,
   though using kernel convolution).  */
int affine_transform_kc(int *dims, OUTTYPE *out_arr, INTYPE *in_arr, double *mat, double *tr, int int_type, double int_param, char *mode, double miss_val);


#endif

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

/*

%%% coord this with aff_tr header types!

 */


#ifndef Z_ROT_EXTN_H
#define Z_ROT_EXTN_H 1



/* Input and output images are Python objects, which know what tyoe
   they are.  You create a C array from the input array casting to
   whatever type you want, but at least in simplest form you must
   create an output array and tell it what your C routine is producing
   (which makes sense as then you don't lose anything). */ 
#define IN_FLOAT32    // if we want float32, default is double/float64
#define OUT_FLOAT32   // if we want float32, default is double/float64

/* Types for data - must be the same as in aff_tr.h */
/* Type of input C array from whatever Python type */
#ifdef IN_FLOAT32
#define INTYPE float            // for C decl
#define PYIN_TYPE NPY_FLOAT32   // for cast
#else
#define INTYPE double
#define PYIN_TYPE NPY_DOUBLE
#endif
/* Type of output C array and of output Python type */
#ifdef OUT_FLOAT32
#define OUTTYPE float            // for C decl
#define PYOUT_TYPE PyArray_FLOAT32  // for NPy decl
#else
#define OUTTYPE double
#define PYOUT_TYPE PyArray_FLOAT64
#endif



#endif

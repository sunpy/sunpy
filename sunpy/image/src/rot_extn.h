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

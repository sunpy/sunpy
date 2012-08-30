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

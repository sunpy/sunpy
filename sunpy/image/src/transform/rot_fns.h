/* Header providing rot_idl function. */

#ifndef Z_ROTFUNS_H
#define Z_ROTFUNS_H

#include "aff_tr.h"  // for INTYPE and OUTTYPE for callers

/* C fn that tries to behave like the IDL ROT fn */
int rot_idl(int *dims, OUTTYPE *dest, INTYPE *src, double angle, double scale, double *x0, double *y0, double cubic, int interp, double missing, int pivot);

/* Rotate/scale/shift as seen from input image. */
int rot(int *dims, OUTTYPE *out_arr, INTYPE *in_arr, double angle, double mag, double *centre, double *shift, int i_type, double i_param, char *mode, double miss_val);

#endif

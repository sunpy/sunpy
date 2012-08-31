/* Collection of C fns providing API's for the rotate/scale/shift.

rot_idl() behaves the same as the IDL ROT fn.
rot() is a more usual fn as seen from the input array 

*/

/*
$$ hmm - types from aff_tr.h - is this the best


 */

#include <stdlib.h>
//#include <stdio.h>  // iff printf etc.
#include <math.h>
#include <assert.h>

#include "aff_tr.h"
#include "rot_fns.h"


/* Calculate the values for rotate and transform vectors for affine
   transform fn. from angle, centre and so forth.

## only need 1st two of mat, so return m[0,0],m[0,1],k[0],k[1]

%% angle units - deg for idl

 */
static double *af_args(double angle, double *centre, double *shift, double scale)
{
  static double res[4];   // the 4 calc params
  double a[2];

  a[0] = centre[0] + shift[0];
  a[1] = centre[1] + shift[1];
  res[0] = cos(angle)/scale;  // m11
  res[1] = sin(angle)/scale;  // m12
  res[2] = centre[0] - res[0]*a[0] - res[1]*a[1];  // k1
  res[3] = centre[1] + res[1]*a[0] - res[0]*a[1];  // k2

  return res;
}

// %% types?
// %% defn dest in caller, so need arg...
/*

Array args and index contents always row,column (e.g. dims[0]==rows),
i.e. y,x.

scale : set to 1 for same size
x0, y0 : shift/pivot points - NULL for centre
cubic : <-1 not bicubic, >0 ->-1.0 (as per IDL)
interp : if not bicubic, says if nearest or bilinear
missing : the value for outside the image (as per IDL, modes not available)

 */
int rot_idl(int *dims, OUTTYPE *dest, INTYPE *src, double angle, double scale, double *x0, double *y0, double cubic, int interp, double missing, int pivot)
{
  int ktype;
  double kparm;
  double ang;
  double centre[2];
  double shift[2];
  double im_c[2];
  double *ap;


  if ((cubic < -1) && (interp == 0))
  {
    ktype = NEAREST;
    kparm = 0;       // dummy
  }
  else if (cubic < -1)
  {
    ktype = BILINEAR;
    kparm = 0;       // dummy
  }
  else
  {
    ktype = BICUBIC;
    kparm = (cubic > 0.0?-1:cubic);
  }
  ang = M_PI*angle/180.0;
  im_c[0] = (dims[0] - 1)/2.0;
  im_c[1] = (dims[1] - 1)/2.0;
  centre[0] = (y0 == NULL)?im_c[0]:*y0;
  centre[1] = (x0 == NULL)?im_c[1]:*x0;
  if (pivot == 1)
  {
    shift[0] = shift[1] = 0.0;
  }
  else
  {
    shift[0] = im_c[0] - centre[0];
    shift[1] = im_c[1] - centre[1];
  }
  ap = af_args(ang, centre, shift, scale);
  affine_transform_kc(dims, dest, src, &ap[0], &ap[2], ktype, kparm, "const", missing);

  return 0;
}


/* Rotate/scale/shift in commonly used terms - mainly referenced to
   input image.

## centre of rotation can be fractional and for even array size will
   often be - e.g. 128 -> 63.5

*/
int rot(int *dims, OUTTYPE *out_arr, INTYPE *in_arr, double angle, double mag, double *centre, double *shift, int int_type, double int_param, char *mode, double miss_val)
{
  double mat[2];
  double tr[2];     // transpose, row col as usual
  double s[2];

  assert(mag > 0.001);

/* Make matrix for reference backwards from output, but keep translate in coord\
   space */
  mat[0] = cos(angle)/mag;   // r11
  mat[1] = sin(angle)/mag;   // r12
  s[0] = centre[0] + shift[0];
  s[1] = centre[1] + shift[1];
  tr[0] = centre[0] - mat[0]*s[0] - mat[1]*s[1];
  tr[1] = centre[1] + mat[1]*s[0] - mat[0]*s[1];

  return affine_transform_kc(dims, out_arr, in_arr, mat, tr, int_type, int_param, "const", miss_val);
}

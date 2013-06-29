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

/* Numpy extension example, for call of a rotate function
This with args of array, number and keywords.

*/

/* TTD :

- data types - float is the same as float64 in Python, viz double

 */

// You always need these
#include <Python.h>
#include <numpy/arrayobject.h>

#include <transform/aff_tr.h>

#include <rot_extn.h>


static char docstring[] =
"sunpy.image.Crotate.affine_transformation(input, matrix, offset=[0.0, 0.0], kernel=Crotate.BICUBIC, cubic=-0.5, mode='constant', cval=0.0)\n\nApply an affine transformation to an image array\n\nParameters\n----------\ninput : ndarray\nmatrix : ndarray\noffset\nkernel\ncubic\nmode\ncval";

/* Function called on Python function call via method definition map
   below.  Takes as args two np arrays and optionally a tuple, an int,
   a string and a double.

   It seems you need all args to have keyword names, which may be
   unused. */
static PyObject *rot_shift_scale_args(PyObject *dummy, PyObject *args, PyObject *kwords) 
{
  static char *kwlist[] = {"image", "rsmat", "offset", "kernel", "cubic", "mode", "cval", NULL}; 
  PyObject *arg1 = NULL;   // image
  PyObject *arg2 = NULL;   // rot & scale
  PyObject *arr1 = NULL;   // image formatted
  PyObject *arr2 = NULL;   // rot & scale formatted

  // ## should these be const?
  double *rotscale;
  double offset[2] = {0.0, 0.0};
  int ktype = BICUBIC;        // default value 
  double interp_param = -0.5; // bicubic iterp param
  char *mode = "constant"; 
  double cval = 0.0;          // optional arg

  int dims[2];
  int ndim;
  int i;
  INTYPE *in_arr;

  int result;
  PyObject *out1;  
  OUTTYPE *out_arr;

/* Map via physical position between C variables and Python variables
   (in kwlist) - unused ones are simply left as previously defined. */
  if (!PyArg_ParseTupleAndKeywords(args, kwords, "OO|(dd)idsd", kwlist, &arg1, &arg2, &offset[0], &offset[1], &ktype, &interp_param, &mode, &cval))
    return NULL;   //NULL says error
  
/* Make a nice numpy array using macro for PyArray_FromAny laid out
   for C access as floats and doubles.  The alignment is a bit tricky
   - sometimes you *don't* need FORCECAST, but often you do. */
// %%  arr1 = PyArray_FROM_OTF(arg1, NPY_FLOAT32, NPY_CARRAY_RO | NPY_FORCECAST);  
  arr1 = PyArray_FROM_OTF(arg1, PYIN_TYPE, NPY_CARRAY_RO | NPY_FORCECAST);  
  arr2 = PyArray_FROM_OTF(arg2, NPY_DOUBLE, NPY_IN_ARRAY); 

  // ##printf("E: type check %d %d\n",PyArray_TYPE(arr1), NPY_FLOAT32);

  // ## printf("E: Constant is %f, offset is %f %f\n", cval, offset[0], offset[1]);
  in_arr = PyArray_DATA(arr1);    // that's where the data for C is   
  ndim = PyArray_NDIM(arr1);
  // ## printf("E: Input array has %d dimensions\n", ndim);

  for (i=0; i<ndim; i++)
  {
    //##  printf("E: Input dim %d is %d with stride %d\n",i,(int)PyArray_DIM(arr1,i),(int)PyArray_STRIDE(arr1,i));
    dims[i] = (int)PyArray_DIM(arr1,i);  // make sure alignment OK
  }
  //##printf("E: Sizeof input elements: %d\n",(int)sizeof(in_arr[0]));
  //##printf("E: Sizeof output elements: %d\n",(int)sizeof(out_arr[0]));

  rotscale = PyArray_DATA(arr2);    // set to the location of the operator matrix
  // ## check dims or size e.g. PyArray_SIZE(arr2))
  /*printf("E: Rotscale in order: "); 
  for (i=0; i<4; i++)
    printf(" %f", rotscale[i]);
  printf("\n");
  */

// free() ##?  - well, the whole story in doc...
// Make a new array the same size as the input array

//%%  out1 = PyArray_SimpleNew(PyArray_NDIM(arr1), PyArray_DIMS(arr1), PyArray_FLOAT32);
  out1 = PyArray_SimpleNew(PyArray_NDIM(arr1), PyArray_DIMS(arr1), PYOUT_TYPE);
  out_arr = PyArray_DATA(out1);           // this is the place for the data from C
   
//##  printf("E: ktype (interp type) %d\n", ktype);
/* Call to function that does the actual work.  This one is external. */ 
  result = affine_transform_kc(
        dims, out_arr, in_arr, 
        rotscale, 
        offset,
        ktype, interp_param, 
        mode, cval
        );

  Py_DECREF(arr1);
  Py_DECREF(arr2);

  if (result != 0)
    return NULL;

  return Py_BuildValue("N", out1);
}


/* This defines the method names and maps to C fns */
static PyMethodDef my_methods[] = {
  { "affine_transform", (PyCFunction)rot_shift_scale_args, METH_VARARGS | METH_KEYWORDS , docstring },
  {NULL, NULL, 0, NULL}        /* Sentinel */
};


PyMODINIT_FUNC initCrotate(void)   // name is "init"+<module name>                 
{
  PyObject *ge_mod;

  ge_mod = Py_InitModule("Crotate", my_methods);
  if (ge_mod == NULL)
    return;
  PyModule_AddIntConstant(ge_mod, "NEAREST", NEAREST);
  PyModule_AddIntConstant(ge_mod, "BILINEAR", BILINEAR);
  PyModule_AddIntConstant(ge_mod, "BICUBIC", BICUBIC);
#ifdef HAVESINC
  PyModule_AddIntConstant(ge_mod, "SINC", SINC);
#endif
  import_array();   // always need - despite the name, gets data for all Numpy
}

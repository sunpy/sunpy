/*!
@brief A C extension for Python to read ana f0 files
@author Tim van Werkhoven <t.i.m.vanwerkhoven@gmail.com>

Based on Michiel van Noort's IDL DLM library 'f0' which contains a cleaned up 
version of the original anarw routines.
*/

// Headers
#include <Python.h>				// For python extension
#include <numpy/arrayobject.h> 	// For numpy
#include <sys/time.h>			// For timestamps
#include <time.h>				// For timestamps
//#include "anadecompress.h"  
//#include "anacompress.h"
#include "types.h"
#include "anarw.h"

// vasprintf() and asprintf() may not be defined, particularly on Windows
#ifndef __USE_GNU
int vasprintf( char **sptr, const char *fmt, va_list argv )
{
    va_list argv2;
    va_copy(argv2, argv);

    int wanted = vsnprintf( *sptr = NULL, 0, fmt, argv );
    if( (wanted < 0) || ((*sptr = malloc( 1 + wanted )) == NULL) )
        return -1;

    int retval;
    retval = vsnprintf( *sptr, 1 + wanted, fmt, argv2 );
    va_end(argv2);
    return retval;
}

int asprintf( char **sptr, const char *fmt, ... )
{
    int retval;
    va_list argv;
    va_start( argv, fmt );
    retval = vasprintf( sptr, fmt, argv );
    va_end( argv );
    return retval;
}
#endif

// Prototypes
static PyObject * pyana_fzread(PyObject *self, PyObject *args);
static PyObject * pyana_fzwrite(PyObject *self, PyObject *args);

// Methods table for this module
static PyMethodDef PyanaMethods[] = {
    {"fzread",  pyana_fzread, METH_VARARGS, "Load an ANA F0 file."},
    {"fzwrite",  pyana_fzwrite, METH_VARARGS, "Save an ANA F0 file."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};


// Init module methods
PyMODINIT_FUNC init_pyana(void) {
    (void) Py_InitModule("_pyana", PyanaMethods);
	// Init numpy usage
	import_array();
}

/*!
@brief load an ANA f0 file data and header
@param [in] filename
@return [out] data, NULL on failure
*/
static PyObject *pyana_fzread(PyObject *self, PyObject *args) {
	// Function arguments
	char *filename;
	int debug=0;
	// Init ANA IO variables
	char *header = NULL;			// ANA header (comments)
	uint8_t *anaraw = NULL;			// Raw data
	int	nd=-1, type=-1, *ds, size=-1, d; // Various properties
	// Data manipulation
	PyArrayObject *anadata;			// Final ndarray
	
	// Parse arguments
	if (!PyArg_ParseTuple(args, "s|i", &filename, &debug)) {
		return NULL;
	}
	
	// Read ANA file
	if (debug == 1)
		printf("pyana_fzread(): Reading in ANA file\n");
	anaraw = ana_fzread(filename, &ds, &nd, &header, &type, &size);
		
	if (NULL == anaraw) {
		PyErr_SetString(PyExc_ValueError, "In pyana_fzread: could not read ana file, data returned is NULL.");
		return NULL;
	}
	if (type == -1) {
		PyErr_SetString(PyExc_ValueError, "In pyana_fzread: could not read ana file, type invalid.");
		return NULL;
	}
	
	// Mold into numpy array
	npy_intp npy_dims[nd];		// Dimensions array 
	int npy_type;					// Numpy datatype
	
	// Calculate total datasize
	if (debug == 1)
		printf("pyana_fzread(): Dimensions: ");
	for (d=0; d<nd; d++) {
		if (debug == 1)
			printf("%d ", ds[d]);
		// ANA stores dimensions the other way around?
		//npy_dims[d] = ds[d];
		npy_dims[nd-1-d] = ds[d];
	}
	if (debug == 1)
		printf("\npyana_fzread(): Datasize: %d\n", size);
	
	// Convert datatype from ANA type to PyArray type
	switch (type) {
		case (INT8): npy_type = PyArray_INT8; break;
		case (INT16): npy_type = PyArray_INT16; break;
		case (INT32): npy_type = PyArray_INT32; break;
		case (FLOAT32): npy_type = PyArray_FLOAT32; break;
		case (FLOAT64): npy_type = PyArray_FLOAT64; break;
		case (INT64): npy_type = PyArray_INT64; break;
		default: 
			PyErr_SetString(PyExc_ValueError, "In pyana_fzread: datatype of ana file unknown/unsupported.");
			return NULL;
	}
	if (debug == 1)
		printf("pyana_fzread(): Read %d bytes, %d dimensions\n", size, nd);
	
	// Create numpy array from the data 
	anadata = (PyArrayObject*) PyArray_SimpleNewFromData(nd, npy_dims,
		npy_type, (void *) anaraw);
	// Make sure Python owns the data, so it will free the data after use
	PyArray_FLAGS(anadata) |= NPY_OWNDATA;
	
	if (!PyArray_CHKFLAGS(anadata, NPY_OWNDATA)) {
		PyErr_SetString(PyExc_RuntimeError, "In pyana_fzread: unable to own the data, will cause memory leak. Aborting");
		return NULL;
	}
	
	// Return the data in a dict with some metainfo attached
	// NB: Use 'N' for PyArrayObject s, because when using 'O' it will create 
	// another reference count such that the memory will never be deallocated.
	// See:
	// http://www.mail-archive.com/numpy-discussion@scipy.org/msg13354.html 
	// ([Numpy-discussion] numpy CAPI questions)
	return Py_BuildValue("{s:N,s:{s:i,s:(ii),s:s}}", 
		"data", anadata,
		"header",
		"size", size,
		"dims", ds[0], ds[1],
		"header", header);
}


/*!
@brief save an ANA format image to disk
@param [in] filename Full path to write data to
@param [in] data Data to write (numpy array)
@param [in] compress Apply (Rice) compression or not
@param [in] header Add a header to the file (or use default)
@return number of bytes read on success, NULL pointer on failure
*/
static PyObject * pyana_fzwrite(PyObject *self, PyObject *args) {
	// Python function arguments
	char *filename = NULL;
	PyArrayObject *anadata;
	int compress = 1, debug=0;
	char *header = NULL;
	// Processed data goes here
	PyObject *anadata_align;
	uint8_t *anadata_bytes;
	// ANA file writing 
	int	type, d;
	
	// Parse arguments from Python function
    if (!PyArg_ParseTuple(args, "sO!|isi", &filename, &PyArray_Type, &anadata, &compress, &header, &debug))
        return NULL;
	
	// Check if filename was parsed correctly (should be, otherwise
	// PyArg_ParseTuple should have raised an error, obsolete?)
	if (NULL == filename) {
		PyErr_SetString(PyExc_ValueError, "In pyana_fzwrite: invalid filename.");
		return NULL;
	}
	// If header is NULL, then set the comment to a default value
	if (NULL == header) {
		if (debug == 1) printf("pyana_fzwrite(): Setting default header\n");
		struct timeval *tv_time=NULL;
		struct tm *tm_time=NULL;
		gettimeofday(tv_time, NULL);
		tm_time = gmtime(&(tv_time->tv_sec));
		asprintf(&header, "#%-42s compress=%d date=%02d:%02d:%02d.%03ld\n",
			filename, 
			compress, 
			tm_time->tm_hour, tm_time->tm_min, tm_time->tm_sec, (long) (tv_time->tv_usec/1000));
	}
	if (debug == 1) printf("pyana_fzwrite(): Header: '%s'\n", header);
	
	// Convert datatype from PyArray type to ANA type, and verify that ANA
	// supports it
	switch (PyArray_TYPE((PyObject *) anadata)) {
		case (PyArray_INT8): 
			type = INT8; 
			if (debug == 1) 
				printf("pyana_fzwrite(): Found type PyArray_INT8\n");
			break;
		case (PyArray_INT16): 
			type = INT16; 
			if (debug == 1) 
				printf("pyana_fzwrite(): Found type PyArray_INT16\n");
			break;
		case (PyArray_FLOAT32): 
			type = FLOAT32; 
			if (debug == 1) 
				printf("pyana_fzwrite(): Found type PyArray_FLOAT32\n");
			break;
		case (PyArray_FLOAT64): 
			type = FLOAT64;
			if (debug == 1)
				printf("pyana_fzwrite(): Found type PyArray_FLOAT64\n");
			break;
		//case (PyArray_INT64): type = INT64; break;
		default:
			PyErr_SetString(PyExc_ValueError, "In pyana_fzwrite: datatype cannot be stored as ANA file.");
			return NULL;
			break;
	}
	// Check if compression flag is sane
	if (compress == 1 && (type == FLOAT32 || type == FLOAT64)) {
		PyErr_SetString(PyExc_RuntimeError, "In pyana_fzwrite: datatype requested cannot be compressed.");
		return NULL;
	}
	if (debug == 1) 
		printf("pyana_fzwrite(): pyarray datatype is %d, ana datatype is %d\n", 
		PyArray_TYPE((PyObject *) anadata), type);


	// Sanitize data, make a new array from the old array and force the
	// NPY_CARRAY_RO requirement which ensures a C-contiguous and aligned
	// array will be made
	anadata_align = PyArray_FromArray(anadata, PyArray_DESCR((PyObject *) anadata), NPY_CARRAY_RO);
	
	// Get a pointer to the aligned data
	anadata_bytes = (uint8_t *) PyArray_BYTES(anadata_align);
	// Get the number of dimensions
	PyArrayObject *arrobj = (PyArrayObject*) anadata_align;
	int nd = arrobj->nd;
	int *dims = malloc(nd*sizeof(int));
	// Get the dimensions and number of elements
	npy_intp *npy_dims = PyArray_DIMS(anadata_align);
	//npy_intp npy_nelem = PyArray_SIZE(anadata_align);
	
	if (debug == 1) printf("pyana_fzwrite(): Dimensions: ");
	for (d=0; d<nd; d++) {
		// ANA stores dimensions the other way around?
		//dims[d] = npy_dims[d];
		dims[d] = npy_dims[nd-1-d];
		if (debug == 1) printf(" %d", dims[d]);
	}
	if (debug == 1) printf("\npyana_fzwrite(): Total is %d-dimensional\n", nd);

	// Write ANA file
	if (debug == 1) printf("pyana_fzwrite(): Compress: %d\n", compress);
	if (compress == 1)
		ana_fcwrite(anadata_bytes, filename, dims, nd, header, type, 5);	
	else
		ana_fzwrite(anadata_bytes, filename, dims, nd, header, type);
	
	free(dims);
	// If we didn't crash up to here, we're probably ok :P
	return Py_BuildValue("i", 1);
}

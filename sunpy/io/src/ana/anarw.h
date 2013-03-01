#ifndef __ANARW_H__      // __ANACOMPRESS_H__
#define __ANARW_H__      // __ANACOMPRESS_H__

#define ANA_VAR_SZ {1,2,4,4,8,8}

#define ANA_LITTLE_ENDIAN    0
#define ANA_BIG_ENDIAN       1

#define INT8	         0
#define INT16	         1
#define INT32	         2
#define FLOAT32	         3
#define FLOAT64	         4
#define INT64	         5

#define ANA2PYTHON_T {PyArray_INT8, PyArray_INT16, PyArray_INT32, PyArray_FLOAT32, PyArray_FLOAT64, PyArray_INT64}

#define M_TM_INPRO       0  
#define M_TM_INFUN      -1


// Helper routines
void bswapi64(int64_t *x, int n);
void bswapi32(int32_t *x, int n);
void bswapi16(int16_t *x, int n);
int ck_synch_hd(FILE *fin, fzhead_t *fh, int t_endian);

// Ana I/O routines
char *ana_fzhead(char *file_name); // fzhead subroutine	
uint8_t *ana_fzread(char *file_name, int **ds, int *nd, char **header, int *type, int *osz); // fzread subroutine
void ana_fzwrite(uint8_t *data, char *file_name, int *ds, int nd, char *header, int py_type);	/* fcwrite subroutine */
void ana_fcwrite(uint8_t *data, char *file_name, int *ds, int nd, char *header, int py_type, int slice);	/* fcwrite subroutine */

#endif				// __ANACOMPRESS_H__

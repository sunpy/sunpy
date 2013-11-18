#ifndef __TYPES_H__  // __TYPES_H__
#define __TYPES_H__

#include <stdint.h>

//! @todo fix include guards to legal versions

typedef float float32_t;
typedef double float64_t;

typedef float64_t fp_t;  //!< The default floating point type

struct complex{
  fp_t re;
  fp_t im;
};
typedef struct complex complex_t;

struct fzhead {                    // first block for fz files
  int synch_pattern;
  uint8_t subf;
  uint8_t source;
  uint8_t nhb,datyp,ndim,file_class;
  uint8_t cbytes[4];	      // can't do as int because of %4 rule
  uint8_t free[178];
  int dim[16];
  char txt[256];
};
typedef struct fzhead fzhead_t;

struct compresshead{
  int tsize,nblocks,bsize;
  uint8_t slice_size,type;
};

void bswapi16(int16_t *x,int n);
void bswapi32(int32_t *x,int n);
void bswapi64(int64_t *x,int n);

#endif               // __TYPES_H__

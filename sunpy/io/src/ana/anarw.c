#include <errno.h>  
#include <stdlib.h>  
#include <string.h>  
#include <sys/stat.h>  
#include <stdio.h>  
#include <stdarg.h>  
#include <stdint.h>
#include "types.h"  

#include "anadecompress.h"  
#include "anacompress.h"

#include "anarw.h"


static __inline int min(int a,int b)
{
  return (a<b)?a:b;
}

static __inline void swap(char *a,char *b)
{
  char c=*a;
  *a=*b;
  *b=c;
}

void bswapi64(int64_t *x,int n)
{
  int sz=sizeof(int64_t),i,j;
  int m=sz/2;
  sz-=1;
  for(i=0;i<=n-1;++i){
    char *p=(char*)(x+i);
    for(j=0;j<=m-1;++j) swap(&p[j],&p[sz-j]);
  }
}

void bswapi32(int32_t *x,int n)
{
  int sz=sizeof(int32_t),i,j;
  int m=sz/2;
  sz-=1;
  for(i=0;i<=n-1;++i){
    char *p=(char*)(x+i);
    for(j=0;j<=m-1;++j) swap(&p[j],&p[sz-j]);
  }
}

void bswapi16(int16_t *x,int n)
{
  int sz=sizeof(int16_t),i,j;
  int m=sz/2;
  sz-=1;
  for(i=0;i<=n-1;++i){
    char *p=(char*)(x+i);
    for(j=0;j<=m-1;++j) swap(&p[j],&p[sz-j]);
  }
}

int ck_synch_hd(FILE *fin,fzhead_t *fh,int t_endian)
{
  int wwflag=0;
  if(fread(fh,1,sizeof(fzhead_t),fin)!=sizeof(fzhead_t)){
    fprintf(stderr,"error in fzread while reading header\n");
    return -1;
  }
  int syncpat=(fh->synch_pattern==0x5555aaaa);
  int revsyncpat=(fh->synch_pattern==0xaaaa5555);
  if(!(syncpat||revsyncpat)){
    fclose(fin);
    fprintf(stderr,"ck_synch_hd: error: file does not have the F0 synch pattern (found 0x%x instead)\n",fh->synch_pattern);
    return -1;
  }
  if(syncpat==t_endian){
    fprintf(stderr,"ck_synch_hd: warning: reversed F0 synch pattern\n");
    wwflag=1; 
  }
  if(fh->nhb>1){                     // if the header is long, read in the rest now
    if(fh->nhb>15){
      fclose(fin);
      fprintf(stderr,"ck_synch_hd: error: annot handle header more than 16 blocks!\n");
      return -1;
    }
    int size=(fh->nhb-1)*sizeof(fzhead_t);
    uint8_t *buf=malloc(size);
    fread(buf,1,size,fin);
    free(buf); // not very useful?
  }
  if(t_endian) bswapi32(fh->dim,fh->ndim); // for big endian machines
  return wwflag; 
}

char *ana_fzhead(char *file_name) // fzhead subroutine	
{
  struct stat stat_buf;
  if(stat(file_name,&stat_buf)<0){
    fprintf(stderr,"ana_fzhead: error: file \"%s\" not found.\n",file_name); 
    return 0;
  }

  int one=1;
  int t_endian=(*(char*)&one==0);                      // an endian detector, taken from SL's tiff library 
//
  FILE *fin=fopen(file_name,"r");
  if(!fin){
    fprintf(stderr,"ana_fzhead: error: could not open file \"%s\": %s!\n",file_name,strerror(errno));
    return 0;
  }
  fzhead_t fh;
  int sef;
  if((sef=ck_synch_hd(fin,&fh,t_endian))<0) return 0;

  char *header=strcpy(malloc(strlen(fh.txt)+1),fh.txt);
  fclose(fin);
  return header; 
}


uint8_t *ana_fzread(char *file_name,int **ds,int *nd,char **header,int *type,int *osz) // fzread subroutine	
{

  struct stat stat_buf;
  if(stat(file_name,&stat_buf)<0){
    fprintf(stderr,"ana_fzread: error: file \"%s\" not found.\n",file_name); 
    return 0;
  }

  int type_sizes[]=ANA_VAR_SZ;
  int one=1;
  int t_endian=(*(char*)&one==0);                      // an endian detector, taken from SL's tiff library 
//
  FILE *fin=fopen(file_name,"r");
  if(!fin){
    fprintf(stderr,"ana_fzread: error: could not open file \"%s\": %s!\n",file_name,strerror(errno));
    return 0;
  }

  fzhead_t fh;
  int sef;
  if((sef=ck_synch_hd(fin,&fh,t_endian))<0) {
    fprintf(stderr,"ana_fzread: error: ck_sync_hd error!\n");
    return 0;
  }

  *header=strcpy(malloc(strlen(fh.txt)+1),fh.txt);

  *ds=((int*)malloc(sizeof(int)*(*nd=fh.ndim)));
  int d;
  for(d=0;d<*nd;++d) (*ds)[d]=fh.dim[d];
  int n_elem=1;
  for(d=0;d<=fh.ndim-1;++d) n_elem*=fh.dim[d]; // compute size of array
  *type=fh.datyp;
  int f_endian=(fh.subf>=128);                     // the top bit of the byte fh->subf denotes endian type, 0 for little and 1 for big
  int swap_endian=(f_endian!=t_endian);            // file has different endianness
  if(sef) swap_endian=(!swap_endian);              // file contains strange data
  int compressed=(fh.subf&1);                      // file is in compressed format
  if(compressed){                                  // compressed format
    struct compresshead ch;
    if(fread(&ch,1,14,fin)<14) fprintf(stderr,"error reading in compression header\n");
// header by default little-endian?
    if(t_endian){ // big endian platform
      bswapi32(&ch.tsize,1);
      bswapi32(&ch.nblocks,1);
      bswapi32(&ch.bsize,1);
    }
// read data
    int size=ch.tsize-14;
	// Tim van Werkhoven, 20090327: 
	// Add 4 to the malloc because anadecrunch() might read some data 
	// beyond the current pixel it's investigating. If the pixel is the last 
	// pixel, it will try to read beyond the malloc'ed area, giving trouble. 
	// This should help prevent that. Problem arises around y.b[1]=x[i+1]; in 
	// anadecrunch(), see Valgrind for more info. Maximum read-ahead is 4 
	// bytes (set by nb), so 4 bytes extra in the malloc should be sufficient.
    uint8_t *buf=malloc(size+4);
    if(fread(buf,1,size,fin)<size) fprintf(stderr,"error reading in compressed data\n");
    fclose(fin);
//
    if(ch.bsize*ch.nblocks>n_elem){ // fix a possible problem with ch.nblocks
      fprintf(stderr,"warning, bad ch.nblocks = %d\ncorrecting to %d, hope this is right!\n",ch.nblocks,n_elem/ch.bsize);
      ch.nblocks=n_elem/ch.bsize;
    }
    if(ch.type%2==*type) fprintf(stderr,"inconsistent compression type\n"); // consistency check
    int rv;
    uint8_t *out=malloc(n_elem*type_sizes[*type]);
    switch(ch.type){
      case(0): rv=anadecrunch(buf,(int16_t*)out,ch.slice_size,ch.bsize,ch.nblocks,t_endian==ANA_LITTLE_ENDIAN); break;
      case(1): rv=anadecrunch8(buf,(int8_t*)out,ch.slice_size,ch.bsize,ch.nblocks,t_endian==ANA_LITTLE_ENDIAN); break;
      case(2): rv=anadecrunchrun(buf,(int16_t*)out,ch.slice_size,ch.bsize,ch.nblocks,t_endian==ANA_LITTLE_ENDIAN); break;
      case(3): rv=anadecrunchrun8(buf,(int8_t*)out,ch.slice_size,ch.bsize,ch.nblocks,t_endian==ANA_LITTLE_ENDIAN); break;
      case(4): rv=anadecrunch32(buf,(int32_t*)out,ch.slice_size,ch.bsize,ch.nblocks,t_endian==ANA_LITTLE_ENDIAN); break;
      default: fprintf(stderr,"error in data type for compressed data, fh.datyp =%d\n",fh.datyp);
    }
    free(buf);
    *osz=n_elem*type_sizes[*type];
    return out;
  }else{                            // uncompressed
    int size=n_elem*type_sizes[*type];
    uint8_t *out=malloc(size);
    if(fread(out,1,size,fin)<size){
      fclose(fin);
      fprintf(stderr,"error: unexpected end of file\n");
    }
    fclose(fin);
    if(swap_endian) // endianness is wrong
      switch(*type){
        case(INT16): bswapi16((int16_t*)out,n_elem); break;
        case(INT32):
        case(FLOAT32): bswapi32((int32_t*)out,n_elem); break;
        case(FLOAT64): bswapi64((int64_t*)out,n_elem); break;
      }
    *osz=size;
    return out; 
  } // end if(compressed)
}

void ana_fzwrite(uint8_t *data,char *file_name,int *ds,int nd,char *header,int type)	/* fcwrite subroutine */
{ // write standard f0 files, compressed format
  FILE *f=fopen(file_name,"w");
  fzhead_t fh;
  memset(&fh,0,sizeof(fzhead_t));
  int one=1;
  int t_endian=(*(char*)&one==0);    // an endian detector, taken from SL's tiff library 
  // switch(idl_type){
  //   case(IDL_TYP_BYTE): type=INT8; break;
  //   case(IDL_TYP_INT):  type=INT16; break;
  //   case(IDL_TYP_LONG): type=INT32; break;
  //   case(IDL_TYP_FLOAT): type=FLOAT32; break;
  //   case(IDL_TYP_DOUBLE): type=FLOAT64; break;
  //   default:{
  //     fprintf(stderr,"ana_fzwrite: error: ana_fzwrite: unknown variable type!\n");
  //     return;
  //   }
  // }
  
  if(t_endian){ // BIG_ENDIAN
    fh.synch_pattern=0xaaaa5555;
  }else{        // LITTLE_ENDIAN
    fh.synch_pattern=0x5555aaaa;
  }
  fh.source=0;
  fh.nhb=1;					// may be changed later
  fh.datyp=type;
  fh.ndim=nd;
//
  int n_elem=1,n;
  for(n=0;n<nd;++n) n_elem*=(fh.dim[n]=ds[n]);


  int type_sizes[]=ANA_VAR_SZ;
  int size=n_elem*type_sizes[type];
  if(t_endian){ // big endian platform
    switch(type){
      case(INT16): bswapi16((int16_t*)data,n_elem); break;
      case(INT32):
      case(FLOAT32): bswapi32((int32_t*)data,n_elem); break;
      case(FLOAT64): bswapi64((int64_t*)data,n_elem); break;
//      case(INT64): fprintf(stderr,"ana_fzwrite: error: ana_fzwrite: unknown variable type!\n");
    }
    bswapi32(fh.dim,nd);
  }
  if(header){
    int len=min(strlen(header),255);
    strncpy(fh.txt,header,len);
    fh.txt[len]=0;
  }
  int res=size;
  if(t_endian) bswapi32(&res,1);
  int i;
  for(i=0;i<4;i++) fh.cbytes[i]=((uint8_t*)&res)[i];
  fwrite(&fh,sizeof(fzhead_t),1,f);  // write header
  fwrite(data,1,size,f);
  fclose(f);
  if(t_endian){ // big endian platform: swap back
    switch(type){
      case(INT16): bswapi16((int16_t*)data,n_elem); break;
      case(INT32):
      case(FLOAT32): bswapi32((int32_t*)data,n_elem); break;
      case(FLOAT64): bswapi64((int64_t*)data,n_elem); break;
    }
  }
}

void ana_fcwrite(uint8_t *data,char *file_name,int *ds,int nd,char *header,int type,int slice)	/* fcwrite subroutine */
{ // write standard f0 files, compressed format
  FILE *f=fopen(file_name,"w");
  fzhead_t fh;
  memset(&fh,0,sizeof(fzhead_t));
  int one=1;
  int t_endian=(*(char*)&one==0);    // an endian detector, taken from SL's tiff library 

  if(t_endian){ // BIG_ENDIAN
    fh.synch_pattern=0xaaaa5555;
    fh.subf=129;
  }else{        // LITTLE_ENDIAN
    fh.synch_pattern=0x5555aaaa;
    fh.subf=1;
  }
  fh.source=0;
  fh.nhb=1;					// may be changed later
  fh.datyp=type;
  fh.ndim=nd;
//
  int n_elem=1,n;
  for(n=0;n<nd;++n) n_elem*=(fh.dim[n]=ds[n]);
  int nx=fh.dim[0];
  int ny=n_elem/nx;
  int type_sizes[]=ANA_VAR_SZ;
  int size=n_elem*type_sizes[type];
  if(t_endian){ // big endian platform
/*
    switch(type){
      case(INT16): swap((int16_t*)data,n_elem); break;
      case(INT32):
      case(FLOAT32): swap((int32_t*)data,n_elem); break;
      case(FLOAT64): swap((int64_t*)data,n_elem); break;
    }
*/
    bswapi32(fh.dim,nd);
  }
  if(header){
    int len=min(strlen(header),255);
    strncpy(fh.txt,header,len);
    fh.txt[len]=0;
  }
/* now compress the array, must be a byte or short */
/* extended to 32 bits 2/4/96 */
  int res,crunch_slice=slice,runlengthflag=0,limit=size+size/2; // reserve a bit extra just in case
  uint8_t *q=malloc(limit);
  switch(type){
    case(0):{
      if(runlengthflag)
        res=anacrunchrun8(q,data,crunch_slice,nx,ny,limit,t_endian);
      else
        res=anacrunch8(q,data,crunch_slice,nx,ny,limit,t_endian);
      break;
    }
    case(1):{
      if(runlengthflag)
        res=anacrunchrun(q,(int16_t*)data,crunch_slice,nx,ny,limit,t_endian);
      else
        res=anacrunch(q,(int16_t*)data,crunch_slice,nx,ny,limit,t_endian);
      break;
    }
    case(2):{
      if(runlengthflag){
        fprintf(stderr,"ana_fcwrite: warning: FCRUNWRITE not supported for I*4 yet\n");
        fclose(f);
        free(q);
        return;
      }else
        res=anacrunch32(q,(int32_t*)data,crunch_slice,nx,ny,limit,t_endian);
      break;
    }
    default:{
      fprintf(stderr,"ana_fcwrite: warning: FCWRITE: unsupported variable type.\n");
      fclose(f);
      free(q);
      return;
    }
  }
  if(res<0){
    fprintf(stderr,"ana_fcwrite: warning: not enough space allocated (%d bytes) for compressed array, trying uncompressed!\n",limit);
    free(q);
    fclose(f);
    ana_fzwrite(data,file_name,ds,nd,header,type);
    return;
  }
  if(res>size){
    fprintf(stderr,"ana_fcwrite: warning: compressed data (%d bytes) larger than raw data (%d bytes), writing uncompressed!\n",limit,size);
    free(q);
    fclose(f);
    ana_fzwrite(data,file_name,ds,nd,header,type);
    return;
  }
  size=res;
  if(t_endian) bswapi32(&res,1);
  int i;
  for(i=0;i<=3;++i) fh.cbytes[i]=((uint8_t*)&res)[i];
  fwrite(&fh,1,sizeof(fzhead_t),f);  // write header
  fwrite(q,1,size,f);
  free(q);
  fclose(f);
}

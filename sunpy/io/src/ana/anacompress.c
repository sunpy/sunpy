#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#include "types.h"
#include "anacompress.h"

int anacrunchrun8(uint8_t *x,uint8_t *array,int slice,int nx,int ny,int limit,int t_endian)
 /* compress 8 bit array into x (a uint8_t array) using ny blocks each of size
 nx, bit slice size slice, returns # of bytes in x */
 {
  uint8_t bits[8]={1,2,4,8,16,32,64,128};
 struct compresshead *ch;
 uint8_t	*p;
 unsigned nb;
 unsigned register i,j,r1;
 int	r0,r2,r3,mask,fac, nrun, lrun, ic;
 int	*dif, *d, nc, zq, yq, *dd;
 int	i2,k,iy;

 union { int i; short w; unsigned char b[4]; } y;
 /* begin execution */
 if (limit<25) { printf("limit (%d) too small in crunchrun8\n", limit); return -1;}
 limit = limit - 24;	/* need 14 for header and some margin since
 			we don't check all times */
 mask=1; for (i=0;i<slice;i++) mask=2*mask;
 fac=mask;       mask=mask-1; /* no inline expon. in C */
 /* determine the # of bytes to transfer to 32 bit int for fixed portion */
 if (slice == 0) nb=0; else { if (slice < 2 ) nb=1;
	 else { if (slice < 10) nb=2; else nb=3;    }};
 y.i=0;
 /* do the compression header */
 ch = (struct compresshead *) x;
 /* important note - can't use the sizeof(struct compresshead) because it
	 is 14 on some machines and rounded up to 16 on others */
 x = x + 14;
 ch->bsize = nx;  ch->nblocks = ny;  ch->slice_size = slice;  ch->type = 3;
 i=0;    r1=0;
 dif = (int *) malloc(nx*4);			/* line buffer */
 for (iy=0;iy<ny;iy++) {                 	/* start of iy (outer) loop */
 /* load the first value */
 x[i++] = array[iy*nx];
 
 /* compute and store the first differences for this line */
 p = (array+nx*iy);	nc=nx-1;
 d=dif; yq=(int) *p++;	zq=(int) *p++;
 while (nc--) { *d++ = zq - yq; yq = zq; zq = (int) *p++; }
 r1=r1+8;
 r2=0;
 p = (array+nx*iy);	nc=nx-1;
 d=dif;
 ic = i++;			/* count position */
 r1=r1+8;		/* to cover first count */
 lrun = 0;		/* literal count) */
 while (1) {
 /* look for a run */
 /* printf("*d, *(d+1) = %d %d, nc = %d\n",*d, *(d+1), nc );*/
 y.i = *d++;
 if (nc > 1) {
 while ( y.i == *d ) {	/* at least a run of 2 */
 dd = d+1;	nrun = 2;
 while ( nc-- > 2 && y.i == *dd) {
 /* printf("run!, y.i, *dd = %d %d, nc = %d\n", y.i, *dd, nc ); */
 nrun++;  dd++; }
 /* short runs are not worth it, make the legal limit 4 */
 /* printf("nrun = %d, nc = %d\n", nrun,nc);*/
 if ( nrun >= 4 ) {	/* code the run */
 /* a previous literal ? */
 if (lrun != 0) {
 /* printf("previous was literal, ic, i = %d %d\n", ic,i);*/
 x[ic] = lrun;	i = (r1+7)/8;	lrun = 0;
 /* printf("now, i = %d\n",i );*/
 } else i=ic;
 while (nrun  > 128 )	{ /* a big one, multiple runs */
 /* printf("big run, nrun = %d\n", nrun); */
 /* need only 2 bytes to represent run, runs can't be 17 bits */
 if (nrun == 129)	/* beware the dreaded 129 */
   { x[i++] = 0x82; nrun -= 127;} else { x[i++] = 0x81;	nrun -= 128; }
  if(t_endian){
    x[i++]=y.b[3];	x[i++]=y.b[2];
  }else{
    x[i++]=y.b[0];	x[i++]=y.b[1];
  }
 }
 /* printf("encoding run, nrun = %d, i=%d, iy = %d\n",nrun,i,iy); */
  if(t_endian){ // big endian
    x[i++] = -(nrun-1);	x[i++]=y.b[3];	x[i++]=y.b[2];
  }else{
    x[i++] = -(nrun-1);	x[i++]=y.b[0];	x[i++]=y.b[1];
  }
 /* prepare for a literal and for next run check */
 nc--;
 if (nc <= 0) goto ended_on_run;
 lrun = 0;	ic = i++;	r1 = 8*i;	d = dd;
 y.i = *d++;
 /* printf("after a run, new y.i = %d, new *d = %d\n", y.i, *d);*/
 } else { nc = nc + nrun -1;	break; }
 }	/* not a run, do next literal, assume setup for literals */
 } else if (nc <= 0) break;
 nc--;
 /* increment the literal count */
 /* printf("literal, lrun = %d, nc = %d, ic,i = %d %d\n", lrun, nc, ic,i);*/
 if (++lrun > 127)	{		/* need a new literal run count */
 x[ic] = 127;
 /* printf("ic = %d, i,r1 = %d %d\n", ic,i,r1); */
 /* bump to next byte boundary */
 i = (r1+7)/8;	ic = i++;     r1 = 8*i;		lrun = 1;
 /* printf("next ic = %d\n", ic); */
 }
 /* first the fixed slice portion */
 /* printf("y.i = %d\n", y.i);*/
 r3=(y.i>>slice);
 i=r1>>3;	/* byte number */
 j=r1 & 7;		/* bit number */
 if ( i > limit ) return -1;			/* bad news, went too far */
 /* now load nb bytes into x */
 /*low order byte of y.i is first in stream */
  if(t_endian){ // big endian
    if (j == 0) {y.i=(y.i & mask); x[i]=y.b[3];}
	 else { y.i=(y.i & mask)<<j; x[i]=x[i] | y.b[3];}
    if (nb>1) { x[i+1]=y.b[2]; if (nb>2) x[i+2]=y.b[1]; }
  }else{
    if (j == 0) {y.i=(y.i & mask); x[i]=y.b[0];}
	 else { y.i=(y.i & mask)<<j; x[i]=x[i] | y.b[0];}
    if (nb>1) { x[i+1]=y.b[1]; if (nb>2) x[i+2]=y.b[2]; }
  }
 
 r1=r1+slice;       			/* bump r1 pass the fixed part */
 i=r1>>3;                j=r1 & 7;
 /* note that r3 is the # of bits required minus 1 */
 if (r3==0) { if (j ==0 ) {x[i]=bits[j];} else {x[i]=x[i]|bits[j];}
 r1+=1;} else    {
 r3=2*r3;        if (r3<0) r3 = -r3-1;
 if (r3<31)  {
 r0=j+r3;        /* this is the bit that needs setting offset from x[i] */
 if (r0 < 8) { if (j == 0) x[i]=bits[r0]; else x[i]=x[i]|bits[r0];}
 else {if (j == 0) x[i]=0;  j=r0%8; if (r0 < 16) x[i+1]=bits[j];
 else { i2=i+r0/8; for (k=i+1;k<i2;k++) x[k]=0;  x[i2]=bits[j]; }
 }
 r1+=1+r3;
 } else {        /* big one exception, should be rare */
 /* does not need to be efficient, used rarely */
 /* printf("big one \n");*/
 if (j == 0) x[i]=0;     /* gotta zero the virgins */
 r0=j+31;        j=r0%8;         i2=i+r0/8;
 for (k=i+1;k<i2;k++) x[k]=0;    x[i2]=bits[j];
 /* recompute the difference and load 9 bits (always 2 bytes) */
 r1=r1+32;
 i=r1/8;
 j=r1%8;
 if (j == 0) x[i]=0;     /* gotta zero the virgins */
 y.i=((*(d-1))& 0x1ff) << j;
  if(t_endian){ // big endian
    x[i]=x[i] | y.b[3]; x[i+1]=y.b[2];
  }else{
    x[i]=x[i] | y.b[0]; x[i+1]=y.b[1];
  }
 r1=r1+9;       } /* end of big one exception */
		 }       /* end of (r3==0) conditional */
 /* printf("end of literal, i,r1 = %d %d\n", i,r1);*/
 /* printf(" x[r1/8] = %d\n",  x[r1/8]);*/
		 }       /* end of ix loop */
 /* some checks here */
 /* bump to next byte boundary */
 /* a final literal ? */
 /* printf("final lrun = %d, ic = %d\n", lrun, ic);*/
 if (lrun != 0) { x[ic] = lrun;	lrun = 0; }
 i=(r1+7)/8;
 ended_on_run:
 r1=8*i;
		  }       /* end of iy loop */
 ch->tsize = i = i + 14;
 /* we have to put these in a form readable by the Vax (these may be used
	 by fcwrite) */
  if(t_endian){ // big endian
   bswapi32(&(ch->tsize),1); bswapi32(&(ch->bsize),1); bswapi32(&(ch->nblocks),1); 
  }
 free(dif);
 return  i;      /*return # of bytes used */
 }


int anacrunch8(uint8_t *x,uint8_t *array,int slice,int nx,int ny,int limit,int t_endian)
 /* compress 8 bit array into x (a byte array) using ny blocks each of size
 nx, bit slice size slice, returns # of bytes in x */
 {
  uint8_t bits[8]={1,2,4,8,16,32,64,128};
 struct compresshead *ch;

 unsigned nb,ixa,ixb;
 unsigned register i,j,r1,in;
 int r0,r2,r3,mask,fac;
 int i2,k,iy;

 union { int i; short w; unsigned char b[4]; } y;
 /* begin execution */
 if (limit<25) { printf("limit (%d) too small in crunch8\n", limit); return -1;}
 limit = limit - 24;	/* need 14 for header and some margin since
 			we don't check all times */
 mask=1; for (i=0;i<slice;i++) mask=2*mask;
 fac=mask;       mask=mask-1; /* no inline expon. in C */
 /* determine the # of bytes to transfer to 32 bit int for fixed portion */
 if (slice > 8) slice = 8;
 if (slice == 0) nb=0; else { if (slice < 2 ) nb=1;
	 else { if (slice < 10) nb=2; else nb=3;    }};
 y.i=0;
 /* do the compression header */
 ch = (struct compresshead *) x;
 /* important note - can't use the sizeof(struct compresshead) because it
	 is 14 on some machines and rounded up to 16 on others */
 /*x = x + sizeof(struct compresshead);*/
 x = x + 14;
 ch->bsize = nx;  ch->nblocks = ny;  ch->slice_size = slice;  ch->type = 1;
 i=0;    r1=0;   in=0;
 for (iy=0;iy<ny;iy++) {                 /* start of iy (outer) loop */
 /* load the first value */
 x[i] = array[in];
 r1=r1+8;
 r2=0;
 ixa=1+iy*nx;    ixb=(iy+1)*nx;
 for (in=ixa; in<ixb; in++)      {               /* start of ix (inner) loop */
 /* first the fixed slice portion */
 y.i= (int) array[in]- (int) array[in-1];
 r3=(y.i>>slice);
 i=r1>>3;
 j=r1%8;
 if ( i > limit ) return -1;		/* bad news, went too far */
 /* now load nb bytes into x */
 /*low order byte of y.i is first in stream */
  if(t_endian){ // big endian
    if (j == 0) {y.i=(y.i & mask); x[i]=y.b[3];}
	 else { y.i=(y.i & mask)<<j; x[i]=x[i] | y.b[3];}
    if (nb>1) { x[i+1]=y.b[2]; }
  }else{
    if (j == 0) {y.i=(y.i & mask); x[i]=y.b[0];}
	 else { y.i=(y.i & mask)<<j; x[i]=x[i] | y.b[0];}
    if (nb>1) { x[i+1]=y.b[1]; }
  }
 r1=r1+slice;       /* bump r1 pass the fixed part */
 i=r1>>3;                j=r1%8;
 /* note that r3 is the # of bits required minus 1 */
 if (r3==0) { if (j ==0 ) {x[i]=bits[j];} else {x[i]=x[i]|bits[j];}
 r1+=1;} else    {
 r3=2*r3;        if (r3<0) r3 = -r3-1;
 if (r3<31)  {
 r0=j+r3;        /* this is the bit that needs setting offset from x[i] */
 if (r0 < 8) { if (j == 0) x[i]=bits[r0]; else x[i]=x[i]|bits[r0];}
 else {if (j == 0) x[i]=0;  j=r0%8; if (r0 < 16) x[i+1]=bits[j];
 else { i2=i+r0/8; for (k=i+1;k<i2;k++) x[k]=0;  x[i2]=bits[j]; }
 }
 r1+=1+r3;
 } else {        /* big one exception, should be rare */
 /* does not need to be efficient, used rarely */
 /* printf("big one \n"); */
 if (j == 0) x[i]=0;     /* gotta zero the virgins */
 r0=j+31;        j=r0%8;         i2=i+r0/8;
 for (k=i+1;k<i2;k++) x[k]=0;    x[i2]=bits[j];
 /* recompute the difference and load 9 bits (always 2 bytes) */
 r1=r1+32;
 i=r1/8;
 j=r1%8;
 if (j == 0) x[i]=0;     /* gotta zero the virgins */
 y.i=((array[in]-array[in-1])& 0x1ff) << j;
  if(t_endian){ // big endian
    x[i]=x[i] | y.b[3]; x[i+1]=y.b[2];
  }else{
    x[i]=x[i] | y.b[0]; x[i+1]=y.b[1];
  }
 r1=r1+9;       } /* end of big one exception */
		 }       /* end of (r3==0) conditional */
		 }       /* end of ix loop */
 /* some checks here */
 /* bump to next byte boundary */
 i=(r1+7)/8;     r1=8*i;                 }       /* end of iy loop */
 ch->tsize = i = i + 14;
 /* we have to put these in a form readable by the Vax (these may be used
	 by fcwrite) */
  if(t_endian){ // big endian
    bswapi32(&(ch->tsize),1); bswapi32(&(ch->bsize),1); bswapi32(&(ch->nblocks),1);
  }
 return  i;      /*return # of bytes used */
 }       /* end of routine */


int anacrunchrun(uint8_t *x,int16_t *array,int slice,int nx,int ny,int limit,int t_endian)
 /* compress 16 bit array into x (a byte array) using ny blocks each of size
 nx, bit slice size slice, returns # of bytes in x */
 {
  uint8_t bits[8]={1,2,4,8,16,32,64,128};
  struct compresshead *ch;
 short *p;
 unsigned nb;
 unsigned register i,j,r1;
 int	r0,r2,r3,mask,fac, nrun, lrun, ic;
 int	*dif, *d, nc, zq, yq, *dd;

 int	i2,k,iy;

 union { int i; short w; unsigned char b[4]; } y;
 /* begin execution */
 if (limit<25) { printf("limit (%d) too small in crunchrun\n", limit); return -1;}
 limit = limit - 24;	/* need 14 for header and some margin since
 			we don't check all times */
 mask=1; for (i=0;i<slice;i++) mask=2*mask;
 fac=mask;       mask=mask-1; /* no inline expon. in C */
 /* determine the # of bytes to transfer to 32 bit int for fixed portion */
 if (slice == 0) nb=0; else { if (slice < 2 ) nb=1;
	 else { if (slice < 10) nb=2; else nb=3;    }};
 y.i=0;
 /* do the compression header */
 ch = (struct compresshead *) x;
 /* important note - can't use the sizeof(struct compresshead) because it
	 is 14 on some machines and rounded up to 16 on others */
 x = x + 14;
 ch->bsize = nx;  ch->nblocks = ny;  ch->slice_size = slice;  ch->type = 2;
 i=0;    r1=0;
 dif = (int *) malloc(nx*4);			/* line buffer */
 for (iy=0;iy<ny;iy++) {                 	/* start of iy (outer) loop */
 /* load the first value, reverse bytes (VAX style)*/
  if(t_endian){ // big endian
    y.w=array[iy*nx]   ;x[i++]=y.b[1]    ;x[i++]=y.b[0];
  }else{ // little endian
    y.w=array[iy*nx]   ;x[i++]=y.b[0]    ;x[i++]=y.b[1];
  }
  /* compute and store the first differences for this line */
 p = (short int*)(array+nx*iy);	nc=nx-1;
 d=dif; yq=(int) *p++;	zq=(int) *p++;
 while (nc--) { *d++ = zq - yq; yq = zq; zq = (int) *p++; }
 r1=r1+16;
 r2=0;
 p = (short int*)(array+nx*iy);	nc=nx-1;
 d=dif;
 ic = i++;			/* count position */
 r1=r1+8;		/* to cover first count */
 lrun = 0;		/* literal count) */
 while (1) {
 /* look for a run */
 /* printf("*d, *(d+1) = %d %d, nc = %d\n",*d, *(d+1), nc );*/
 y.i = *d++;
 if (nc > 1) {
 while ( y.i == *d ) {	/* at least a run of 2 */
 dd = d+1;	nrun = 2;
 while ( nc-- > 2 && y.i == *dd) {
 /* printf("run!, y.i, *dd = %d %d, nc = %d\n", y.i, *dd, nc ); */
 nrun++;  dd++; }
 /* short runs are not worth it, make the legal limit 4 */
 /* printf("nrun = %d, nc = %d\n", nrun,nc);*/
 if ( nrun >= 4 ) {	/* code the run */
 /* a previous literal ? */
 if (lrun != 0) {
 /* printf("previous was literal, ic, i = %d %d\n", ic,i);*/
 x[ic] = lrun;	i = (r1+7)/8;	lrun = 0;
 /* printf("now, i = %d\n",i );*/
 } else i=ic;
 while (nrun  > 128 )	{ /* a big one, multiple runs */
 /* need only 2 bytes to represent run, runs can't be 17 bits */
 if (nrun == 129)	/* beware the dreaded 129 */
   { x[i++] = 0x82; nrun -= 127;} else { x[i++] = 0x81;	nrun -= 128; }
  if(t_endian){ // big endian
    x[i++]=y.b[3];	x[i++]=y.b[2];
  }else{
    x[i++]=y.b[0];	x[i++]=y.b[1];
  }
 }
 /* printf("encoding run, nrun = %d, i=%d, iy = %d\n",nrun,i,iy); */
  if(t_endian){ // big endian
    x[i++] = -(nrun-1);	x[i++]=y.b[3];	x[i++]=y.b[2];
  }else{
    x[i++] = -(nrun-1);	x[i++]=y.b[0];	x[i++]=y.b[1];
  }
 /* prepare for a literal and for next run check */
 nc--;
 if (nc <= 0) goto ended_on_run;
 lrun = 0;	ic = i++;	r1 = 8*i;	d = dd;
 y.i = *d++;
 /* printf("after a run, new y.i = %d, new *d = %d\n", y.i, *d);*/
 } else { nc = nc + nrun -1;	break; }
 }	/* not a run, do next literal, assume setup for literals */
 } else if (nc <= 0) break;
 nc--;
 /* increment the literal count */
 /* printf("literal, lrun = %d, nc = %d, ic,i = %d %d\n", lrun, nc, ic,i);*/
 if (++lrun > 127)	{		/* need a new literal run count */
 x[ic] = 127;
 /* printf("ic = %d, i,r1 = %d %d\n", ic,i,r1); */
 /* bump to next byte boundary */
 i = (r1+7)/8;	ic = i++;     r1 = 8*i;		lrun = 1;
 /* printf("next ic = %d\n", ic); */
 }
 /* first the fixed slice portion */
 /* printf("y.i = %d\n", y.i);*/
 r3=(y.i>>slice);
 i=r1>>3;	/* byte number */
 j=r1 & 7;		/* bit number */
 if ( i > limit ) return -1;			/* bad news, went too far */
 /* now load nb bytes into x */
 /*low order byte of y.i is first in stream */
  if(t_endian){ // big endian
    if (j == 0) {y.i=(y.i & mask); x[i]=y.b[3];}
	 else { y.i=(y.i & mask)<<j; x[i]=x[i] | y.b[3];}
    if (nb>1) { x[i+1]=y.b[2]; if (nb>2) x[i+2]=y.b[1]; }
  }else{
    if (j == 0) {y.i=(y.i & mask); x[i]=y.b[0];}
	 else { y.i=(y.i & mask)<<j; x[i]=x[i] | y.b[0];}
    if (nb>1) { x[i+1]=y.b[1]; if (nb>2) x[i+2]=y.b[2]; }
  }
 r1=r1+slice;       			/* bump r1 pass the fixed part */
 i=r1>>3;                j=r1 & 7;
 /* note that r3 is the # of bits required minus 1 */
 /* printf("r3 = %d\n", r3);*/
 if (r3==0) { if (j ==0 ) {x[i]=bits[j];} else {x[i]=x[i]|bits[j];}
 r1+=1;} else    {
 r3=2*r3;        if (r3<0) r3 = -r3-1;
 if (r3<31)  {
 r0=j+r3;        /* this is the bit that needs setting offset from x[i] */
 if (r0 < 8) { if (j == 0) x[i]=bits[r0]; else x[i]=x[i]|bits[r0];}
 else {if (j == 0) x[i]=0;  j=r0%8; if (r0 < 16) x[i+1]=bits[j];
 else { i2=i+r0/8; for (k=i+1;k<i2;k++) x[k]=0;  x[i2]=bits[j]; }
 }
 r1+=1+r3;
 } else {        /* big one exception, should be rare */
 /* does not need to be efficient, used rarely */
 /* printf("big one \n");*/
 if (j == 0) x[i]=0;     /* gotta zero the virgins */
 r0=j+31;        j=r0%8;         i2=i+r0/8;
 for (k=i+1;k<i2;k++) x[k]=0;    x[i2]=bits[j];
 /* recompute the difference and load 17 bits (always 3 bytes) */
 r1=r1+32;
 i=r1/8;
 j=r1%8;
 if (j == 0) x[i]=0;     /* gotta zero the virgins */
 y.i=((*(d-1))& 0x1ffff) << j;
  if(t_endian){ // big endian
    x[i]=x[i] | y.b[3]; x[i+1]=y.b[2];      x[i+2]=y.b[1];
  }else{
    x[i]=x[i] | y.b[0]; x[i+1]=y.b[1];      x[i+2]=y.b[2];
  }
 r1=r1+17;       } /* end of big one exception */
		 }       /* end of (r3==0) conditional */
 /* printf("end of literal, i,r1 = %d %d\n", i,r1);*/
 /* printf(" x[r1/8] = %d\n",  x[r1/8]);*/
		 }       /* end of ix loop */
 /* some checks here */
 /* bump to next byte boundary */
 /* a final literal ? */
 /* printf("final lrun = %d, ic = %d\n", lrun, ic);*/
 if (lrun != 0) { x[ic] = lrun;	lrun = 0; }
 i=(r1+7)/8;
 ended_on_run:
 r1=8*i;
		  }       /* end of iy loop */
 ch->tsize = i = i + 14;
 /* we have to put these in a form readable by the Vax (these may be used
	 by fcwrite) */
 if(t_endian){ // big endian
   bswapi32(&(ch->tsize),1); bswapi32(&(ch->bsize),1); bswapi32(&(ch->nblocks),1); 
 }
 free(dif);
 return  i;      /*return # of bytes used */
 }       /* end of routine */



int anacrunch(uint8_t *x,int16_t *array,int slice,int nx,int ny,int limit,int t_endian)
// compress 16 bit array into x (a byte array) using ny blocks each of size
// nx, bit slice size slice, returns # of bytes in x 
 {
  uint8_t bits[8]={1,2,4,8,16,32,64,128};
  unsigned register i,j,r1,in;
  int r0,r2,r3,mask,fac;

  union{
    int i;
    short w;
    unsigned char b[4];
  } y;

  if(limit<25){
    printf("limit (%d) too small in crunch\n", limit);
    return -1;
  }
  limit-=24;                              // need 14 for header and some margin since we don't check all times
  mask=1;
  for(i=0;i<slice;i++) mask*=2;
  fac=mask;
  mask-=1;                                // no inline expon. in C
  unsigned nb;                            // determine the # of bytes to transfer to 32 bit int for fixed portion
  if(slice==0){
    nb=0; 
  }else{
    if(slice<2){ 
      nb=1;
    }else{
      if(slice<10) nb=2; else nb=3;
    }
  }
  y.i=0;
  struct compresshead *ch=(struct compresshead*)x;  // the compression header
// important note: can't use the sizeof(struct compresshead) because it
//                 is 14 on some machines and rounded up to 16 on others
//  x+=sizeof(struct compresshead);
  x+=14;
  ch->bsize=nx;
  ch->nblocks=ny;
  ch->slice_size=slice;
  ch->type=0;
  r1=0;                                   // r1 is the bit index in the stream...?
  in=0;                                   // in is the byte index in the uncompressed stream...?
  i=0;                                    // i is the byte index in the compressed stream...?
  int iy;
  for(iy=0;iy<ny;++iy){               // start of iy (outer) loop
    y.w=array[in];
    if(t_endian){                         // load the first value, reverse bytes (VAX style)
      x[i]=y.b[1];
      x[i+1]=y.b[0];
    }else{
      x[i+1]=y.b[1];
      x[i]=y.b[0];
    }
    r1+=16;
    r2=0;
    unsigned int iynx=iy*nx;
    for(in=iynx+1;in<iynx+nx;++in){        // start of ix (inner) loop
      y.i=array[in]-array[in-1];           // first the fixed slice portion
      r3=(y.i>>slice);
      i=r1>>3;                             // compressed data size (number of bits/8)
      j=r1%8;
      if(i>limit) return -1;               // bad news: compressed data too big...
      if(j==0){                            // now load nb bytes into x, low order byte of y.i is first in stream
        y.i=(y.i&mask);
        x[i]=(uint8_t)y.i;
        if(slice>8) x[i+1]=(uint8_t)(y.i>>8); // since we started at bit 0, spillover to the next byte is determined as follows (and is unlikely since slice gt 8 is unusual
      }else{
        y.i=(y.i&mask)<<j;
        x[i]=x[i]|(uint8_t)y.i;
	if(nb>1){                          // spillover more likely here
          x[i+1]=(uint8_t)(y.i>>8);
	  if(nb>2) x[i+2]=(uint8_t)(y.i>>16);
        }
      }
      r1+=slice;                           // bump r1 pass the fixed part
      i=r1>>3;
      j=r1%8;
      if(r3==0){                           // note that r3 is the # of bits required minus 1
        if(j==0){
          x[i]=*bits;
        }else{
          x[i]=x[i]|bits[j];
        }
        ++r1;
      }else{
        r3*=2;
        if(r3<0) r3=-r3-1;
        if(r3<31){
          r0=j+r3;                         // this is the bit that needs setting offset from x[i]
          if(r0<8){
            if(j==0) x[i]=bits[r0]; else x[i]=x[i]|bits[r0];
          }else{                           // note, discovered on 2/3/96 that the j==0 case not done above for the sunbow version, was OK in the umbra version, may have happened while cleaning up code?, caused extra bits to be set if x[i] wasn't zero
            if(j==0) x[i]=0;
            j=r0%8;
            if(r0<16){
              x[i+1]=bits[j];
            }else{
              int i2=i+r0/8;
              int k;
              for(k=i+1;k<i2;++k) x[k]=0;
              x[i2]=bits[j];
            }
          }
          r1+=r3+1;
        }else{             // big one exception, should be rare so does not need to be efficient
          if(j==0) x[i]=0; // gotta zero the virgins
          r0=j+31;
          j=r0%8;
          int i2=i+r0/8;
          int k;
          for(k=i+1;k<i2;++k) x[k]=0;
          x[i2]=bits[j];
          r1+=32;          // recompute the difference and load 17 bits (always 3 bytes)
          i=r1/8;
          j=r1%8;
          if(j==0) x[i]=0; // gotta zero the virgins
          y.i=((array[in]-array[in-1])&0x1ffff)<<j;
          if(t_endian){    // big endian
            x[i]=x[i]|y.b[3];
            x[i+1]=y.b[2];
            x[i+2]=y.b[1];
          }else{
            x[i]=x[i]|y.b[0];
            x[i+1]=y.b[1];
            x[i+2]=y.b[2];
          }
          r1+=17;
        }        // end of big one exception
      }          // end of (r3==0) conditional
    }            // end of ix loop
    i=(r1+7)/8;  // bump to next byte boundary
    r1=8*i;
  }              // end of iy loop
  ch->tsize=(i+=14);
  if(t_endian){  // we have to put these in a form readable by the Vax (these may be used by fcwrite)
    bswapi32(&(ch->tsize),1);
    bswapi32(&(ch->bsize),1);
    bswapi32(&(ch->nblocks),1); 
  }
  return i;     // return # of bytes used
}


int anacrunch32(uint8_t *x,int32_t *array,int slice,int nx,int ny,int limit,int t_endian)
 /* compress 32 bit array into x (a byte array) using ny blocks each of size
 nx, bit slice size slice, returns # of bytes in x */
 {
  uint8_t bits[8]={1,2,4,8,16,32,64,128};
  struct compresshead *ch;

 unsigned int nb,ixa,ixb,big=0;
 unsigned register i,j,r1,in;
 int r0,r2,fac;
 long long	r3, mask, y64;
 int i2,k,iy;

 union { int i; short w; unsigned char b[4]; } y;
 union { long long l64; unsigned char b[8];  } yy;
 /* begin execution */
 if (limit<25) { printf("limit (%d) too small in crunch32\n", limit); return -1;}
 limit = limit - 24;	/* need 14 for header and some margin since
 			we don't check all times */
 mask=1; for (i=0;i<slice;i++) mask=2*mask;
 fac=mask;       mask=mask-1; /* no inline expon. in C */
 /* determine the # of bytes to transfer to 32 bit int for fixed portion */
 nb = (slice + 14)/8;	/* range 1 to 5 */
 if (slice == 0) nb=0;	/* but slice = 0 a special case */

 y.i=0;
 /* do the compression header */
 ch = (struct compresshead *) x;
 /* important note - can't use the sizeof(struct compresshead) because it
	 is 14 on some machines and rounded up to 16 on others */
 /*x = x + sizeof(struct compresshead);*/
 x = x + 14;
 ch->bsize = nx;  ch->nblocks = ny;  ch->slice_size = slice;  ch->type = 4;
 i=0;    r1=0;   in=0;
 for (iy=0;iy<ny;iy++) {                 /* start of iy (outer) loop */
 /* load the first value, reverse bytes (VAX style)*/
 if(t_endian){ // big endian
   y.i=array[in]; x[i]=y.b[3]; x[i+1]=y.b[2]; x[i+2]=y.b[1]; x[i+3]=y.b[0];
 }else{
   y.i=array[in]; x[i]=y.b[0]; x[i+1]=y.b[1]; x[i+2]=y.b[2]; x[i+3]=y.b[3];
 }
 r1=r1+32;
 r2=0;
 ixa=1+iy*nx;    ixb=(iy+1)*nx;
 for (in=ixa; in<ixb; in++)      {               /* start of ix (inner) loop */
 /* first the fixed slice portion */
 y64 = (long long) array[in] - (long long) array[in-1];
 r3 = (y64>>slice);
 i=r1>>3;
 j=r1%8;
 if ( i > limit ) return -1;		/* bad news, went too far */
 /* now load nb bytes into x */
 /*low order byte of y.i is first in stream */
 if (j == 0) {
	y64=(y64 & mask);	x[i]= (uint8_t) y64;
	/* since we started at bit 0, spillover to the next byte is
	determined as follows (and is unlikely since slice gt 8 is unusual */
	if (slice > 8) { x[i+1]= (uint8_t) (y64 >> 8);
	 if (slice > 16) { x[i+2]= (uint8_t) (y64 >> 16);
	  if (slice > 24) { x[i+3]= (uint8_t) (y64 >> 24);}}}
} else {
	y64=(y64 & mask)<<j;	x[i]=x[i] | (uint8_t) y64;
	/* spillover more likely here */
	if (nb>1) { x[i+1] = (uint8_t) (y64 >> 8);
	 if (nb>2) { x[i+2] = (uint8_t) (y64 >> 16);
	  if (nb>3) { x[i+3] = (uint8_t) (y64 >> 24);
	   if (nb>4) { x[i+4] = (uint8_t) (y64 >> 32);}}}}
}

 r1=r1+slice;       /* bump r1 pass the fixed part */
 i=r1>>3;                j=r1%8;
 /* note that r3 is the # of bits required minus 1 */
 if (r3==0) { if (j ==0 ) {x[i]= *bits;} else {x[i]=x[i]| bits[j];}
 r1++;} else    {
 r3=2*r3;        if (r3<0) r3 = -r3-1;
 if (r3<31)  {
 r0=j+r3;        /* this is the bit that needs setting offset from x[i] */
 if (r0 < 8) { if (j == 0) x[i]=bits[r0]; else x[i]=x[i]|bits[r0];}
   else  {
    if (j == 0) x[i]=0;
     j=r0%8;
     if (r0 < 16) x[i+1]=bits[j]; else {
 	i2=i+r0/8;
	for (k=i+1;k<i2;k++) x[k]=0;
	x[i2]=bits[j]; }
   }
   r1+=1+r3;
 } else {        /* big one exception, should be rare */
 /* does not need to be efficient, used rarely */
 /* printf("big one for I*4\n"); */
 big++;
 if (j == 0) x[i]=0;     /* gotta zero the virgins */
 r0=j+31;        j=r0%8;         i2=i+r0/8;
 for (k=i+1;k<i2;k++) x[k]=0;    x[i2]=bits[j];
 /* recompute the difference and load 33 bits (always 5 bytes) */
 r1=r1+32;
 i=r1/8;
 j=r1%8;
 if (j == 0) x[i]=0;     /* gotta zero the virgins */
 yy.l64=(((long long int) array[in]- (long long int) array[in-1])& 0x1ffffffffLL) << j;
  if(t_endian){ // big endian
   x[i]=x[i] | yy.b[7]; x[i+1]=yy.b[6]; x[i+2]=yy.b[5];
   x[i+3]=yy.b[4]; x[i+4]=yy.b[3];
  }else{
   x[i]=x[i] | yy.b[0]; x[i+1]=yy.b[1]; x[i+2]=yy.b[2];
   x[i+3]=yy.b[3]; x[i+4]=yy.b[4];
  }
 r1=r1+33;       } /* end of big one exception */
		 }       /* end of (r3==0) conditional */
 /*in=in+1; */                   }       /* end of ix loop */
 /* some checks here */
 /* bump to next byte boundary */
 i=(r1+7)/8;     r1=8*i;                 }       /* end of iy loop */
 ch->tsize = i = i + 14;
 /* we have to put these in a form readable by the Vax (these may be used
	 by fcwrite) */
 if(t_endian){ // big endian
   bswapi32(&(ch->tsize),1); bswapi32(&(ch->bsize),1); bswapi32(&(ch->nblocks),1); 
 }
 /* printf("number of big ones for this I*4 = %d\n", big); */
 return  i;      /*return # of bytes used */
 }       /* end of routine */


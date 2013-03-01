#ifndef __ANACOMPRESS_H__      // __ANACOMPRESS_H__
#define __ANACOMPRESS_H__      // __ANACOMPRESS_H__

int anacrunchrun8(uint8_t *x,uint8_t *array,int slice,int nx,int ny,int limit,int t_endian);
int anacrunch8(uint8_t *x,uint8_t *array,int slice,int nx,int ny,int limit,int t_endian);
int anacrunchrun(uint8_t *x,int16_t *array,int slice,int nx,int ny,int limit,int t_endian);
int anacrunch(uint8_t *x,int16_t *array,int slice,int nx,int ny,int limit,int t_endian);
int anacrunch32(uint8_t *x,int32_t *array,int slice,int nx,int ny,int limit,int t_endian);

#endif                         // __ANACOMPRESS_H__

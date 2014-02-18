#ifndef __ANADECOMPRESS_H__      // __ANADECOMPRESS_H__
#define __ANADECOMPRESS_H__      // __ANADECOMPRESS_H__

int anadecrunch32(unsigned char *x,int32_t *array,int r9,int nx,int ny,int);
int anadecrunchrun8(unsigned char *x,int8_t *array,int r9,int nx,int ny,int);
int anadecrunchrun(unsigned char *x,int16_t *array,int r9,int nx,int ny,int);
int anadecrunch8(unsigned char *x,int8_t *array,int r9,int nx,int ny,int);
int anadecrunch(unsigned char *x,int16_t *array,int r9,int nx,int ny,int);

#endif                           // __ANADECOMPRESS_H__

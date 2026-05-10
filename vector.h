/*
Parallelised by Samhain Ackerman and Henry Leap
*/

#ifndef __TYPES_H__
#define __TYPES_H__

typedef double vector3[3];
#define FILL_VECTOR(vector,a,b,c) {vector[0]=a;vector[1]=b;vector[2]=c;}
#define COPY_VECTOR(dest, src) {dest[0] = src[0]; dest[1] = src[1]; dest[2] = src[2];}
#define ADD_VECTORS(dest, src) {dest[0] += src[0]; dest[1] += src[1]; dest[2] += src[2];}

#endif

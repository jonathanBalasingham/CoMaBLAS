//
// Created by jon on 2/20/21.
//

#ifndef COMABLAS_COMA_ROT_H
#define COMABLAS_COMA_ROT_H
#include <complex.h>

void srot(unsigned int n, float* sx, int incx, float* sy, int incy, float c, float s);
void drot(unsigned int n, double* sx, int incx, double* sy, int incy, double c, double s);
void csrot(unsigned int n,complex float* sx, int incx,complex float* sy, int incy, float c, float s);
void zdrotf(unsigned int n,complex double* sx, int incx,complex double* sy, int incy, double c, double s);


#endif //COMABLAS_COMA_ROT_H

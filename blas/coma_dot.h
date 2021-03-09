//
// Created by jon on 2/20/21.
//

#ifndef COMABLAS_COMA_DOT_H
#define COMABLAS_COMA_DOT_H

#include <complex.h>

float sdot(unsigned int n, float* sx, int incx, float* sy, int incy);
double ddot(unsigned int n, const double* dx, int incx, const double* dy, int incy);
complex float cdotu(unsigned int n, const complex float* cx, int incx, const complex float* cy, int incy);
complex double zdotu(unsigned int n, const complex double* zx, int incx, const complex double* zy, int incy);
complex float cdotc(unsigned int n, const complex float* cx, int incx, const complex float* cy, int incy);
complex double zdotc(unsigned int n, const complex double* zx, int incx, const complex double* zy, int incy);


#endif //COMABLAS_COMA_DOT_H

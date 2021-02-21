//
// Created by jon on 2/20/21.
//

#ifndef COMABLAS_COMA_DOT_H
#define COMABLAS_COMA_DOT_H

#include <complex.h>

float sdot(unsigned int n, float* sx, int incx, float* sy, int incy);
double ddot(int n, const double* dx, int incx, const double* dy, int incy);
complex float cdotu(int n, complex float* cx, int incx, complex float* cy, int incy);
complex double zdotu(int n, complex double* zx, int incx, complex double* zy, int incy);

#endif //COMABLAS_COMA_DOT_H

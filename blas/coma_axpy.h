//
// Created by jon on 2/20/21.
//

#ifndef COMABLAS_COMA_AXPY_H
#define COMABLAS_COMA_AXPY_H
#include <complex.h>

void saxpy(unsigned int n, float a, const float* x, int incx, float* y, int incy);
void daxpy(unsigned int n, double a, const double* x, int incx, double* y, int incy);
void caxpy(unsigned int n, complex float a, const complex float* x, int incx, complex float* y, int incy);
void zaxpy(unsigned int n, complex double a, const complex double* x, int incx, complex double* y, int incy);

#endif //COMABLAS_COMA_AXPY_H

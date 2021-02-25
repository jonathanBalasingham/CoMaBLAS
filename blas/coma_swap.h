//
// Created by jon on 2/20/21.
//

#ifndef COMABLAS_COMA_SWAP_H
#define COMABLAS_COMA_SWAP_H

#include <complex.h>

void sswap(unsigned int n, float* x, int incx, float* y, int incy);
void dswap(unsigned int n, double* x, int incx, double* y, int incy);
void cswap(unsigned int n, complex float* x, int incx, complex float* y, int incy);
void zswap(unsigned int n, complex double* x, int incx, complex double* y, int incy);

#endif //COMABLAS_COMA_SWAP_H

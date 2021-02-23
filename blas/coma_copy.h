//
// Created by jon on 2/20/21.
//

#ifndef COMABLAS_COMA_COPY_H
#define COMABLAS_COMA_COPY_H

#include <complex.h>

void scopy(unsigned int n, float *x, int incx, float *y, int incy);
void dcopy(unsigned int n, double *x, int incx, double *y, int incy);
void ccopy(unsigned int n,const complex float *x, int incx,complex float *y, int incy);
void zcopy(unsigned int n,const complex double *x, int incx,complex double *y, int incy);


#endif //COMABLAS_COMA_COPY_H

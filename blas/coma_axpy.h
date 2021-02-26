//
// Created by jon on 2/20/21.
//

#ifndef COMABLAS_COMA_AXPY_H
#define COMABLAS_COMA_AXPY_H

void saxpy(unsigned int n, float a, const float* x, int incx, float* y, int incy);
void daxpy(unsigned int n, double a, const double* x, int incx, double* y, int incy);

#endif //COMABLAS_COMA_AXPY_H

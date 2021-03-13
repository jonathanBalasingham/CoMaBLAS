//
// Created by jon on 3/12/21.
//

#ifndef COMABLAS_COMA_GEMV_H
#define COMABLAS_COMA_GEMV_H
#include <complex.h>

void sgemv(char trans, unsigned int m, unsigned int n, float alpha, float** A, int LDA, float* x, unsigned int incx, float beta, float* y, unsigned int incy);
void dgemv(char trans, unsigned int m, unsigned int n, double alpha, double** A, int LDA, double* x, unsigned int incx, double beta, double* y, unsigned int incy);
void cgemv(char trans, unsigned int m, unsigned int n, complex float alpha, complex float** A, int LDA, complex float* x, unsigned int incx, complex float beta, complex float* y, unsigned int incy);
void zgemv(char trans, unsigned int m, unsigned int n, complex double alpha,complex double** A, int LDA,complex double* x, unsigned int incx,complex double beta,complex double* y, unsigned int incy);



#endif //COMABLAS_COMA_GEMV_H

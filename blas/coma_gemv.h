//
// Created by jon on 3/12/21.
//

#ifndef COMABLAS_COMA_GEMV_H
#define COMABLAS_COMA_GEMV_H
#include <complex.h>

void sgemv(char trans, int m, int n, float alpha, float** A, int lda, float* x, int incx, float beta, float* y, int incy);
void dgemv(char trans, int m, int n, double alpha, double** A, int lda, double* x, int incx, double beta, double* y, int incy);
void cgemv(char trans, int m, int n, complex float alpha, complex float** A, int lda,
           complex float* x, int incx, complex float beta, complex float* y, int incy);
void zgemv(char trans, int m, int n, complex double alpha,complex double** A, int lda,
           complex double* x, int incx,complex double beta,complex double* y, int incy);



#endif //COMABLAS_COMA_GEMV_H

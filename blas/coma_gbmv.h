//
// Created by jon on 3/13/21.
//

#ifndef COMABLAS_COMA_GBMV_H
#define COMABLAS_COMA_GBMV_H
#include <complex.h>

void sgbmv(char trans, int m, int n, int kl, int ku, float alpha, float** A,
           int lda, float* x, int incx, float beta, float* y, int incy);
void dgbmv(char trans, int m, int n, int kl, int ku, double alpha, double** A,
           int lda, double* x, int incx, double beta, double* y, int incy);
void cgbmv(char trans, int m, int n, int kl, int ku, complex float alpha, complex float** A,
           int lda, complex float* x, int incx, complex float beta, complex float* y, int incy);
void zgbmv(char trans, int m, int n, int kl, int ku, complex double alpha, complex double** A,
           int lda,complex double* x, int incx,complex double beta,complex double* y, int incy);


#endif //COMABLAS_COMA_GBMV_H

//
// Created by jon on 4/4/21.
//

#ifndef COMABLAS_COMA_GEMM_H
#define COMABLAS_COMA_GEMM_H

#include <complex.h>

void sgemm(char transa, char transb, int m, int n, int k, float alpha, float** A, int lda,
           float** B, int ldb, float beta, float** C, int ldc);
void dgemm(char transa, char transb, int m, int n, int k, double alpha, double** A, int lda,
           double** B, int ldb, double beta, double** C, int ldc);
void cgemm(char transa, char transb, int m, int n, int k,complex float alpha, complex float** A, int lda,
           complex float** B, int ldb, complex float beta, complex float** C, int ldc);
void zgemm(char transa, char transb, int m, int n, int k, complex double alpha, complex double** A, int lda,
           complex double** B, int ldb, double beta, complex double** C, int ldc);



#endif //COMABLAS_COMA_GEMM_H

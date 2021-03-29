//
// Created by jon on 3/28/21.
//

#include <complex.h>

#ifndef COMABLAS_COMA_TBMV_H
#define COMABLAS_COMA_TBMV_H

void stbmv(char uplo, char trans, char diag, int n, int k, float** A, int lda, float* x, int incx);
void dtbmv(char uplo, char trans, char diag, int n, int k, double** A, int lda, double* x, int incx);
void ctbmv(char uplo, char trans, char diag, int n, int k, complex float** A, int lda, complex float* x, int incx);
void ztbmv(char uplo, char trans, char diag, int n, int k, complex double** A, int lda, complex double* x, int incx);

#endif //COMABLAS_COMA_TBMV_H

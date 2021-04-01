//
// Created by jon on 4/1/21.
//

#ifndef COMABLAS_COMA_TRSV_H
#define COMABLAS_COMA_TRSV_H

#include <complex.h>

void strsv(char uplo, char trans, char diag, int n, float** A, int lda, float* x, int incx);
void dtrsv(char uplo, char trans, char diag, int n, double** A, int lda, double* x, int incx);
void ctrsv(char uplo, char trans, char diag, int n, complex float** A, int lda, complex float* x, int incx);
void ztrsv(char uplo, char trans, char diag, int n, complex double** A, int lda, complex double* x, int incx);

#endif //COMABLAS_COMA_TRSV_H

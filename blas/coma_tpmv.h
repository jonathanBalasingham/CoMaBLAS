//
// Created by jon on 3/31/21.
//

#ifndef COMABLAS_COMA_TPMV_H
#define COMABLAS_COMA_TPMV_H

#include <complex.h>

void stpmv(char uplo, char trans, char diag, int n, int k, float** A, float* x, int incx);
void dtpmv(char uplo, char trans, char diag, int n, int k, double** A, double* x, int incx);
void ctpmv(char uplo, char trans, char diag, int n, int k, complex float** A, complex float* x, int incx);
void ztpmv(char uplo, char trans, char diag, int n, int k, complex double** A, complex double* x, int incx);


#endif //COMABLAS_COMA_TPMV_H

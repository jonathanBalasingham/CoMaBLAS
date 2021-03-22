//
// Created by jon on 3/22/21.
//

#ifndef COMABLAS_COMA_SBMV_H
#define COMABLAS_COMA_SBMV_H

void ssbmv(char uplo, int n, int k, float alpha, float** A, int lda, float* x, int incx, float beta, float* y, int incy);
void dsbmv(char uplo, int n, int k, double alpha, double** A, int lda, double* x, int incx, double beta, double* y, int incy);


#endif //COMABLAS_COMA_SBMV_H

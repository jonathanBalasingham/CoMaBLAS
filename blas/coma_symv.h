//
// Created by jon on 3/19/21.
//

#ifndef COMABLAS_COMA_SYMV_H
#define COMABLAS_COMA_SYMV_H

void ssymv(char uplo, int n, float alpha, float** A, int lda, float* x, int incx, float beta, float* y, int incy);
void dsymv(char uplo, int n, double alpha, double** A, int lda, double* x, int incx, double beta, double* y, int incy);

#endif //COMABLAS_COMA_SYMV_H

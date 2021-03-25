#ifndef COMABLAS_COMA_SPMV_H
#define COMABLAS_COMA_SPMV_H

void sspmv(char uplo, int n, float alpha, float* AP, float* x, int incx, float beta, float* y, int incy);
void dspmv(char uplo, int n, double alpha, double* AP, double* x, int incx, double beta, double* y, int incy);


#endif //COMABLAS_COMA_SPMV_H
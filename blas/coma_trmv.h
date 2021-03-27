//
// Created by jon on 3/27/21.
//

#ifndef COMABLAS_COMA_TRMV_H
#define COMABLAS_COMA_TRMV_H

void strmv(char uplo, char trans, char diag, int n, float** A, int lda, float* x, int incx);
void dtrmv(char uplo, char trans, char diag, int n, double** A, int lda, double* x, int incx);


#endif //COMABLAS_COMA_TRMV_H

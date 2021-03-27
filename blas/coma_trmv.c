//
// Created by jon on 3/27/21.
//

#include <ctype.h>
#include <stdbool.h>
#include "coma_trmv.h"

int _validate_trmv_inputs(char uplo, char trans, char diag, int n, int lda, int incx){
    int info = 0;
    if (uplo != 'U' && uplo != 'L')
            info = 1;
    else if (trans != 'N' && trans != 'T' && trans != 'C')
            info = 2;
    else if (diag != 'U' && diag != 'N')
            info = 3;
    else if (n < 0) 
            info = 4;
    else if (lda < (1 > n ? 1: n))
            info = 6;
    else if (incx  == 0)  
            info = 8;

    return info;
}

void strmv(char uplo, char trans, char diag, int n, float **A, int lda, float *x, int incx) {
    uplo = (char) toupper(uplo);
    trans = (char) toupper(trans);
    diag = (char) toupper(diag);

    int info = _validate_trmv_inputs(uplo, trans, diag, n, lda, incx);
    if (info != 0){
        return;
    }

    if (n == 0)
        return;

    bool nounit = diag == 'N';
    float temp;
    if (trans == 'N') {
        if (uplo == 'U') {
            if (incx == 1) {
                for (int j = 0; j < n; ++j) {
                    if (x[j] != 0) {
                        temp = x[j];
                        for (int i = 0; i < j - 1; ++i) {
                            x[i] += temp * A[i][j];
                        }
                        if (nounit)
                            x[j] *= A[j][j];
                    }
                }
            } else {

            }
        } else {
            if (incx == 1) {

            } else {

            }
        }
    } else {
        if (uplo == 'U') {
            if (incx == 1) {

            } else {

            }
        } else {
            if (incx == 1) {

            } else {

            }
        }
    }

}

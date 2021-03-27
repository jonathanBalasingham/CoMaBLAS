//
// Created by jon on 3/27/21.
//

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
    int info = _validate_trmv_inputs(uplo, trans, diag, n, lda, incx);
    if (info != 0){
        return;
    }


}

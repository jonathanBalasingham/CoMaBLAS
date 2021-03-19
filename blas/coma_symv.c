//
// Created by jon on 3/19/21.
//

#include <ctype.h>
#include "coma_symv.h"

int _validate_symv_inputs(char uplo, int n, int lda, int incx, int incy){
    if (uplo != 'U' && uplo != 'L')
        return 1;
    else if (n < 0)
        return 2;
    else if (lda < (1 < n ? n : 1))
        return 5;
    else if (incx == 0)
        return 7;
    else if (incy == 0)
        return 10;

    return 0;
}

void ssymv(char uplo, int n, float alpha, float **A, int lda, int incx, float beta, float *y, int incy) {
    int info = _validate_symv_inputs(uplo, n, lda, incx, incy);
    if (info != 0)
        return;

    if (n == 0 || (alpha == 0 && beta == 1))
        return;

    int kx, ky;
    if (incx > 0)
        kx = 0;
    else
        kx = 1 - (n - 1) * incx;

    if (incy > 0)
        ky = 0;
    else
        ky = 1 - (n - 1) * incy;

    if (beta != 1) {
        if (incy == 1) {
            if (beta == 0) {
                for (int i = 0; i < n; ++i) {
                    y[i] = 0;
                }
            } else {
                for (int i = 0; i < n; ++i) {
                    y[i] = beta * y[i];
                }
            }
        } else {
            int iy = ky;
            if (beta == 0) {
                for (int i = 0; i < n; ++i) {
                    y[iy] = 0;
                    iy += incy;
                }
            } else {
                for (int i = 0; i < n; ++i) {
                    y[iy] = beta * y[iy];
                    iy += incy;
                }
            }
        }
    }

    if (alpha == 0)
        return;

    if (toupper(uplo) == 'U') {

    } else {

    }

}

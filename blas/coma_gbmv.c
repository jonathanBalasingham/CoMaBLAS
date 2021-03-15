//
// Created by jon on 3/13/21.
//

#include <ctype.h>
#include "coma_gbmv.h"

int _validate_gbmv_inputs(char trans, int n, int m, int kl, int ku, int lda, int incx, int incy){
    int info = 0;
    if (trans != 'N' && trans != 'T' && trans != 'C')
        info = 1;
    else if (m < 0)
        info = 2;
    else if (n < 0)
        info = 3;
    else if (kl < 0)
        info = 4;
    else if (ku < 0)
        info = 5;
    else if (lda < kl + ku + 1)
        info = 8;
    else if (lda < 1)
        info = 6;
    else if (incx == 0)
        info = 8;
    else if (incy == 0)
        info = 11;

    return info;
}

void
sgbmv(char trans, int m, int n, int kl, int ku, float alpha, float **A, int lda,
      float *x, int incx, float beta, float *y, int incy) {

    int info = _validate_gbmv_inputs((char)toupper(trans), m,n,kl,ku,lda,incx,incy);

    if (info != 0){
        // error
        return;
    }

    if (m == 0 || n == 0 || (alpha == 0 && beta == 1))
        return;

    int lenx, leny, kx, ky;
    if (trans == 'N') {
        lenx = n;
        leny = m;
    } else {
        lenx = m;
        leny = n;
    }

    if (incx > 0)
        kx = 0;
    else
        kx = 1 - (lenx - 1) * incx;

    if (incy > 0)
        ky = 0;
    else
        ky = 1 - (leny - 1) * incy;

    if (beta != 1) {
        if (incy == 1) {
            if (beta == 0) {
                for (int i = 0; i < leny; ++i) {
                    y[i] = 0;
                }
            } else {
                for (int i = 0; i < leny; ++i) {
                    y[i] = beta * y[i];
                }
            }
        } else {
            int iy = ky;
            if (beta == 0) {
                for (int i = 0; i < leny; ++i) {
                    y[iy] = 0;
                    iy += incy;
                }
            } else {
                for (int i = 0; i < leny; ++i) {
                    y[iy] = beta * y[iy];
                    iy += incy;
                }
            }
        }
    }

    if (alpha == 0)
        return;

    int kup1 = ku + 1;
    float temp;
    if (trans == 'N') {
        int jx = kx;
        if (incy == 1) {
            for (int i = 0; i < n; ++i) {
                temp = alpha * x[jx];
                k = kup1 - i;
                for (int j = max(0, i - ku); j < min(m, i + kl); ++j) {

                }
            }
        } else {

        }
    } else {

    }

}

//
// Created by jon on 3/12/21.
//

#include <ctype.h>
#include "coma_gemv.h"
#include <stdbool.h>

void sgemv(char trans, unsigned int m, unsigned int n, float alpha, float** A, int lda, float *x, unsigned int incx,
           float beta, float *y, unsigned int incy) {

    int info = 0;
    trans = toupper(trans);

    if (trans != 'N' && trans != 'T' && trans != 'C')
        info = 1;
    else if (lda < 1)
        info = 6;
    else if (incx == 0)
        info = 8;
    else if (incy == 0)
        info = 11;

    if (info != 0){
        // throw error
    }

    bool noconj = trans == 'T';

    unsigned int lenx, leny, kx, ky;
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
            unsigned int iy = ky;
            if (beta == 0) {
                for (int i = 0; i < leny; ++i) {
                    y[iy] = 0;
                    iy += incy;
                }
            } else {
                for (int i = 0; i < leny; ++i) {
                    y[iy] = beta * y[i];
                    iy += incy;
                }
            }
        }
    }

    if (alpha == 0)
        return;

    if (trans == 'N') {
        unsigned int jx = kx;
        float temp;
        if (incy == 1) {
            for (int i = 0; i < n; ++i) {
                temp = alpha * x[jx];
                for (int j = 0; j < m; ++j) {
                    y[j] += temp * A[j][i];
                }
                jx += incx;
            }
        } else {
            unsigned int iy;
            for (int i = 0; i < n; ++i) {
                temp = alpha * x[jx];
                iy = ky;
                for (int j = 0; j < m; ++j) {
                    y[iy] += temp * A[j][i];
                    iy += incy;
                }
                jx += incx;
            }
        }
    } else {
        unsigned int jy = ky;
        if (incx == 1) {
            float temp;
            for (int i = 0; i < n; ++i) {
                temp = 0;
                if (noconj) {
                    for (int j = 0; j < m; ++j) {
                        temp += A[j][i] * x[j];
                    }
                } else {
                    for (int j = 0; j < m; ++j) {
                        temp += conjg(A[j][i]) * x[j];
                    }
                }
            }
        } else {

        }
    }
}

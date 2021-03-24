//
// Created by jon on 3/22/21.
//

#include <ctype.h>
#include "coma_sbmv.h"

int _validate_sbmv_inputs(char uplo, int n, int k, int lda, int incx, int incy){
    if (uplo != 'U' && uplo != 'L')
        return 1;
    else if (n < 0)
        return 2;
    else if (k < 0)
        return 3;
    else if (lda < k + 1)
        return 6;
    else if (incx == 0)
        return 8;
    else if (incy == 0)
        return 11;

    return 0;
}

void
ssbmv(char uplo, int n, int k, float alpha, float **A, int lda, float *x, int incx, float beta, float *y, int incy) {
    uplo = (char) toupper(uplo);

    int info = _validate_sbmv_inputs(uplo, n, k, lda, incx, incy);
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

    int l, jx, jy, ix, iy;
    if (uplo == 'U') {
        int kplus1 = k + 1;
        if (incx == 1 && incy == 1) {
            float temp1, temp2;
            for (int i = 0; i < n; ++i) {
                temp1 = alpha * x[i];
                temp2 = 0;
                l = kplus1 - i;
                for (int j = (1 > i-k ? 1 : i-k); j < i - 1; ++j) {
                    y[j] += temp1 * A[l+j][i];
                    temp2 += A[l+j][i] * x[j];
                }
                y[i] += temp1 * A[kplus1][i] + alpha * temp2;
            }
        } else {
            jx = kx;
            jy = ky;
            float temp1, temp2;
            for (int i = 0; i < n; ++i) {
                temp1 = alpha * x[i];
                temp2 = 0;
                ix = kx;
                iy = ky;
                l = kplus1 - i;
                for (int j = (1 > i-k ? 1 : i-k); j < i - 1; ++j) {
                    y[iy] += temp1 * A[l+j][i];
                    temp2 += A[l+j][i] * x[ix];
                    ix += incx;
                    ix += incy;
                }
                y[jy] += temp1 * A[kplus1][i] + alpha * temp2;
                jy += incy;
                jx += incx;
                if (i > k) {
                    kx += incx;
                    ky += incy;
                }
            }
        }
    } else {
        int kplus1 = k + 1;
        if (incx == 1 && incy == 1) {
            float temp1, temp2;
            for (int i = 0; i < n; ++i) {
                temp1 = alpha * x[i];
                temp2 = 0;
                l = 1 - i;
                for (int j = i + 1; j < (n < i+k ? n : i+k); ++j) {
                    y[j] += temp1 * A[l+j][i];
                    temp2 += A[l+j][i] * x[j];
                }
                y[i] += alpha * temp2;
            }
        } else {
            jx = kx;
            jy = ky;
            float temp1, temp2;
            for (int i = 0; i < n; ++i) {
                temp1 = alpha * x[i];
                temp2 = 0;
                ix = kx;
                iy = ky;
                l = 1 - i;
                for (int j = i + 1; j < (n < i+k ? n : i+k); ++j) {
                    ix += incx;
                    ix += incy;
                    y[iy] += temp1 * A[l+j][i];
                    temp2 += A[l+j][i] * x[ix];
                }
                y[jy] += alpha * temp2;
                jy += incy;
                jx += incx;
            }
        }
    }
}

void dsbmv(char uplo, int n, int k, double alpha, double **A, int lda, double *x, int incx, double beta, double *y,
           int incy) {
    uplo = (char) toupper(uplo);

    int info = _validate_sbmv_inputs(uplo, n, k, lda, incx, incy);
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

    int l, jx, jy, ix, iy;
    if (uplo == 'U') {
        int kplus1 = k + 1;
        if (incx == 1 && incy == 1) {
            double temp1, temp2;
            for (int i = 0; i < n; ++i) {
                temp1 = alpha * x[i];
                temp2 = 0;
                l = kplus1 - i;
                for (int j = (1 > i-k ? 1 : i-k); j < i - 1; ++j) {
                    y[j] += temp1 * A[l+j][i];
                    temp2 += A[l+j][i] * x[j];
                }
                y[i] += temp1 * A[kplus1][i] + alpha * temp2;
            }
        } else {
            jx = kx;
            jy = ky;
            double temp1, temp2;
            for (int i = 0; i < n; ++i) {
                temp1 = alpha * x[i];
                temp2 = 0;
                ix = kx;
                iy = ky;
                l = kplus1 - i;
                for (int j = (1 > i-k ? 1 : i-k); j < i - 1; ++j) {
                    y[iy] += temp1 * A[l+j][i];
                    temp2 += A[l+j][i] * x[ix];
                    ix += incx;
                    ix += incy;
                }
                y[jy] += temp1 * A[kplus1][i] + alpha * temp2;
                jy += incy;
                jx += incx;
                if (i > k) {
                    kx += incx;
                    ky += incy;
                }
            }
        }
    } else {
        if (incx == 1 && incy == 1) {
            float temp1, temp2;
            for (int i = 0; i < n; ++i) {
                temp1 = alpha * x[i];
                temp2 = 0;
                l = 1 - i;
                for (int j = i + 1; j < (n < i+k ? n : i+k); ++j) {
                    y[j] += temp1 * A[l+j][i];
                    temp2 += A[l+j][i] * x[j];
                }
                y[i] += alpha * temp2;
            }
        } else {
            jx = kx;
            jy = ky;
            float temp1, temp2;
            for (int i = 0; i < n; ++i) {
                temp1 = alpha * x[i];
                temp2 = 0;
                ix = kx;
                iy = ky;
                l = 1 - i;
                for (int j = i + 1; j < (n < i+k ? n : i+k); ++j) {
                    ix += incx;
                    ix += incy;
                    y[iy] += temp1 * A[l+j][i];
                    temp2 += A[l+j][i] * x[ix];
                }
                y[jy] += alpha * temp2;
                jy += incy;
                jx += incx;
            }
        }
    }
}

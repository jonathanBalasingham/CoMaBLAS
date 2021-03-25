
#include <ctype.h>
#include "coma_spmv.h"

int _validate_spmv_inputs(char uplo, int n, int incx, int incy) {
    if (uplo != 'U' && uplo != 'L')
        return 1;
    else if (n < 0)
        return 2;
    else if (incx == 0)
        return 6;
    else if (incy == 0)
        return 9;

    return 0;
}

void sspmv(char uplo, int n, float alpha, float *AP, float *x, int incx, float beta, float *y, int incy) {
    uplo = (char) toupper(uplo);
    int info = _validate_spmv_inputs(uplo, n, incx, incy);

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

    int k,kk = 0;
    float temp1, temp2;
    int ix,iy;
    if (uplo == 'U') {
        if (incx == 1 && incy == 1) {
            for (int i = 0; i < n; ++i) {
                temp1 = alpha * x[i];
                temp2 = 0;
                k = kk;
                for (int j = 0; j < i - 1; ++j) {
                    y[j] += temp1 * AP[k];
                    temp2 += AP[k] * x[j];
                    k += 1;
                }
                y[i] += temp1 * AP[kk + i - 1] + alpha * temp2;
                kk += i;
            }
        } else {
            int jx = kx;
            int jy = ky;
            for (int i = 0; i < n; ++i) {
                temp1 = alpha * x[i];
                temp2 = 0;
                ix = kx;
                iy = ky;
                for (int k = kk; k < kk + i -2; ++k) {
                    y[iy] += temp1 * AP[k];
                    temp2 += AP[k] * x[ix];
                    ix += incx;
                    iy += incy;
                }
                y[jy] += temp1 * AP[kk+i-1] + alpha * temp2;
                jx += incx;
                jy += incy;
                kk += i;
            }
        }
    } else {
        if (incx == 1 && incy == 1) {
            for (int i = 0; i < n; ++i) {
                temp1 = alpha * x[i];
                temp2 = 0;
                y[i] += temp1 * AP[k];
                k = kk + 1;
                for (int j = i + 1; j < n; ++j) {
                    y[j] += temp1 * AP[k];
                    temp2 += AP[k] * x[j];
                    k += 1;
                }
                y[i] += alpha * temp2;
                kk += n - i + 1;
            }
        } else {
            int jx = kx;
            int jy = ky;
            for (int i = 0; i < n; ++i) {
                temp1 = alpha * x[i];
                temp2 = 0;
                y[jy] += temp1 * AP[kk];
                ix = kx;
                iy = ky;
                for (int k = kk + 1; k < kk + i - i; ++k) {
                    ix += incx;
                    iy += incy;
                    y[iy] += temp1 * AP[k];
                    temp2 += AP[k] * x[ix];
                }
                y[jy] += temp1 * AP[kk+i-1] + alpha * temp2;
                jx += incx;
                jy += incy;
                kk += n - i + 1;
            }
        }
    }

}

void dspmv(char uplo, int n, double alpha, double *AP, double *x, int incx, double beta, double *y, int incy) {
    uplo = (char) toupper(uplo);
    int info = _validate_spmv_inputs(uplo, n, incx, incy);

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

    int k,kk = 0;
    double temp1, temp2;
    int ix,iy;
    if (uplo == 'U') {
        if (incx == 1 && incy == 1) {
            for (int i = 0; i < n; ++i) {
                temp1 = alpha * x[i];
                temp2 = 0;
                k = kk;
                for (int j = 0; j < i - 1; ++j) {
                    y[j] += temp1 * AP[k];
                    temp2 += AP[k] * x[j];
                    k += 1;
                }
                y[i] += temp1 * AP[kk + i - 1] + alpha * temp2;
                kk += i;
            }
        } else {
            int jx = kx;
            int jy = ky;
            for (int i = 0; i < n; ++i) {
                temp1 = alpha * x[i];
                temp2 = 0;
                ix = kx;
                iy = ky;
                for (int k = kk; k < kk + i -2; ++k) {
                    y[iy] += temp1 * AP[k];
                    temp2 += AP[k] * x[ix];
                    ix += incx;
                    iy += incy;
                }
                y[jy] += temp1 * AP[kk+i-1] + alpha * temp2;
                jx += incx;
                jy += incy;
                kk += i;
            }
        }
    } else {
        if (incx == 1 && incy == 1) {
            for (int i = 0; i < n; ++i) {
                temp1 = alpha * x[i];
                temp2 = 0;
                y[i] += temp1 * AP[k];
                k = kk + 1;
                for (int j = i + 1; j < n; ++j) {
                    y[j] += temp1 * AP[k];
                    temp2 += AP[k] * x[j];
                    k += 1;
                }
                y[i] += alpha * temp2;
                kk += n - i + 1;
            }
        } else {
            int jx = kx;
            int jy = ky;
            for (int i = 0; i < n; ++i) {
                temp1 = alpha * x[i];
                temp2 = 0;
                y[jy] += temp1 * AP[kk];
                ix = kx;
                iy = ky;
                for (int k = kk + 1; k < kk + i - i; ++k) {
                    ix += incx;
                    iy += incy;
                    y[iy] += temp1 * AP[k];
                    temp2 += AP[k] * x[ix];
                }
                y[jy] += temp1 * AP[kk+i-1] + alpha * temp2;
                jx += incx;
                jy += incy;
                kk += n - i + 1;
            }
        }
    }
}

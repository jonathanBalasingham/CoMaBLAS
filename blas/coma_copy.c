//
// Created by jon on 2/20/21.
//

#include "coma_copy.h"

void scopy(unsigned int n, float *x, int incx, float *y, int incy) {
    int m = n % 7;
    if (n == 0)
        return;

    if (incx == 1 && incy == 1) {
        if (m != 0){
            for (int i = 0; i < m; ++i) {
                y[i] = x[i];
            }
        }

        for (int i = m; i < n; i += 7) {
            y[i] = x[i];
            y[i+1] = x[i+1];
            y[i+2] = x[i+2];
            y[i+3] = x[i+3];
            y[i+4] = x[i+4];
            y[i+5] = x[i+5];
            y[i+6] = x[i+6];
        }
    } else {
        int ix = 1, iy = 1;
        if (incx < 0)
            ix = (-1*n+1) * incx;
        if (incx < 0)
            iy = (-1*n+1) * incy;

        for (int i = 0; i < n; ++i) {
            y[iy] = x[ix];
            ix += incx;
            iy += incy;
        }
    }
}

void dcopy(unsigned int n, double *x, int incx, double *y, int incy) {
    int m = n % 7;
    if (n == 0)
        return;

    if (incx == 1 && incy == 1) {
        if (m != 0){
            for (int i = 0; i < m; ++i) {
                y[i] = x[i];
            }
        }

        for (int i = m; i < n; i += 7) {
            y[i] = x[i];
            y[i+1] = x[i+1];
            y[i+2] = x[i+2];
            y[i+3] = x[i+3];
            y[i+4] = x[i+4];
            y[i+5] = x[i+5];
            y[i+6] = x[i+6];
        }
    } else {
        int ix = 1, iy = 1;
        if (incx < 0)
            ix = (-1*n+1) * incx;
        if (incx < 0)
            iy = (-1*n+1) * incy;

        for (int i = 0; i < n; ++i) {
            y[iy] = x[ix];
            ix += incx;
            iy += incy;
        }
    }
}

void ccopy(unsigned int n, const complex float *x, int incx, complex float *y, int incy) {
    if (n == 0)
        return;

    if (incx == 1 && incy == 1) {
        for (int i = 0; i < n; ++i) {
            y[i] = x[i];
        }
    } else {
        int ix = 1, iy = 1;
        if (incx < 0)
            ix = (-1*n+1) * incx;
        if (incx < 0)
            iy = (-1*n+1) * incy;

        for (int i = 0; i < n; ++i) {
            y[iy] = x[ix];
            ix += incx;
            iy += incy;
        }
    }
}

void zcopy(unsigned int n, const complex double *x, int incx, complex double *y, int incy) {
    if (n == 0)
        return;

    if (incx == 1 && incy == 1) {
        for (int i = 0; i < n; ++i) {
            y[i] = x[i];
        }
    } else {
        int ix = 1, iy = 1;
        if (incx < 0)
            ix = (-1*n+1) * incx;
        if (incx < 0)
            iy = (-1*n+1) * incy;

        for (int i = 0; i < n; ++i) {
            y[iy] = x[ix];
            ix += incx;
            iy += incy;
        }
    }
}

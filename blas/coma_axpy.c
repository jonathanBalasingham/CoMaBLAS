//
// Created by jon on 2/20/21.
//

#include "coma_axpy.h"
#include <complex.h>

void saxpy(unsigned int n, float a, const float *x, int incx, float *y, int incy) {
    if (n == 0 || a == 0)
        return;

    if (incx == 1 && incy == 1) {
        int m = n % 4;
        if (m != 0) {
            for (int i = 0; i < m; ++i) {
                y[i] = y[i] + a * x[i];
            }
        }
        if (n < 4)
            return;

        for (int i = m; i < n; i += 4) {
            y[i] = y[i] + a * x[i];
            y[i+1] = y[i+1] + a * x[i+1];
            y[i+2] = y[i+2] + a * x[i+2];
            y[i+3] = y[i+3] + a * x[i+3];
        }
    } else {
        int ix = 0, iy = 0;
        if (incx < 0)
            ix = (-1*n+1) * incx + 1;
        if (incy < 0)
            iy = (-1*n+1) * incy + 1;

        for (int i = 0; i < n; ++i) {
            y[iy] = y[iy] + a * x[ix];
            ix += incx;
            iy += incy;
        }
    }
}

void daxpy(unsigned int n, double a, const double *x, int incx, double *y, int incy) {
    if (n == 0 || a == 0)
        return;

    if (incx == 1 && incy == 1) {
        int m = n % 4;
        if (m != 0) {
            for (int i = 0; i < m; ++i) {
                y[i] = y[i] + a * x[i];
            }
        }
        if (n < 4)
            return;

        for (int i = m; i < n; i += 4) {
            y[i] = y[i] + a * x[i];
            y[i+1] = y[i+1] + a * x[i+1];
            y[i+2] = y[i+2] + a * x[i+2];
            y[i+3] = y[i+3] + a * x[i+3];
        }
    } else {
        int ix = 0, iy = 0;
        if (incx < 0)
            ix = (-1*n+1) * incx + 1;
        if (incy < 0)
            iy = (-1*n+1) * incy + 1;

        for (int i = 0; i < n; ++i) {
            y[iy] = y[iy] + a * x[ix];
            ix += incx;
            iy += incy;
        }
    }
}

void caxpy(unsigned int n, float _Complex a, const float _Complex *x, int incx, float _Complex *y, int incy) {
    if (n == 0 || cabsf(a) == 0)
        return;

    if (incx == 1 && incy == 1) {
        for (int i = 0; i < n; i += 4) {
            y[i] = y[i] + a * x[i];
        }
    } else {
        int ix = 0, iy = 0;
        if (incx < 0)
            ix = (-1*n+1) * incx + 1;
        if (incy < 0)
            iy = (-1*n+1) * incy + 1;

        for (int i = 0; i < n; ++i) {
            y[iy] = y[iy] + a * x[ix];
            ix += incx;
            iy += incy;
        }
    }
}

void zaxpy(unsigned int n, double _Complex a, const double _Complex *x, int incx, double _Complex *y, int incy) {
    if (n == 0 || cabs(a) == 0)
        return;

    if (incx == 1 && incy == 1) {
        for (int i = 0; i < n; i += 4) {
            y[i] = y[i] + a * x[i];
        }
    } else {
        int ix = 0, iy = 0;
        if (incx < 0)
            ix = (-1*n+1) * incx + 1;
        if (incy < 0)
            iy = (-1*n+1) * incy + 1;

        for (int i = 0; i < n; ++i) {
            y[iy] = y[iy] + a * x[ix];
            ix += incx;
            iy += incy;
        }
    }
}

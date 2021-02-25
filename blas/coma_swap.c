//
// Created by jon on 2/20/21.
//

#include "coma_swap.h"


void sswap(unsigned int n, float *x, int incx, float *y, int incy) {
    float temp;

    if (n == 0)
        return;

    if (incx == 1 && incy == 1) {
        int m = n % 3;
        if (m != 0) {
            for (int i = 0; i < m; ++i) {
                temp = x[i];
                x[i] = y[i];
                y[i] = temp;
            }
            if (n < 3)
                return;

        }
        for (int i = m; i < n; i += 3) {
            temp = x[i];
            x[i] = y[i];
            y[i] = temp;
            temp = x[i+1];
            x[i+1] = y[i+1];
            y[i+1] = temp;
            temp = x[i+2];
            x[i+2] = y[i+2];
            y[i+2] = temp;
        }
    } else {
        int ix = 1, iy = 1;
        if (incx < 0)
            ix = (-1*n+1) * incx + 1;
        if (incy < 0)
            iy = (-1*n+1) * incy + 1;

        for (int i = 0; i < n; ++i) {
            temp = x[ix];
            x[ix] = y[iy];
            y[iy] = temp;
            ix += incx;
            iy += incy;
        }
    }
}

void dswap(unsigned int n, double *x, int incx, double *y, int incy) {
    double temp;

    if (n == 0)
        return;

    if (incx == 1 && incy == 1) {
        int m = n % 3;
        if (m != 0) {
            for (int i = 0; i < m; ++i) {
                temp = x[i];
                x[i] = y[i];
                y[i] = temp;
            }
            if (n < 3)
                return;

        }
        for (int i = m; i < n; i += 3) {
            temp = x[i];
            x[i] = y[i];
            y[i] = temp;
            temp = x[i+1];
            x[i+1] = y[i+1];
            y[i+1] = temp;
            temp = x[i+2];
            x[i+2] = y[i+2];
            y[i+2] = temp;
        }
    } else {
        int ix = 1, iy = 1;
        if (incx < 0)
            ix = (-1*n+1) * incx + 1;
        if (incy < 0)
            iy = (-1*n+1) * incy + 1;

        for (int i = 0; i < n; ++i) {
            temp = x[ix];
            x[ix] = y[iy];
            y[iy] = temp;
            ix += incx;
            iy += incy;
        }
    }
}

void cswap(unsigned int n, complex float *x, int incx, complex float *y, int incy) {
    complex float temp;

    if (n == 0)
        return;

    if (incx == 1 && incy == 1) {
        for (int i = 0; i < n; ++i) {
            temp = x[i];
            x[i] = y[i];
            y[i] = temp;
        }
    } else {
        int ix = 0, iy = 0;
        if (incx < 0)
            ix = (-1*n+1) * incx + 1;
        if (incy < 0)
            iy = (-1*n+1) * incy + 1;

        for (int i = 0; i < n; ++i) {
            temp = x[ix];
            x[ix] = y[iy];
            y[iy] = temp;
            ix += incx;
            iy += incy;
        }
    }
}

void zswap(unsigned int n, complex double *x, int incx, complex double *y, int incy) {
    complex double temp;

    if (n == 0)
        return;

    if (incx == 1 && incy == 1) {
        for (int i = 0; i < n; ++i) {
            temp = x[i];
            x[i] = y[i];
            y[i] = temp;
        }
    } else {
        int ix = 1, iy = 1;
        if (incx < 0)
            ix = (-1*n+1) * incx + 1;
        if (incy < 0)
            iy = (-1*n+1) * incy + 1;

        for (int i = 0; i < n; ++i) {
            temp = x[ix];
            x[ix] = y[iy];
            y[iy] = temp;
            ix += incx;
            iy += incy;
        }
    }
}

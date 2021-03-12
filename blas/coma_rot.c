//
// Created by jon on 2/20/21.
//

#include "coma_rot.h"

void srot(unsigned int n, float *x, int incx, float *y, int incy, float c, float s) {
    float temp;

    if (n == 0)
        return;

    if (incx == 1 && incy == 1) {
        for (int i = 0; i < n; ++i) {
            temp = c*x[i] + s * y[i];
            y[i] = c * y[i] - s * x[i];
            x[i] = temp;
        }
    } else {
        int ix = 0, iy = 0;
        if (incx < 0)
            ix = (-1*n+1) * incx + 1;
        if (incy < 0)
            iy = (-1*n+1) * incy + 1;

        for (int i = 0; i < n; ++i) {
            temp = c * x[ix] + s * y[iy];
            y[iy] = c * y[iy] - s * x[ix];
            x[ix] = temp;
            ix += incx;
            iy += incy;
        }

    }
}

void drot(unsigned int n, double *x, int incx, double *y, int incy, double c, double s) {
    double temp;

    if (n == 0)
        return;

    if (incx == 1 && incy == 1) {
        for (int i = 0; i < n; ++i) {
            temp = c*x[i] + s * y[i];
            y[i] = c * y[i] - s * x[i];
            x[i] = temp;
        }
    } else {
        int ix = 0, iy = 0;
        if (incx < 0)
            ix = (-1*n+1) * incx + 1;
        if (incy < 0)
            iy = (-1*n+1) * incy + 1;

        for (int i = 0; i < n; ++i) {
            temp = c * x[ix] + s * y[iy];
            y[iy] = c * y[iy] - s * x[ix];
            x[ix] = temp;
            ix += incx;
            iy += incy;
        }
    }
}

void csrot(unsigned int n, complex float *x, int incx, complex float *y, int incy, float c, float s) {
    complex float temp;

    if (n == 0)
        return;

    if (incx == 1 && incy == 1) {
        for (int i = 0; i < n; ++i) {
            temp = c*x[i] + s * y[i];
            y[i] = c * y[i] - s * x[i];
            x[i] = temp;
        }
    } else {
        int ix = 0, iy = 0;
        if (incx < 0)
            ix = (-1*n+1) * incx + 1;
        if (incy < 0)
            iy = (-1*n+1) * incy + 1;

        for (int i = 0; i < n; ++i) {
            temp = c * x[ix] + s * y[iy];
            y[iy] = c * y[iy] - s * x[ix];
            x[ix] = temp;
            ix += incx;
            iy += incy;
        }
    }
}

void zdrotf(unsigned int n, complex double *x, int incx, complex double *y, int incy, double c, double s) {
    complex double temp;

    if (n == 0)
        return;

    if (incx == 1 && incy == 1) {
        for (int i = 0; i < n; ++i) {
            temp = c*x[i] + s * y[i];
            y[i] = c * y[i] - s * x[i];
            x[i] = temp;
        }
    } else {
        int ix = 0, iy = 0;
        if (incx < 0)
            ix = (-1*n+1) * incx + 1;
        if (incy < 0)
            iy = (-1*n+1) * incy + 1;

        for (int i = 0; i < n; ++i) {
            temp = c * x[ix] + s * y[iy];
            y[iy] = c * y[iy] - s * x[ix];
            x[ix] = temp;
            ix += incx;
            iy += incy;
        }
    }
}

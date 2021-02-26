//
// Created by jon on 2/25/21.
//

#include "coma_dsdot.h"

double sdsdot(unsigned int n, float b, float *x, int incx, float *y, int incy) {
    double dsdot = b;
    if (n == 0)
        return dsdot;

    if (incx == incy && incx > 0) {
        int ns = n * incx;
        for (int i = 0; i < ns; i+=incx) {
            dsdot += (double) x[i] * (double) y[i];
        }
    } else {
        int ix = 0, iy = 0;
        if (incx < 0)
            ix = (-1*n+1) * incx;
        if (incy < 0)
            iy = (-1*n+1) * incy;

        for (int j = 0; j < n; ++j) {
            dsdot += (double) x[ix] * (double) y[iy];
            ix += incx;
            iy += incy;
        }
    }
    return dsdot;
}

double dsdot(unsigned int n, double *x, int incx, double *y, int incy) {
    double dsdot = 0;
    if (n == 0)
        return dsdot;

    if (incx == incy && incx > 0) {
        int ns = n * incx;
        for (int i = 0; i < ns; i+=incx) {
            dsdot += (double) x[i] * (double) y[i];
        }
    } else {
        int ix = 0, iy = 0;
        if (incx < 0)
            ix = (-1*n+1) * incx;
        if (incy < 0)
            iy = (-1*n+1) * incy;

        for (int j = 0; j < n; ++j) {
            dsdot += (double) x[ix] * (double) y[iy];
            ix += incx;
            iy += incy;
        }
    }
    return dsdot;}

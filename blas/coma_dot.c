//
// Created by jon on 2/20/21.
//

#include "coma_dot.h"
#include <stdio.h>

float sdot(unsigned int n, float *sx, int incx, float *sy, int incy) {
    float stemp = 0.0, sdot = 0.0;
    int i,ix,iy,m,mp1;

    if (n == 0)
        return sdot;

    if (incx == 1 && incy == 1) {
        m = n % 5;
        if (m != 0){
            for (int j = 0; j < m; ++j) {
                stemp += sx[j] * sy[j];
            }

            if (n < 5){
                sdot = stemp;
                return sdot;
            }
        }
        //mp1 = m + 1; unecessary since indexed at 0
        for (int j = m; j < n; j += 5) {
            stemp += sx[j]*sy[j] + sx[j+1]*sy[j+1] + sx[j+2]*sy[j+2] + sx[j+3]*sy[j+3] + sx[j+4]*sy[j+4];
        }

    } else {
        ix = 1, iy = 1;
        if (incx < 0)
            ix = (-1*n+1) * incx;
        if (incy < 0)
            iy = (-1*n+1) * incy;

        for (int j = 0; j < n-1; ++j) {
            stemp += sx[ix] * sy[iy];
            ix += incx;
            iy += incy;
        }

    }

    sdot = stemp;
    return sdot;
}

double ddot(unsigned int n, const double *dx, int incx, const double *dy, int incy) {
    double temp = 0, ddot = 0;
    int i,ix,iy,m,mp1;

    if (n == 0)
        return ddot;

    if (incx == 1 && incy == 1) {
        m = n % 5;
        if (m != 0){
            for (int j = 0; j < m; ++j) {
                temp += dx[j] * dy[j];
            }

            if (n < 5){
                ddot = temp;
                return ddot;
            }
        }
        //mp1 = m + 1; unecessary since indexed at 0
        for (int j = m; j < n; j += 5) {
            temp += dx[j]*dy[j] + dx[j+1]*dy[j+1] + dx[j+2]*dy[j+2] + dx[j+3]*dy[j+3] + dx[j+4]*dy[j+4];
        }

    } else {
        ix = 1, iy = 1;
        if (incx < 0)
            ix = (-1*n+1) * incx;
        if (incy < 0)
            iy = (-1*n+1) * incy;

        for (int j = 0; j < n-1; ++j) {
            temp += dx[ix] * dy[iy];
            ix += incx;
            iy += incy;
        }

    }

    ddot = temp;
    return ddot;
}

complex float cdotu(unsigned int n, complex float *cx, int incx, complex float *cy, int incy) {
    complex float temp, dotu;
    int ix, iy;

    temp = 0.0 + 0.0i, dotu = 0.0 + 0.0i;
    if (n == 0)
        return dotu;

    if (incx == 1 && incy == 1){
        for (int i = 0; i < n; ++i) {
            temp += cx[i] * cy[i];
        }
    } else {
        ix = 1, iy = 1;
        if (incx < 0)
            ix = (-1*n+1) * incx;
        if (incy < 0)
            iy = (-1*n+1) * incy;

        for (int j = 0; j < n; ++j) {
            temp += cx[ix] * cy[iy];
            ix += incx;
            iy += incy;
        }
    }

    dotu = temp;
    return dotu;
}

complex double zdotu(unsigned int n, complex double *zx, int incx, complex double *zy, int incy) {
    complex double temp, dotu;
    int ix, iy;

    temp = 0.0 + 0.0i, dotu = 0.0 + 0.0i;
    if (n == 0)
        return dotu;

    if (incx == 1 && incy == 1){
        for (int i = 0; i < n; ++i) {
            temp += zx[i] * zy[i];
        }
    } else {
        ix = 1, iy = 1;
        if (incx < 0)
            ix = (-1*n+1) * incx;
        if (incy < 0)
            iy = (-1*n+1) * incy;

        for (int j = 0; j < n; ++j) {
            temp += zx[ix] * zy[iy];
            ix += incx;
            iy += incy;
        }
    }

    dotu = temp;
    return dotu;}

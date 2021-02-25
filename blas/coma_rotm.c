//
// Created by jon on 2/20/21.
//

#include "coma_rotm.h"

void srotm(int n, float *x, int incx, float *y, int incy, float *param) {
    float flag,sh11,sh12,sh21,sh22,w,z;

    flag = param[0];
    if (n == 0 || flag + 2 == 0)
        return;

    if (incx == incy && incx > 0) {
        int nsteps = n * incx;
        if (flag < 0) {
            sh11 = param[1];
            sh12 = param[2];
            sh21 = param[3];
            sh22 = param[4];
            for (int i = 0; i < nsteps; i += incx) {
                w = x[i];
                z = y[i];
                x[i] = w * sh11 + z * sh12;
                y[i] = w * sh21 + z * sh22;
            }
        } else if (flag == 0) {
            sh12 = param[3];
            sh21 = param[2];
            for (int i = 0; i < nsteps; i += incx) {
                w = x[i];
                z = y[i];
                x[i] = w + z * sh12;
                y[i] = w * sh21 + z;
            }
        } else {
            sh11 = param[1];
            sh22 = param[4];
            for (int i = 0; i < nsteps; i += incx) {
                w = x[i];
                z = y[i];
                x[i] = w * sh11 + z;
                y[i] = -1* w + sh22 * z;
            }
        }
    } else {
        int kx = 0, ky = 0;
        if (incx < 0)
            kx = 1 + (1-n) * incx;
        if (incy < 0)
            ky = 1 + (1-n) * incy;

        if (flag < 0) {
            sh11 = param[1];
            sh12 = param[2];
            sh21 = param[3];
            sh22 = param[4];
            for (int i = 0; i < n; ++i) {
                w = x[kx];
                z = y[ky];
                x[kx] = w * sh11 + z * sh12;
                y[ky] = w * sh21 + z * sh22;
                kx += incx;
                ky += incy;
            }
        } else if (flag == 0) {
            sh12 = param[3];
            sh21 = param[2];
            for (int i = 0; i < n; ++i) {
                w = x[kx];
                z = y[ky];
                x[kx] = w + z * sh12;
                y[ky] = w * sh21 + z;
                kx += incx;
                ky += incy;
            }
        } else {
            sh11 = param[1];
            sh22 = param[4];
            for (int i = 0; i < n; ++i) {
                w = x[kx];
                z = y[ky];
                x[kx] = w * sh11 + z;
                y[ky] = -w + sh22 * z;
                kx += incx;
                ky += incy;
            }
        }
    }
}

void drotm(int n, double *x, int incx, double *y, int incy, double *param) {
    double flag,sh11,sh12,sh21,sh22,w,z;

    flag = param[0];
    if (n == 0 || flag + 2 == 0)
        return;

    if (incx == incy && incx > 0) {
        int nsteps = n * incx;
        if (flag < 0) {
            sh11 = param[1];
            sh12 = param[2];
            sh21 = param[3];
            sh22 = param[4];
            for (int i = 0; i < nsteps; i += incx) {
                w = x[i];
                z = y[i];
                x[i] = w * sh11 + z * sh12;
                y[i] = w * sh21 + z * sh22;
            }
        } else if (flag == 0) {
            sh12 = param[3];
            sh21 = param[2];
            for (int i = 0; i < nsteps; i += incx) {
                w = x[i];
                z = y[i];
                x[i] = w + z * sh12;
                y[i] = w * sh21 + z;
            }
        } else {
            sh11 = param[1];
            sh22 = param[4];
            for (int i = 0; i < nsteps; i += incx) {
                w = x[i];
                z = y[i];
                x[i] = w * sh11 + z;
                y[i] = -1* w + sh22 * z;
            }
        }
    } else {
        int kx = 0, ky = 0;
        if (incx < 0)
            kx = 1 + (1-n) * incx;
        if (incy < 0)
            ky = 1 + (1-n) * incy;

        if (flag < 0) {
            sh11 = param[1];
            sh12 = param[2];
            sh21 = param[3];
            sh22 = param[4];
            for (int i = 0; i < n; ++i) {
                w = x[kx];
                z = y[ky];
                x[kx] = w * sh11 + z * sh12;
                y[ky] = w * sh21 + z * sh22;
                kx += incx;
                ky += incy;
            }
        } else if (flag == 0) {
            sh12 = param[3];
            sh21 = param[2];
            for (int i = 0; i < n; ++i) {
                w = x[kx];
                z = y[ky];
                x[kx] = w + z * sh12;
                y[ky] = w * sh21 + z;
                kx += incx;
                ky += incy;
            }
        } else {
            sh11 = param[1];
            sh22 = param[4];
            for (int i = 0; i < n; ++i) {
                w = x[kx];
                z = y[ky];
                x[kx] = w * sh11 + z;
                y[ky] = -w + sh22 * z;
                kx += incx;
                ky += incy;
            }
        }
    }
}


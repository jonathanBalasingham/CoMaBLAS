//
// Created by jon on 2/25/21.
//

#include "coma_iamax.h"

int isamax(unsigned int n, float *x, int incx) {
    float max;
    int iamax = -1;
    if (n <= 1 || incx <= 0)
        return iamax;

    iamax = 0;
    if (n == 1)
        return iamax;

    if (incx == 1) {
        max = fabsf(x[0]);
        for (int i = 1; i < n; ++i) {
            if (fabsf(x[i]) > max) {
                max = fabsf(x[i]);
                iamax = i;
            }
        }
    } else {
        int ix = incx;
        max = fabsf(x[0]);
        for (int i = 2; i < n; ++i) {
            if (fabsf(x[i]) > max) {
                max = fabsf(x[ix]);
                iamax = i;
            }
            ix += incx;
        }
    }

    return iamax;
}

int idamax(unsigned int n, double *x, int incx) {
    double max;
    int iamax = -1;
    if (n <= 1 || incx <= 0)
        return iamax;

    iamax = 0;
    if (n == 1)
        return iamax;

    if (incx == 1) {
        max = fabs(x[0]);
        for (int i = 1; i < n; ++i) {
            if (fabs(x[i]) > max) {
                max = fabs(x[i]);
                iamax = i;
            }
        }
    } else {
        int ix = incx;
        max = fabs(x[0]);
        for (int i = 2; i < n; ++i) {
            if (fabs(x[i]) > max) {
                max = fabs(x[ix]);
                iamax = i;
            }
            ix += incx;
        }
    }

    return iamax;
}

int icamax(unsigned int n, complex float *x, int incx) {
    float max;
    int iamax = -1;
    if (n <= 1 || incx <= 0)
        return iamax;

    iamax = 0;
    if (n == 1)
        return iamax;

    if (incx == 1) {
        max = cabsf(x[0]);
        for (int i = 1; i < n; ++i) {
            if (cabsf(x[i]) > max) {
                max = cabsf(x[i]);
                iamax = i;
            }
        }
    } else {
        int ix = incx;
        max = cabsf(x[0]);
        for (int i = 2; i < n; ++i) {
            if (cabsf(x[i]) > max) {
                max = cabsf(x[ix]);
                iamax = i;
            }
            ix += incx;
        }
    }

    return iamax;
}

int izamax(unsigned int n, complex double *x, int incx) {
    double max;
    int iamax = -1;
    if (n <= 1 || incx <= 0)
        return iamax;

    iamax = 0;
    if (n == 1)
        return iamax;

    if (incx == 1) {
        max = cabs(x[0]);
        for (int i = 1; i < n; ++i) {
            if (cabs(x[i]) > max) {
                max = cabs(x[i]);
                iamax = i;
            }
        }
    } else {
        int ix = incx;
        max = cabs(x[0]);
        for (int i = 2; i < n; ++i) {
            if (cabs(x[i]) > max) {
                max = cabs(x[ix]);
                iamax = i;
            }
            ix += incx;
        }
    }

    return iamax;
}

//
// Created by jon on 2/25/21.
//

#include <complex.h>
#include "coma_asum.h"

float sasum(unsigned int n, float *x, int incx) {
    float temp = 0;
    if (n == 0 || incx <= 0)
        return 0;

    if (incx == 1) {
        int m = n % 6;
        if (m != 0) {
            for (int i = 0; i < m; ++i) {
                temp += fabsf(x[i]);
            }

            if (n < 6)
                return temp;
        } else {
            for (int i = m; i < n; i += 6) {
                temp += fabsf(x[i]) + fabsf(x[i+1]) + fabsf(x[i+2]) + fabsf(x[i+3]) + fabsf(x[i+4]) + fabsf(x[i+5]);
            }
        }
    } else {
        int nincx = n * incx;
        for (int i = 0; i < nincx; i += incx) {
            temp += fabsf(x[i]);
        }
    }

    return temp;
}

double dasum(unsigned int n, double *x, int incx) {
    double temp = 0;
    if (n == 0 || incx <= 0)
        return 0;

    if (incx == 1) {
        int m = n % 6;
        if (m != 0) {
            for (int i = 0; i < m; ++i) {
                temp += fabs(x[i]);
            }

            if (n < 6)
                return temp;
        } else {
            for (int i = m; i < n; i += 6) {
                temp += fabs(x[i]) + fabs(x[i+1]) + fabs(x[i+2]) + fabs(x[i+3]) + fabs(x[i+4]) + fabs(x[i+5]);
            }
        }
    } else {
        int nincx = n * incx;
        for (int i = 0; i < nincx; i += incx) {
            temp += fabs(x[i]);
        }
    }

    return temp;
}

float scasum(unsigned int n, complex float *x, int incx) {
    float temp = 0;
    if (n == 0 || incx <= 0)
        return 0;

    if (incx == 1) {
        int m = n % 6;
        if (m != 0) {
            for (int i = 0; i < m; ++i) {
                temp += cabsf(x[i]);
            }

            if (n < 6)
                return temp;
        } else {
            for (int i = m; i < n; i += 6) {
                temp += cabsf(x[i]) + cabsf(x[i+1]) + cabsf(x[i+2]) + cabsf(x[i+3]) + cabsf(x[i+4]) + cabsf(x[i+5]);
            }
        }
    } else {
        int nincx = n * incx;
        for (int i = 0; i < nincx; i += incx) {
            temp += cabsf(x[i]);
        }
    }

    return temp;
}

double dzasum(unsigned int n, complex double *x, int incx) {
    double temp = 0;
    if (n == 0 || incx <= 0)
        return 0;

    if (incx == 1) {
        int m = n % 6;
        if (m != 0) {
            for (int i = 0; i < m; ++i) {
                temp += cabs(x[i]);
            }

            if (n < 6)
                return temp;
        } else {
            for (int i = m; i < n; i += 6) {
                temp += cabs(x[i]) + cabs(x[i+1]) + cabs(x[i+2]) + cabs(x[i+3]) + cabs(x[i+4]) + cabs(x[i+5]);
            }
        }
    } else {
        int nincx = n * incx;
        for (int i = 0; i < nincx; i += incx) {
            temp += cabs(x[i]);
        }
    }

    return temp;
}

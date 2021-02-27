//
// Created by jon on 2/27/21.
//

#include <math.h>
#include "coma_nrm2.h"

float snrm2(unsigned int n, float *x, int incx) {
    float norm;
    float absxi;
    if (n == 0 || incx < 1)
        return 0;
    else if (n == 1){
        return fabsf(x[0]);
    } else {
        float scale = 0;
        float ssq = 1;
        for (int i = 0; i < 1 + (n-1) * incx; i += incx) {
            if (x[i] != 0) {
                absxi = fabsf(x[i]);
                if (scale < absxi) {
                    ssq = 1 + ssq * pow((scale / absxi), 2);
                    scale = absxi;
                } else {
                    ssq += pow(absxi / scale, 2);
                }
            }
        }
        norm = scale * sqrt(ssq);
    }
    return norm;
}

double dnrm2(unsigned int n, double *x, int incx) {
    double norm;
    double absxi;
    if (n == 0 || incx < 1)
        return 0;
    else if (n == 1){
        return fabs(x[0]);
    } else {
        double scale = 0;
        double ssq = 1;
        for (int i = 0; i < 1 + (n-1) * incx; i += incx) {
            if (x[i] != 0) {
                absxi = fabs(x[i]);
                if (scale < absxi) {
                    ssq = 1 + ssq * pow((scale / absxi), 2);
                    scale = absxi;
                } else {
                    ssq += pow(absxi / scale, 2);
                }
            }
        }
        norm = scale * sqrt(ssq);
    }
    return norm;
}

float scnrm2(unsigned int n, complex float *x, int incx) {
    float norm;
    float temp;
    if (n == 0 || incx < 1)
        return 0;
    else {
        float scale = 0;
        float ssq = 1;
        for (int i = 0; i < 1 + (n-1) * incx; i += incx) {
            if (crealf(x[i]) != 0) {
                temp = fabsf(crealf(x[i]));
                if (scale < temp) {
                    ssq = 1 + ssq * pow((scale / temp), 2);
                    scale = temp;
                } else {
                    ssq += pow(temp / scale, 2);
                }
            }
            if (cimagf(x[i]) != 0) {
                temp = fabsf(cimagf(x[i]));
                if (scale < temp) {
                    ssq = 1 + ssq * pow((scale / temp), 2);
                    scale = temp;
                } else {
                    ssq += pow(temp / scale, 2);
                }
            }
        }
        norm = scale * sqrt(ssq);
    }
    return norm;
}

double dznrm2(unsigned int n, complex double *x, int incx) {
    double norm;
    double temp;
    if (n == 0 || incx < 1)
        return 0;
    else {
        double scale = 0;
        double ssq = 1;
        for (int i = 0; i < 1 + (n-1) * incx; i += incx) {
            if (creal(x[i]) != 0) {
                temp = fabs(creal(x[i]));
                if (scale < temp) {
                    ssq = 1 + ssq * pow((scale / temp), 2);
                    scale = temp;
                } else {
                    ssq += pow(temp / scale, 2);
                }
            }
            if (cimagf(x[i]) != 0) {
                temp = fabs(cimag(x[i]));
                if (scale < temp) {
                    ssq = 1 + ssq * pow((scale / temp), 2);
                    scale = temp;
                } else {
                    ssq += pow(temp / scale, 2);
                }
            }
        }
        norm = scale * sqrt(ssq);
    }
    return norm;
}

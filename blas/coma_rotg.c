//
// Created by jon on 2/20/21.
//

#include "coma_rotg.h"
#include <math.h>

void srotg(float* a, float* b, float* c, float* s) {
    float r,roe,scale,z;
    roe = *b;
    if (fabsf(*a) > fabsf(*b))
        roe = *a;

    scale = fabsf(*a) + fabsf(*b);

    if (scale == 0) {
        *c = 1, *s = 0, r = 0, z = 0;
    } else {
        r = scale * sqrt(pow((*a / scale), 2) + pow((*b / scale), 2));
        r = roe / fabsf(roe) * r;
        *c = *a / r;
        *s = *b / r;
        z = 1;

        if (fabsf(*a) > fabsf(*b))
            z = *s;
        if(fabsf(*b) >= fabsf(*a) && *c != 0)
            z = 1 / *c;
    }
    *a = r;
    *b = z;
}

void drotg(double* a, double* b, double* c, double* s) {
    double r,roe,scale,z;
    roe = *b;
    if (fabs(*a) > fabs(*b))
        roe = *a;

    scale = fabs(*a) + fabs(*b);

    if (scale == 0) {
        *c = 1, *s = 0, r = 0, z = 0;
    } else {
        r = scale * sqrt(pow((*a / scale), 2) + pow((*b / scale), 2));
        r = roe / fabs(roe) * r;
        *c = *a / r;
        *s = *b / r;
        z = 1;

        if (fabs(*a) > fabs(*b))
            z = *s;
        if(fabs(*b) >= fabs(*a) && *c != 0)
            z = 1 / *c;
    }
    *a = r;
    *b = z;
}

void crotg(complex float *a, complex float *b, float *c, float *s) {
    if (cabsf(*a) == 0){
        *c = 0;
        *s = 1+0*I;
        *a = *b;
    } else {
        float scale = cabsf(*a) + cabsf(*b);
        float norm = scale * sqrtf(pow((cabs(*a/scale)),2)+ pow(cabs(*b / scale),2));
        complex float alpha = *a / cabsf(*a);
        *c = cabsf(*a) / norm;
        *s = alpha * conjf(*b) / norm;
        *a = alpha * norm;
    }
}

void zrotg(complex double *a, complex double *b, double *c, double *s) {
    if (cabs(*a) == 0){
        *c = 0;
        *s = 1+0*I;
        *a = *b;
    } else {
        double scale = cabs(*a) + cabs(*b);
        double norm = scale * sqrt(pow((cabs(*a/scale)),2)+ pow(cabs(*b / scale),2));
        complex double alpha = *a / cabs(*a);
        *c = cabs(*a) / norm;
        *s = alpha * conjf(*b) / norm;
        *a = alpha * norm;
    }
}

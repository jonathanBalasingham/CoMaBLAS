//
// Created by jon on 2/25/21.
//

#ifndef COMABLAS_COMA_ASUM_H
#define COMABLAS_COMA_ASUM_H

#include <math.h>

float sasum(unsigned int n, float* x, int incx);
double dasum(unsigned int n, double* x, int incx);
float scasum(unsigned int n, complex float *x, int incx);
double dzasum(unsigned int n, complex double *x, int incx);

#endif //COMABLAS_COMA_ASUM_H

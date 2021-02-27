//
// Created by jon on 2/27/21.
//

#ifndef COMABLAS_COMA_NRM2_H
#define COMABLAS_COMA_NRM2_H

#include <complex.h>

float snrm2(unsigned int n, float* x, int incx);
double dnrm2(unsigned int n, double* x, int incx);
float scnrm2(unsigned int n, complex float* x, int incx);
double dznrm2(unsigned int n, complex double* x, int incx);

#endif //COMABLAS_COMA_NRM2_H

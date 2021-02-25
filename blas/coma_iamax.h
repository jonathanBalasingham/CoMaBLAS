//
// Created by jon on 2/25/21.
//

#ifndef COMABLAS_COMA_IAMAX_H
#define COMABLAS_COMA_IAMAX_H

#include <complex.h>
#include <math.h>

int isamax(unsigned int n, float* x, int incx);
int idamax(unsigned int n, double* x, int incx);
int icamax(unsigned int n, complex float* x, int incx);
int izamax(unsigned int n, complex double* x, int incx);

#endif //COMABLAS_COMA_IAMAX_H

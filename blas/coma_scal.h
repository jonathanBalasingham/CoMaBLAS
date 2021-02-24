//
// Created by jon on 2/20/21.
//

#ifndef COMABLAS_COMA_SCAL_H
#define COMABLAS_COMA_SCAL_H

#include <complex.h>

void sscal(unsigned int n, float sa, float* sx, int incx);
void dscal(unsigned int n, double sa, double* sx, int incx);
void cscal(int n, complex float sa,complex float* sx, int incx);
void zscal(int n,complex double sa,complex double* sx, int incx);



#endif //COMABLAS_COMA_SCAL_H

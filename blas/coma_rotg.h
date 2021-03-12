//
// Created by jon on 2/20/21.
//

#ifndef COMABLAS_COMA_ROTG_H
#define COMABLAS_COMA_ROTG_H

#include <stdlib.h>
#include <math.h>
#include <complex.h>

void srotg(float* sa, float* sb, float* c, float* s);
void drotg(double* sa, double* sb, double* c, double* s);
void crotg(complex float* sa, complex float* sb, float* c, float* s);
void zrotg(complex double* sa, complex double* sb, double* c, double* s);


#endif //COMABLAS_COMA_ROTG_H

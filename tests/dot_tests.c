//
// Created by jon on 2/21/21.
//
#include <stdio.h>
#include "../blas/coma_dot.h"
#include "../munit/munit.h"

void t_mod(int* a, int* b, int* c) {
    *a = 1;
    *b = 2;
    *c = 3;
}

void main(){
    float t1[3] = {1,2,3};
    float t2[3] = {1,2,3};
    int n = 3;

    float answer = sdot(n, (float *) &t1, 1, (float *) &t2, 1);
    munit_assert_float(14.0, ==, answer);

    float t3[7] = {1,1,1,1,1,1,1};
    float t4[7] = {1,1,1,1,1,1,1};
    int n2 = 7;

    float answer2 = sdot(n2, (float *) &t3, 1, (float *) &t4, 1);
    munit_assert_float(7.0, ==, answer2);

    double t5[3] = {1,2,3};
    double t6[3] = {1,2,3};

    double answer3 = ddot(n, (double *) &t5, 1, (double *) &t6, 1);
    munit_assert_double(14.0, ==, answer3);

    double t7[7] = {1,1,1,1,1,1,1};
    double t8[7] = {1,1,1,1,1,1,1};

    double answer4 = ddot(n2, (double *) &t7, 1, (double *) &t8, 1);
    munit_assert_double(7.0, ==, answer4);


    int a,b,c;
    a=0,b=0,c=0;
    t_mod(&a,&b,&c);
    printf("%d, %d, %d", a, b, c);
}
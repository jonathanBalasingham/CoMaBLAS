//
// Created by jon on 2/21/21.
//
#include <stdio.h>
#include "../blas/coma_dot.h"

int main(){
    float t1[3] = {1,2,3};
    float t2[3] = {1,2,3};
    int n = 3;

    float answer = sdot(n, (float *) &t1, 1, (float *) &t2, 1);
    printf("answer: %f\n", answer);

    float t3[7] = {1,1,1,1,1,1,1};
    float t4[7] = {1,1,1,1,1,1,1};
    int n2 = 7;

    float answer2 = sdot(n2, (float *) &t3, 1, (float *) &t4, 1);
    printf("answer: %f\n", answer2);
}
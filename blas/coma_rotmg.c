//
// Created by jon on 2/20/21.
//

#include "coma_rotmg.h"
#include <math.h>

float srotmg(float* d1, float* d2, float* x1, float* y1, float *param) {
    float gam,gamsq,rgamsq,sflag,sh11,sh12,sh21,sh22,sp1,sp2,sq1,sq2,stemp,su;
    gam = 4096, gamsq = 1.67772e7, rgamsq =  5.96046e-8;

    if (d1 < 0 ){
        sflag = -1;
        sh11 = 0, sh12 = 0, sh21 = 0, sh22 = 0;
        *d1 = 0, *d2 = 0, *x1 = 0;
    } else {
        sp2 = *d2 * *y1;
        if (sp2 == 0) {
            sflag = -2;
            param[0] = sflag;
        }
        sp1 = *d1 * *x1;
        sq2 = sp2 * *y1;
        sq1 = sp1 * *x1;

        if (fabsf(sq1) > fabsf(sq2)) {
            sh21 = -1* *y1 / *x1;
            sh12 = sp2 / sp1;
            su = 1 - sh12 * sh21;

            if (su > 0){
                sflag = 0;
                *d1 = *d1 / su;
                *d2 = *d2 / su;
                *x1 = *x1 * su;
            }
        } else {
            if (sq2 < 0) {
                sflag = -1;
                sh11 = 0, sh12 = 0, sh21 = 0, sh22 = 0;
                *d1 = 0, *d2 = 0, *x1 = 0;
            } else {
                sflag = 1;
                sh11 = sp1 / sp2;
                sh22 = *x1 / *y1;
                su = 1 + sh11 * sh22;
                stemp = *d2 / su;
                *d2 = *d1 / su;
                *d1 = stemp;
                *x1 = *y1 * su;
            }
        }

        if (*d2 != 0) {
            while (*d1 <= rgamsq || *d1 >= gamsq) {
                if (sflag == 0){
                    sh11 = 1, sh22 = 1, sflag = -1;
                } else {
                    sh21 = -1, sh12 = 1, sflag = -1;
                }

                if (*d1 <= rgamsq) {
                    *d1 = *d1 * pow(gam,2);
                    *x1 = *x1 / gam;
                    sh11 = sh11 / gam;
                    sh12 = sh12 / gam;
                } else {
                    *d1 = *d1 / pow(gam,2);
                    *x1 = *x1 * gam;
                    sh11 = sh11 * gam;
                    sh12 = sh12 * gam;
                }
            }
        }

        if (*d2 != 0) {
            while (fabsf(*d2) <= rgamsq || fabsf(*d2) >= gamsq) {
                if (sflag == 0){
                    sh11 = 1;
                    sh22 = 1;
                    sflag = -1;
                } else {
                    sh21 = -1;
                    sh12 = 1;
                    sflag = -1;
                }

                if (fabsf(*d2) <= rgamsq) {
                    *d2 = *d2 * pow(gam,2);
                    sh21 = sh21 / gam;
                    sh22 = sh22 / gam;
                } else {
                    *d2 = *d2 / pow(gam,2);
                    sh21 = sh21 * gam;
                    sh22 = sh22 * gam;
                }
            }
        }
    }

    if (sflag < 0) {
        param[1] = sh11;
        param[2] = sh21;
        param[3] = sh12;
        param[4] = sh22;
    } else if (sflag == 0) {
        param[2] = sh21;
        param[3] = sh12;
    } else {
        param[1] = sh11;
        param[4] = sh22;
    }

    param[0] = sflag;
}

float drotmg(double* d1, double* d2, double* x1, double* y1, double *param) {
    double gam,gamsq,rgamsq,sflag,sh11,sh12,sh21,sh22,sp1,sp2,sq1,sq2,stemp,su;
    gam = 4096, gamsq = 1.67772e7, rgamsq =  5.96046e-8;

    if (d1 < 0 ){
        sflag = -1;
        sh11 = 0, sh12 = 0, sh21 = 0, sh22 = 0;
        *d1 = 0, *d2 = 0, *x1 = 0;
    } else {
        sp2 = *d2 * *y1;
        if (sp2 == 0) {
            sflag = -2;
            param[0] = sflag;
        }
        sp1 = *d1 * *x1;
        sq2 = sp2 * *y1;
        sq1 = sp1 * *x1;

        if (fabsf(sq1) > fabsf(sq2)) {
            sh21 = -1* *y1 / *x1;
            sh12 = sp2 / sp1;
            su = 1 - sh12 * sh21;

            if (su > 0){
                sflag = 0;
                *d1 = *d1 / su;
                *d2 = *d2 / su;
                *x1 = *x1 * su;
            }
        } else {
            if (sq2 < 0) {
                sflag = -1;
                sh11 = 0, sh12 = 0, sh21 = 0, sh22 = 0;
                *d1 = 0, *d2 = 0, *x1 = 0;
            } else {
                sflag = 1;
                sh11 = sp1 / sp2;
                sh22 = *x1 / *y1;
                su = 1 + sh11 * sh22;
                stemp = *d2 / su;
                *d2 = *d1 / su;
                *d1 = stemp;
                *x1 = *y1 * su;
            }
        }

        if (*d2 != 0) {
            while (*d1 <= rgamsq || *d1 >= gamsq) {
                if (sflag == 0){
                    sh11 = 1, sh22 = 1, sflag = -1;
                } else {
                    sh21 = -1, sh12 = 1, sflag = -1;
                }

                if (*d1 <= rgamsq) {
                    *d1 = *d1 * pow(gam,2);
                    *x1 = *x1 / gam;
                    sh11 = sh11 / gam;
                    sh12 = sh12 / gam;
                } else {
                    *d1 = *d1 / pow(gam,2);
                    *x1 = *x1 * gam;
                    sh11 = sh11 * gam;
                    sh12 = sh12 * gam;
                }
            }
        }

        if (*d2 != 0) {
            while (fabsf(*d2) <= rgamsq || fabsf(*d2) >= gamsq) {
                if (sflag == 0){
                    sh11 = 1;
                    sh22 = 1;
                    sflag = -1;
                } else {
                    sh21 = -1;
                    sh12 = 1;
                    sflag = -1;
                }

                if (fabsf(*d2) <= rgamsq) {
                    *d2 = *d2 * pow(gam,2);
                    sh21 = sh21 / gam;
                    sh22 = sh22 / gam;
                } else {
                    *d2 = *d2 / pow(gam,2);
                    sh21 = sh21 * gam;
                    sh22 = sh22 * gam;
                }
            }
        }
    }

    if (sflag < 0) {
        param[1] = sh11;
        param[2] = sh21;
        param[3] = sh12;
        param[4] = sh22;
    } else if (sflag == 0) {
        param[2] = sh21;
        param[3] = sh12;
    } else {
        param[1] = sh11;
        param[4] = sh22;
    }

    param[0] = sflag;
}

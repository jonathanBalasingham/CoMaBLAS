//
// Created by jon on 2/20/21.
//

#include "coma_rotg.h"
#include <math.h>

void srotg(float* sa, float* sb, float* c, float* s) {
    float r,roe,scale,z;
    roe = *sb;
    if (fabsf(*sa) > fabsf(*sb))
        roe = *sa;

    scale = fabsf(*sa) + fabsf(*sb);

    if (scale == 0) {
        *c = 1, *s = 0, r = 0, z = 0;
    } else {
        r = scale * sqrt(pow((*sa / scale), 2) + pow((*sb / scale), 2));
        r = roe / fabsf(roe) * r;
        *c = *sa / r;
        *s = *sb / r;
        z = 1;

        if (fabsf(*sa) > fabsf(*sb))
            z = *s;
        if(fabsf(*sb) >= fabsf(*sa) && *c != 0)
            z = 1 / *c;
    }
    *sa = r;
    *sb = z;
}

void drotg(double* sa, double* sb, double* c, double* s) {
    double r,roe,scale,z;
    roe = *sb;
    if (fabs(*sa) > fabs(*sb))
        roe = *sa;

    scale = fabs(*sa) + fabs(*sb);

    if (scale == 0) {
        *c = 1, *s = 0, r = 0, z = 0;
    } else {
        r = scale * sqrt(pow((*sa / scale), 2) + pow((*sb / scale), 2));
        r = roe / fabs(roe) * r;
        *c = *sa / r;
        *s = *sb / r;
        z = 1;

        if (fabs(*sa) > fabs(*sb))
            z = *s;
        if(fabs(*sb) >= fabs(*sa) && *c != 0)
            z = 1 / *c;
    }
    *sa = r;
    *sb = z;
}

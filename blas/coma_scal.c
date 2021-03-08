//
// Created by jon on 2/20/21.
//

#include "coma_scal.h"

void sscal(unsigned int n, float sa, float *sx, int incx) {
    int m, nincx;

    if (n == 0)
        return;

    if (incx == 1) {
        m = n % 5;
        if (m != 0){
            for (int j = 0; j < m; ++j) {
                sx[j] = sx[j] * sa;
            }

            if (n < 5)
                return;
        }

        for (int j = m; j < n; j += 5) {
            sx[j] = sx[j] * sa;
            sx[j+1] = sx[j+1] * sa;
            sx[j+2] = sx[j+2] * sa;
            sx[j+3] = sx[j+3] * sa;
            sx[j+4] = sx[j+4] * sa;
        }

    } else {
        nincx = n * incx;
        for (int i = 0; i < nincx; i += incx) {
            sx[i] = sx[i] * sa;
        }
    }

}

void dscal(unsigned int n, double sa, double *sx, int incx) {
    int m, nincx;

    if (n == 0)
        return;

    if (incx == 1) {
        m = n % 5;
        if (m != 0){
            for (int j = 0; j < m; ++j) {
                sx[j] = sx[j] * sa;
            }

            if (n < 5)
                return;
        }

        for (int j = m; j < n; j += 5) {
            sx[j] = sx[j] * sa;
            sx[j+1] = sx[j+1] * sa;
            sx[j+2] = sx[j+2] * sa;
            sx[j+3] = sx[j+3] * sa;
            sx[j+4] = sx[j+4] * sa;
        }

    } else {
        nincx = n * incx;
        for (int i = 0; i < nincx; i += incx) {
            sx[i] = sx[i] * sa;
        }
    }
}

void cscal(int n, complex float sa, complex float *sx, int incx) {
    if (n == 0)
        return;

    if (incx == 1){
        for (int i = 0; i < n; ++i) {
            sx[i] = sx[i] * sa;
        }
    } else {
        int nincx = n * incx;
        for (int i = 0; i < nincx; i += incx) {
            sx[i] = sx[i] * sa;
        }
    }
}

void csscal(int n, float sa, complex float *sx, int incx) {
    if (n == 0)
        return;

    if (incx == 1){
        for (int i = 0; i < n; ++i) {
            sx[i] = sx[i] * sa;
        }
    } else {
        int nincx = n * incx;
        for (int i = 0; i < nincx; i += incx) {
            sx[i] = sx[i] * sa;
        }
    }
}

void zscal(int n, complex double sa, complex double *sx, int incx) {
    if (n == 0)
        return;

    if (incx == 1){
        for (int i = 0; i < n; ++i) {
            sx[i] = sx[i] * sa;
        }
    } else {
        int nincx = n * incx;
        for (int i = 0; i < nincx; i += incx) {
            sx[i] = sx[i] * sa;
        }
    }
}

void zdscal(int n, double sa, complex double *sx, int incx) {
    if (n == 0)
        return;

    if (incx == 1){
        for (int i = 0; i < n; ++i) {
            sx[i] = sx[i] * sa;
        }
    } else {
        int nincx = n * incx;
        for (int i = 0; i < nincx; i += incx) {
            sx[i] = sx[i] * sa;
        }
    }
}

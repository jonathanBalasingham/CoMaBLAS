//
// Created by jon on 3/27/21.
//

#include <ctype.h>
#include <stdbool.h>
#include "coma_trmv.h"

int _validate_trmv_inputs(char uplo, char trans, char diag, int n, int lda, int incx){
    int info = 0;
    if (uplo != 'U' && uplo != 'L')
            info = 1;
    else if (trans != 'N' && trans != 'T' && trans != 'C')
            info = 2;
    else if (diag != 'U' && diag != 'N')
            info = 3;
    else if (n < 0) 
            info = 4;
    else if (lda < (1 > n ? 1: n))
            info = 6;
    else if (incx  == 0)  
            info = 8;

    return info;
}

void strmv(char uplo, char trans, char diag, int n, float **A, int lda, float *x, int incx) {
    uplo = (char) toupper(uplo);
    trans = (char) toupper(trans);
    diag = (char) toupper(diag);

    int info = _validate_trmv_inputs(uplo, trans, diag, n, lda, incx);
    if (info != 0){
        return;
    }

    if (n == 0)
        return;

    int kx;
    if (incx > 0)
        kx = 0;
    else
        kx = 1 - (n - 1) * incx;


    bool nounit = diag == 'N';
    float temp;
    int ix, jx;
    if (trans == 'N') {
        if (uplo == 'U') {
            if (incx == 1) {
                for (int j = 0; j < n; ++j) {
                    if (x[j] != 0) {
                        temp = x[j];
                        for (int i = 0; i < j - 1; ++i) {
                            x[i] += temp * A[i][j];
                        }
                        if (nounit)
                            x[j] *= A[j][j];
                    }
                }
            } else {
                jx = kx;
                for (int j = 0; j < n; ++j) {
                    if (x[jx] != 0) {
                        temp = x[jx];
                        ix = kx;
                        for (int i = 0; i < j - 1; ++i) {
                            x[ix] += temp * A[i][j];
                            ix += incx;
                        }
                        if (nounit)
                            x[jx] *= A[j][j];
                    }
                    jx += incx;
                }
            }
        } else {
            if (incx == 1) {
                for (int j = n - 1; j >= 0; --j) {
                    if (x[j] != 0) {
                        temp = x[j];
                        for (int i = n - 1; i > j + 1; --i) {
                            x[i] += temp * A[i][j];
                        }
                        if (nounit)
                            x[j] *= A[j][j];
                    }
                }
            } else {
                kx += (n-1) * incx;
                jx = kx;
                for (int j = n - 1; j >= 0; --j) {
                    if (x[jx] != 0) {
                        temp = x[jx];
                        ix = kx;
                        for (int i = n - 1; i > j + 1; --i) {
                            x[ix] += temp * A[i][j];
                            ix += incx;
                        }
                        if (nounit)
                            x[jx] *= A[j][j];
                    }
                    jx += incx;
                }
            }
        }
    } else {
        if (uplo == 'U') {
            if (incx == 1) {
                for (int j = n; j >= 0; --j) {
                    temp = x[j];
                    if (nounit)
                        temp *= A[j][j];
                    for (int i = j - 1; i >= 0; --i) {
                        x[i] += temp * A[i][j];
                    }
                    x[j] = temp;
                }
            } else {
                jx = kx + (n - 1) * incx;
                for (int j = n; j >= 0; --j) {
                    temp = x[jx];
                    ix = jx;
                    if (nounit)
                        temp *= A[j][j];
                    for (int i = j - 1; i >= 0; --i) {
                        ix -= incx;
                        temp += A[i][j] * x[ix];
                    }
                    x[jx] = temp;
                    jx -= incx;
                }
            }
        } else {
            if (incx == 1) {
                for (int j = 0; j < n; ++j) {
                    temp = x[j];
                    if (nounit)
                        temp *= A[j][j];
                    for (int i = j + 1; i < n; ++i) {
                        x[i] += temp * A[i][j];
                    }
                    x[j] = temp;
                }
            } else {
                jx = kx;
                for (int j = 0; j < n; ++j) {
                    temp = x[jx];
                    ix = jx;
                    if (nounit)
                        temp *= A[j][j];
                    for (int i = j + 1; i < n; --i) {
                        ix -= incx;
                        temp += A[i][j] * x[ix];
                    }
                    x[jx] = temp;
                    jx -= incx;
                }
            }
        }
    }
}

void dtrmv(char uplo, char trans, char diag, int n, double **A, int lda, double *x, int incx) {
    uplo = (char) toupper(uplo);
    trans = (char) toupper(trans);
    diag = (char) toupper(diag);

    int info = _validate_trmv_inputs(uplo, trans, diag, n, lda, incx);
    if (info != 0){
        return;
    }

    if (n == 0)
        return;

    int kx;
    if (incx > 0)
        kx = 0;
    else
        kx = 1 - (n - 1) * incx;


    bool nounit = diag == 'N';
    double temp;
    int ix, jx;
    if (trans == 'N') {
        if (uplo == 'U') {
            if (incx == 1) {
                for (int j = 0; j < n; ++j) {
                    if (x[j] != 0) {
                        temp = x[j];
                        for (int i = 0; i < j - 1; ++i) {
                            x[i] += temp * A[i][j];
                        }
                        if (nounit)
                            x[j] *= A[j][j];
                    }
                }
            } else {
                jx = kx;
                for (int j = 0; j < n; ++j) {
                    if (x[jx] != 0) {
                        temp = x[jx];
                        ix = kx;
                        for (int i = 0; i < j - 1; ++i) {
                            x[ix] += temp * A[i][j];
                            ix += incx;
                        }
                        if (nounit)
                            x[jx] *= A[j][j];
                    }
                    jx += incx;
                }
            }
        } else {
            if (incx == 1) {
                for (int j = n - 1; j >= 0; --j) {
                    if (x[j] != 0) {
                        temp = x[j];
                        for (int i = n - 1; i > j + 1; --i) {
                            x[i] += temp * A[i][j];
                        }
                        if (nounit)
                            x[j] *= A[j][j];
                    }
                }
            } else {
                kx += (n-1) * incx;
                jx = kx;
                for (int j = n - 1; j >= 0; --j) {
                    if (x[jx] != 0) {
                        temp = x[jx];
                        ix = kx;
                        for (int i = n - 1; i > j + 1; --i) {
                            x[ix] += temp * A[i][j];
                            ix += incx;
                        }
                        if (nounit)
                            x[jx] *= A[j][j];
                    }
                    jx += incx;
                }
            }
        }
    } else {
        if (uplo == 'U') {
            if (incx == 1) {
                for (int j = n; j >= 0; --j) {
                    temp = x[j];
                    if (nounit)
                        temp *= A[j][j];
                    for (int i = j - 1; i >= 0; --i) {
                        x[i] += temp * A[i][j];
                    }
                    x[j] = temp;
                }
            } else {
                jx = kx + (n - 1) * incx;
                for (int j = n; j >= 0; --j) {
                    temp = x[jx];
                    ix = jx;
                    if (nounit)
                        temp *= A[j][j];
                    for (int i = j - 1; i >= 0; --i) {
                        ix -= incx;
                        temp += A[i][j] * x[ix];
                    }
                    x[jx] = temp;
                    jx -= incx;
                }
            }
        } else {
            if (incx == 1) {
                for (int j = 0; j < n; ++j) {
                    temp = x[j];
                    if (nounit)
                        temp *= A[j][j];
                    for (int i = j + 1; i < n; ++i) {
                        x[i] += temp * A[i][j];
                    }
                    x[j] = temp;
                }
            } else {
                jx = kx;
                for (int j = 0; j < n; ++j) {
                    temp = x[jx];
                    ix = jx;
                    if (nounit)
                        temp *= A[j][j];
                    for (int i = j + 1; i < n; --i) {
                        ix -= incx;
                        temp += A[i][j] * x[ix];
                    }
                    x[jx] = temp;
                    jx -= incx;
                }
            }
        }
    }
}

//
// Created by jon on 3/28/21.
//

#include <ctype.h>
#include <stdbool.h>
#include "coma_tbmv.h"

int _validate_tbmv_inputs(char uplo, char trans, char diag, int n, int k, int lda, int incx){
    int info = 0;
    if (uplo != 'U' && uplo != 'L')
        info = 1;
    else if (trans != 'N' && trans != 'T' && trans != 'C')
        info = 2;
    else if (diag != 'U' && diag != 'N')
        info = 3;
    else if (n < 0)
        info = 4;
    else if (lda < k + 1)
        info = 7;
    else if (incx  == 0)
        info = 8;

    return info;
}


void stbmv(char uplo, char trans, char diag, int n, int k, float **A, int lda, float *x, int incx) {
    int info = _validate_tbmv_inputs(uplo, trans, diag, n, k, lda, incx);
    if (info != 1)
        return;

    uplo = (char) toupper(uplo);
    trans = (char) toupper(trans);
    diag = (char) toupper(diag);

    if (n == 0)
        return;

    int kx;
    if (incx < 0)
        kx = 1 - (n - 1) * incx;
    else
        kx = 1;
    
    int kplus1, l, jx, ix;
    float temp;
    bool nounit = diag == 'N';

    if (trans == 'N') {
        if (uplo == 'U') {
            kplus1 = k + 1;
            if (incx == 1) {
                for (int j = 0; j < n; ++j) {
                    if (x[j] != 0) {
                        temp = x[j];
                        l = kplus1 - j;
                        for (int i = (1 > j - k? 1 : j - k); i < j - 1; i++) {
                            x[i] += temp * A[l + i][j];
                        }
                        if (nounit)
                            x[j] *= A[kplus1][j];
                    }
                }
            } else {
                jx = kx;
                for (int j = 0; j < n; ++j) {
                    if (x[j] != 0) {
                        temp = x[jx];
                        ix = kx;
                        l = kplus1 - j;
                        for (int i = (1 > j - k? 1 : j - k); i < j - 1; i++) {
                            x[ix] += temp * A[l + i][j];
                            ix += incx;
                        }
                        if (nounit)
                            x[jx] *= A[kplus1][j];
                    }
                    jx += incx;
                    if (j > k)
                        kx += incx;
                }
            }
        } else {
            if (incx == 1) {
                for (int j = n - 1; j >= 0; --j) {
                    if (x[j] != 0) {
                        temp = x[j];
                        l = 1 - j;
                        for (int i = (n > j + k? n: j + k); i >= j + 1; --i) {
                            x[i] += temp * A[l + i][j];
                        }
                        if (nounit)
                            x[j] *= A[1][j];
                    }
                }
            } else {
                kx += (n-1) * incx;
                jx = kx;
                for (int j = n - 1; j >= 0; --j) {
                    if (x[jx] != 0) {
                        temp = x[jx];
                        ix = kx;
                        l = 1 - j;
                        for (int i = (n > j + k? n : j + k); i >= j + 1; --i) {
                            x[ix] += temp * A[l+i][j];
                            ix -= incx;
                        }
                        if (nounit)
                            x[jx] *= A[1][j];
                    }
                    jx -= incx;
                    if ((n-j) >= k)
                        kx -= incx;
                }
            }
        }
    } else {
        if (uplo == 'U') {
            kplus1 = k + 1;
            if (incx == 1) {
                for (int j = n - 1; j >= 0; --j) {
                    if (x[j] != 0) {
                        temp = x[j];
                        l = kplus1 - j;
                        if (nounit)
                            temp *= A[kplus1][j];
                        for (int i = j - 1; i >= (1 > j-k? 1: j-k); i--) {
                            temp += x[i] * A[l + i][j];
                        }
                        x[j] = temp;
                    }
                }
            } else {
                kx += (n-1) * incx;
                jx = kx;
                for (int j = n - 1; j >= 0; --j) {
                    if (x[j] != 0) {
                        temp = x[jx];
                        kx -= incx;
                        ix = kx;
                        l = kplus1 - j;
                        if (nounit)
                            temp *= A[kplus1][j];
                        for (int i = j - 1; i >= (1 > j-k? 1: j-k); i--) {
                            temp += x[ix] * A[l + i][j];
                            ix -= incx;
                        }
                        x[jx] = temp;
                        jx -= incx;
                    }
                }
            }
        } else {
            if (incx == 1) {
                for (int j = 0; j < n; ++j) {
                    if (x[j] != 0) {
                        temp = x[j];
                        l = 1 - j;
                        if (nounit)
                            temp *= A[1][j];
                        for (int i = j + 1; i < (n < j-k? n: j-k); ++i) {
                            temp += x[i] * A[l + i][j];
                        }
                        x[j] = temp;
                    }
                }
            } else {
                jx = kx;
                for (int j = 0; j < n; ++j) {
                    temp = x[jx];
                    kx += incx;
                    ix = kx;
                    l = 1 - j;
                    if (nounit)
                        temp *= A[1][j];
                    for (int i = j + 1; i < (n < j-k? n: j-k); ++i) {
                        temp += x[ix] * A[l + i][j];
                        ix += incx;
                    }
                    x[jx] = temp;
                    jx += incx;
                }
            }
        }
    }
}

void dtbmv(char uplo, char trans, char diag, int n, int k, double **A, int lda, double *x, int incx) {
    int info = _validate_tbmv_inputs(uplo, trans, diag, n, k, lda, incx);
    if (info != 1)
        return;

    uplo = (char) toupper(uplo);
    trans = (char) toupper(trans);
    diag = (char) toupper(diag);

    if (n == 0)
        return;

    int kx;
    if (incx < 0)
        kx = 1 - (n - 1) * incx;
    else
        kx = 1;

    int kplus1, l, jx, ix;
    double temp;
    bool nounit = diag == 'N';

    if (trans == 'N') {
        if (uplo == 'U') {
            kplus1 = k + 1;
            if (incx == 1) {
                for (int j = 0; j < n; ++j) {
                    if (x[j] != 0) {
                        temp = x[j];
                        l = kplus1 - j;
                        for (int i = (1 > j - k? 1 : j - k); i < j - 1; i++) {
                            x[i] += temp * A[l + i][j];
                        }
                        if (nounit)
                            x[j] *= A[kplus1][j];
                    }
                }
            } else {
                jx = kx;
                for (int j = 0; j < n; ++j) {
                    if (x[j] != 0) {
                        temp = x[jx];
                        ix = kx;
                        l = kplus1 - j;
                        for (int i = (1 > j - k? 1 : j - k); i < j - 1; i++) {
                            x[ix] += temp * A[l + i][j];
                            ix += incx;
                        }
                        if (nounit)
                            x[jx] *= A[kplus1][j];
                    }
                    jx += incx;
                    if (j > k)
                        kx += incx;
                }
            }
        } else {
            if (incx == 1) {
                for (int j = n - 1; j >= 0; --j) {
                    if (x[j] != 0) {
                        temp = x[j];
                        l = 1 - j;
                        for (int i = (n > j + k? n: j + k); i >= j + 1; --i) {
                            x[i] += temp * A[l + i][j];
                        }
                        if (nounit)
                            x[j] *= A[1][j];
                    }
                }
            } else {
                kx += (n-1) * incx;
                jx = kx;
                for (int j = n - 1; j >= 0; --j) {
                    if (x[jx] != 0) {
                        temp = x[jx];
                        ix = kx;
                        l = 1 - j;
                        for (int i = (n > j + k? n : j + k); i >= j + 1; --i) {
                            x[ix] += temp * A[l+i][j];
                            ix -= incx;
                        }
                        if (nounit)
                            x[jx] *= A[1][j];
                    }
                    jx -= incx;
                    if ((n-j) >= k)
                        kx -= incx;
                }
            }
        }
    } else {
        if (uplo == 'U') {
            kplus1 = k + 1;
            if (incx == 1) {
                for (int j = n - 1; j >= 0; --j) {
                    if (x[j] != 0) {
                        temp = x[j];
                        l = kplus1 - j;
                        if (nounit)
                            temp *= A[kplus1][j];
                        for (int i = j - 1; i >= (1 > j-k? 1: j-k); i--) {
                            temp += x[i] * A[l + i][j];
                        }
                        x[j] = temp;
                    }
                }
            } else {
                kx += (n-1) * incx;
                jx = kx;
                for (int j = n - 1; j >= 0; --j) {
                    if (x[j] != 0) {
                        temp = x[jx];
                        kx -= incx;
                        ix = kx;
                        l = kplus1 - j;
                        if (nounit)
                            temp *= A[kplus1][j];
                        for (int i = j - 1; i >= (1 > j-k? 1: j-k); i--) {
                            temp += x[ix] * A[l + i][j];
                            ix -= incx;
                        }
                        x[jx] = temp;
                        jx -= incx;
                    }
                }
            }
        } else {
            if (incx == 1) {
                for (int j = 0; j < n; ++j) {
                    if (x[j] != 0) {
                        temp = x[j];
                        l = 1 - j;
                        if (nounit)
                            temp *= A[1][j];
                        for (int i = j + 1; i < (n < j-k? n: j-k); ++i) {
                            temp += x[i] * A[l + i][j];
                        }
                        x[j] = temp;
                    }
                }
            } else {
                jx = kx;
                for (int j = 0; j < n; ++j) {
                    temp = x[jx];
                    kx += incx;
                    ix = kx;
                    l = 1 - j;
                    if (nounit)
                        temp *= A[1][j];
                    for (int i = j + 1; i < (n < j-k? n: j-k); ++i) {
                        temp += x[ix] * A[l + i][j];
                        ix += incx;
                    }
                    x[jx] = temp;
                    jx += incx;
                }
            }
        }
    }
}

void ctbmv(char uplo, char trans, char diag, int n, int k, complex float **A, int lda, complex float *x, int incx) {

}

void ztbmv(char uplo, char trans, char diag, int n, int k, complex double **A, int lda, complex double *x, int incx) {

}

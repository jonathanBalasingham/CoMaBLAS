//
// Created by jon on 3/12/21.
//

#include <ctype.h>
#include "coma_gemv.h"
#include <stdbool.h>

int _validate_inputs(char trans, int n, int m, int lda, int incx, int incy) {
    int info = 0;
    if (trans != 'N' && trans != 'T' && trans != 'C')
        info = 1;
    else if (m < 0)
        info = 2;
    else if (n < 0)
        info = 3;
    else if (lda < 1)
        info = 6;
    else if (incx == 0)
        info = 8;
    else if (incy == 0)
        info = 11;

    return info;
}

void sgemv(char trans, int m, int n, float alpha, float** A, int lda, float *x, int incx,
           float beta, float *y, int incy) {

    int info = _validate_inputs((char)toupper(trans), m,n,lda,incx,incy);

    if (info != 0){
        // throw error
    }

    int lenx, leny, kx, ky;
    if (trans == 'N') {
        lenx = n;
        leny = m;
    } else {
        lenx = m;
        leny = n;
    }

    if (incx > 0)
        kx = 0;
    else
        kx = 1 - (lenx - 1) * incx;

    if (incy > 0)
        ky = 0;
    else
        ky = 1 - (leny - 1) * incy;

    if (beta != 1) {
        if (incy == 1) {
            if (beta == 0) {
                for (int i = 0; i < leny; ++i) {
                    y[i] = 0;
                }
            } else {
                for (int i = 0; i < leny; ++i) {
                    y[i] = beta * y[i];
                }
            }
        } else {
            int iy = ky;
            if (beta == 0) {
                for (int i = 0; i < leny; ++i) {
                    y[iy] = 0;
                    iy += incy;
                }
            } else {
                for (int i = 0; i < leny; ++i) {
                    y[iy] = beta * y[i];
                    iy += incy;
                }
            }
        }
    }

    if (alpha == 0)
        return;

    if (trans == 'N') {
        int jx = kx;
        float temp;
        if (incy == 1) {
            for (int i = 0; i < n; ++i) {
                temp = alpha * x[jx];
                for (int j = 0; j < m; ++j) {
                    y[j] += temp * A[j][i];
                }
                jx += incx;
            }
        } else {
            int iy;
            for (int i = 0; i < n; ++i) {
                temp = alpha * x[jx];
                iy = ky;
                for (int j = 0; j < m; ++j) {
                    y[iy] += temp * A[j][i];
                    iy += incy;
                }
                jx += incx;
            }
        }
    } else {
        int jy = ky;
        if (incx == 1) {
            float temp;
            for (int i = 0; i < n; ++i) {
                temp = 0;
                for (int j = 0; j < m; ++j) {
                    temp += A[j][i] * x[j];
                }
                y[jy] += alpha * temp;
                jy += incy;
            }
        } else {
            float temp;
            int ix;
            for (int i = 0; i < n; ++i) {
                temp = 0;
                ix = kx;
                for (int j = 0; j < m; ++j) {
                    temp += A[j][i] * x[ix];
                    ix += incx;
                }
                y[jy] += alpha * temp;
                jy += incy;
            }
        }
    }
}

void dgemv(char trans, int m, int n, double alpha, double **A, int lda, double *x, int incx,
           double beta, double *y, int incy) {

    int info = 0;
    trans = (char) toupper(trans);

    if (trans != 'N' && trans != 'T' && trans != 'C')
        info = 1;
    else if (lda < 1)
        info = 6;
    else if (incx == 0)
        info = 8;
    else if (incy == 0)
        info = 11;

    if (info != 0){
        // throw error
    }

    int lenx, leny, kx, ky;
    if (trans == 'N') {
        lenx = n;
        leny = m;
    } else {
        lenx = m;
        leny = n;
    }

    if (incx > 0)
        kx = 0;
    else
        kx = 1 - (lenx - 1) * incx;

    if (incy > 0)
        ky = 0;
    else
        ky = 1 - (leny - 1) * incy;

    if (beta != 1) {
        if (incy == 1) {
            if (beta == 0) {
                for (int i = 0; i < leny; ++i) {
                    y[i] = 0;
                }
            } else {
                for (int i = 0; i < leny; ++i) {
                    y[i] = beta * y[i];
                }
            }
        } else {
            int iy = ky;
            if (beta == 0) {
                for (int i = 0; i < leny; ++i) {
                    y[iy] = 0;
                    iy += incy;
                }
            } else {
                for (int i = 0; i < leny; ++i) {
                    y[iy] = beta * y[i];
                    iy += incy;
                }
            }
        }
    }

    if (alpha == 0)
        return;

    if (trans == 'N') {
        int jx = kx;
        double temp;
        if (incy == 1) {
            for (int i = 0; i < n; ++i) {
                temp = alpha * x[jx];
                for (int j = 0; j < m; ++j) {
                    y[j] += temp * A[j][i];
                }
                jx += incx;
            }
        } else {
            int iy;
            for (int i = 0; i < n; ++i) {
                temp = alpha * x[jx];
                iy = ky;
                for (int j = 0; j < m; ++j) {
                    y[iy] += temp * A[j][i];
                    iy += incy;
                }
                jx += incx;
            }
        }
    } else {
        int jy = ky;
        if (incx == 1) {
            double temp;
            for (int i = 0; i < n; ++i) {
                temp = 0;
                for (int j = 0; j < m; ++j) {
                    temp += A[j][i] * x[j];
                }
                y[jy] += alpha * temp;
                jy += incy;
            }
        } else {
            double temp;
            int ix;
            for (int i = 0; i < n; ++i) {
                temp = 0;
                ix = kx;
                for (int j = 0; j < m; ++j) {
                    temp += A[j][i] * x[ix];
                    ix += incx;
                }
                y[jy] += alpha * temp;
                jy += incy;
            }
        }
    }
}

void
cgemv(char trans, int m, int n, complex float alpha, complex float **A, int lda, complex float *x,
      int incx, complex float beta, complex float *y, int incy) {

    int info = 0;
    trans = (char) toupper(trans);

    if (trans != 'N' && trans != 'T' && trans != 'C')
        info = 1;
    else if (lda < 1)
        info = 6;
    else if (incx == 0)
        info = 8;
    else if (incy == 0)
        info = 11;

    if (info != 0){
        // throw error
        return;
    }

    int lenx, leny, kx, ky;
    if (trans == 'N') {
        lenx = n;
        leny = m;
    } else {
        lenx = m;
        leny = n;
    }

    if (incx > 0)
        kx = 0;
    else
        kx = 1 - (lenx - 1) * incx;

    if (incy > 0)
        ky = 0;
    else
        ky = 1 - (leny - 1) * incy;

    int noconj = trans == 'T';

    if (beta != 1) {
        if (incy == 1) {
            if (beta == 0) {
                for (int i = 0; i < leny; ++i) {
                    y[i] = 0;
                }
            } else {
                for (int i = 0; i < leny; ++i) {
                    y[i] = beta * y[i];
                }
            }
        } else {
            int iy = ky;
            if (beta == 0) {
                for (int i = 0; i < leny; ++i) {
                    y[iy] = 0;
                    iy += incy;
                }
            } else {
                for (int i = 0; i < leny; ++i) {
                    y[iy] = beta * y[i];
                    iy += incy;
                }
            }
        }
    }

    if (alpha == 0)
        return;

    if (trans == 'N') {
        int jx = kx;
        complex float temp;
        if (incy == 1) {
            for (int i = 0; i < n; ++i) {
                temp = alpha * x[jx];
                for (int j = 0; j < m; ++j) {
                    y[j] += temp * A[j][i];
                }
                jx += incx;
            }
        } else {
            int iy;
            for (int i = 0; i < n; ++i) {
                temp = alpha * x[jx];
                iy = ky;
                for (int j = 0; j < m; ++j) {
                    y[iy] += temp * A[j][i];
                    iy += incy;
                }
                jx += incx;
            }
        }
    } else {
        int jy = ky;
        if (incx == 1) {
            complex float temp;
            for (int i = 0; i < n; ++i) {
                temp = 0;
                if (noconj) {
                    for (int j = 0; j < m; ++j) {
                        temp += A[j][i] * x[j];
                    }
                } else {
                    for (int j = 0; j < m; ++j) {
                        temp += conjf(A[j][i]) * x[j];
                    }
                }

                y[jy] += alpha * temp;
                jy += incy;
            }
        } else {
            complex float temp;
            int ix;
            for (int i = 0; i < n; ++i) {
                temp = 0;
                ix = kx;
                if (noconj) {
                    for (int j = 0; j < m; ++j) {
                        temp += A[j][i] * x[j];
                        ix += incx;
                    }
                } else {
                    for (int j = 0; j < m; ++j) {
                        temp += conjf(A[j][i]) * x[j];
                        ix += incx;
                    }
                }
                y[jy] += alpha * temp;
                jy += incy;
            }
        }
    }
}

void
zgemv(char trans, int m, int n, complex double alpha, complex double **A, int lda, complex double *x,
      int incx, complex double beta, complex double *y, int incy) {

    int info = 0;
    trans = (char) toupper(trans);

    if (trans != 'N' && trans != 'T' && trans != 'C')
        info = 1;
    else if (lda < 1)
        info = 6;
    else if (incx == 0)
        info = 8;
    else if (incy == 0)
        info = 11;

    if (info != 0){
        // throw error
        return;
    }

    int lenx, leny, kx, ky;
    if (trans == 'N') {
        lenx = n;
        leny = m;
    } else {
        lenx = m;
        leny = n;
    }

    if (incx > 0)
        kx = 0;
    else
        kx = 1 - (lenx - 1) * incx;

    if (incy > 0)
        ky = 0;
    else
        ky = 1 - (leny - 1) * incy;

    int noconj = trans == 'T';

    if (beta != 1) {
        if (incy == 1) {
            if (beta == 0) {
                for (int i = 0; i < leny; ++i) {
                    y[i] = 0;
                }
            } else {
                for (int i = 0; i < leny; ++i) {
                    y[i] = beta * y[i];
                }
            }
        } else {
            int iy = ky;
            if (beta == 0) {
                for (int i = 0; i < leny; ++i) {
                    y[iy] = 0;
                    iy += incy;
                }
            } else {
                for (int i = 0; i < leny; ++i) {
                    y[iy] = beta * y[i];
                    iy += incy;
                }
            }
        }
    }

    if (alpha == 0)
        return;

    if (trans == 'N') {
        int jx = kx;
        complex double temp;
        if (incy == 1) {
            for (int i = 0; i < n; ++i) {
                temp = alpha * x[jx];
                for (int j = 0; j < m; ++j) {
                    y[j] += temp * A[j][i];
                }
                jx += incx;
            }
        } else {
            int iy;
            for (int i = 0; i < n; ++i) {
                temp = alpha * x[jx];
                iy = ky;
                for (int j = 0; j < m; ++j) {
                    y[iy] += temp * A[j][i];
                    iy += incy;
                }
                jx += incx;
            }
        }
    } else {
        int jy = ky;
        if (incx == 1) {
            complex double temp;
            for (int i = 0; i < n; ++i) {
                temp = 0;
                if (noconj) {
                    for (int j = 0; j < m; ++j) {
                        temp += A[j][i] * x[j];
                    }
                } else {
                    for (int j = 0; j < m; ++j) {
                        temp += conjf(A[j][i]) * x[j];
                    }
                }

                y[jy] += alpha * temp;
                jy += incy;
            }
        } else {
            complex double temp;
            int ix;
            for (int i = 0; i < n; ++i) {
                temp = 0;
                ix = kx;
                if (noconj) {
                    for (int j = 0; j < m; ++j) {
                        temp += A[j][i] * x[j];
                        ix += incx;
                    }
                } else {
                    for (int j = 0; j < m; ++j) {
                        temp += conj(A[j][i]) * x[j];
                        ix += incx;
                    }
                }
                y[jy] += alpha * temp;
                jy += incy;
            }
        }
    }
}


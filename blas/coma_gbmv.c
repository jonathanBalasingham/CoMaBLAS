//
// Created by jon on 3/13/21.
//

#include <ctype.h>
#include "coma_gbmv.h"

int _validate_gbmv_inputs(char trans, int n, int m, int kl, int ku, int lda, int incx, int incy){
    int info = 0;
    if (trans != 'N' && trans != 'T' && trans != 'C')
        info = 1;
    else if (m < 0)
        info = 2;
    else if (n < 0)
        info = 3;
    else if (kl < 0)
        info = 4;
    else if (ku < 0)
        info = 5;
    else if (lda < kl + ku + 1)
        info = 8;
    else if (lda < 1)
        info = 6;
    else if (incx == 0)
        info = 8;
    else if (incy == 0)
        info = 11;

    return info;
}

void
sgbmv(char trans, int m, int n, int kl, int ku, float alpha, float **A, int lda,
      float *x, int incx, float beta, float *y, int incy) {

    int info = _validate_gbmv_inputs((char)toupper(trans), m,n,kl,ku,lda,incx,incy);

    if (info != 0){
        // error
        return;
    }

    if (m == 0 || n == 0 || (alpha == 0 && beta == 1))
        return;

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
                    y[iy] = beta * y[iy];
                    iy += incy;
                }
            }
        }
    }

    if (alpha == 0)
        return;

    int kup1 = ku + 1;
    float temp;
    int k, start, stop;
    if (trans == 'N') {
        int jx = kx;
        if (incy == 1) {
            for (int i = 0; i < n; ++i) {
                temp = alpha * x[jx];
                k = kup1 - i;
                start = 0 > i - ku ? 0 : i - ku;
                stop = m < i + kl ? m : i + kl;
                for (int j = start; j < stop; ++j) {
                    y[j] += temp * A[k+j][i];
                }
                jx += incx;
            }
        } else {
            int iy;
            for (int i = 0; i < n; ++i) {
                temp = alpha * x[jx];
                iy = ky;
                k = kup1 - i;
                start = 0 > i - ku ? 0 : i - ku;
                stop = m < i + kl ? m : i + kl;
                for (int j = start; j < stop; ++j) {
                    y[iy] += temp * A[k+j][i];
                    iy += incy;
                }
                jx += incx;
                if (i > ku)
                    ky += incy;
            }
        }
    } else {
        int jy = ky;
        if (incx == 1) {
            for (int i = 0; i < n; ++i) {
                temp = 0;
                k = kup1 - i;
                start = 0 > i - ku ? 0 : i - ku;
                stop = m < i + kl ? m : i + kl;
                for (int j = start; j < stop; ++j) {
                    temp += A[k+j][i]*x[j];
                }
                y[jy] += alpha * temp;
                jy += incy;
            }
        } else {
            int ix;
            for (int i = 0; i < n; ++i) {
                temp = 0;
                ix = kx;
                k = kup1 - i;
                start = 0 > i - ku ? 0 : i - ku;
                stop = m < i + kl ? m : i + kl;
                for (int j = start; j < stop; ++j) {
                    temp += A[k+j][i] * x[ix];
                    ix += incx;
                }
                y[jy] += alpha * temp;
                jy += incy;
                if (i > ku)
                    kx += incx;
            }
        }
    }

}



void
dgbmv(char trans, int m, int n, int kl, int ku, double alpha, double **A, int lda,
      double *x, int incx, double beta, double *y, int incy) {

    int info = _validate_gbmv_inputs((char)toupper(trans), m,n,kl,ku,lda,incx,incy);

    if (info != 0){
        // error
        return;
    }

    if (m == 0 || n == 0 || (alpha == 0 && beta == 1))
        return;

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
                    y[iy] = beta * y[iy];
                    iy += incy;
                }
            }
        }
    }

    if (alpha == 0)
        return;

    int kup1 = ku + 1;
    double temp;
    int k, start, stop;
    if (trans == 'N') {
        int jx = kx;
        if (incy == 1) {
            for (int i = 0; i < n; ++i) {
                temp = alpha * x[jx];
                k = kup1 - i;
                start = 0 > i - ku ? 0 : i - ku;
                stop = m < i + kl ? m : i + kl;
                for (int j = start; j < stop; ++j) {
                    y[j] += temp * A[k+j][i];
                }
                jx += incx;
            }
        } else {
            int iy;
            for (int i = 0; i < n; ++i) {
                temp = alpha * x[jx];
                iy = ky;
                k = kup1 - i;
                start = 0 > i - ku ? 0 : i - ku;
                stop = m < i + kl ? m : i + kl;
                for (int j = start; j < stop; ++j) {
                    y[iy] += temp * A[k+j][i];
                    iy += incy;
                }
                jx += incx;
                if (i > ku)
                    ky += incy;
            }
        }
    } else {
        int jy = ky;
        if (incx == 1) {
            for (int i = 0; i < n; ++i) {
                temp = 0;
                k = kup1 - i;
                start = 0 > i - ku ? 0 : i - ku;
                stop = m < i + kl ? m : i + kl;
                for (int j = start; j < stop; ++j) {
                    temp += A[k+j][i]*x[j];
                }
                y[jy] += alpha * temp;
                jy += incy;
            }
        } else {
            int ix;
            for (int i = 0; i < n; ++i) {
                temp = 0;
                ix = kx;
                k = kup1 - i;
                start = 0 > i - ku ? 0 : i - ku;
                stop = m < i + kl ? m : i + kl;
                for (int j = start; j < stop; ++j) {
                    temp += A[k+j][i] * x[ix];
                    ix += incx;
                }
                y[jy] += alpha * temp;
                jy += incy;
                if (i > ku)
                    kx += incx;
            }
        }
    }
}

void cgbmv(char trans, int m, int n, int kl, int ku, complex float alpha, complex float **A, int lda, complex float *x,
           int incx, complex float beta, complex float *y, int incy) {
    int info = _validate_gbmv_inputs((char)toupper(trans), m,n,kl,ku,lda,incx,incy);

    if (info != 0){
        // error
        return;
    }

    if (m == 0 || n == 0 || (alpha == 0 && beta == 1))
        return;


    int noconj = (char) toupper(trans) == 'T';
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
                    y[iy] = beta * y[iy];
                    iy += incy;
                }
            }
        }
    }

    if (alpha == 0)
        return;

    int kup1 = ku + 1;
    complex float temp;
    int k, start, stop;
    if (trans == 'N') {
        int jx = kx;
        if (incy == 1) {
            for (int i = 0; i < n; ++i) {
                temp = alpha * x[jx];
                k = kup1 - i;
                start = 0 > i - ku ? 0 : i - ku;
                stop = m < i + kl ? m : i + kl;
                for (int j = start; j < stop; ++j) {
                    y[j] += temp * A[k+j][i];
                }
                jx += incx;
            }
        } else {
            int iy;
            for (int i = 0; i < n; ++i) {
                temp = alpha * x[jx];
                iy = ky;
                k = kup1 - i;
                start = 0 > i - ku ? 0 : i - ku;
                stop = m < i + kl ? m : i + kl;
                for (int j = start; j < stop; ++j) {
                    y[iy] += temp * A[k+j][i];
                    iy += incy;
                }
                jx += incx;
                if (i > ku)
                    ky += incy;
            }
        }
    } else {
        int jy = ky;
        if (incx == 1) {
            for (int i = 0; i < n; ++i) {
                temp = 0;
                k = kup1 - i;
                start = 0 > i - ku ? 0 : i - ku;
                stop = m < i + kl ? m : i + kl;
                if (noconj) {
                    for (int j = start; j < stop; ++j) {
                        temp += A[k+j][i]*x[j];
                    }
                } else {
                    for (int j = start; j < stop; ++j) {
                        temp += conjf(A[k+j][i])*x[j];
                    }
                }
                y[jy] += alpha * temp;
                jy += incy;
            }
        } else {
            int ix;
            for (int i = 0; i < n; ++i) {
                temp = 0;
                ix = kx;
                k = kup1 - i;
                start = 0 > i - ku ? 0 : i - ku;
                stop = m < i + kl ? m : i + kl;
                if (noconj) {
                    for (int j = start; j < stop; ++j) {
                        temp += A[k+j][i] * x[ix];
                        ix += incx;
                    }
                } else {
                    for (int j = start; j < stop; ++j) {
                        temp += conjf(A[k+j][i]) * x[ix];
                        ix += incx;
                    }
                }
                y[jy] += alpha * temp;
                jy += incy;
                if (i > ku)
                    kx += incx;
            }
        }
    }
}

void
zgbmv(char trans, int m, int n, int kl, int ku, complex double alpha, complex double **A, int lda, complex double *x,
      int incx, complex double beta, complex double *y, int incy) {
    int info = _validate_gbmv_inputs((char)toupper(trans), m,n,kl,ku,lda,incx,incy);

    if (info != 0){
        // error
        return;
    }

    if (m == 0 || n == 0 || (alpha == 0 && beta == 1))
        return;


    int noconj = (char) toupper(trans) == 'T';
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
                    y[iy] = beta * y[iy];
                    iy += incy;
                }
            }
        }
    }

    if (alpha == 0)
        return;

    int kup1 = ku + 1;
    complex double temp;
    int k, start, stop;
    if (trans == 'N') {
        int jx = kx;
        if (incy == 1) {
            for (int i = 0; i < n; ++i) {
                temp = alpha * x[jx];
                k = kup1 - i;
                start = 0 > i - ku ? 0 : i - ku;
                stop = m < i + kl ? m : i + kl;
                for (int j = start; j < stop; ++j) {
                    y[j] += temp * A[k+j][i];
                }
                jx += incx;
            }
        } else {
            int iy;
            for (int i = 0; i < n; ++i) {
                temp = alpha * x[jx];
                iy = ky;
                k = kup1 - i;
                start = 0 > i - ku ? 0 : i - ku;
                stop = m < i + kl ? m : i + kl;
                for (int j = start; j < stop; ++j) {
                    y[iy] += temp * A[k+j][i];
                    iy += incy;
                }
                jx += incx;
                if (i > ku)
                    ky += incy;
            }
        }
    } else {
        int jy = ky;
        if (incx == 1) {
            for (int i = 0; i < n; ++i) {
                temp = 0;
                k = kup1 - i;
                start = 0 > i - ku ? 0 : i - ku;
                stop = m < i + kl ? m : i + kl;
                if (noconj) {
                    for (int j = start; j < stop; ++j) {
                        temp += A[k+j][i]*x[j];
                    }
                } else {
                    for (int j = start; j < stop; ++j) {
                        temp += conj(A[k+j][i])*x[j];
                    }
                }
                y[jy] += alpha * temp;
                jy += incy;
            }
        } else {
            int ix;
            for (int i = 0; i < n; ++i) {
                temp = 0;
                ix = kx;
                k = kup1 - i;
                start = 0 > i - ku ? 0 : i - ku;
                stop = m < i + kl ? m : i + kl;
                if (noconj) {
                    for (int j = start; j < stop; ++j) {
                        temp += A[k+j][i] * x[ix];
                        ix += incx;
                    }
                } else {
                    for (int j = start; j < stop; ++j) {
                        temp += conj(A[k+j][i]) * x[ix];
                        ix += incx;
                    }
                }
                y[jy] += alpha * temp;
                jy += incy;
                if (i > ku)
                    kx += incx;
            }
        }
    }
}

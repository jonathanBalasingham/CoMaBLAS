//
// Created by jon on 4/4/21.
//

#include <stdbool.h>
#include <ctype.h>
#include "coma_gemm.h"

int _validate_gemm_inputs(char transa, char transb, int m, int n, int k, int lda, int ldb, int ldc, bool nota,
                          bool notb, int nrowa, int nrowb) {
    int info = 0;
    if (!nota && transa != 'C' && !(transa != 'T'))
        info = 1;
    else if (!notb && transb != 'C' && !(transb != 'T'))
        info = 2;
    else if (m < 0)
        info = 3;
    else if (n < 0)
        info = 4;
    else if (k < 0)
        info = 5;
    else if (lda < (1 < nrowa ? nrowa : 1))
        info = 8;
    else if (ldb < (1 < nrowb ? nrowb : 1))
        info = 10;
    else if (ldc < (1 < m ? m : 1))
        info = 13;

    return info;
}

void
sgemm(char transa, char transb, int m, int n, int k, float alpha, float **A, int lda, float **B, int ldb, float beta,
      float **C, int ldc) {
    transa = (char) toupper(transa);
    transb = (char) toupper(transb);
    bool nota = transa == 'N';
    bool notb = transb == 'N';

    int nrowa, nrowb;
    if (nota)
        nrowa = m;
    else
        nrowa = k;
    if (notb)
        nrowb = k;
    else
        nrowb = n;

    int info = _validate_gemm_inputs(transa, transb, m, n, k, lda, ldb, ldc, nota, notb, nrowa, nrowb);
    if (info != 0)
        return;

    if (m == 0 || n == 0 || (alpha == 0 || (k == 0 && beta == 1)))
        return;

    if (alpha == 0) {
        if (beta == 0) {
            for (int j = 0; j < n; ++j) {
                for (int i = 0; i < m; ++i) {
                    C[i][j] = 0;
                }
            }
        } else {
            for (int j = 0; j < n; ++j) {
                for (int i = 0; i < m; ++i) {
                    C[i][j] = beta * C[i][j];
                }
            }
        }
    }
    float temp;
    if (notb) {
        if (nota) {
            for (int j = 0; j < n; ++j) {
                if (beta == 0) {
                    for (int i = 0; i < m; ++i) {
                        C[i][j] = 0;
                    }
                } else if (beta != 1) {
                    for (int i = 0; i < m; ++i) {
                        C[i][j] = beta * C[i][j];
                    }
                }
                for (int l = 0; l < k; ++l) {
                    temp = alpha * B[l][j];
                    for (int i = 0; i < m; ++i) {
                        C[i][j] += temp * A[i][l];
                    }
                }
            }
        } else {
            for (int j = 0; j < n; ++j) {
                for (int i = 0; i < m; ++i) {
                    temp = 0;
                    for (int l = 0; l < k; ++l) {
                        temp += A[l][i] * B[l][j];
                    }
                    if (beta == 0) {
                        C[i][j] = alpha * temp;
                    } else {
                        C[i][j] = alpha * temp + beta * C[i][j];
                    }
                }
            }
        }
    } else {
        if (nota) {
            for (int j = 0; j < n; ++j) {
                if (beta == 0) {
                    for (int i = 0; i < m; ++i) {
                        C[i][j] = 0;
                    }
                } else if (beta != 1) {
                    for (int i = 0; i < m; ++i) {
                        C[i][j] = beta * C[i][j];
                    }
                }
                for (int l = 0; l < k; ++l) {
                    temp = alpha * B[j][l];
                    for (int i = 0; i < m; ++i) {
                        C[i][j] += temp * A[i][l];
                    }
                }
            }
        } else {
            for (int j = 0; j < n; ++j) {
                for (int i = 0; i < m; ++i) {
                    temp = 0;
                    for (int l = 0; l < k; ++l) {
                        temp += A[l][i] * B[j][l];
                    }
                    if (beta == 0) {
                        C[i][j] = alpha * temp;
                    } else {
                        C[i][j] = alpha * temp + beta * C[i][j];
                    }
                }
            }
        }
    }
}

void dgemm(char transa, char transb, int m, int n, int k, double alpha, double **A, int lda, double **B, int ldb,
           double beta, double **C, int ldc) {
    transa = (char) toupper(transa);
    transb = (char) toupper(transb);
    bool nota = transa == 'N';
    bool notb = transb == 'N';

    int nrowa, nrowb;
    if (nota)
        nrowa = m;
    else
        nrowa = k;
    if (notb)
        nrowb = k;
    else
        nrowb = n;

    int info = _validate_gemm_inputs(transa, transb, m, n, k, lda, ldb, ldc, nota, notb, nrowa, nrowb);
    if (info != 0)
        return;

    if (m == 0 || n == 0 || (alpha == 0 || (k == 0 && beta == 1)))
        return;

    if (alpha == 0) {
        if (beta == 0) {
            for (int j = 0; j < n; ++j) {
                for (int i = 0; i < m; ++i) {
                    C[i][j] = 0;
                }
            }
        } else {
            for (int j = 0; j < n; ++j) {
                for (int i = 0; i < m; ++i) {
                    C[i][j] = beta * C[i][j];
                }
            }
        }
    }
    double temp;
    if (notb) {
        if (nota) {
            for (int j = 0; j < n; ++j) {
                if (beta == 0) {
                    for (int i = 0; i < m; ++i) {
                        C[i][j] = 0;
                    }
                } else if (beta != 1) {
                    for (int i = 0; i < m; ++i) {
                        C[i][j] = beta * C[i][j];
                    }
                }
                for (int l = 0; l < k; ++l) {
                    temp = alpha * B[l][j];
                    for (int i = 0; i < m; ++i) {
                        C[i][j] += temp * A[i][l];
                    }
                }
            }
        } else {
            for (int j = 0; j < n; ++j) {
                for (int i = 0; i < m; ++i) {
                    temp = 0;
                    for (int l = 0; l < k; ++l) {
                        temp += A[l][i] * B[l][j];
                    }
                    if (beta == 0) {
                        C[i][j] = alpha * temp;
                    } else {
                        C[i][j] = alpha * temp + beta * C[i][j];
                    }
                }
            }
        }
    } else {
        if (nota) {
            for (int j = 0; j < n; ++j) {
                if (beta == 0) {
                    for (int i = 0; i < m; ++i) {
                        C[i][j] = 0;
                    }
                } else if (beta != 1) {
                    for (int i = 0; i < m; ++i) {
                        C[i][j] = beta * C[i][j];
                    }
                }
                for (int l = 0; l < k; ++l) {
                    temp = alpha * B[j][l];
                    for (int i = 0; i < m; ++i) {
                        C[i][j] += temp * A[i][l];
                    }
                }
            }
        } else {
            for (int j = 0; j < n; ++j) {
                for (int i = 0; i < m; ++i) {
                    temp = 0;
                    for (int l = 0; l < k; ++l) {
                        temp += A[l][i] * B[j][l];
                    }
                    if (beta == 0) {
                        C[i][j] = alpha * temp;
                    } else {
                        C[i][j] = alpha * temp + beta * C[i][j];
                    }
                }
            }
        }
    }
}

void
cgemm(char transa, char transb, int m, int n, int k, complex float alpha, complex float **A, int lda, complex float **B,
      int ldb, complex float beta, complex float **C, int ldc) {
    transa = (char) toupper(transa);
    transb = (char) toupper(transb);
    bool nota = transa == 'N';
    bool notb = transb == 'N';

    int nrowa, nrowb;
    if (nota)
        nrowa = m;
    else
        nrowa = k;
    if (notb)
        nrowb = k;
    else
        nrowb = n;

    int info = _validate_gemm_inputs(transa, transb, m, n, k, lda, ldb, ldc, nota, notb, nrowa, nrowb);
    if (info != 0)
        return;

    if (m == 0 || n == 0 || (alpha == 0 || (k == 0 && beta == 1)))
        return;

    if (alpha == 0) {
        if (beta == 0) {
            for (int j = 0; j < n; ++j) {
                for (int i = 0; i < m; ++i) {
                    C[i][j] = 0;
                }
            }
        } else {
            for (int j = 0; j < n; ++j) {
                for (int i = 0; i < m; ++i) {
                    C[i][j] = beta * C[i][j];
                }
            }
        }
    }
    complex float temp;
    if (notb) {
        if (nota) {
            for (int j = 0; j < n; ++j) {
                if (beta == 0) {
                    for (int i = 0; i < m; ++i) {
                        C[i][j] = 0;
                    }
                } else if (beta != 1) {
                    for (int i = 0; i < m; ++i) {
                        C[i][j] = beta * C[i][j];
                    }
                }
                for (int l = 0; l < k; ++l) {
                    temp = alpha * B[l][j];
                    for (int i = 0; i < m; ++i) {
                        C[i][j] += temp * A[i][l];
                    }
                }
            }
        } else {
            for (int j = 0; j < n; ++j) {
                for (int i = 0; i < m; ++i) {
                    temp = 0;
                    for (int l = 0; l < k; ++l) {
                        temp += A[l][i] * B[l][j];
                    }
                    if (beta == 0) {
                        C[i][j] = alpha * temp;
                    } else {
                        C[i][j] = alpha * temp + beta * C[i][j];
                    }
                }
            }
        }
    } else {
        if (nota) {
            for (int j = 0; j < n; ++j) {
                if (beta == 0) {
                    for (int i = 0; i < m; ++i) {
                        C[i][j] = 0;
                    }
                } else if (beta != 1) {
                    for (int i = 0; i < m; ++i) {
                        C[i][j] = beta * C[i][j];
                    }
                }
                for (int l = 0; l < k; ++l) {
                    temp = alpha * B[j][l];
                    for (int i = 0; i < m; ++i) {
                        C[i][j] += temp * A[i][l];
                    }
                }
            }
        } else {
            for (int j = 0; j < n; ++j) {
                for (int i = 0; i < m; ++i) {
                    temp = 0;
                    for (int l = 0; l < k; ++l) {
                        temp += A[l][i] * B[j][l];
                    }
                    if (beta == 0) {
                        C[i][j] = alpha * temp;
                    } else {
                        C[i][j] = alpha * temp + beta * C[i][j];
                    }
                }
            }
        }
    }
}

void zgemm(char transa, char transb, int m, int n, int k, complex double alpha, complex double **A, int lda,
           complex double **B, int ldb, double beta, complex double **C, int ldc) {
    transa = (char) toupper(transa);
    transb = (char) toupper(transb);
    bool nota = transa == 'N';
    bool notb = transb == 'N';

    int nrowa, nrowb;
    if (nota)
        nrowa = m;
    else
        nrowa = k;
    if (notb)
        nrowb = k;
    else
        nrowb = n;

    int info = _validate_gemm_inputs(transa, transb, m, n, k, lda, ldb, ldc, nota, notb, nrowa, nrowb);
    if (info != 0)
        return;

    if (m == 0 || n == 0 || (alpha == 0 || (k == 0 && beta == 1)))
        return;

    if (alpha == 0) {
        if (beta == 0) {
            for (int j = 0; j < n; ++j) {
                for (int i = 0; i < m; ++i) {
                    C[i][j] = 0;
                }
            }
        } else {
            for (int j = 0; j < n; ++j) {
                for (int i = 0; i < m; ++i) {
                    C[i][j] = beta * C[i][j];
                }
            }
        }
    }
    complex double temp;
    if (notb) {
        if (nota) {
            for (int j = 0; j < n; ++j) {
                if (beta == 0) {
                    for (int i = 0; i < m; ++i) {
                        C[i][j] = 0;
                    }
                } else if (beta != 1) {
                    for (int i = 0; i < m; ++i) {
                        C[i][j] = beta * C[i][j];
                    }
                }
                for (int l = 0; l < k; ++l) {
                    temp = alpha * B[l][j];
                    for (int i = 0; i < m; ++i) {
                        C[i][j] += temp * A[i][l];
                    }
                }
            }
        } else {
            for (int j = 0; j < n; ++j) {
                for (int i = 0; i < m; ++i) {
                    temp = 0;
                    for (int l = 0; l < k; ++l) {
                        temp += A[l][i] * B[l][j];
                    }
                    if (beta == 0) {
                        C[i][j] = alpha * temp;
                    } else {
                        C[i][j] = alpha * temp + beta * C[i][j];
                    }
                }
            }
        }
    } else {
        if (nota) {
            for (int j = 0; j < n; ++j) {
                if (beta == 0) {
                    for (int i = 0; i < m; ++i) {
                        C[i][j] = 0;
                    }
                } else if (beta != 1) {
                    for (int i = 0; i < m; ++i) {
                        C[i][j] = beta * C[i][j];
                    }
                }
                for (int l = 0; l < k; ++l) {
                    temp = alpha * B[j][l];
                    for (int i = 0; i < m; ++i) {
                        C[i][j] += temp * A[i][l];
                    }
                }
            }
        } else {
            for (int j = 0; j < n; ++j) {
                for (int i = 0; i < m; ++i) {
                    temp = 0;
                    for (int l = 0; l < k; ++l) {
                        temp += A[l][i] * B[j][l];
                    }
                    if (beta == 0) {
                        C[i][j] = alpha * temp;
                    } else {
                        C[i][j] = alpha * temp + beta * C[i][j];
                    }
                }
            }
        }
    }
}

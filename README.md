# CoMaBLAS
A pure C implementation of BLAS, based on the LAPACK Fortran implementation.

## Development 

### Level 1
- SROTG - setup Givens rotation [IMPLEMENTED]
- SROTMG - setup modified Givens rotation [IMPLEMENTED]
- SROT - apply Givens rotation [IMPLEMENTED]
- SROTM - apply modified Givens rotation [IMPLEMENTED]
- SSWAP - swap x and y [TESTED]
- SSCAL - x = a*x [TESTED]
- SCOPY - copy x into y [TESTED]
- SAXPY - y = a*x + y [TESTED]
- SDOT - dot product [TESTED]
- SDSDOT - dot product with extended precision accumulation [TESTED]
- SNRM2 - Euclidean norm [TESTED]
- SCNRM2- Euclidean norm [TESTED]
- SASUM - sum of absolute values [TESTED]
- ISAMAX - index of max abs value [TESTED]
- DROTG - setup Givens rotation [IMPLEMENTED]
- DROTMG - setup modified Givens rotation [IMPLEMENTED]
- DROT - apply Givens rotation [IMPLEMENTED]
- DROTM - apply modified Givens rotation [IMPLEMENTED]
- DSWAP - swap x and y [TESTED]
- DSCAL - x = a*x [TESTED]
- DCOPY - copy x into y [TESTED]
- DAXPY - y = a*x + y [TESTED]
- DDOT - dot product [TESTED]
- DSDOT - dot product with extended precision accumulation [TESTED]
- DNRM2 - Euclidean norm [TESTED]
- DZNRM2 - Euclidean norm [TESTED]
- DASUM - sum of absolute values [TESTED]
- IDAMAX - index of max abs value [TESTED]
- CROTG - setup Givens rotation [IMPLEMENTED]
- CSROT - apply Givens rotation [IMPLEMENTED]
- CSWAP - swap x and y [TESTED]
- CSCAL - x = a*x [TESTED]
- CSSCAL - x = a*x [TESTED]
- CCOPY - copy x into y [TESTED]
- CAXPY - y = a*x + y [TESTED]
- CDOTU - dot product [TESTED]
- CDOTC - dot product, conjugating the first vector [TESTED]
- SCASUM - sum of absolute values [TESTED]
- ICAMAX - index of max abs value [TESTED]
- ZROTG - setup Givens rotation [IMPLEMENTED]
- ZDROTF - apply Givens rotation [IMPLEMENTED]
- ZSWAP - swap x and y [TESTED]
- ZSCAL - x = a*x [TESTED]
- ZDSCAL - x = a*x [TESTED]
- ZCOPY - copy x into y [TESTED]
- ZAXPY - y = a*x + y [TESTED]
- ZDOTU - dot product [TESTED]
- ZDOTC - dot product, conjugating the first vector [TESTED]
- DZASUM - sum of absolute values [TESTED]
- IZAMAX - index of max abs value [TESTED]

### Level 2
- SGEMV - matrix vector multiply [IMPLEMENTED]
- SGBMV - banded matrix vector multiply [IMPLEMENTED]
- SSYMV - symmetric matrix vector multiply [IMPLEMENTED]
- SSBMV - symmetric banded matrix vector multiply
- SSPMV - symmetric packed matrix vector multiply
- STRMV - triangular matrix vector multiply
- STBMV - triangular banded matrix vector multiply
- STPMV - triangular packed matrix vector multiply
- STRSV - solving triangular matrix problems
- STBSV - solving triangular banded matrix problems
- STPSV - solving triangular packed matrix problems
- SGER - performs the rank 1 operation A := alpha*x*y' + A
- SSYR - performs the symmetric rank 1 operation A := alpha*x*x' + A
- SSPR - symmetric packed rank 1 operation A := alpha*x*x' + A
- SSYR2 - performs the symmetric rank 2 operation, A := alpha*x*y' + alpha*y*x' + A
- SSPR2 - performs the symmetric packed rank 2 operation, A := alpha*x*y' + alpha*y*x' + A
- DGEMV - matrix vector multiply [IMPLEMENTED]
- DGBMV - banded matrix vector multiply [IMPLEMENTED]
- DSYMV - symmetric matrix vector multiply
- DSBMV - symmetric banded matrix vector multiply
- DSPMV - symmetric packed matrix vector multiply
- DTRMV - triangular matrix vector multiply
- DTBMV - triangular banded matrix vector multiply
- DTPMV - triangular packed matrix vector multiply
- DTRSV - solving triangular matrix problems
- DTBSV - solving triangular banded matrix problems
- DTPSV - solving triangular packed matrix problems
- DGER - performs the rank 1 operation A := alpha*x*y' + A
- DSYR - performs the symmetric rank 1 operation A := alpha*x*x' + A
- DSPR - symmetric packed rank 1 operation A := alpha*x*x' + A
- DSYR2 - performs the symmetric rank 2 operation, A := alpha*x*y' + alpha*y*x' + A
- DSPR2 - performs the symmetric packed rank 2 operation, A := alpha*x*y' + alpha*y*x' + A
- CGEMV - matrix vector multiply [IMPLEMENTED]
- CGBMV - banded matrix vector multiply [IMPLEMENTED]
- CHEMV - hermitian matrix vector multiply
- CHBMV - hermitian banded matrix vector multiply
- CHPMV - hermitian packed matrix vector multiply
- CTRMV - triangular matrix vector multiply
- CTBMV - triangular banded matrix vector multiply
- CTPMV - triangular packed matrix vector multiply
- CTRSV - solving triangular matrix problems
- CTBSV - solving triangular banded matrix problems
- CTPSV - solving triangular packed matrix problems
- CGERU - performs the rank 1 operation A := alpha*x*y' + A
- CGERC - performs the rank 1 operation A := alpha*x*conjg( y' ) + A
- CHER - hermitian rank 1 operation A := alpha*x*conjg(x') + A
- CHPR - hermitian packed rank 1 operation A := alpha*x*conjg( x' ) + A
- CHER2 - hermitian rank 2 operation
- CHPR2 - hermitian packed rank 2 operation
- ZGEMV - matrix vector multiply [IMPLEMENTED]
- ZGBMV - banded matrix vector multiply [IMPLEMENTED]
- ZHEMV - hermitian matrix vector multiply
- ZHBMV - hermitian banded matrix vector multiply
- ZHPMV - hermitian packed matrix vector multiply
- ZTRMV - triangular matrix vector multiply
- ZTBMV - triangular banded matrix vector multiply
- ZTPMV - triangular packed matrix vector multiply
- ZTRSV - solving triangular matrix problems
- ZTBSV - solving triangular banded matrix problems
- ZTPSV - solving triangular packed matrix problems
- ZGERU - performs the rank 1 operation A := alpha*x*y' + A
- ZGERC - performs the rank 1 operation A := alpha*x*conjg( y' ) + A
- ZHER - hermitian rank 1 operation A := alpha*x*conjg(x') + A
- ZHPR - hermitian packed rank 1 operation A := alpha*x*conjg( x' ) + A
- ZHER2 - hermitian rank 2 operation
- ZHPR2 - hermitian packed rank 2 operation
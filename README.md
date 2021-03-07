# CoMaBLAS
A pure C implementation of BLAS, based on the LAPACK Fortran implementation.

## Development 

### Level 1
- SROTG - setup Givens rotation
- SROTMG - setup modified Givens rotation
- SROT - apply Givens rotation
- SROTM - apply modified Givens rotation
- SSWAP - swap x and y [TESTED]
- SSCAL - x = a*x [TESTED]
- SCOPY - copy x into y [TESTED]
- SAXPY - y = a*x + y [TESTED]
- SDOT - dot product [TESTED]
- SDSDOT - dot product with extended precision accumulation
- SNRM2 - Euclidean norm [TESTED]
- SCNRM2- Euclidean norm [TESTED]
- SASUM - sum of absolute values
- ISAMAX - index of max abs value [TESTED]



- DROTG - setup Givens rotation
- DROTMG - setup modified Givens rotation
- DROT - apply Givens rotation
- DROTM - apply modified Givens rotation
- DSWAP - swap x and y [TESTED]
- DSCAL - x = a*x [TESTED]
- DCOPY - copy x into y [TESTED]
- DAXPY - y = a*x + y [TESTED]
- DDOT - dot product [TESTED]
- DSDOT - dot product with extended precision accumulation
- DNRM2 - Euclidean norm [TESTED]
- DZNRM2 - Euclidean norm [TESTED]
- DASUM - sum of absolute values
- IDAMAX - index of max abs value [TESTED]
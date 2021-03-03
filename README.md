# CoMaBLAS
A pure C implementation of BLAS, based on the LAPACK Fortran implementation.

## Development 

### Level 1
- SROTG - setup Givens rotation
- SROTMG - setup modified Givens rotation
- SROT - apply Givens rotation
- SROTM - apply modified Givens rotation
- SSWAP - swap x and y
- SSCAL - x = a*x
- SCOPY - copy x into y
- SAXPY - y = a*x + y
- SDOT - dot product [TESTED]
- SDSDOT - dot product with extended precision accumulation
- SNRM2 - Euclidean norm [TESTED]
- SCNRM2- Euclidean norm [TESTED]
- SASUM - sum of absolute values
- ISAMAX - index of max abs value



- DROTG - setup Givens rotation
- DROTMG - setup modified Givens rotation
- DROT - apply Givens rotation
- DROTM - apply modified Givens rotation
- DSWAP - swap x and y
- DSCAL - x = a*x
- DCOPY - copy x into y [TESTED]
- DAXPY - y = a*x + y
- DDOT - dot product [TESTED]
- DSDOT - dot product with extended precision accumulation
- DNRM2 - Euclidean norm [TESTED]
- DZNRM2 - Euclidean norm [TESTED]
- DASUM - sum of absolute values
- IDAMAX - index of max abs value
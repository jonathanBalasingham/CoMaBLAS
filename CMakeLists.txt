cmake_minimum_required(VERSION 3.17)
project(CoMaBLAS C)

set(CMAKE_C_STANDARD 11)

add_library(CoMaBLAS SHARED library.c library.h blas/coma_rotg.c blas/coma_rotg.h blas/coma_rotmg.c blas/coma_rotmg.h blas/coma_rot.c blas/coma_rot.h blas/coma_rotm.c blas/coma_rotm.h blas/coma_swap.c blas/coma_swap.h blas/coma_scal.c blas/coma_scal.h blas/coma_copy.c blas/coma_copy.h blas/coma_axpy.c blas/coma_axpy.h blas/coma_dot.c blas/coma_dot.h blas/coma_asum.c blas/coma_asum.h blas/coma_iamax.c blas/coma_iamax.h blas/coma_dsdot.c blas/coma_dsdot.h blas/coma_nrm2.c blas/coma_nrm2.h)
add_executable(test tests/dot_tests.c munit/munit.c munit/munit.h tests/copy_tests.c blas/coma_asum.c blas/coma_asum.h)
target_link_libraries(test CoMaBLAS m)
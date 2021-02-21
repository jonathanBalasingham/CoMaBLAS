cmake_minimum_required(VERSION 3.17)
project(CoMaBLAS C)

set(CMAKE_C_STANDARD 11)

add_library(CoMaBLAS SHARED library.c library.h blas/coma_srotg.c blas/coma_srotg.h blas/coma_srotmg.c blas/coma_srotmg.h blas/coma_rot.c blas/coma_rot.h blas/coma_srotm.c blas/coma_srotm.h blas/coma_swap.c blas/coma_swap.h blas/coma_scal.c blas/coma_scal.h blas/coma_copy.c blas/coma_copy.h blas/coma_axpy.c blas/coma_axpy.h blas/coma_dot.c blas/coma_dot.h)
add_executable(test tests/dot_tests.c)
target_link_libraries(test CoMaBLAS)
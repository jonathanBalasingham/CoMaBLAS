//
// Created by jon on 3/1/21.
//
#define MUNIT_ENABLE_ASSERT_ALIASES
#include "../munit/munit.h"
#include "../blas/coma_dot.h"
#include "../blas/coma_nrm2.h"
#include "../blas/coma_copy.h"
#include "../blas/coma_scal.h"
#include "../blas/coma_axpy.h"
#include <math.h>
/*
 * Dot Products
 */
MunitResult sdot_under5(const MunitParameter params[], void* user_data_or_fixture) {
    float t1[3] = {1,2,3};
    float t2[3] = {1,2,3};
    int n = 3;

    float answer = sdot(n, (float *) &t1, 1, (float *) &t2, 1);
    munit_assert_float(14.0, ==, answer);
    return MUNIT_OK;
}

MunitResult sdot_over5(const MunitParameter params[], void* user_data_or_fixture) {
    float t3[7] = {1,1,1,1,1,1,1};
    float t4[7] = {1,1,1,1,1,1,1};
    int n2 = 7;

    float answer2 = sdot(n2, (float *) &t3, 1, (float *) &t4, 1);
    munit_assert_float(7.0, ==, answer2);
    return MUNIT_OK;
}

MunitResult ddot_under5(const MunitParameter params[], void* user_data_or_fixture) {
    int n = 3;
    double t5[3] = {1,2,3};
    double t6[3] = {1,2,3};

    double answer3 = ddot(n, (double *) &t5, 1, (double *) &t6, 1);
    munit_assert_double(14.0, ==, answer3);
    return MUNIT_OK;
}

MunitResult ddot_over5(const MunitParameter params[], void* user_data_or_fixture) {
    int n = 7;
    double t7[7] = {1,1,1,1,1,1,1};
    double t8[7] = {1,1,1,1,1,1,1};

    double answer4 = ddot(n, (double *) &t7, 1, (double *) &t8, 1);
    munit_assert_double(7.0, ==, answer4);
    return MUNIT_OK;
}

/*
 * Euclidean Norm
 */
MunitResult test_snrm2(const MunitParameter params[], void* user_data_or_fixture) {
    unsigned int n = 3;
    float t1[3] = {1,2,3};
    munit_assert_float(3.74165, ==, snrm2(n, (float *) &t1, 1));
}

// TODO: What precision?
MunitResult test_dnrm2(const MunitParameter params[], void* user_data_or_fixture) {
    unsigned int n = 3;
    float t1[3] = {1,2,3};
    munit_assert_double_equal(3.7416573867739, dnrm2(n, (double *) &t1, 1), 10);
}

MunitResult test_scnrm2(const MunitParameter params[], void* user_data_or_fixture) {
    unsigned int n = 3;
    complex float t1[3] = {1 + 1*I, 2 + 2*I, 3 + 3*I};
    assert_float((float) 2*sqrt(7), ==, scnrm2(n, (complex float *) &t1, 1));
}

// TODO: precision needed?
MunitResult test_dznrm2(const MunitParameter params[], void* user_data_or_fixture) {
    unsigned int n = 3;
    complex float t1[3] = {1 + 1*I, 2 + 2*I, 3 + 3*I};
    assert_double_equal((double) 2*sqrt(7), scnrm2(n, (complex float *) &t1, 1), 10);
}

MunitResult test_scopy(const MunitParameter params[], void* user_data_or_fixture) {
    unsigned int n = 6;
    float t1[6] = {1,2,3,4,5,6};
    float* t2;
    t2 = malloc(sizeof(t1));
    scopy(n, &t1, 1, t2, 1);
    for (int i = 0; i < n; ++i) {
        assert_float(t1[i], ==, t2[i]);
    }

    free(t2);
}

MunitResult test_dcopy(const MunitParameter params[], void* user_data_or_fixture) {
    unsigned int n = 6;
    double t1[6] = {1,2,3,4,5,6};
    double* t2;
    t2 = malloc(sizeof(t1));
    dcopy(n, (double *) &t1, 1, t2, 1);
    for (int i = 0; i < n; ++i) {
        assert_double(t1[i], ==, t2[i]);
    }

    free(t2);
}

// TODO: Munit doesnt have complex number comparisons
MunitResult test_ccopy(const MunitParameter params[], void* user_data_or_fixture) {
    unsigned int n = 6;
    complex float t1[6] = {1,2,3,4,5,6};
    complex float* t2;
    t2 = malloc(sizeof(t1));
    scopy(n, &t1, 1, (float *) t2, 1);
    for (int i = 0; i < n; ++i) {
        assert_float(t1[i], ==, t2[i]);
    }

    free(t2);
}

MunitResult test_zcopy(const MunitParameter params[], void* user_data_or_fixture) {
    unsigned int n = 6;
    const complex double t1[6] = {1+1*I,2+2*I,3+3*I,4+4*I,5+5*I,6+6*I};
    complex double* t2;
    t2 = malloc(sizeof(t1));
    zcopy(n, (const complex double *) &t1, 1, t2, 1);
    for (int i = 0; i < n; ++i) {
        assert_double(t1[i], ==, t2[i]);
    }

    free(t2);
}

MunitResult test_sscal(const MunitParameter params[], void* user_data_or_fixture) {
    unsigned int n = 6;
    float t1[6] = {1,2,3,4,5,6};
    float t2[6] = {2,4,6,8,10,12};
    float sa = 2;
    sscal(n, sa, (float *) &t1, 1);

    for (int i = 0; i < n; ++i) {
        assert_float(t1[i], ==, t2[i]);
    }
}

MunitResult test_dscal(const MunitParameter params[], void* user_data_or_fixture) {
    unsigned int n = 6;
    double t1[6] = {1,2,3,4,5,6};
    double t2[6] = {2,4,6,8,10,12};
    double sa = 2;
    dscal(n, sa, (double *) &t1, 1);

    for (int i = 0; i < n; ++i) {
        assert_double(t1[i], ==, t2[i]);
    }
}

// TODO: adjust assertion
MunitResult test_cscal(const MunitParameter params[], void* user_data_or_fixture) {
    unsigned int n = 6;
    complex float t1[6] = {1+1*I,2+2*I,3+3*I,4+4*I,5+5*I,6+6*I};
    const complex float t2[6] = {2+2*I,4+4*I,6+6*I,8+8*I,10+10*I,12+12*I};
    complex float sa = 2 + 1*I;
    cscal(n, sa, (complex float *) &t1, 1);

    for (int i = 0; i < n; ++i) {
        assert_float(t1[i], ==, t2[i]);
    }
}

// TODO: adjust assertion
MunitResult test_zscal(const MunitParameter params[], void* user_data_or_fixture) {
    unsigned int n = 6;
    complex double t1[6] = {1+1*I,2+2*I,3+3*I,4+4*I,5+5*I,6+6*I};
    const complex double t2[6] = {2+2*I,4+4*I,6+6*I,8+8*I,10+10*I,12+12*I};
    complex double sa = 2 + 1*I;
    cscal(n, sa, (complex float *) &t1, 1);

    for (int i = 0; i < n; ++i) {
        assert_double(t1[i], ==, t2[i]);
    }
}

MunitResult test_saxpy(const MunitParameter params[], void* user_data_or_fixture) {
    unsigned int n = 6;
    float t1[6] = {1,2,3,4,5,6};
    float t2[6] = {2,4,6,8,10,12};
    float t3[6] = {4,8,12,16,20,24};

    float sa = 2;
    saxpy(n, sa, (float *) &t1, 1, (float *) &t2, 1);

    for (int i = 0; i < n; ++i) {
        assert_float(t1[i], ==, t3[i]);
    }
}

MunitResult test_daxpy(const MunitParameter params[], void* user_data_or_fixture) {
    unsigned int n = 6;
    double t1[6] = {1,2,3,4,5,6};
    double t2[6] = {2,4,6,8,10,12};
    float t3[6] = {4,8,12,16,20,24};

    double sa = 2;
    daxpy(n, sa, (double *) &t1, 1, (double *) &t2, 1);

    for (int i = 0; i < n; ++i) {
        assert_double(t1[i], ==, t3[i]);
    }
}

typedef MunitResult (*munit_test)(const MunitParameter[], void*);

static const munit_test level1_tests[16] = {&sdot_under5, &sdot_over5, &ddot_under5, &ddot_under5,
                                           &test_snrm2, &test_scnrm2, &test_dnrm2, &test_dznrm2,
                                           &test_scopy, &test_dcopy, &test_ccopy, &test_zcopy,
                                           &test_sscal, &test_dscal, &test_cscal, &test_zscal};

static const MunitSuite suite = {
        "/my-tests", /* name */
        &level1_tests, /* tests */
        NULL, /* suites */
        1, /* iterations */
        MUNIT_SUITE_OPTION_NONE /* options */
};

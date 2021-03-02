//
// Created by jon on 3/1/21.
//
#define MUNIT_ENABLE_ASSERT_ALIASES
#include "../munit/munit.h"
#include "../blas/coma_dot.h"
#include "../blas/coma_nrm2.h"
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


typedef MunitResult (*munit_test)(const MunitParameter[], void*);

static const munit_test level1_tests[8] = {sdot_under5, sdot_over5, ddot_under5, ddot_under5,
                                           test_snrm2, test_scnrm2, test_dnrm2, test_dznrm2};

static const MunitSuite suite = {
        "/my-tests", /* name */
        level1_tests, /* tests */
        NULL, /* suites */
        1, /* iterations */
        MUNIT_SUITE_OPTION_NONE /* options */
};

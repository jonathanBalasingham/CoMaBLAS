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
#include "../blas/coma_swap.h"
#include "../blas/coma_iamax.h"
#include "../blas/coma_asum.h"
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
    munit_assert_double_equal(3.741657, snrm2(n, (float *) &t1, 1), 6);
    return MUNIT_OK;
}

// TODO: What precision?
MunitResult test_dnrm2(const MunitParameter params[], void* user_data_or_fixture) {
    unsigned int n = 3;
    double t1[3] = {1,2,3};
    munit_assert_double_equal(3.7416573867739, dnrm2(n, (double *) &t1, 1), 10);
    return MUNIT_OK;
}

MunitResult test_scnrm2(const MunitParameter params[], void* user_data_or_fixture) {
    unsigned int n = 3;
    complex float t1[3] = {1 + 1*I, 2 + 2*I, 3 + 3*I};
    assert_double_equal((float) 2*sqrt(7), scnrm2(n, (complex float *) &t1, 1), 5);
    return MUNIT_OK;
}

// TODO: precision needed?
MunitResult test_dznrm2(const MunitParameter params[], void* user_data_or_fixture) {
    unsigned int n = 3;
    complex double t1[3] = {1 + 1*I, 2 + 2*I, 3 + 3*I};
    assert_double_equal((double) 2*sqrt(7), dznrm2(n, (complex double *) &t1, 1), 5);
    return MUNIT_OK;
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
    return MUNIT_OK;
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
    return MUNIT_OK;
}

MunitResult test_ccopy(const MunitParameter params[], void* user_data_or_fixture) {
    unsigned int n = 6;
    complex float t1[6] = {1 + I,2 + I,3 + I,4 + I,5 + I,6 + I};
    complex float* t2;
    t2 = malloc(n*sizeof(t1[0]));
    ccopy(n, (const complex float *) &t1, 1, t2, 1);
    for (int i = 0; i < n; ++i) {
        assert_float(crealf(t1[i]), ==, crealf(t2[i]));
        assert_float(cimagf(t1[i]), ==, cimagf(t2[i]));
    }

    free(t2);
    return MUNIT_OK;
}

MunitResult test_zcopy(const MunitParameter params[], void* user_data_or_fixture) {
    unsigned int n = 6;
    const complex double t1[6] = {1+1*I,2+2*I,3+3*I,4+4*I,5+5*I,6+6*I};
    complex double* t2;
    t2 = malloc(sizeof(t1));
    zcopy(n, (const complex double *) &t1, 1, t2, 1);
    for (int i = 0; i < n; ++i) {
        assert_double(creal(t1[i]), ==, creal(t2[i]));
        assert_double(cimag(t1[i]), ==, cimag(t2[i]));
    }

    free(t2);
    return MUNIT_OK;
}

MunitResult test_sscal(const MunitParameter params[], void* user_data_or_fixture) {
    unsigned int n = 6;
    float t1[6] = {1,2,3,4,5,6};
    float t2[6] = {2,4,6,8,10,12};
    float sa = 2;
    sscal(n, sa, (float *) &t1, 1);
    assert_float(t1[0], ==, t2[0]);

    return MUNIT_OK;
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

    return MUNIT_OK;
}

MunitResult test_cscal(const MunitParameter params[], void* user_data_or_fixture) {
    unsigned int n = 6;
    complex float t1[6] = {1+1*I,2+2*I,3+3*I,4+4*I,5+5*I,6+6*I};
    const complex float t2[6] = {1+1*I,2+2*I,3+3*I,4+4*I,5+5*I,6+6*I};
    complex float sa = 2 + 1*I;
    cscal(n, sa, (complex float *) &t1, 1);

    for (int i = 0; i < n; ++i) {
        assert_float(crealf(t1[i]), ==, crealf(sa) * crealf(t2[i]) - cimagf(sa) * cimagf(t2[i]));
        assert_float(cimagf(t1[i]), ==, crealf(sa) * cimagf(t2[i]) + cimagf(sa) * crealf(t2[i]));
    }

    return MUNIT_OK;
}

MunitResult test_zscal(const MunitParameter params[], void* user_data_or_fixture) {
    unsigned int n = 6;
    complex double t1[6] = {1+1*I,2+2*I,3+3*I,4+4*I,5+5*I,6+6*I};
    const complex double t2[6] = {1+1*I,2+2*I,3+3*I,4+4*I,5+5*I,6+6*I};
    complex double sa = 2 + 1*I;
    zscal(n, sa, (complex double *) &t1, 1);

    for (int i = 0; i < n; ++i) {
        assert_double(creal(t1[i]), ==, creal(sa) * creal(t2[i]) - cimag(sa) * cimag(t2[i]));
        assert_double(cimag(t1[i]), ==, creal(sa) * cimag(t2[i]) + cimag(sa) * creal(t2[i]));
    }

    return MUNIT_OK;
}

MunitResult test_saxpy(const MunitParameter params[], void* user_data_or_fixture) {
    unsigned int n = 6;
    float t1[6] = {1,2,3,4,5,6};
    float t2[6] = {2,4,6,8,10,12};
    float t3[6] = {4,8,12,16,20,24};

    float sa = 2;
    saxpy(n, sa, (float *) &t1, 1, (float *) &t2, 1);

    for (int i = 0; i < n; ++i) {
        assert_float(t2[i], ==, t3[i]);
    }

    return MUNIT_OK;
}

MunitResult test_daxpy(const MunitParameter params[], void* user_data_or_fixture) {
    unsigned int n = 6;
    double t1[6] = {1,2,3,4,5,6};
    double t2[6] = {2,4,6,8,10,12};
    double t3[6] = {4,8,12,16,20,24};

    double sa = 2;
    daxpy(n, sa, (double *) &t1, 1, (double *) &t2, 1);

    for (int i = 0; i < n; ++i) {
        assert_double(t2[i], ==, t3[i]);
    }
    return MUNIT_OK;
}

MunitResult test_sswap(const MunitParameter params[], void* user_data_or_fixture) {
    unsigned int n = 6;
    float t1[6] = {1,2,3,4,5,6};
    float t2[6] = {2,4,6,8,10,12};

    float t3[6] = {1,2,3,4,5,6};
    float t4[6] = {2,4,6,8,10,12};

    sswap(n, (float *) &t1, 1, (float *) &t2, 1);

    for (int i = 0; i < n; ++i) {
        assert_float(t2[i], ==, t3[i]);
        assert_float(t1[i], ==, t4[i]);
    }

    return MUNIT_OK;
}

MunitResult test_dswap(const MunitParameter params[], void* user_data_or_fixture) {
    unsigned int n = 6;
    double t1[6] = {1,2,3,4,5,6};
    double t2[6] = {2,4,6,8,10,12};

    double t3[6] = {1,2,3,4,5,6};
    double t4[6] = {2,4,6,8,10,12};

    dswap(n, (double *) &t1, 1, (double *) &t2, 1);

    for (int i = 0; i < n; ++i) {
        assert_double(t2[i], ==, t3[i]);
        assert_double(t1[i], ==, t4[i]);
    }

    return MUNIT_OK;
}

MunitResult test_cswap(const MunitParameter params[], void* user_data_or_fixture) {
    unsigned int n = 6;
    complex float t1[6] = {1 + I,2 + I,3 + I,4 + I,5 + I,6 + I};
    complex float t2[6] = {2 + I,4 + I,6 + I,8 + I,10 + I,12 + I};

    complex float t3[6] = {1 + I,2 + I,3 + I,4 + I,5 + I,6 + I};
    complex float t4[6] = {2 + I,4 + I,6 + I,8 + I,10 + I,12 + I};

    cswap(n, (complex float *) &t1, 1, (complex float *) &t2, 1);

    for (int i = 0; i < n; ++i) {
        assert_float(crealf(t2[i]), ==, crealf(t3[i]));
        assert_float(crealf(t1[i]), ==, crealf(t4[i]));
        assert_float(cimagf(t2[i]), ==, cimagf(t3[i]));
        assert_float(cimagf(t1[i]), ==, cimagf(t4[i]));
    }

    return MUNIT_OK;
}

MunitResult test_zswap(const MunitParameter params[], void* user_data_or_fixture) {
    unsigned int n = 6;
    double t1[6] = {1 + I,2 + I,3 + I,4 + I,5 + I,6 + I};
    double t2[6] = {2 + I,4 + I,6 + I,8 + I,10 + I,12 + I};

    double t3[6] = {1 + I,2 + I,3 + I,4 + I,5 + I,6 + I};
    double t4[6] = {2 + I,4 + I,6 + I,8 + I,10 + I,12 + I};

    zswap(n, (complex double *) &t1, 1, (complex double *) &t2, 1);

    for (int i = 0; i < n; ++i) {
        assert_double_equal(creal(t2[i]), creal(t3[i]), 8);
        assert_double_equal(creal(t1[i]), creal(t4[i]), 8);
        assert_double_equal(cimag(t2[i]), cimag(t3[i]), 8);
        assert_double_equal(cimag(t1[i]), cimag(t4[i]), 8);
    }

    return MUNIT_OK;
}

MunitResult test_isamax(const MunitParameter params[], void* user_data_or_fixture) {
    unsigned int n = 6;
    unsigned int n2 = 3;
    float t1[6] = {1,2,3,4,5,6};
    float t2[3] = {3,2,1};
    int ind  = isamax(n, (float *) &t1, 1);
    int ind2  = isamax(n2, (float *) &t2, 1);
    assert_int(ind, ==, 5);
    assert_int(ind2, ==, 0);

    return MUNIT_OK;
}

MunitResult test_idamax(const MunitParameter params[], void* user_data_or_fixture) {
    unsigned int n = 6;
    unsigned int n2 = 2;
    double t1[6] = {1,2,3,4,5,6};
    double t2[3] = {3,2,1};
    int ind  = idamax(n, (double *) &t1, 1);
    int ind2  = idamax(n2, (double *) &t2, 1);
    assert_int(ind, ==, 5);
    assert_int(ind2, ==, 0);

    return MUNIT_OK;
}

MunitResult test_icamax(const MunitParameter params[], void* user_data_or_fixture) {
    unsigned int n = 6;
    unsigned int n2 = 2;
    complex float t1[6] = {1+I,2+I,3+I,4+I,5+I,6+I};
    complex float t2[3] = {3+I,2+I,1+I};
    int ind  = icamax(n, (complex float *) &t1, 1);
    int ind2  = icamax(n2, (complex float *) &t2, 1);
    assert_int(ind, ==, 5);
    assert_int(ind2, ==, 0);

    return MUNIT_OK;
}
MunitResult test_izamax(const MunitParameter params[], void* user_data_or_fixture) {
    unsigned int n = 6;
    unsigned int n2 = 2;
    complex double t1[6] = {1+I,2+I,3+I,4+I,5+I,6+I};
    complex double t2[3] = {3+I,2+I,1+I};
    int ind  = izamax(n, (complex double *) &t1, 1);
    int ind2  = izamax(n2, (complex double *) &t2, 1);
    assert_int(ind, ==, 5);
    assert_int(ind2, ==, 0);

    return MUNIT_OK;
}

MunitResult test_sasum(const MunitParameter params[], void* user_data_or_fixture) {
    unsigned int n = 6;
    unsigned int n2 = 3;
    float t1[6] = {1,-1,1,-1,-1,1.5};
    float t2[3] = {-3,-2,-1};
    float sum1  = sasum(n, (float *) &t1, 1);
    float sum2  = sasum(n2, (float *) &t2, 1);
    assert_float(sum1, ==, 6.5);
    assert_float(sum2, ==, 6);

    return MUNIT_OK;
}

MunitResult test_dasum(const MunitParameter params[], void* user_data_or_fixture) {
    unsigned int n = 6;
    unsigned int n2 = 3;
    double t1[6] = {1,-1,1,-1,-1,1.5};
    double t2[3] = {-3,-2,-1};
    double sum1  = dasum(n, (double *) &t1, 1);
    double sum2  = dasum(n2, (double *) &t2, 1);
    assert_double_equal(sum1, 6.5, 8);
    assert_double_equal(sum2, 6, 8);

    return MUNIT_OK;
}

MunitResult test_scasum(const MunitParameter params[], void* user_data_or_fixture) {
    unsigned int n = 6;
    unsigned int n2 = 3;
    complex float t1[6] = {1 + 2*I,-1 + 2*I,1 + 2*I,-1 + 2 * I,-1 + 2* I,1.5 + 2*I};
    complex float t2[3] = {-3 + I,-2 + I,-1 + I};
    float sum1  = scasum(n, (complex float *) &t1, 1);
    float sum2  = scasum(n2, (complex float *) &t2, 1);
    assert_double_equal(sum1, 13.680339887498949, 6);
    assert_double_equal(sum2, 6.812559200041265, 6);

    return MUNIT_OK;
}

MunitResult test_dzasum(const MunitParameter params[], void* user_data_or_fixture) {
    unsigned int n = 6;
    unsigned int n2 = 3;
    complex double t1[6] = {1 + 2*I,-1 + 2*I,1 + 2*I,-1 + 2 * I,-1 + 2* I,1.5 + 2*I};
    complex double t2[3] = {-3 + I,-2 + I,-1 + I};
    double sum1  = dzasum(n, (complex double *) &t1, 1);
    double sum2  = dzasum(n2, (complex double *) &t2, 1);
    assert_double_equal(sum1, 13.680339887498949, 16);
    assert_double_equal(sum2, 6.812559200041265, 16);
    return MUNIT_OK;
}


static const MunitTest level1_tests[] = {
        {
                (char*) "/dot/sdot1",
                sdot_over5,
                MUNIT_TEST_OPTION_NONE,
                NULL
        },
        {
                (char*) "/dot/sdot2",
                sdot_under5,
                MUNIT_TEST_OPTION_NONE,
                NULL
        },
        {
                (char*) "/dot/ddot1",
                ddot_over5,
                MUNIT_TEST_OPTION_NONE,
                NULL
        },
        {
                (char*) "/dot/ddot2",
                ddot_under5,
                MUNIT_TEST_OPTION_NONE,
                NULL
        },
        {
                (char*) "/nrm/snrm2",
                test_snrm2,
                MUNIT_TEST_OPTION_NONE,
                NULL
        },
        {
                (char*) "/nrm/dnrm2",
                test_dnrm2,
                MUNIT_TEST_OPTION_NONE,
                NULL
        },
        {
                (char*) "/nrm/scnrm2",
                test_scnrm2,
                MUNIT_TEST_OPTION_NONE,
                NULL
        },
        {
                (char*) "/nrm/dznrm2",
                test_dznrm2,
                MUNIT_TEST_OPTION_NONE,
                NULL
        },
        {
                (char*) "/copy/scopy",
                test_scopy,
                MUNIT_TEST_OPTION_NONE,
                NULL
        },
        {
                (char*) "/copy/dcopy",
                test_dcopy,
                MUNIT_TEST_OPTION_NONE,
                NULL
        },
        {
                (char*) "/copy/ccopy",
                test_ccopy,
                MUNIT_TEST_OPTION_NONE,
                NULL
        },
        {
                (char*) "/copy/zcopy",
                test_zcopy,
                MUNIT_TEST_OPTION_NONE,
                NULL
        },
        {
                (char*) "/scal/sscal",
                test_sscal,
                MUNIT_TEST_OPTION_NONE,
                NULL
        },
        {
                (char*) "/scal/dscal",
                test_dscal,
                MUNIT_TEST_OPTION_NONE,
                NULL
        },
        {
                (char*) "/scal/cscal",
                test_cscal,
                MUNIT_TEST_OPTION_NONE,
                NULL
        },
        {
                (char*) "/scal/zscal",
                test_zscal,
                MUNIT_TEST_OPTION_NONE,
                NULL
        },
        {
                (char*) "/axpy/saxpy",
                test_saxpy,
                MUNIT_TEST_OPTION_NONE,
                NULL
        },
        {
                (char*) "/axpy/daxpy",
                test_daxpy,
                MUNIT_TEST_OPTION_NONE,
                NULL
        },
        {
                (char*) "/swap/sswap",
                test_sswap,
                MUNIT_TEST_OPTION_NONE,
                NULL
        },
        {
                (char*) "/swap/dswap",
                test_dswap,
                MUNIT_TEST_OPTION_NONE,
                NULL
        },
        {
                (char*) "/swap/cswap",
                test_cswap,
                MUNIT_TEST_OPTION_NONE,
                NULL
        },
        {
                (char*) "/swap/zswap",
                test_zswap,
                MUNIT_TEST_OPTION_NONE,
                NULL
        },
        {
                (char*) "/amax/isamax",
                test_isamax,
                MUNIT_TEST_OPTION_NONE,
                NULL
        },
        {
                (char*) "/amax/idamax",
                test_idamax,
                MUNIT_TEST_OPTION_NONE,
                NULL
        },
        {
                (char*) "/amax/icamax",
                test_icamax,
                MUNIT_TEST_OPTION_NONE,
                NULL
        },
        {
                (char*) "/amax/izamax",
                test_izamax,
                MUNIT_TEST_OPTION_NONE,
                NULL
        },
        {
                (char*) "/asum/sasum",
                test_sasum,
                MUNIT_TEST_OPTION_NONE,
                NULL
        },
        {
                (char*) "/asum/dasum",
                test_dasum,
                MUNIT_TEST_OPTION_NONE,
                NULL
        },
        {
                (char*) "/asum/scasum",
                test_scasum,
                MUNIT_TEST_OPTION_NONE,
                NULL
        },
        {
                (char*) "/asum/dzasum",
                test_dzasum,
                MUNIT_TEST_OPTION_NONE,
                NULL
        },
        { NULL, NULL, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL }
};

static const MunitSuite suite = {
        "/level1", /* name */
        level1_tests, /* tests */
        NULL, /* suites */
        1, /* iterations */
        MUNIT_SUITE_OPTION_NONE /* options */
};

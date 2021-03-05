//
// Created by jon on 3/1/21.
//
#include <stdio.h>
#include "level1_tests.h"

#include "../munit/munit.h"

int main(int argc, char* argv[MUNIT_ARRAY_PARAM(argc + 1)]) {
    return munit_suite_main(&suite, NULL, argc, argv);
}
/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file ReduceSerialTests.cpp
 * @brief Unit tests for the global reductions in the serial (no-MPI) build
 *****************************************************************************/

// ********************************** Headers *********************************

// External library headers
#include <catch2/catch_test_macros.hpp>

// Project headers
#include "Reduce.h"
#include "Vector.h"
#include "Integer.h"
#include "Scalar.h"

// ***************************** Serial Reductions ***************************

TEST_CASE("Serial reductions are the identity", "[parallel]")
{
    // Each reduction returns its single-rank input unchanged.
    REQUIRE(globalSum(S(2.5)) == S(2.5));
    REQUIRE(globalSum(Count{42}) == Count{42});

    REQUIRE
    (
        globalSum(Vector(S(1.0), S(2.0), S(3.0)))
     == Vector(S(1.0), S(2.0), S(3.0))
    );

    REQUIRE(globalMax(S(3.0)) == S(3.0));
    REQUIRE(globalMin(S(-1.0)) == S(-1.0));
}

// ****************************** Serial globalOr ****************************

TEST_CASE("Serial globalOr is the identity", "[parallel]")
{
    REQUIRE(globalOr(true));
    REQUIRE_FALSE(globalOr(false));
}
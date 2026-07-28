/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file ReduceParallelTests.cpp
 * @brief Multi-rank tests for the global reductions (run at np = 1/2/4)
 *****************************************************************************/

// ********************************** Headers *********************************

// External library headers
#include <catch2/catch_test_macros.hpp>

// Project headers
#include "Reduce.h"
#include "Comm.h"
#include "Vector.h"
#include "Integer.h"
#include "Scalar.h"

// ****************************** globalSum Count *****************************

TEST_CASE("globalSum of rank-dependent counts", "[mpi][parallel]")
{
    const Count numRanks = Comm::numProcessors();

    // local = rank + 1. The sum over ranks is the triangular number
    const Count local = Comm::myProcessorNum() + 1;

    const Count total = globalSum(local);

    REQUIRE(total == numRanks * (numRanks + 1) / 2);
}

// ************************** globalMax and globalMin *************************

TEST_CASE("globalMax and globalMin across ranks", "[mpi][parallel]")
{
    const Count numRanks = Comm::numProcessors();
    const Scalar local = S(Comm::myProcessorNum());

    const Scalar maximum = globalMax(local);
    const Scalar minimum = globalMin(local);

    REQUIRE(maximum == S(numRanks - 1));
    REQUIRE(minimum == S(0.0));
}

// ****************************** globalSum Vector ****************************

TEST_CASE("globalSum of a rank-dependent vector", "[mpi][parallel]")
{
    const Count numRanks = Comm::numProcessors();
    const Scalar rank = S(Comm::myProcessorNum());

    const Vector local(rank, S(2.0) * rank, S(3.0) * rank);

    // sum of the ranks 0 .. numRanks-1 is another triangular number
    const Scalar sumOfRanks = S(numRanks * (numRanks - 1) / 2);

    const Vector total = globalSum(local);

    REQUIRE
    (
        total
     == Vector(sumOfRanks, S(2.0) * sumOfRanks, S(3.0) * sumOfRanks)
    );
}

// ********************************* globalOr *********************************

TEST_CASE("globalOr across ranks", "[mpi][parallel]")
{
    // Only the last rank contributes true, every rank must see true
    const bool trueOnLastRank =
        Comm::myProcessorNum() == Comm::numProcessors() - 1;

    const bool anyLastRank = globalOr(trueOnLastRank);

    // All-false stays false everywhere
    const bool anyFalse = globalOr(false);

    REQUIRE(anyLastRank);
    REQUIRE_FALSE(anyFalse);
}
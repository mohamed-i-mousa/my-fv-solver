/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file ReduceConsistencyTests.cpp
 * @brief Reductions of a fixed global set are independent of the rank count
 *****************************************************************************/

// ********************************** Headers *********************************

// Standard library headers
#include <limits>

// External library headers
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

// Project headers
#include "Reduce.h"
#include "Comm.h"
#include "Vector.h"
#include "Scalar.h"
#include "Integer.h"
#include "TestTolerances.h"

using Catch::Matchers::WithinRel;

// ***************************** Internal Helpers *****************************

namespace
{

/// Length of the fixed global sequence (divisible by 1, 2, and 4)
constexpr Count sequenceLength = 24;

} // namespace

// *********************** Rank-Count-Invariant Sums ************************

TEST_CASE
(
    "Scalar reductions match their closed forms at any np",
    "[mpi][parallel]"
)
{
    const Index rank = Comm::myProcessorNum();
    const Count size = Comm::numProcessors();

    Scalar localSum = S(0.0);
    Scalar localMax = std::numeric_limits<Scalar>::lowest();
    Scalar localMin = std::numeric_limits<Scalar>::max();
    Count localCount = 0;

    for (Count g = rank; g < sequenceLength; g += size)
    {
        const Scalar value = S(g + 1);
        localSum += value;
        localMax = std::max(localMax, value);
        localMin = std::min(localMin, value);
        ++localCount;
    }

    // Every collective completes before the first assertion: a failing
    // REQUIRE unwinds one rank only, and the rest would hang in the next one
    const Scalar totalSum = globalSum(localSum);
    const Scalar totalMax = globalMax(localMax);
    const Scalar totalMin = globalMin(localMin);
    const Count totalCount = globalSum(localCount);
    const Count contributingRanks = globalSum(Count{1});

    // sum 1..N = N(N+1)/2, max = N, min = 1, regardless of the partition
    REQUIRE_THAT
    (
        totalSum,
        WithinRel
        (
            S(sequenceLength * (sequenceLength + 1) / 2),
            TestTolerances::relTight
        )
    );
    REQUIRE_THAT
    (
        totalMax,
        WithinRel(S(sequenceLength),
        TestTolerances::relTight)
    );
    REQUIRE_THAT
    (
        totalMin,
        WithinRel(S(1.0),
        TestTolerances::relTight)
    );

    // Every value is accounted for exactly once across the ranks
    REQUIRE(totalCount == sequenceLength);
    REQUIRE(contributingRanks == size);
}

// ********************** Rank-Count-Invariant Vector **********************

TEST_CASE
(
    "Vector reduction matches its closed form at any np",
    "[mpi][parallel]"
)
{
    const Index rank = Comm::myProcessorNum();
    const Count size = Comm::numProcessors();

    Vector localSum(S(0.0), S(0.0), S(0.0));
    for (Count g = rank; g < sequenceLength; g += size)
    {
        localSum += Vector(S(g + 1), S(0.0), S(0.0));
    }

    const Vector total = globalSum(localSum);
    REQUIRE_THAT
    (
        total.x(),
        WithinRel
        (
            S(sequenceLength * (sequenceLength + 1) / 2),
            TestTolerances::relTight
        )
    );
    REQUIRE(total.y() == S(0.0));
    REQUIRE(total.z() == S(0.0));
}

// ***************************** Logical Or ********************************

TEST_CASE
(
    "globalOr is true if any rank contributes true",
    "[mpi][parallel]"
)
{
    const Index rank = Comm::myProcessorNum();
    const Count size = Comm::numProcessors();

    // Exactly one rank (the last) contributes true, then no rank does
    const bool anyLastRank = globalOr(rank + 1 == size);
    const bool anyFalse = globalOr(false);

    REQUIRE(anyLastRank);
    REQUIRE_FALSE(anyFalse);
}
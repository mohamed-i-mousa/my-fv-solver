/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file GlobalIndexParallelTests.cpp
 * @brief Multi-rank tests for the rank-major GlobalIndex (run at np = 1/2/4)
 *****************************************************************************/

// ********************************** Headers *********************************

// External library headers
#include <catch2/catch_test_macros.hpp>

// Project headers
#include "GlobalIndex.h"
#include "Comm.h"
#include "Reduce.h"
#include "Integer.h"

// **************************** Contiguous Offsets ****************************

TEST_CASE("Contiguous rank offsets", "[mpi][parallel]")
{
    // Each rank owns a rank-dependent count. The GlobalIndex constructor is
    // collective at np > 1, so it is built unconditionally and every rank
    // reaches it in lockstep before any assertion.
    const Count local = Count{3} + Comm::myProcessorNum();
    const GlobalIndex g(local);

    // totalCount is identical on every rank: the global sum of local counts
    const Count total = globalSum(local);

    REQUIRE(g.totalCount() == total);

    // offset is this rank's exclusive prefix. Replicate it with pure integer
    // arithmetic that needs no further collective, keeping lockstep intact.
    const Index myRank = Comm::myProcessorNum();
    Index expectedOffset = 0;
    for (Index r = 0; r < myRank; ++r)
    {
        expectedOffset += 3 + r;
    }

    REQUIRE(g.offset() == expectedOffset);

    // Local item 0 maps to the offset itself
    REQUIRE(g.toGlobal(0) == g.offset());
}
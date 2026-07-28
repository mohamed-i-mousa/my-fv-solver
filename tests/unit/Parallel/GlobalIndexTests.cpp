/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file GlobalIndexTests.cpp
 * @brief Unit tests for the serial GlobalIndex rank-major numbering
 *****************************************************************************/

// ********************************** Headers *********************************

// External library headers
#include <catch2/catch_test_macros.hpp>

// Project headers
#include "GlobalIndex.h"
#include "Integer.h"

// ************************** Serial Identity Mapping *************************

TEST_CASE("Serial global numbering is the identity", "[parallel]")
{
    // In the serial binary the constructor takes the no-MPI branch: this
    // rank owns every item, so its slice starts at global 0 and toGlobal()
    // is the identity on owned local indices.
    const GlobalIndex g(Count{7});

    REQUIRE(g.offset() == 0);
    REQUIRE(g.totalCount() == 7);

    for (Index i = 0; i <= 6; ++i)
    {
        REQUIRE(g.toGlobal(i) == i);
    }
}
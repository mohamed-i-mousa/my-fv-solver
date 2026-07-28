/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file CommTests.cpp
 * @brief Unit tests for the Comm rank-identity queries
 *****************************************************************************/

// ********************************** Headers *********************************

// External library headers
#include <catch2/catch_test_macros.hpp>

// Project headers
#include "Comm.h"

// *************************** Uninitialized Identity **************************

TEST_CASE("Uninitialized MPI reports a single master rank", "[parallel]")
{
    // The serial test binary never calls MPI_Init, so Comm answers from its
    // uninitialized single-rank identity: not parallel, one rank, rank zero,
    // and that rank is the master.
    REQUIRE(Comm::parallelRun() == false);
    REQUIRE(Comm::numProcessors() == 1);
    REQUIRE(Comm::myProcessorNum() == 0);
    REQUIRE(Comm::master() == true);
}
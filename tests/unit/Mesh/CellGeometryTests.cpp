/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file CellGeometryTests.cpp
 * @brief Unit tests for Cell volume and centroid on the hex-box fixture
 *****************************************************************************/

// ********************************** Headers *********************************

// External library headers
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

// Project headers
#include "MeshFixtures.h"
#include "TestTolerances.h"
#include "Mesh.h"
#include "Cell.h"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

// ***************************** Unit Cube Volume *****************************

TEST_CASE("Single unit cube volume and centroid", "[mesh]")
{
    const TestMesh box(1, 1, 1);
    const Cell& cell = box.mesh().cells()[0];

    REQUIRE_THAT
    (
        cell.volume(),
        WithinRel(S(1.0), TestTolerances::relTight)
    );

    REQUIRE_THAT(cell.centroid().x(), WithinAbs(S(0.5), TestTolerances::absTight));
    REQUIRE_THAT(cell.centroid().y(), WithinAbs(S(0.5), TestTolerances::absTight));
    REQUIRE_THAT(cell.centroid().z(), WithinAbs(S(0.5), TestTolerances::absTight));
}

// **************************** Scaled Cube Volume ****************************

TEST_CASE("Scaled cube volume", "[mesh]")
{
    // A single cell of edge 0.25 has volume 0.25^3 = 0.015625
    const TestMesh box(1, 1, 1, S(0.25));
    const Cell& cell = box.mesh().cells()[0];

    REQUIRE_THAT
    (
        cell.volume(),
        WithinRel(S(0.015625), TestTolerances::relTight)
    );

    REQUIRE_THAT
    (
        cell.centroid().x(),
        WithinAbs(S(0.125), TestTolerances::absTight)
    );
}

// **************************** Interior Cell Count ***************************

TEST_CASE("Hex box reports the expected cell and face counts", "[mesh]")
{
    const TestMesh box(2, 2, 2);

    // 2x2x2 = 8 cells; 12 internal + 24 boundary faces = 36
    REQUIRE(box.mesh().numCells() == 8);
    REQUIRE(box.mesh().numFaces() == 36);
    REQUIRE(box.mesh().numOwnedCells() == box.mesh().numCells());
    REQUIRE(box.mesh().numGhostCells() == 0);

    // Every interior cell of the box carries a positive volume
    for (Index cellIdx = 0; cellIdx < box.mesh().numCells(); ++cellIdx)
    {
        REQUIRE(box.mesh().cells()[cellIdx].volume() > S(0.0));
    }
}

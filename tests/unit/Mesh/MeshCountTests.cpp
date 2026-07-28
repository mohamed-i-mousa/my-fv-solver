/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file MeshCountTests.cpp
 * @brief Unit tests for the static cell/face counts that size field containers
 *****************************************************************************/

// ********************************** Headers *********************************

// External library headers
#include <catch2/catch_test_macros.hpp>

// Project headers
#include "MeshFixtures.h"
#include "Mesh.h"
#include "CellData.h"
#include "Vector.h"

// ****************************** Static Counts *******************************

TEST_CASE("Static counts follow the constructed mesh", "[mesh]")
{
    const TestMesh box(2, 2, 2);

    // 2x2x2 = 8 cells; 12 internal + 24 boundary faces = 36
    REQUIRE(Mesh::cellCount() == 8);
    REQUIRE(Mesh::faceCount() == 36);
    REQUIRE(box.mesh().numOwnedCells() == box.mesh().numCells());
    REQUIRE(box.mesh().numGhostCells() == 0);
}

// ***************************** Count Slot Reuse *****************************

TEST_CASE("resetCounts releases the slot for the next mesh", "[mesh]")
{
    // The first mesh claims the single global count slot; its destructor
    // calls Mesh::resetCounts() on scope exit. Only one mesh may be alive at
    // a time.
    {
        const TestMesh box(2, 2, 2);
        REQUIRE(Mesh::cellCount() == 8);
    }

    // The slot is free again: a differently sized mesh reclaims it
    const TestMesh box2(3, 1, 1);
    REQUIRE(Mesh::cellCount() == 3);
}

// ******************************* Field Sizing *******************************

TEST_CASE("Field containers size from the counts", "[mesh]")
{
    const TestMesh box(2, 1, 1);

    // Fields size themselves from Mesh::cellCount(), fixed at construction
    // of the mesh above, so they must be built after the TestMesh exists.
    const ScalarField s;
    REQUIRE(s.size() == box.mesh().numCells());

    const VectorField v(Vector(S(1.0), S(0.0), S(0.0)));
    REQUIRE(v.size() == box.mesh().numCells());
    REQUIRE(v[0] == Vector(S(1.0), S(0.0), S(0.0)));
}
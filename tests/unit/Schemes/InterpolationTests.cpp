/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file InterpolationTests.cpp
 * @brief Unit tests for linear cell-to-face interpolation
 *****************************************************************************/

// ********************************** Headers *********************************

// External library headers
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

// Project headers
#include "MeshFixtures.h"
#include "TestTolerances.h"
#include "LinearInterpolation.h"
#include "CellData.h"
#include "Mesh.h"
#include "Cell.h"
#include "Face.h"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

// ***************************** Internal Helpers ****************************

namespace
{

// Locate the single internal (non-boundary) face of the box(2,1,1) fixture
[[nodiscard]] const Face& internalFace(const Mesh& mesh)
{
    const FaceList& faces = mesh.faces();

    for (Index faceIdx = 0; faceIdx < faces.size(); ++faceIdx)
    {
        if (!faces[faceIdx].isBoundary())
        {
            return faces[faceIdx];
        }
    }

    FatalError("box(2,1,1) fixture must contain exactly one internal face");
}

} // namespace

// ****************************** Constant Field *****************************

TEST_CASE("Constant field interpolates to itself", "[schemes]")
{
    const TestMesh box(2, 1, 1);

    // ScalarField sizes itself from Mesh::cellCount(); build it after the mesh
    ScalarField phi;
    phi[0] = S(7.5);
    phi[1] = S(7.5);

    const Face& face = internalFace(box.mesh());

    REQUIRE_THAT
    (
        interpolateToFace(face, phi),
        WithinAbs(S(7.5), TestTolerances::absTight)
    );
}

// **************************** Midpoint Weighting **************************

TEST_CASE("Equal distances give the midpoint", "[schemes]")
{
    const TestMesh box(2, 1, 1);

    // Uniform spacing makes dPf == dNf, so the face value is the midpoint
    ScalarField phi;
    phi[0] = S(1.0);
    phi[1] = S(3.0);

    const Face& face = internalFace(box.mesh());

    REQUIRE_THAT
    (
        interpolateToFace(face, phi),
        WithinRel(S(2.0), TestTolerances::relTight)
    );
}

// ************************** Linear Field Exactness ************************

TEST_CASE("Linear field is exact at the face", "[schemes]")
{
    const TestMesh box(2, 1, 1);

    // phi(x) = x sampled at the cell centroids (0.5 and 1.5); the internal
    // face sits at x = 1, so the interpolated value must be exactly 1.0
    ScalarField phi;

    const CellList& cells = box.mesh().cells();
    for (Index cellIdx = 0; cellIdx < cells.size(); ++cellIdx)
    {
        phi[cellIdx] = cells[cellIdx].centroid().x();
    }

    const Face& face = internalFace(box.mesh());

    REQUIRE_THAT
    (
        interpolateToFace(face, phi),
        WithinRel(S(1.0), TestTolerances::relTight)
    );
}
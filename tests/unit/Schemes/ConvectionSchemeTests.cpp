/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file ConvectionSchemeTests.cpp
 * @brief Unit tests for the convection-scheme deferred-correction terms
 *****************************************************************************/

// ********************************** Headers *********************************

// External library headers
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

// Standard library headers
#include <algorithm>

// Project headers
#include "MeshFixtures.h"
#include "TestTolerances.h"
#include "ConvectionScheme.h"
#include "Upwind.h"
#include "SecondOrderUpwind.h"
#include "CentralDifference.h"
#include "Mesh.h"
#include "Face.h"
#include "Vector.h"
#include "CellData.h"
#include "StringTypes.h"
#include "Integer.h"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

// ***************************** Internal Helpers *****************************

namespace
{

/// The single internal face of a 2x1x1 hex box (emitted before boundaries)
[[nodiscard]] const Face& internalFace(const Mesh& mesh)
{
    for (Index faceIdx = 0; faceIdx < mesh.numFaces(); ++faceIdx)
    {
        if (!mesh.faces()[faceIdx].isBoundary())
        {
            return mesh.faces()[faceIdx];
        }
    }

    // A 2x1x1 box always has exactly one internal face, so this is unreached
    return mesh.faces()[0];
}

/// Whether a scheme-name list carries the given entry
[[nodiscard]] bool contains(const NameList& names, const Name& target)
{
    return std::find(names.begin(), names.end(), target) != names.end();
}

} // namespace

// ***************************** Upwind Correction ***************************

TEST_CASE("Upwind correction is always zero", "[schemes]")
{
    const TestMesh box(2, 1, 1);
    const Face& face = internalFace(box.mesh());
    const ScalarField phi;

    const Upwind up;

    const Vector gradA(S(2.0), S(-1.0), S(3.0));
    const Vector gradB(S(-4.0), S(5.0), S(0.5));
    const Vector zero(S(0.0), S(0.0), S(0.0));

    // First-order upwind never adds a deferred correction
    REQUIRE_THAT
    (
        up.correction(face, phi, gradA, gradB, S(3.0)),
        WithinAbs(S(0.0), TestTolerances::absTight)
    );

    REQUIRE_THAT
    (
        up.correction(face, phi, gradA, gradB, S(-3.0)),
        WithinAbs(S(0.0), TestTolerances::absTight)
    );

    REQUIRE_THAT
    (
        up.correction(face, phi, gradB, gradA, S(0.0)),
        WithinAbs(S(0.0), TestTolerances::absTight)
    );

    REQUIRE_THAT
    (
        up.correction(face, phi, zero, zero, S(100.0)),
        WithinAbs(S(0.0), TestTolerances::absTight)
    );
}

// *************************** Second Order Upwind ***************************

TEST_CASE("SecondOrderUpwind uses the upwind gradient", "[schemes]")
{
    const TestMesh box(2, 1, 1);
    const Face& face = internalFace(box.mesh());
    const ScalarField phi;

    const SecondOrderUpwind so;

    const Vector gradP(S(2.0), S(0.0), S(0.0));
    const Vector gradN(S(2.0), S(0.0), S(0.0));

    // Positive flux upwinds to the owner: the correction projects the owner
    // gradient onto the owner-to-face vector
    const Scalar flowPos = S(3.0);
    const Scalar expectedPos = flowPos * dot(gradP, face.dPf());

    REQUIRE_THAT
    (
        so.correction(face, phi, gradP, gradN, flowPos),
        WithinRel(expectedPos, TestTolerances::relTight)
    );

    // Negative flux upwinds to the neighbor: the correction projects the
    // neighbor gradient onto the neighbor-to-face vector
    const Scalar flowNeg = S(-3.0);
    const Scalar expectedNeg = flowNeg * dot(gradN, face.dNf().value());

    REQUIRE_THAT
    (
        so.correction(face, phi, gradP, gradN, flowNeg),
        WithinRel(expectedNeg, TestTolerances::relTight)
    );

    // Zero flux gives zero correction
    REQUIRE_THAT
    (
        so.correction(face, phi, gradP, gradN, S(0.0)),
        WithinAbs(S(0.0), TestTolerances::absTight)
    );
}

// **************************** Central Difference ***************************

TEST_CASE("CentralDifference vanishes for constant phi", "[schemes]")
{
    const TestMesh box(2, 1, 1);
    const Face& face = internalFace(box.mesh());

    // A uniform field: the linearly interpolated face value equals the upwind
    // cell value, so the deferred correction cancels for any non-zero flux
    ScalarField phi;
    phi[0] = S(4.0);
    phi[1] = S(4.0);

    const CentralDifference cd;

    const Vector zero(S(0.0), S(0.0), S(0.0));

    REQUIRE_THAT
    (
        cd.correction(face, phi, zero, zero, S(3.0)),
        WithinAbs(S(0.0), TestTolerances::absTight)
    );
}

// ***************************** Runtime Selection ***************************

TEST_CASE("Convection scheme factory", "[schemes]")
{
    REQUIRE(ConvectionScheme::create("Upwind") != nullptr);
    REQUIRE(ConvectionScheme::create("CentralDifference") != nullptr);
    REQUIRE(ConvectionScheme::create("SecondOrderUpwind") != nullptr);
    REQUIRE(ConvectionScheme::create("LUST") != nullptr);

    const NameList names = ConvectionScheme::availableSchemes();

    REQUIRE(contains(names, "Upwind"));
    REQUIRE(contains(names, "CentralDifference"));
    REQUIRE(contains(names, "SecondOrderUpwind"));
}
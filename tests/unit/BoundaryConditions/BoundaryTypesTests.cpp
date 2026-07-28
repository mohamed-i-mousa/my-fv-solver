/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file BoundaryTypesTests.cpp
 * @brief Unit tests for the BoundaryType hierarchy linearization formulas
 *
 * @details Every BoundaryType method is a pure function of hand-built scalar
 * and vector arguments, so the whole hierarchy is testable with no mesh and no
 * manager. Each type's diagonal, source, and face-value contributions are
 * pinned against the closed forms documented in its header.
 *****************************************************************************/

// ********************************** Headers *********************************

// Standard library headers
#include <algorithm>

// External library headers
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

// Project headers
#include "BoundaryType.h"
#include "FixedValue.h"
#include "ZeroGradient.h"
#include "FixedGradient.h"
#include "NoSlip.h"
#include "Field.h"
#include "Vector.h"
#include "TestTolerances.h"

using Catch::Matchers::WithinRel;

// ***************************** Internal Helpers *****************************

namespace
{

/// Representative per-face data handed to a boundary type by the manager
const Scalar testFlux = S(2.0);
const Scalar testGammaSf = S(3.0);
const Scalar testDiffMetric = S(4.0);
const Scalar testNormalDistance = S(0.5);
const Scalar testOwnerValue = S(1.5);
const Vector testNormal(S(1.0), S(0.0), S(0.0));
const Vector testOwnerVelocity(S(0.0), S(0.0), S(0.0));

/// Whether a name list contains a given token
[[nodiscard]] bool contains(const NameList& names, const Name& target)
{
    return std::find(names.begin(), names.end(), target) != names.end();
}

} // namespace

// ******************************** FixedValue ********************************

TEST_CASE("FixedValue linearization", "[bc]")
{
    const FixedValue bc(Field::p, S(2.0));

    REQUIRE(bc.typeName() == "fixedValue");
    REQUIRE(bc.fixesValue());

    // diag = GammaSf * diffMetric = 3 * 4 = 12
    REQUIRE_THAT
    (
        bc.addToDiagonal(testFlux, testGammaSf, testDiffMetric, testNormal),
        WithinRel(S(12.0), TestTolerances::relTight)
    );

    // source = value * (GammaSf * diffMetric - flux) = 2 * (12 - 2) = 20
    REQUIRE_THAT
    (
        bc.addToSource
        (
            testFlux,
            testGammaSf,
            testDiffMetric,
            testNormalDistance,
            testNormal,
            testOwnerVelocity
        ),
        WithinRel(S(20.0), TestTolerances::relTight)
    );

    // The prescribed value is the face value regardless of the owner value
    REQUIRE_THAT
    (
        bc.faceValue
        (
            testOwnerValue,
            testNormalDistance,
            testNormal,
            testOwnerVelocity
        ),
        WithinRel(S(2.0), TestTolerances::relTight)
    );
}

// ******************************* ZeroGradient *******************************

TEST_CASE("ZeroGradient linearization", "[bc]")
{
    const ZeroGradient bc(Field::p);

    REQUIRE(bc.typeName() == "zeroGradient");

    // diag = flux
    REQUIRE_THAT
    (
        bc.addToDiagonal(testFlux, testGammaSf, testDiffMetric, testNormal),
        WithinRel(testFlux, TestTolerances::relTight)
    );

    // source = 0
    REQUIRE_THAT
    (
        bc.addToSource
        (
            testFlux,
            testGammaSf,
            testDiffMetric,
            testNormalDistance,
            testNormal,
            testOwnerVelocity
        ),
        Catch::Matchers::WithinAbs(S(0.0), TestTolerances::absTight)
    );

    // face value is the owner value
    REQUIRE_THAT
    (
        bc.faceValue
        (
            testOwnerValue,
            testNormalDistance,
            testNormal,
            testOwnerVelocity
        ),
        WithinRel(testOwnerValue, TestTolerances::relTight)
    );
}

// ******************************* FixedGradient ******************************

TEST_CASE("FixedGradient linearization", "[bc]")
{
    const FixedGradient bc(Field::p, S(0.5));

    REQUIRE(bc.typeName() == "fixedGradient");
    REQUIRE(bc.correctsBoundaryFlux());

    // diag = flux
    REQUIRE_THAT
    (
        bc.addToDiagonal(testFlux, testGammaSf, testDiffMetric, testNormal),
        WithinRel(testFlux, TestTolerances::relTight)
    );

    // source = GammaSf * g - flux * g * dn = 3*0.5 - 2*0.5*0.5 = 1.0
    REQUIRE_THAT
    (
        bc.addToSource
        (
            testFlux, testGammaSf, testDiffMetric,
            testNormalDistance, testNormal, testOwnerVelocity
        ),
        WithinRel(S(1.0), TestTolerances::relTight)
    );

    // phi_f = phi_P + g * dn = 1.5 + 0.5 * 0.5 = 1.75
    REQUIRE_THAT
    (
        bc.faceValue
        (
            testOwnerValue,
            testNormalDistance,
            testNormal,
            testOwnerVelocity
        ),
        WithinRel(S(1.75), TestTolerances::relTight)
    );
}

// ********************************** NoSlip **********************************

TEST_CASE("NoSlip is a zero-valued Dirichlet", "[bc]")
{
    const NoSlip bc(Field::Ux);

    REQUIRE(bc.typeName() == "noSlip");
    REQUIRE(bc.fixesValue());

    // A zero fixed value: the face value is zero for any owner value
    REQUIRE_THAT
    (
        bc.faceValue
        (
            testOwnerValue,
            testNormalDistance,
            testNormal,
            testOwnerVelocity
        ),
        Catch::Matchers::WithinAbs(S(0.0), TestTolerances::absTight)
    );
}

// *************************** Selectable Type Names **************************

TEST_CASE("availableTypes lists the case-file-selectable BCs", "[bc]")
{
    const NameList velocityTypes = BoundaryType::availableTypes(Field::Ux);

    REQUIRE(contains(velocityTypes, "fixedValue"));
    REQUIRE(contains(velocityTypes, "fixedGradient"));
    REQUIRE(contains(velocityTypes, "noSlip"));
    REQUIRE(contains(velocityTypes, "zeroGradient"));

    // The pressure correction field is not user-selectable
    REQUIRE(BoundaryType::availableTypes(Field::pCorr).empty());
}
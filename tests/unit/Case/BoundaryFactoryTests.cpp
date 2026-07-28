/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file BoundaryFactoryTests.cpp
 * @brief BoundaryType::create from case sections and write() round-trips
 *****************************************************************************/

// ********************************** Headers *********************************

// Standard library headers
#include <memory>
#include <sstream>
#include <string>

// External library headers
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

// Project headers
#include "BoundaryType.h"
#include "Symmetry.h"
#include "CaseReader.h"
#include "Field.h"
#include "Vector.h"
#include "StringTypes.h"
#include "TestTolerances.h"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

// ***************************** Internal Helpers *****************************

namespace
{

/// Per-face geometry a boundary type ignores when it returns a fixed value
const Scalar anyNormalDistance = S(0.5);
const Vector anyNormal(S(1.0), S(0.0), S(0.0));
const Vector anyVelocity(S(0.0), S(0.0), S(0.0));

/// The boundaryTypes section of the committed parser fixture
[[nodiscard]] CaseReader boundaryTypesFixture()
{
    const CaseReader reader
    (
        FilePath(TURBLYZE_TEST_FIXTURE_DIR) + "/cases/parserCase"
    );
    return reader.section("boundaryTypes");
}

/// Whether write() output begins with the type name
[[nodiscard]] bool writeStartsWithTypeName(const BoundaryType& bc)
{
    std::ostringstream oss;
    bc.write(oss);
    return oss.str().rfind(bc.typeName(), 0) == 0;
}

} // namespace

// ************************* Fixed-Value Velocity ****************************

TEST_CASE("Factory builds a fixedValue velocity component", "[bc][selection]")
{
    const CaseReader bcSection = boundaryTypesFixture();
    const CaseReader& patch = bcSection.section("fixedVelocity");

    // Each velocity component draws its own scalar from the (1 2 3) vector
    const auto bcX = BoundaryType::create("fixedValue", Field::Ux, patch);
    const auto bcZ = BoundaryType::create("fixedValue", Field::Uz, patch);

    REQUIRE(bcX->typeName() == "fixedValue");
    REQUIRE(bcX->field() == Field::Ux);
    REQUIRE(bcX->fixesValue());

    REQUIRE_THAT
    (
        bcX->faceValue(S(0.0), anyNormalDistance, anyNormal, anyVelocity),
        WithinRel(S(1.0), TestTolerances::relTight)
    );
    REQUIRE_THAT
    (
        bcZ->faceValue(S(0.0), anyNormalDistance, anyNormal, anyVelocity),
        WithinRel(S(3.0), TestTolerances::relTight)
    );

    REQUIRE(writeStartsWithTypeName(*bcX));
}

// ************************** Fixed-Value Scalar *****************************

TEST_CASE("Factory builds a fixedValue scalar", "[bc][selection]")
{
    const CaseReader bcSection = boundaryTypesFixture();
    const auto bc = BoundaryType::create
    (
        "fixedValue", Field::p, bcSection.section("fixedScalar")
    );

    REQUIRE(bc->typeName() == "fixedValue");
    REQUIRE_THAT
    (
        bc->faceValue(S(1.0), anyNormalDistance, anyNormal, anyVelocity),
        WithinRel(S(2.5), TestTolerances::relTight)
    );
    REQUIRE(writeStartsWithTypeName(*bc));
}

// ************************** Fixed-Gradient Scalar **************************

TEST_CASE("Factory builds a fixedGradient scalar", "[bc][selection]")
{
    const CaseReader bcSection = boundaryTypesFixture();
    const auto bc = BoundaryType::create
    (
        "fixedGradient", Field::p, bcSection.section("gradientScalar")
    );

    REQUIRE(bc->typeName() == "fixedGradient");
    REQUIRE(bc->correctsBoundaryFlux());

    // phi_f = phi_P + gradient * normalDistance = 1.0 + 0.5 * 0.5 = 1.25
    REQUIRE_THAT
    (
        bc->faceValue(S(1.0), anyNormalDistance, anyNormal, anyVelocity),
        WithinRel(S(1.25), TestTolerances::relTight)
    );
    REQUIRE(writeStartsWithTypeName(*bc));
}

// ****************************** Trait Types ********************************

TEST_CASE("Factory builds the zero-parameter trait types", "[bc][selection]")
{
    const CaseReader bcSection = boundaryTypesFixture();

    const auto slip = BoundaryType::create
    (
        "zeroGradient", Field::p, bcSection.section("slipWall")
    );
    REQUIRE(slip->typeName() == "zeroGradient");
    REQUIRE_THAT
    (
        slip->faceValue(S(4.0), anyNormalDistance, anyNormal, anyVelocity),
        WithinRel(S(4.0), TestTolerances::relTight)
    );

    const auto wall = BoundaryType::create
    (
        "noSlip", Field::Ux, bcSection.section("solidWall")
    );
    REQUIRE(wall->typeName() == "noSlip");
    REQUIRE(wall->fixesValue());
    REQUIRE_THAT
    (
        wall->faceValue(S(4.0), anyNormalDistance, anyNormal, anyVelocity),
        WithinAbs(S(0.0), TestTolerances::absTight)
    );

    REQUIRE(writeStartsWithTypeName(*slip));
    REQUIRE(writeStartsWithTypeName(*wall));
}

// ******************************* Symmetry **********************************

TEST_CASE("Symmetry mirrors velocity and passes scalars through", "[bc]")
{
    // Mesh-derived, never case-file selectable, so it is constructed directly
    const Symmetry symmetryScalar(Field::p);
    const Symmetry symmetryVelocity(Field::Ux);

    REQUIRE(symmetryScalar.typeName() == "symmetry");
    REQUIRE(symmetryVelocity.constrainsZeroFlux());
    REQUIRE_FALSE(symmetryVelocity.contributesToLimiterHull());

    // A scalar sees the plane as zero-gradient: the face value is the owner
    REQUIRE_THAT
    (
        symmetryScalar.faceValue(S(3.3), anyNormalDistance, anyNormal, anyVelocity),
        WithinRel(S(3.3), TestTolerances::relTight)
    );

    // For an x-normal plane the mirrored Ux is (1 - 1) * U_P - 1 * 0 = 0
    const Vector xNormal(S(1.0), S(0.0), S(0.0));
    const Vector ownerVel(S(5.0), S(2.0), S(3.0));
    REQUIRE_THAT
    (
        symmetryVelocity.faceValue(S(5.0), anyNormalDistance, xNormal, ownerVel),
        WithinAbs(S(0.0), TestTolerances::absTight)
    );
}
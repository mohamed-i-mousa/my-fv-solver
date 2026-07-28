/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file GradientSchemeTests.cpp
 * @brief Least-squares gradient exactness on a linear field
 *****************************************************************************/

// ********************************** Headers *********************************

// External library headers
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

// Project headers
#include "MeshFixtures.h"
#include "TestTolerances.h"
#include "LeastSquares.h"
#include "GradientScheme.h"
#include "BoundaryConditions.h"
#include "CellData.h"
#include "Vector.h"
#include "Field.h"

using Catch::Matchers::WithinAbs;

// ***************************** Internal Helpers *****************************

namespace
{

/// The linear field phi = 1.5 x + 2.5 y - 3.5 z + 4.2
/// its exact gradient is (1.5, 2.5, -3.5)
[[nodiscard]] Scalar linearField(const Vector& x) noexcept
{
    return S(1.5) * x.x() + S(2.5) * x.y() - S(3.5) * x.z() + S(4.2);
}

/// Linear index of the interior centre cell of a 3x3x3 box
constexpr Index centreCell = 13;

} // namespace

// ************************** Linear-Field Exactness **************************

TEST_CASE("Least-squares gradient of a linear field is exact", "[schemes]")
{
    const TestMesh box(3, 3, 3);
    const BoundaryConditions bc;
    const LeastSquares leastSquares(box.mesh(), bc);

    ScalarField phi;

    for (Index cellIdx = 0; cellIdx < box.mesh().numCells(); ++cellIdx)
    {
        phi[cellIdx] = linearField(box.mesh().cells()[cellIdx].centroid());
    }

    const Vector gradient =
        leastSquares.cellGradient(Field::p, phi, centreCell);

    REQUIRE_THAT
    (
        gradient.x(), WithinAbs(S(1.5),
        TestTolerances::absOperator)
    );
    REQUIRE_THAT
    (
        gradient.y(), WithinAbs(S(2.5),
        TestTolerances::absOperator)
    );
    REQUIRE_THAT
    (
        gradient.z(), WithinAbs(S(-3.5),
        TestTolerances::absOperator)
    );
}

// **************************** Constant-Field Zero ***************************

TEST_CASE("Least-squares gradient of a constant field is zero", "[schemes]")
{
    const TestMesh box(3, 3, 3);
    const BoundaryConditions bc;
    const LeastSquares leastSquares(box.mesh(), bc);

    ScalarField phi(S(4.2));

    const Vector gradient =
        leastSquares.cellGradient(Field::p, phi, centreCell);

    REQUIRE_THAT(gradient.x(), WithinAbs(S(0.0), TestTolerances::absTight));
    REQUIRE_THAT(gradient.y(), WithinAbs(S(0.0), TestTolerances::absTight));
    REQUIRE_THAT(gradient.z(), WithinAbs(S(0.0), TestTolerances::absTight));
}

// ******************************** Factory ***********************************

TEST_CASE("Gradient-scheme factory builds leastSquares", "[schemes]")
{
    const TestMesh box(2, 1, 1);
    const BoundaryConditions bc;

    const auto scheme =
        GradientScheme::create("leastSquares", box.mesh(), bc);

    REQUIRE(scheme != nullptr);

    const NameList available = GradientScheme::availableSchemes();

    REQUIRE_FALSE(available.empty());
    REQUIRE(available.front() == "leastSquares");
}
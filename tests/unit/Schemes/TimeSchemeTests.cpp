/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file TimeSchemeTests.cpp
 * @brief Unit tests for the time-derivative discretization schemes
 *****************************************************************************/

// ********************************** Headers *********************************

// External library headers
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

// Standard library headers
#include <algorithm>
#include <memory>

// Project headers
#include "TimeScheme.h"
#include "SteadyState.h"
#include "ImplicitEuler.h"
#include "CrankNicolson.h"
#include "StringTypes.h"
#include "TestTolerances.h"

using Catch::Matchers::WithinRel;

// ***************************** Internal Helpers *****************************

namespace
{

/// Whether a name list carries the target scheme name
[[nodiscard]] bool contains(const NameList& names, const Name& target)
{
    return std::find(names.begin(), names.end(), target) != names.end();
}

/// Assert two schemes report identical diagonal and source for one tuple
void requireSameContribution
(
    const TimeScheme& a,
    const TimeScheme& b,
    Scalar volume,
    Scalar deltaT,
    Scalar phiPrevStep,
    Scalar ddtPrevStep
)
{
    const TimeContribution ca =
        a.coefficients(volume, deltaT, phiPrevStep, ddtPrevStep);
    const TimeContribution cb =
        b.coefficients(volume, deltaT, phiPrevStep, ddtPrevStep);

    REQUIRE_THAT(ca.diag, WithinRel(cb.diag, TestTolerances::relTight));
    REQUIRE_THAT(ca.source, WithinRel(cb.source, TestTolerances::relTight));
}

} // namespace

// ****************************** Steady State ********************************

TEST_CASE("steadyState contributes nothing", "[schemes]")
{
    const SteadyState scheme;

    const TimeContribution c =
        scheme.coefficients(S(2.0), S(0.5), S(3.0), S(5.0));

    // The null-object scheme adds nothing to either the matrix or the source
    REQUIRE(c.diag == S(0.0));
    REQUIRE(c.source == S(0.0));
    REQUIRE(scheme.isTransient() == false);
}

// ***************************** Implicit Euler ******************************

TEST_CASE("implicitEuler diagonal and source", "[schemes]")
{
    const ImplicitEuler scheme;

    // V/dt = 2/0.5 = 4, and 4 * phi^n = 4 * 3 = 12, both exact
    const TimeContribution c =
        scheme.coefficients(S(2.0), S(0.5), S(3.0), S(0.0));

    REQUIRE(c.diag == S(4.0));
    REQUIRE(c.source == S(12.0));
    REQUIRE(scheme.isTransient() == true);
}

// ***************************** Crank Nicolson ******************************

TEST_CASE("crankNicolson with coeff 0 matches implicitEuler", "[schemes]")
{
    const CrankNicolson crank(S(0.0));
    const ImplicitEuler euler;

    // With CrankNicolsonCoeff = 0 the stored derivative drops out and the
    // scheme degenerates to backward Euler for any input tuple
    requireSameContribution(crank, euler, S(2.0), S(0.5), S(3.0), S(5.0));
    requireSameContribution(crank, euler, S(1.0), S(0.1), S(-2.0), S(7.0));
    requireSameContribution(crank, euler, S(0.25), S(0.2), S(4.0), S(-1.0));

    // Header formula with CrankNicolsonCoeff = 1: coefft = 2, so
    //   diag   = coefft * V/dt = 2 * (2/0.5) = 8
    //   source = diag * phi^n + coeff * ddtPrevStep = 8*3 + 1*5 = 29
    const CrankNicolson crankFull(S(1.0));

    const TimeContribution c =
        crankFull.coefficients(S(2.0), S(0.5), S(3.0), S(5.0));

    const Scalar coefft = S(1.0) + S(1.0);
    const Scalar rDeltaT = coefft * S(2.0) / S(0.5);
    const Scalar expectedDiag = rDeltaT;
    const Scalar expectedSource = rDeltaT * S(3.0) + S(1.0) * S(5.0);

    REQUIRE_THAT(c.diag, WithinRel(expectedDiag, TestTolerances::relTight));
    REQUIRE_THAT
    (
        c.source,
        WithinRel(expectedSource, TestTolerances::relTight)
    );
    REQUIRE(crankFull.isTransient() == true);
}

// ****************************** Runtime Selection ***************************

TEST_CASE("factory lists the three schemes", "[schemes]")
{
    const std::unique_ptr<TimeScheme> scheme =
        TimeScheme::create("implicitEuler");

    REQUIRE(scheme->isTransient() == true);

    const NameList schemes = TimeScheme::availableSchemes();

    REQUIRE(contains(schemes, "steadyState"));
    REQUIRE(contains(schemes, "implicitEuler"));
    REQUIRE(contains(schemes, "CrankNicolson"));
}
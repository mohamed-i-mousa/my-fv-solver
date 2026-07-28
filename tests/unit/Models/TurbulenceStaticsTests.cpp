/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file TurbulenceStaticsTests.cpp
 * @brief Unit tests for the RANS/k-omega SST static inlet calculators
 *****************************************************************************/

// ********************************** Headers *********************************

// External library headers
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

// Project headers
#include "RANS.h"
#include "kOmegaSST.h"
#include "Vector.h"
#include "Scalar.h"
#include "TestTolerances.h"

using Catch::Matchers::WithinRel;

// ********************** Inlet Turbulent Kinetic Energy **********************

TEST_CASE("inletK from turbulence intensity", "[models]")
{
    // uPrime = I * |U| = 0.05 * 10 = 0.5; k = 1.5 * uPrime^2 = 1.5 * 0.25
    const Vector velocity(S(10.0), S(0.0), S(0.0));

    REQUIRE_THAT
    (
        RANS::inletK(velocity, S(0.05)),
        WithinRel(S(0.375), TestTolerances::relTight)
    );

    // A zero velocity drives k to zero, so the max() clamps it to the floor
    // (smallValue = machine epsilon, from Scalar.h)
    const Vector zeroVelocity(S(0.0), S(0.0), S(0.0));

    REQUIRE(RANS::inletK(zeroVelocity, S(0.05)) == smallValue);
}

// ********************** Inlet Specific Dissipation Rate *********************

TEST_CASE("inletOmega from k and diameter", "[models]")
{
    // omega = sqrt(k) / (betaStar^0.25 * lengthScale),
    // lengthScale = 0.07 * hydraulicDiameter
    //
    // Derivation for inletOmega(0.375, 0.1):
    //   lengthScale       = 0.07 * 0.1            = 0.007
    //   betaStar^0.25     = 0.09^0.25             = 0.5477225575051661
    //   denominator       = 0.54772255 * 0.007    = 0.003834057902536163
    //   sqrt(0.375)                               = 0.6123724356957945
    //   omega = 0.6123724356 / 0.003834057902     = 159.71914124998494
    const Scalar expectedOmega = S(159.71914124998494);

    REQUIRE_THAT
    (
        kOmegaSST::inletOmega(S(0.375), S(0.1)),
        WithinRel(expectedOmega, S(1.0e-10))
    );
}
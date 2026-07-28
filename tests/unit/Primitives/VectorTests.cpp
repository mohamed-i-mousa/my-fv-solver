/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file VectorTests.cpp
 * @brief Unit tests for the Vector primitive
 *****************************************************************************/

// ********************************** Headers *********************************

// External library headers
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

// Project headers
#include "Vector.h"
#include "TestTolerances.h"

using Catch::Matchers::WithinRel;

// ******************************* Dot Product ********************************

TEST_CASE("Vector dot product", "[primitives]")
{
    const Vector a(S(1.0), S(2.0), S(3.0));
    const Vector b(S(4.0), S(5.0), S(6.0));

    // 1*4 + 2*5 + 3*6 = 32, exact in floating point
    REQUIRE(dot(a, b) == S(32.0));

    // Commutativity is exact
    REQUIRE(dot(a, b) == dot(b, a));

    // Dot with self is the squared magnitude
    REQUIRE(dot(a, a) == magnitudeSquared(a));
}

// ****************************** Cross Product *******************************

TEST_CASE("Vector cross product", "[primitives]")
{
    const Vector xHat(S(1.0), S(0.0), S(0.0));
    const Vector yHat(S(0.0), S(1.0), S(0.0));
    const Vector zHat(S(0.0), S(0.0), S(1.0));

    REQUIRE(cross(xHat, yHat) == zHat);

    const Vector a(S(1.0), S(2.0), S(3.0));

    // A vector crossed with itself is the zero vector
    REQUIRE(cross(a, a) == Vector(S(0.0), S(0.0), S(0.0)));

    // The cross product is orthogonal to both operands (integer-exact here:
    // cross(a, b) = (-3, 6, -3), dot with a = -3 + 12 - 9 = 0)
    const Vector b(S(4.0), S(5.0), S(6.0));
    const Vector c = cross(a, b);

    REQUIRE(dot(c, a) == S(0.0));
    REQUIRE(dot(c, b) == S(0.0));
}

// ************************ Magnitude and Normalization ***********************

TEST_CASE("Vector magnitude and normalization", "[primitives]")
{
    const Vector v(S(3.0), S(4.0), S(0.0));

    // 3-4-5 triangle: magnitude is exactly 5
    REQUIRE(magnitude(v) == S(5.0));
    REQUIRE(magnitudeSquared(v) == S(25.0));

    const Vector ones(S(1.0), S(1.0), S(1.0));

    REQUIRE_THAT
    (
        magnitude(normalized(ones)),
        WithinRel(S(1.0), TestTolerances::relTight)
    );
}

// ***************************** Arithmetic Operators *************************

TEST_CASE("Vector arithmetic operators", "[primitives]")
{
    const Vector a(S(1.0), S(2.0), S(3.0));
    const Vector b(S(4.0), S(5.0), S(6.0));

    REQUIRE(a + b == Vector(S(5.0), S(7.0), S(9.0)));
    REQUIRE(b - a == Vector(S(3.0), S(3.0), S(3.0)));
    REQUIRE(S(2.0) * a == Vector(S(2.0), S(4.0), S(6.0)));
    REQUIRE(a * S(2.0) == Vector(S(2.0), S(4.0), S(6.0)));
    REQUIRE(b / S(2.0) == Vector(S(2.0), S(2.5), S(3.0)));

    Vector accumulator = a;
    accumulator += b;
    REQUIRE(accumulator == Vector(S(5.0), S(7.0), S(9.0)));

    accumulator -= b;
    REQUIRE(accumulator == a);
}

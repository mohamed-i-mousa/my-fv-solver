/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file TensorTests.cpp
 * @brief Unit tests for the Tensor primitive
 *****************************************************************************/

// ********************************** Headers *********************************

// External library headers
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

// Standard library headers
#include <cmath>

// Project headers
#include "Tensor.h"
#include "Vector.h"
#include "TestTolerances.h"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

// ***************************** Internal Helpers *****************************

namespace
{

/// Assert two tensors are componentwise bit-identical
void requireExactEqual(const Tensor& a, const Tensor& b)
{
    REQUIRE(a.xx() == b.xx());
    REQUIRE(a.xy() == b.xy());
    REQUIRE(a.xz() == b.xz());
    REQUIRE(a.yx() == b.yx());
    REQUIRE(a.yy() == b.yy());
    REQUIRE(a.yz() == b.yz());
    REQUIRE(a.zx() == b.zx());
    REQUIRE(a.zy() == b.zy());
    REQUIRE(a.zz() == b.zz());
}

/// Assert two tensors match componentwise to an absolute tolerance
void requireCloseEqual(const Tensor& a, const Tensor& b, Scalar tol)
{
    REQUIRE_THAT(a.xx(), WithinAbs(b.xx(), tol));
    REQUIRE_THAT(a.xy(), WithinAbs(b.xy(), tol));
    REQUIRE_THAT(a.xz(), WithinAbs(b.xz(), tol));
    REQUIRE_THAT(a.yx(), WithinAbs(b.yx(), tol));
    REQUIRE_THAT(a.yy(), WithinAbs(b.yy(), tol));
    REQUIRE_THAT(a.yz(), WithinAbs(b.yz(), tol));
    REQUIRE_THAT(a.zx(), WithinAbs(b.zx(), tol));
    REQUIRE_THAT(a.zy(), WithinAbs(b.zy(), tol));
    REQUIRE_THAT(a.zz(), WithinAbs(b.zz(), tol));
}

} // namespace

// ************************** Transpose and Symmetry **************************

TEST_CASE("Tensor transpose and symmetric decomposition", "[primitives]")
{
    // An asymmetric integer-valued tensor exercises every off-diagonal term
    const Tensor T
    (
        S(1.0), S(2.0), S(3.0),
        S(4.0), S(5.0), S(6.0),
        S(7.0), S(8.0), S(9.0)
    );

    // Transposing twice returns the original tensor exactly
    requireExactEqual(T.transpose().transpose(), T);

    // The symmetric and antisymmetric parts sum back to the original tensor
    const Tensor reconstructed = T.symm() + T.skew();

    requireCloseEqual(reconstructed, T, TestTolerances::absTight);

    // The skew part is trace-free by construction
    REQUIRE_THAT
    (
        T.skew().trace(),
        WithinAbs(S(0.0), TestTolerances::absTight)
    );
}

// *************************** Outer and Double-Dot ***************************

TEST_CASE("Tensor outer and double-dot identities", "[primitives]")
{
    const Vector a(S(1.0), S(2.0), S(3.0));
    const Vector b(S(4.0), S(5.0), S(6.0));

    // outer(a, b)^T == outer(b, a): swapping the operands transposes result
    requireExactEqual(outer(a, b).transpose(), outer(b, a));

    const Tensor A
    (
        S(1.0), S(2.0), S(3.0),
        S(4.0), S(5.0), S(6.0),
        S(7.0), S(8.0), S(9.0)
    );

    const Tensor B
    (
        S(9.0), S(8.0), S(7.0),
        S(6.0), S(5.0), S(4.0),
        S(3.0), S(2.0), S(1.0)
    );

    // The double-dot product is symmetric in its arguments
    REQUIRE(doubleDot(A, B) == doubleDot(B, A));

    // Frobenius squared magnitude is the tensor double-dotted with itself
    REQUIRE(A.magnitudeSquared() == doubleDot(A, A));
}

// *************************** Strain-Rate Magnitude **************************

TEST_CASE("Strain-rate magnitude of pure shear", "[primitives]")
{
    // A pure-shear velocity gradient: only the xy entry is non-zero
    const Tensor gradU
    (
        S(0.0), S(1.0), S(0.0),
        S(0.0), S(0.0), S(0.0),
        S(0.0), S(0.0), S(0.0)
    );

    // symm(gradU) has xy = yx = 0.5, so magnitudeSquared = 0.5; the
    // strain-rate magnitude sqrt(2 S:S) then evaluates to exactly 1
    const Scalar strainRate =
        std::sqrt(S(2.0) * gradU.symm().magnitudeSquared());

    REQUIRE_THAT(strainRate, WithinRel(S(1.0), TestTolerances::relTight));
}
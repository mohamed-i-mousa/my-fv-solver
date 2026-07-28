/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file TestTolerances.h
 * @brief Named floating-point tolerances shared across the test suite
 *
 * @details Catch2's Catch::Matchers::WithinRel / WithinAbs take a tolerance;
 * these named constants keep the policy in one place. Policy: exact small-
 * integer arithmetic asserts with plain ==; one/two-op FP kernels use
 * relTight; quantities that are legitimately zero use absTight (never a
 * relative tolerance near zero); a discrete operator that must reproduce an
 * analytic field exactly (least-squares gradient of a linear field) uses
 * absOperator; an iteratively solved value uses absSolve against a solver
 * tolerance an order of magnitude tighter. Approx is deliberately unused.
 * Every value is precision-dependent: the FP64 set sits far below FLT_EPSILON,
 * so the FP32 build carries its own.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include "Scalar.h"

// ******************************* Tolerances ********************************

namespace TestTolerances
{

#ifdef PROJECT_USE_DOUBLE_PRECISION

/// One-or-two-op FP kernels (sqrt, divide): relative
constexpr Scalar relTight = S(1.0e-12);

/// Quantities that are legitimately zero: absolute
constexpr Scalar absTight = S(1.0e-12);

/// Discrete-operator exactness on a well-conditioned box
constexpr Scalar absOperator = S(1.0e-10);

/// Iteratively solved value against its analytic answer
constexpr Scalar absSolve = S(1.0e-8);

/// Krylov solver relative tolerance for the tiny assembled systems
constexpr Scalar solverTolerance = S(1.0e-12);

#else

/// One-or-two-op FP kernels (sqrt, divide): relative
constexpr Scalar relTight = S(1.0e-5);

/// Quantities that are legitimately zero: absolute
constexpr Scalar absTight = S(1.0e-6);

/// Discrete-operator exactness on a well-conditioned box
constexpr Scalar absOperator = S(1.0e-4);

/// Iteratively solved value against its analytic answer
constexpr Scalar absSolve = S(1.0e-4);

/// Krylov solver relative tolerance for the tiny assembled systems
constexpr Scalar solverTolerance = S(1.0e-6);

#endif

} // namespace TestTolerances
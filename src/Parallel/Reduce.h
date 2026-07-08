/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file Reduce.h
 * @brief Global reductions: combine per-rank partial values across the run
 *
 * @details Every residual, statistic, and convergence quantity is the
 * reduction of per-cell (or per-face) values. In a parallel run each rank
 * holds only its part of the mesh, so the local partial must be combined
 * across ranks before it means anything.
 *
 * All reductions return the combined value on EVERY rank (MPI_Allreduce),
 * because reduced values feed control flow, convergence tests, Courant
 * limits, and every rank must take the same branch, or the next
 * collective operation deadlocks.
 *
 * On a single rank each function returns its argument unchanged, so the
 * serial solver is bit-identical with and without these calls.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <span>

// Project headers
#include "Integer.h"
#include "Scalar.h"
#include "Vector.h"

// **************************** Global Reductions *****************************

/// Sum a value across every rank
[[nodiscard]] Scalar globalSum(Scalar value);

/// Sum a count across every rank
[[nodiscard]] Count globalSum(Count value);

/// Sum a vector across every rank
[[nodiscard]] Vector globalSum(const Vector& value);

/// Maximum of a value across every rank
[[nodiscard]] Scalar globalMax(Scalar value);

/// Minimum of a value across every rank
[[nodiscard]] Scalar globalMin(Scalar value);

// Batched variants: element-wise reduction of several values through ONE
// collective call, replacing a sequence of latency-bound single-value
// round-trips (results overwrite the input in place)

/// Sum each element across every rank, in place
void globalSum(std::span<Scalar> values);

/// Maximum of each element across every rank, in place
void globalMax(std::span<Scalar> values);

/// Minimum of each element across every rank, in place
void globalMin(std::span<Scalar> values);
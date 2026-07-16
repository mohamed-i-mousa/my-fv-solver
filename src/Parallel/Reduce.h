/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file Reduce.h
 * @brief Global reductions: combine per-rank partial values across the run
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

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

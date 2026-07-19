/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file Logger.h
 * @brief Formatting helpers for solver console output
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <span>
#include <utility>

// Project headers
#include "Scalar.h"
#include "StringTypes.h"
#include "Integer.h"

// ***************************** namespace Logger *****************************

namespace Logger
{
    using Residuals = std::pair<Name, Scalar>;

    /// Print a generic 80-char framed banner with the given title
    void sectionHeader(const Message& title);

    /// Print the per-iteration banner
    void iterationHeader(Count n);

    /// Print the closing rule that terminates a framed block
    void iterationFooter();

    /// Print the column header row for the linear solver configuration table
    void linearSolverConfigHeader();

    /// Print one row of the linear solver configuration table
    void linearSolverConfigRow
    (
        const Name& equation,
        const Name& solver,
        Scalar tolerance,
        Count maxIters
    );

    /// Print one indented label-value row with a Scalar value
    void keyValue(const Message& label, Scalar value);

    /// Print one indented label-value row with a fixed-point value and unit
    void keyValue(const Message& label, Scalar value, const Message& unit);

    /// Print one indented label-value row with an int value
    void keyValue(const Message& label, int value);

    /// Print one indented label-value row with a Count value
    void keyValue(const Message& label, Count value);

    /// Print one indented label-value row with a string value
    void keyValue(const Message& label, const Message& value);

    /// Print the column header row for the table
    void residualTableHeader();

    /// Print one row of the per-iteration residual table
    void residualRow
    (
        const Name& equation,
        const Name& solver,
        int iterations,
        Scalar linearSolverResidual
    );

    /// Print a sub-section title line with two-space indentation
    void subsection(const Message& title);

    /// Print the header for a pressure/friction/total breakdown table
    void breakdownHeader(const Message& cornerLabel);

    /// Print one labelled pressure/friction/total breakdown row
    void breakdownRow
    (
        const Message& label,
        Scalar pressure,
        Scalar friction,
        Scalar total
    );

    /// Print one min/max/mean statistics line for a scalar field
    void scalarStat
    (
        const Name& name,
        Scalar minVal,
        Scalar maxVal,
        Scalar meanVal
    );

    /// Print one labelled scaled-residual line
    void scaledResidual(const Name& name, Scalar value);

    /// Print the non-debug one-line per-iteration residual summary,
    void residualSummary
    (
        Scalar mass,
        Scalar velocity,
        Scalar pressure,
        std::span<const Residuals> residuals = {}
    );

} // namespace Logger

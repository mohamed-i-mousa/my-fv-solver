/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file Logger.cpp
 * @brief Implementations of solver-output formatting helpers
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "Logger.h"

// Standard library headers
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

// ***************************** namespace Logger *****************************

void Logger::sectionHeader(const Message& title)
{
    std::cout
        << "========================================"
        << "========================================"
        << '\n' << " " << title << '\n'
        << "----------------------------------------"
        << "----------------------------------------"
        << '\n';
}


void Logger::iterationHeader(Count n)
{
    sectionHeader("Iteration " + std::to_string(n));
}


void Logger::iterationFooter()
{
    std::cout
        << "========================================"
        << "========================================"
        << '\n';
}


void Logger::residualTableHeader()
{
    std::ostringstream table;

    table
        << "  "
        << std::left  << std::setw(11) << "Equation"
        << std::left  << std::setw(11) << "Solver"
        << std::right << std::setw(5)  << "Iters"
        << "    " << "Linear Solver Residual" << '\n';

    table
        << "  "
        << std::left  << std::setw(11) << "--------"
        << std::left  << std::setw(11) << "--------"
        << std::right << std::setw(5)  << "-----"
        << "    " << "-----------" << '\n';

    std::cout << table.str();
}


void Logger::residualRow
(
    const Name& equation,
    const Name& solver,
    int iterations,
    Scalar linearSolverResidual
)
{
    std::ostringstream row;

    row
        << "  "
        << std::left  << std::setw(11) << equation
        << std::left  << std::setw(11) << solver
        << std::right << std::setw(5)  << iterations
        << "    "
        << std::scientific << std::setprecision(6)
        << linearSolverResidual;

    std::cout << row.str() << '\n';
}


void Logger::subsection(const Message& title)
{
    std::cout << '\n' << "  " << title << '\n';
}


void Logger::breakdownHeader(const Message& cornerLabel)
{
    std::ostringstream header;

    header
        << '\n' << "  "
        << std::left  << std::setw(14) << cornerLabel
        << std::right
        << std::setw(16) << "Pressure"
        << std::setw(16) << "Friction"
        << std::setw(16) << "Total" << '\n';

    header
        << "  "
        << std::left  << std::setw(14) << "----------"
        << std::right
        << std::setw(16) << "--------"
        << std::setw(16) << "--------"
        << std::setw(16) << "-----" << '\n';

    std::cout << header.str();
}


void Logger::breakdownRow
(
    const Message& label,
    Scalar pressure,
    Scalar friction,
    Scalar total
)
{
    std::ostringstream row;

    row
        << "  "
        << std::left  << std::setw(14) << label
        << std::scientific << std::setprecision(6)
        << std::right
        << std::setw(16) << pressure
        << std::setw(16) << friction
        << std::setw(16) << total;

    std::cout << row.str() << '\n';
}


void Logger::scalarStat
(
    const Name& name,
    Scalar minVal,
    Scalar maxVal,
    Scalar meanVal
)
{
    std::ostringstream stat;

    stat
        << "    "
        << std::left << std::setw(7) << name
        << std::scientific << std::setprecision(2)
        << "min="  << minVal
        << "  max="  << maxVal
        << "  mean=" << meanVal;

    std::cout << stat.str() << '\n';
}


void Logger::scaledResidual(const Name& name, Scalar value)
{
    std::ostringstream residual;

    residual
        << "    "
        << std::left << std::setw(10) << name
        << std::scientific << std::setprecision(6)
        << value;

    std::cout << residual.str() << '\n';
}


void Logger::residualSummary
(
    Scalar mass,
    Scalar velocity,
    Scalar pressure,
    std::span<const Residuals> residuals
)
{
    std::ostringstream summary;

    summary
        << " - Mass: " << std::scientific << mass
        << ", Velocity: " << velocity
        << ", Pressure: " << pressure;

    for (const Residuals& residual : residuals)
    {
        summary
            << ", " << residual.first << ": " << residual.second;
    }

    std::cout << summary.str() << '\n';
}


void Logger::linearSolverConfigHeader()
{
    std::ostringstream header;

    header
        << "    "
        << std::left  << std::setw(11) << "Equation"
        << std::left  << std::setw(11) << "Solver"
        << std::left  << std::setw(12) << "Tolerance"
        << std::right << std::setw(13) << "Max Iters"
        << '\n';

    header
        << "    "
        << std::left  << std::setw(11) << "--------"
        << std::left  << std::setw(11) << "--------"
        << std::left  << std::setw(12) << "----------"
        << std::right << std::setw(13) << "---------"
        << '\n';

    std::cout << header.str();
}


void Logger::linearSolverConfigRow
(
    const Name& equation,
    const Name& solver,
    Scalar tolerance,
    Count maxIters
)
{
    std::ostringstream row;

    row
        << "    " << std::left  << std::setw(11) << equation
        << std::left  << std::setw(11) << solver
        << std::scientific << std::setprecision(6) << tolerance
        << std::right << std::setw(13) << maxIters;

    std::cout << row.str() << '\n';
}


void Logger::keyValue(const Message& label, Scalar value)
{
    std::ostringstream row;

    row
        << "    " << std::left << std::setw(24) << label
        << "  " << std::scientific << std::setprecision(6) << value;

    std::cout << row.str() << '\n';
}


void Logger::keyValue(const Message& label, Scalar value, const Message& unit)
{
    std::ostringstream row;

    row
        << "    " << std::left << std::setw(24) << label
        << "  " << std::fixed << std::setprecision(6) << value
        << ' ' << unit;

    std::cout << row.str() << '\n';
}


void Logger::keyValue(const Message& label, int value)
{
    std::ostringstream row;

    row
        << "    " << std::left  << std::setw(24) << label
        << "  " << std::right << std::setw(12) << value;

    std::cout << row.str() << '\n';
}


void Logger::keyValue(const Message& label, Count value)
{
    std::ostringstream row;

    row
        << "    " << std::left  << std::setw(24) << label
        << "  " << std::right << std::setw(12) << value;

    std::cout << row.str() << '\n';
}


void Logger::keyValue(const Message& label, const Message& value)
{
    std::ostringstream row;
    
    row
        << "    " << std::left << std::setw(24) << label
        << "  " << value;

    std::cout << row.str() << '\n';
}

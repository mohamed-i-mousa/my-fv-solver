/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file ErrorHandler.h
 * @brief Fatal error and warning functions for program diagnostics
 *
 * @note FatalError is declared `[[noreturn]] noexcept` and calls
 * std::abort() unconditionally. It can therefore be invoked from any
 * function marked `noexcept` (including accessors and destructors) without
 * risking a silent std::terminate from a propagated exception.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include <cstdlib>
#include <iostream>
#include <source_location>

#include "StringTypes.h"

// *********************************** Alias **********************************

using Location = std::source_location;

// ****************************** Abort Handler *******************************

/// Handler FatalError terminates through: std::abort until the parallel
/// runtime installs an MPI-wide abort at startup, so a fatal error on one
/// rank cannot leave the other ranks deadlocked in a collective
using AbortHandler = void (*)();

/// Access the active abort handler (nullptr selects std::abort)
[[nodiscard]] inline AbortHandler& activeAbortHandler() noexcept
{
    static AbortHandler handler = nullptr;
    return handler;
}

/// Install the abort handler FatalError terminates through
inline void setAbortHandler(AbortHandler handler) noexcept
{
    activeAbortHandler() = handler;
}

// ************************* Error Handling Functions *************************

/// Print a fatal error message and abort the program
[[noreturn]] inline void FatalError
(
    const Message& errorMessage,
    const Location errorLocation = Location::current()
) noexcept
{
    std::cerr
        << '\n' << '\n' << "FATAL ERROR"
        << '\n' << "    " << errorLocation.file_name() << ':'
        << errorLocation.line()
        << '\n' << "    " << errorMessage << '\n' << '\n' << std::endl;

    if (activeAbortHandler() != nullptr)
    {
        activeAbortHandler()();
    }

    std::abort();
}


/// Print a warning message and continue execution
inline void Warning
(
    const Message& warningMessage,
    const Location warningLocation = Location::current()
) noexcept
{
    std::cerr
        << '\n' << "[WARNING] (" << warningLocation.file_name() << ':'
        << warningLocation.line() << ") " << warningMessage << std::endl;
}

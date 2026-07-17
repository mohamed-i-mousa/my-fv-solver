/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file ErrorHandler.h
 * @brief Fatal error and warning functions for program diagnostics
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <cstdlib>
#include <iostream>
#include <source_location>

// External library headers
#include <mpi.h>

// Project headers
#include "StringTypes.h"

// *********************************** Alias **********************************

using Location = std::source_location;

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

    // Abort parallel run by aborting all ranks
    int MPIInitialized = 0;
    int MPIFinalized = 0;
    MPI_Initialized(&MPIInitialized);
    MPI_Finalized(&MPIFinalized);

    if (MPIInitialized != 0 && MPIFinalized == 0)
    {
        MPI_Abort(MPI_COMM_WORLD, 1);
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

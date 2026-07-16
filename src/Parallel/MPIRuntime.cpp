/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file MPIRuntime.cpp
 * @brief RAII ownership of the MPI runtime
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "MPIRuntime.h"

// External library headers
#include <mpi.h>

// Project headers
#include "Comm.h"
#include "ErrorHandler.h"

// ************************* Special Member Functions *************************

MPIRuntime::MPIRuntime()
{
    int providedThreadSupport = 0;

    MPI_Init_thread
    (
        nullptr,
        nullptr,
        MPI_THREAD_FUNNELED,
        &providedThreadSupport
    );

    Comm::init();

    // OpenMP threads compute within each rank, but only the main thread
    // makes MPI calls - anything below funneled support cannot honor that
    if (providedThreadSupport < MPI_THREAD_FUNNELED && Comm::master())
    {
        Warning("MPI runtime provides no funneled-thread support");
    }
}


MPIRuntime::~MPIRuntime() noexcept
{
    MPI_Finalize();
}

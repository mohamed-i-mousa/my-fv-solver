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

// ************************* Special Member Functions *************************

MPIRuntime::MPIRuntime()
{
    MPI_Init(nullptr, nullptr);

    Comm::init();
}


MPIRuntime::~MPIRuntime() noexcept
{
    MPI_Finalize();
}

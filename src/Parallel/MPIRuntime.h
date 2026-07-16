/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file MPIRuntime.h
 * @brief RAII ownership of the MPI runtime
 *
 * @details One MPIRuntime is created at the very top of main(), before
 * PETScRuntime and every other object: MPI is the process-wide parallel
 * substrate that mesh decomposition, halo exchanges, reductions, and
 * collective HDF5 output all stand on, so its lifetime belongs to
 * src/Parallel/, not to any one consumer. PetscInitialize finds MPI
 * already initialized and attaches to it as a guest; PetscFinalize then
 * leaves MPI running, and the reverse destruction order in main()
 * guarantees this destructor's MPI_Finalize runs last.
 *
 * The constructor requests MPI_THREAD_FUNNELED: OpenMP threads compute
 * within each rank, but only the main thread makes MPI calls.
 *
 * @class MPIRuntime
 * - Constructor runs MPI_Init_thread, then Comm::init()
 * - Destructor runs MPI_Finalize
 *****************************************************************************/

#pragma once

// **************************** class MPIRuntime ******************************

class MPIRuntime
{
public:

// ************************* Special Member Functions *************************

    /// Initialize the MPI runtime and cache this rank's identity
    MPIRuntime();

    /// Copy constructor and assignment - Not copyable (owns the runtime)
    MPIRuntime(const MPIRuntime&) = delete;
    MPIRuntime& operator=(const MPIRuntime&) = delete;

    /// Move constructor and assignment - Not movable
    MPIRuntime(MPIRuntime&&) = delete;
    MPIRuntime& operator=(MPIRuntime&&) = delete;

    /// Finalize the MPI runtime
    ~MPIRuntime() noexcept;
};
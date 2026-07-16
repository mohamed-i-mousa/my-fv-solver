/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file Comm.h
 * @brief Rank identity of this process within the MPI run
 *
 * @details Comm answers the four questions parallel-aware code may ask —
 * am I in a parallel run, how many ranks, which one am I, am I the
 * master, without exposing <mpi.h> to its callers. The MPI runtime
 * itself is owned by MPIRuntime, whose constructor calls init() so the
 * cached identity is valid before any other Comm query.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include "Integer.h"

// ***************************** namespace Comm *******************************

namespace Comm
{
    /// Initialize the cached MPI_COMM_WORLD
    void init();

    /// True when the run spans more than one MPI rank
    [[nodiscard]] bool parallelRun();

    /// Number of MPI ranks in the run
    [[nodiscard]] Count numProcessors();

    /// This processor's rank
    [[nodiscard]] Index myProcessorNum();

    /// True on the master rank (rank 0)
    [[nodiscard]] bool master();

    /// Abort every rank of the run
    void abortAllRanks();

} // namespace Comm

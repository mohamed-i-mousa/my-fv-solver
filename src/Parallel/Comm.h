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
 * itself is owned by PETScRuntime (PetscInitialize starts MPI), which is
 * constructed in main() before any Comm query.
 *
 * abortAllRanks() is installed as FatalError's abort handler at startup:
 * a rank that dies alone would leave the other ranks deadlocked inside
 * their next collective, so a fatal error on any rank must bring down
 * the whole run.
 *  
 * Rank and size never change after MPI_Init, so both are queried once on
 * first use (after PETScRuntime has initialized MPI) and cached
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include "Integer.h"

// ***************************** namespace Comm *******************************

namespace Comm
{
    /// True when the run spans more than one MPI rank
    [[nodiscard]] bool parallelRun();

    /// Number of MPI ranks in the run
    [[nodiscard]] Count nProcs();

    /// This process's rank
    [[nodiscard]] Index myProcNo();

    /// True on the master rank (rank 0)
    [[nodiscard]] bool master();

    /// Abort every rank of the run
    void abortAllRanks();

} // namespace Comm
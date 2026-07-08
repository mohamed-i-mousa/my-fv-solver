/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file GlobalIndex.cpp
 * @brief Rank-major global numbering built from the per-rank counts
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "GlobalIndex.h"

// External library headers
#include <mpi.h>

// Project headers
#include "Comm.h"

// **************************** Build-Time Invariants **************************

// Count buffers are handed to MPI as MPI_UINT64_T
static_assert
(
    sizeof(Count) == 8,
    "GlobalIndex reductions assume a 64-bit Count"
);

// ************************* Special Member Functions *************************

GlobalIndex::GlobalIndex(Count localCount)
{
    if (!Comm::parallelRun())
    {
        offset_ = 0;
        totalCount_ = localCount;
        return;
    }

    Count exclusivePrefix = 0;

    MPI_Exscan
    (
        &localCount,
        &exclusivePrefix,
        1,
        MPI_UINT64_T,
        MPI_SUM,
        MPI_COMM_WORLD
    );

    offset_ = Comm::master() ? 0 : static_cast<Index>(exclusivePrefix);

    MPI_Allreduce
    (
        &localCount,
        &totalCount_,
        1,
        MPI_UINT64_T,
        MPI_SUM,
        MPI_COMM_WORLD
    );
}
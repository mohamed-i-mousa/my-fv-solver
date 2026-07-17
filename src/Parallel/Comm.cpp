/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file Comm.cpp
 * @brief Rank identity queries over the world communicator
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "Comm.h"

// External library headers
#include <mpi.h>

// ****************************** Internal Helpers ****************************

namespace
{

// Cached MPI_COMM_WORLD identity
int size = 1;
int rank = 0;

} // namespace

// ***************************** namespace Comm *******************************

void Comm::init()
{
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
}


bool Comm::parallelRun()
{
    return size > 1;
}


Count Comm::numProcessors()
{
    return static_cast<Count>(size);
}


Index Comm::myProcessorNum()
{
    return static_cast<Index>(rank);
}


bool Comm::master()
{
    return rank == 0;
}

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

int cachedNProcs()
{
    static int nProcs = 0;

    if (nProcs == 0)
    {
        MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    }

    return nProcs;
}


int cachedMyProcNo()
{
    static int myProcNo = -1;

    if (myProcNo == -1)
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &myProcNo);
    }

    return myProcNo;
}

} // namespace

// ***************************** namespace Comm *******************************

bool Comm::parallelRun()
{
    return cachedNProcs() > 1;
}


Count Comm::nProcs()
{
    return static_cast<Count>(cachedNProcs());
}


Index Comm::myProcNo()
{
    return static_cast<Index>(cachedMyProcNo());
}


bool Comm::master()
{
    return cachedMyProcNo() == 0;
}


void Comm::abortAllRanks()
{
    MPI_Abort(MPI_COMM_WORLD, 1);
}
/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file Reduce.cpp
 * @brief Global reductions over the world communicator
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "Reduce.h"

// External library headers
#include <mpi.h>

// Project headers
#include "MPIScalarType.h"
#include "Comm.h"

// ****************************** Internal Helpers ****************************

namespace
{

static_assert
(
    sizeof(Count) == 8,
    "Count reductions assume a 64-bit Count"
);

} // namespace

// **************************** Global Reductions *****************************

Scalar globalSum(Scalar value)
{
    if (!Comm::parallelRun())
    {
        return value;
    }

    Scalar result = S(0.0);
    MPI_Allreduce
    (
        &value,
        &result,
        1,
        MPIScalarType(),
        MPI_SUM,
        MPI_COMM_WORLD
    );

    return result;
}


Count globalSum(Count value)
{
    if (!Comm::parallelRun())
    {
        return value;
    }

    Count result = 0;
    MPI_Allreduce(&value, &result, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);

    return result;
}


Vector globalSum(const Vector& value)
{
    if (!Comm::parallelRun())
    {
        return value;
    }

    const Scalar partial[3] = {value.x(), value.y(), value.z()};
    Scalar result[3] = {S(0.0), S(0.0), S(0.0)};
    MPI_Allreduce
    (
        partial,
        result,
        3,
        MPIScalarType(),
        MPI_SUM,
        MPI_COMM_WORLD
    );

    return {result[0], result[1], result[2]};
}


Scalar globalMax(Scalar value)
{
    if (!Comm::parallelRun())
    {
        return value;
    }

    Scalar result = S(0.0);
    MPI_Allreduce
    (
        &value,
        &result,
        1,
        MPIScalarType(),
        MPI_MAX,
        MPI_COMM_WORLD
    );

    return result;
}


Scalar globalMin(Scalar value)
{
    if (!Comm::parallelRun())
    {
        return value;
    }

    Scalar result = S(0.0);
    MPI_Allreduce
    (
        &value,
        &result,
        1,
        MPIScalarType(),
        MPI_MIN,
        MPI_COMM_WORLD
    );

    return result;
}


bool globalOr(bool value)
{
    if (!Comm::parallelRun())
    {
        return value;
    }

    const int localValue = value ? 1 : 0;
    int result = 0;
    MPI_Allreduce
    (
        &localValue,
        &result,
        1,
        MPI_INT,
        MPI_LOR,
        MPI_COMM_WORLD
    );

    return result != 0;
}

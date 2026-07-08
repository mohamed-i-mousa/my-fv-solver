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
#include "Comm.h"

// ****************************** Internal Helpers ****************************

namespace
{

MPI_Datatype mpiScalarType()
{
    return sizeof(Scalar) == 8 ? MPI_DOUBLE : MPI_FLOAT;
}

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
        mpiScalarType(),
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
        mpiScalarType(),
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
        mpiScalarType(),
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
        mpiScalarType(),
        MPI_MIN,
        MPI_COMM_WORLD
    );

    return result;
}

// ************************* Batched Global Reductions ************************

void globalSum(std::span<Scalar> values)
{
    if (!Comm::parallelRun())
    {
        return;
    }

    MPI_Allreduce
    (
        MPI_IN_PLACE,
        values.data(),
        static_cast<int>(values.size()),
        mpiScalarType(),
        MPI_SUM,
        MPI_COMM_WORLD
    );
}


void globalMax(std::span<Scalar> values)
{
    if (!Comm::parallelRun())
    {
        return;
    }

    MPI_Allreduce
    (
        MPI_IN_PLACE,
        values.data(),
        static_cast<int>(values.size()),
        mpiScalarType(),
        MPI_MAX,
        MPI_COMM_WORLD
    );
}


void globalMin(std::span<Scalar> values)
{
    if (!Comm::parallelRun())
    {
        return;
    }

    MPI_Allreduce
    (
        MPI_IN_PLACE,
        values.data(),
        static_cast<int>(values.size()),
        mpiScalarType(),
        MPI_MIN,
        MPI_COMM_WORLD
    );
}
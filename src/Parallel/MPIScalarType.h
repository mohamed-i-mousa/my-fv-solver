/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file MPIScalarType.h
 * @brief The MPI datatypes matching project primitive types
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// External library headers
#include <mpi.h>

// Project headers
#include "Scalar.h"

// **************************** Build-Time Invariants *************************

static_assert
(
    sizeof(Count) == 8,
    "Count operations assume a 64-bit Count (MPI_UINT64_T)"
);

// ****************************** MPI Data Types ******************************

/// MPI datatype of one Scalar
[[nodiscard]] inline MPI_Datatype MPIScalarType()
{
    return sizeof(Scalar) == 8 ? MPI_DOUBLE : MPI_FLOAT;
}

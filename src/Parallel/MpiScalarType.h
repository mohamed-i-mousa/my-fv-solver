/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file MpiScalarType.h
 * @brief The MPI datatype matching Scalar's build-time precision
 *
 * @details One definition shared by every src/Parallel implementation
 * file: a precision change edited into only one hand-written copy would
 * silently corrupt MPI buffers. Includes <mpi.h> — for use inside
 * src/Parallel .cpp files only.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// External library headers
#include <mpi.h>

// Project headers
#include "Scalar.h"

// ****************************** MPI Scalar Type *****************************

/// MPI datatype of one Scalar
[[nodiscard]] inline MPI_Datatype mpiScalarType()
{
    return sizeof(Scalar) == 8 ? MPI_DOUBLE : MPI_FLOAT;
}

/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file HaloExchange.h
 * @brief Update ghost cells from their owning ranks
 *
 * @details The field-level halo: after a rank changes owned-cell values
 * that face kernels or stencils will read across an inter-rank cut, the
 * ghost copies on the neighboring ranks are stale until exchanged. One
 * exchangeHalos call packs this rank's send cells per neighbor, ships
 * one message per neighbor (non-blocking, all cuts in flight at once),
 * and fills the ghost tail with the received values.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <initializer_list>

// Project headers
#include "CellData.h"
#include "Mesh.h"

// ***************************** Halo Exchange ********************************

/// Update ghost cells of several scalar fields
void exchangeHalos
(
    const Mesh& mesh,
    std::initializer_list<ScalarField*> fields
);

/// Update ghost cells of one vector field
void exchangeHalos
(
    const Mesh& mesh,
    VectorField& field
);

/// Update ghost cells of several vector fields (one message per neighbor)
void exchangeHalos
(
    const Mesh& mesh,
    std::initializer_list<VectorField*> fields
);

/// Update ghost cells of one tensor field
void exchangeHalos
(
    const Mesh& mesh,
    TensorField& field
);

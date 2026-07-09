/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file DecompositionChecker.h
 * @brief Collective validation of a freshly distributed mesh
 *
 * @details The geometry gate of the decomposition: every rank exchanges
 * its side of each inter-rank cut with the neighbor and asserts, face by
 * face and ghost by ghost, that both sides describe the same geometry —
 * BIT-identical, since cut faces are shipped with unchanged node data and
 * ghost geometry originates from the same complete-mesh computation the
 * owning rank reproduces locally. Also verifies the owned-cell counts sum
 * to the complete mesh and that a least-squares gradient of an analytic
 * linear field is recovered exactly across the cuts (ghost stencils).
 *
 * Every check FatalErrors on failure (aborting all ranks); success prints
 * a one-block summary on the master. Collective: all ranks must call.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Project headers
#include "Integer.h"
#include "Mesh.h"

// ********************** namespace DecompositionChecker **********************

namespace DecompositionChecker
{
    /// Run every decomposition check (collective)
    void check
    (
        const Mesh& mesh,
        Count totalCellCount
    );

} // namespace DecompositionChecker

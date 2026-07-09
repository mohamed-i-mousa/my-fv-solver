/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file MeshDistributor.h
 * @brief Ship each rank its submesh block from the partitioning rank
 *
 * @details distribute() is collective: the master calls it holding one
 * SubmeshData per rank (from MeshDecomposer::decompose), every other
 * rank calls it with an empty list, and each rank returns owning its own
 * block. Blocks travel as one contiguous byte message each (a length
 * message followed by the payload), serialized field-by-field — same-
 * architecture ranks are assumed, as everywhere in this solver.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <vector>

// Project headers
#include "SubmeshData.h"

// ************************* namespace MeshDistributor ************************

namespace MeshDistributor
{
    /// Collective: scatter the blocks; every rank returns its own
    [[nodiscard]] SubmeshData distribute(std::vector<SubmeshData> blocks);

} // namespace MeshDistributor

/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file MeshDecomposer.h
 * @brief METIS partitioning of the complete mesh into per-rank submeshes
 *
 * @details Runs on the partitioning rank only, on the complete, fully
 * prepared mesh (geometry computed). partition() labels every cell with
 * a rank through METIS_PartGraphKway on the cell dual graph (cells are
 * vertices, internal faces are edges); decompose() extracts one flat
 * SubmeshData per rank: owned cells, every face they touch, the node
 * subset, physical patch fragments, processor patches for the inter-rank
 * cuts, and a one-cell ghost layer carrying centroid/volume/global id.
 *
 * Cut faces keep their ORIGINAL orientation and owner/neighbor roles on
 * both ranks — one side maps the remote cell to a ghost index. Both
 * copies are therefore built from identical node data, making their
 * locally recomputed geometry bit-identical across the cut (the
 * decomposition checker asserts exactly that). Global cell numbering is
 * rank-major: rank r's cells follow all cells of ranks 0..r-1, in
 * original-mesh order within the rank. This class is serial and MPI-free
 * by design — MeshDistributor does the shipping.
 *
 * @class MeshDecomposer
 * - partition(): cell -> rank labels from METIS (all zeros for 1 part)
 * - decompose(): per-rank SubmeshData blocks, ready to distribute
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <vector>

// Project headers
#include "Integer.h"
#include "Mesh.h"
#include "SubmeshData.h"

// *************************** class MeshDecomposer ***************************

class MeshDecomposer
{
public:

// ************************* Special Member Functions *************************

    /// Constructor
    MeshDecomposer
    (
        const Mesh& mesh,
        Count numParts
    );

    /// Copy constructor and assignment - Not copyable (const T& members)
    MeshDecomposer(const MeshDecomposer&) = delete;
    MeshDecomposer& operator=(const MeshDecomposer&) = delete;

    /// Move constructor and assignment - Not movable (const T& members)
    MeshDecomposer(MeshDecomposer&&) = delete;
    MeshDecomposer& operator=(MeshDecomposer&&) = delete;

    /// Destructor
    ~MeshDecomposer() noexcept = default;

// ****************************** Public Methods ******************************

    /// Label every cell with its rank
    [[nodiscard]] IndexList partition() const;

    /// Extract a submesh block per rank
    [[nodiscard]] std::vector<SubmeshData> decompose() const;

// ****************************** Private Members *****************************

private:

    /// Complete mesh reference
    const Mesh& mesh_;

    /// Number of parts to split into
    Count numParts_;

// ****************************** Private Methods *****************************

    /// Extract the submesh block of one rank
    [[nodiscard]] SubmeshData extractSubmesh
    (
        Index rank,
        const IndexList& cellToRank,
        const IndexList& localCellIdx,
        const IndexList& rankOffsets,
        const IndexList& faceToPatch
    ) const;
};
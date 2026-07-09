/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file ProcessorPatch.h
 * @brief Per-neighbor metadata of one inter-rank cut in a decomposed mesh
 *
 * @details A ProcessorPatch describes this rank's side of the cut shared
 * with one neighboring rank: the contiguous range of local cut faces, the
 * owned cells whose values must be sent to fill the neighbor's ghosts,
 * and the contiguous slice of this rank's ghost-cell tail that the
 * neighbor's sends fill in return.
 *
 * Ordering contract: sendCellIndices() is ordered so that packing it and
 * sending it verbatim lands element i in the neighbor's ghost cell i of
 * its matching patch — both sides derive the order from the cut faces
 * sorted by their original (undecomposed) face index, so no index maps
 * ever cross the wire. Pure data, no MPI.
 *
 * @class ProcessorPatch
 * - Identity: this rank's neighbor across the cut
 * - Topology: local cut-face range, send-cell list, ghost-tail slice
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

// Standard library headers
#include <utility>
#include <vector>

// Project headers
#include "Integer.h"

// *************************** class ProcessorPatch ***************************

class ProcessorPatch
{
public:

// ************************* Special Member Functions *************************

    /// Construct one side of an inter-rank cut
    ProcessorPatch
    (
        Index neighborRank,
        Index firstFaceIdx,
        Index lastFaceIdx,
        IndexList sendCellIndices,
        Index ghostFirstCell,
        Count ghostCellCount
    )
    :
        neighborRank_(neighborRank),
        firstFaceIdx_(firstFaceIdx),
        lastFaceIdx_(lastFaceIdx),
        sendCellIndices_(std::move(sendCellIndices)),
        ghostFirstCell_(ghostFirstCell),
        ghostCellCount_(ghostCellCount)
    {}

// ***************************** Accessor Methods *****************************

    /// Rank on the other side of this cut
    [[nodiscard]] Index neighborRank() const noexcept
    {
        return neighborRank_;
    }

    /// First local cut-face index (cut faces are contiguous)
    [[nodiscard]] Index firstFaceIdx() const noexcept
    {
        return firstFaceIdx_;
    }

    /// Last local cut-face index
    [[nodiscard]] Index lastFaceIdx() const noexcept
    {
        return lastFaceIdx_;
    }

    /// Number of cut faces shared with the neighbor
    [[nodiscard]] Count numFaces() const noexcept
    {
        return lastFaceIdx_ - firstFaceIdx_ + 1;
    }

    /// Owned cells to send, in the neighbor's ghost order
    [[nodiscard]] IndexListRef sendCellIndices() const noexcept
    {
        return sendCellIndices_;
    }

    /// Local index of the first ghost cell this patch fills
    [[nodiscard]] Index ghostFirstCell() const noexcept
    {
        return ghostFirstCell_;
    }

    /// Number of ghost cells this patch fills
    [[nodiscard]] Count ghostCellCount() const noexcept
    {
        return ghostCellCount_;
    }

// ****************************** Private Members *****************************

private:

    /// Rank on the other side of this cut
    Index neighborRank_;

    /// First local cut-face index
    Index firstFaceIdx_;

    /// Last local cut-face index
    Index lastFaceIdx_;

    /// Owned cells to send, in the neighbor's ghost order
    IndexList sendCellIndices_;

    /// Local index of the first ghost cell this patch fills
    Index ghostFirstCell_;

    /// Number of ghost cells this patch fills
    Count ghostCellCount_;
};

// ********************************** Aliases *********************************

/// An ordered list of processor patches, one per neighboring rank
using ProcessorPatchList = std::vector<ProcessorPatch>;

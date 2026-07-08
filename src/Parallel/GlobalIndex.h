/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file GlobalIndex.h
 * @brief Rank-major mapping from locally-owned indices to global indices
 *
 * @details In a distributed mesh every rank numbers its owned items
 * 0..localCount-1, while globally the items of rank r follow all items of
 * ranks 0..r-1 (rank-major ordering). A GlobalIndex captures this rank's
 * slice of that ordering: offset() is the global index of local item 0,
 * and toGlobal() shifts any owned local index into global numbering.
 *
 * @class GlobalIndex
 * - Constructor combines the per-rank counts (collective at np > 1)
 * - toGlobal() maps owned local indices only; ghost indices belong to
 *   other ranks' slices and get their own lookup path (Phase 5)
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include "Integer.h"

// **************************** class GlobalIndex *****************************

class GlobalIndex
{
public:

// ************************* Special Member Functions *************************

    /// Build this rank's slice of the rank-major global numbering
    explicit GlobalIndex(Count localCount);

// ***************************** Accessor Methods *****************************

    /// Global index of this rank's local item 0
    [[nodiscard]] Index offset() const noexcept
    {
        return offset_;
    }

    /// Total item count across every rank
    [[nodiscard]] Count totalCount() const noexcept
    {
        return totalCount_;
    }

    /// Map a locally-owned index to its global index
    [[nodiscard]] Index toGlobal(Index localIdx) const noexcept
    {
        return offset_ + localIdx;
    }

// ****************************** Private Members *****************************

private:

    /// Global index of this rank's first owned item
    Index offset_ = 0;

    /// Total item count across every rank
    Count totalCount_ = 0;
};
/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file LinearInterpolation.h
 * @brief Linear interpolation of a cell-centered field to internal faces
 *
 * @details Distance-weighted linear interpolation of CellData<T> from owner
 * and neighbour cell centres to the shared face. For boundary-face values,
 * use the patch's BoundaryType::faceValue before calling this function.
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include "Face.h"
#include "CellData.h"
#include "ErrorHandler.h"

// ************************** Interpolation Functions *************************

/// Linear interpolation of a cell-centered field to an internal face
template<CellFieldType T>
[[nodiscard]] T interpolateToFace
(
    const Face& targetFace,
    const CellData<T>& field
)
{
    if (targetFace.isBoundary())
    {
        FatalError
        (
            "interpolateToFace must not be called on boundary "
            "faces; resolve BC values at the call site"
        );
    }

    const Index P = targetFace.ownerCell();
    const Index N = targetFace.neighborCell().value();

    // Distance weight for the neighbour cell: wN = dP / (dP + dN)
    const Scalar dP = targetFace.dPfMag();
    const Scalar dN = targetFace.dNfMag().value();
    const Scalar wN = dP / (dP + dN);

    return (S(1.0) - wN) * field[P] + wN * field[N];
}

/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file Cell.cpp
 * @brief Implementation of cell geometric properties and operations
 *****************************************************************************/

// ********************************** Headers *********************************

// Implementation header
#include "Cell.h"

// Standard library headers
#include <cmath>

// Project headers
#include "ErrorHandler.h"

// ************************ Geometric Property Methods ************************

void Cell::geometricProperties
(
    const std::vector<FaceIntegrals>& allFaceIntegrals
)
{
    volume_ = S(0.0);
    centroid_ = Vector{};
    Vector centroidSum;

    for (Index faceIdx = 0; faceIdx < faceIndices_.size(); ++faceIdx)
    {
        const Index faceIndex = faceIndices_[faceIdx];
        const Scalar faceSign = S(faceSigns_[faceIdx]);
        const FaceIntegrals& integrals = allFaceIntegrals[faceIndex];

        volume_ += faceSign * integrals.volume;

        centroidSum += Vector
        (
            faceSign * integrals.x2,
            faceSign * integrals.y2,
            faceSign * integrals.z2
        );
    }

    volume_ /= S(3.0);

    if (std::abs(volume_) < vSmallValue)
    {
        FatalError
        (
            "Cell " + std::to_string(idx_) + " has zero volume"
        );
    }

    centroid_ = centroidSum / (S(2.0) * volume_);
}

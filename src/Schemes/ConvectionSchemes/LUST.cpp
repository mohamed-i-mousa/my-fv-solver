/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file LUST.cpp
 * @brief Implementation of the LUST convection scheme
 *****************************************************************************/

// ********************************** Headers *********************************

#include "LUST.h"
#include "LinearInterpolation.h"

// ****************************** Public Methods ******************************

Scalar LUST::correction
(
    const Face& face,
    const ScalarField& phi,
    const Vector& gradPhiP,
    const Vector& gradPhiN,
    Scalar flowRate
) const
{
    // 1. Central Difference deferred correction:
    const Scalar phiFaceCentral = interpolateToFace(face, phi);

    const Index upwindCell =
        (flowRate >= S(0.0)) ? face.ownerCell() : face.neighborCell().value();

    const Scalar phiFaceUDS = phi[upwindCell];
    const Scalar corrCDS = flowRate * (phiFaceCentral - phiFaceUDS);

    // 2. Second-Order Linear Upwind deferred correction:
    const Scalar gradientProjection =
        (flowRate >= S(0.0))
      ? dot(gradPhiP, face.dPf())
      : dot(gradPhiN, face.dNf().value());

    const Scalar corrLU = flowRate * gradientProjection;

    // 3. Blend: alpha * CDS + (1 - alpha) * LUD
    return alpha_ * corrCDS + (S(1.0) - alpha_) * corrLU;
}

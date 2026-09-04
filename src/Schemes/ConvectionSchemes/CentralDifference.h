/******************************************************************************

                                     Turblyze
                           3D incompressible CFD solver
                       Copyright (C) 2025-2026 Mohamed Mousa
                        SPDX-License-Identifier: Apache-2.0

 ------------------------------------------------------------------------------
 * @file CentralDifference.h
 * @brief Central difference convection scheme
 *
 * @details Central difference adds a deferred correction from the implicit
 * upwind face value to the linearly interpolated face value.
 *
 * @class CentralDifference
 * - Second-order central difference convection discretization
 * - Deferred correction from implicit upwind to linearly interpolated face
 *   values
 *****************************************************************************/

#pragma once

// ********************************** Headers *********************************

#include "ConvectionScheme.h"

// ************************** class CentralDifference *************************

class CentralDifference final : public ConvectionScheme
{
public:

// ****************************** Public Methods ******************************

    [[nodiscard]] Scalar correction
    (
        const Face& face,
        const ScalarField& phi,
        const Vector& gradPhiP,
        const Vector& gradPhiN,
        Scalar flowRate
    ) const override;
};
